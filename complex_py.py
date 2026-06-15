#!/usr/bin/env python3
"""
ComPlEx re-implementation in Python/numpy.

Algorithm (mirrors ComPlEx.R exactly):
  1. Load expression matrices (VST, featureSelect-filtered)
  2. Filter to genes with orthologs in both species
  3. Pearson correlation → N×N matrix per species
  4. Mutual Rank normalisation: MR[i,j] = sqrt(rank_i(j) × rank_j(i))
  5. Density threshold: keep top `density_thr` fraction of MR scores
  6. For each ortholog pair (S1→S2 and S2→S1):
       - neighbourhood = genes above threshold
       - map S2 neighbourhood back to S1 via orthologs
       - hypergeometric p-value on overlap
  7. BH FDR correction, max p-value per pair
  8. Filter to FDR < 0.05 and save as TSV

Speed vs R:
  - float32 BLAS matmul for correlation: ~50× faster than R cor()
  - scipy.stats.rankdata: vectorised row-wise ranking
  - Avoids repeated R memory allocations on N×N matrices

Usage:
  python3 complex_py.py \
    --s1-expr  SC_expression.txt \
    --s2-expr  PC_expression.txt \
    --orthologs doc/genes_ortholog_categories.tsv \
    --s1-name  spruce --s2-name pine \
    --out-dir  /path/to/output \
    [--density 0.03] [--workers 8] [--cor-method pearson]
"""

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, hypergeom
from statsmodels.stats.multitest import multipletests

# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="ComPlEx co-expressolog analysis")
    p.add_argument("--s1-expr",    required=True)
    p.add_argument("--s2-expr",    required=True)
    p.add_argument("--orthologs",  required=True,
                   help="TSV with columns: gene, species, Ortholog_Group")
    p.add_argument("--s1-name",    default="species1")
    p.add_argument("--s2-name",    default="species2")
    p.add_argument("--out-dir",    required=True)
    p.add_argument("--density",    type=float, default=0.03,
                   help="Fraction of top MR edges to keep (default 0.03)")
    p.add_argument("--cor-method", default="pearson",
                   choices=["pearson", "spearman"])
    p.add_argument("--min-expr",   type=float, default=1.0,
                   help="VST threshold for featureSelect filter")
    p.add_argument("--min-samples",type=int,   default=2,
                   help="Min samples >= min-expr for featureSelect")
    # ── Co-expressolog candidate / FDR-pool definition ────────────────────────
    # The set of ortholog pairs entering Benjamini-Hochberg correction (the FDR
    # denominator). This does NOT change which pairs can be CALLED co-expressologs
    # (a call requires Max.p < threshold, which already implies overlap >= 2 in
    # both directions, since overlap <= 1 yields p = 1). It only sets how many
    # tests the FDR correction divides over, and therefore how conservative the
    # calls are.
    #   --min-overlap 1 --overlap-mode and  -> reproduces the published R ComPlEx
    #        (Netotea et al. 2014): pool = pairs with overlap > 0 in BOTH
    #        directions. Larger pool, more conservative (the default/headline).
    #   --min-overlap 2 --overlap-mode or   -> the earlier optimisation that drops
    #        x=1 pairs before FDR: smaller pool, less conservative (sensitivity).
    p.add_argument("--min-overlap", type=int, default=1,
                   help="Minimum neighbourhood overlap (x) for a pair to enter "
                        "the FDR pool (default 1, = published R behaviour)")
    p.add_argument("--overlap-mode", default="and", choices=["and", "or"],
                   help="Require min-overlap in BOTH directions ('and', "
                        "= published R) or EITHER direction ('or'). Default 'and'.")
    return p.parse_args()


# ── Load expression ───────────────────────────────────────────────────────────
def load_expr(path, min_expr, min_samples):
    df = pd.read_csv(path, sep="\t", index_col=0)
    # featureSelect: keep genes expressed in >= min_samples samples
    sel = (df >= min_expr).sum(axis=1) >= min_samples
    df  = df.loc[sel]
    # remove zero-variance genes
    df  = df.loc[df.std(axis=1) > 0]
    print(f"  {Path(path).name}: {len(df):,} genes × {df.shape[1]} samples "
          f"after filtering", flush=True)
    return df


# ── Build orthogroup table ────────────────────────────────────────────────────
def build_ortho(ortho_path, s1_name, s2_name, s1_genes, s2_genes):
    og = pd.read_csv(ortho_path, sep="\t")
    # Drop unassigned-gene placeholders: merging on the group label would otherwise
    # join every unassigned gene to every other, producing a cartesian explosion of
    # spurious "ortholog pairs". The published R input (parquet) contains no such rows.
    before = len(og)
    og = og[~og["Ortholog_Group"].isin(["No orthogroup", "No_orthogroup", ""])]
    og = og[og["Ortholog_Group"].notna()]
    if len(og) < before:
        print(f"  dropped {before - len(og):,} rows without an orthogroup", flush=True)
    rename = {"Picea_abies": s1_name, "Pinus_sylvestris": s2_name}
    og["species"] = og["species"].replace(rename)

    s1 = og[og["species"] == s1_name][["gene","Ortholog_Group"]].rename(
            columns={"gene": "Species1"})
    s2 = og[og["species"] == s2_name][["gene","Ortholog_Group"]].rename(
            columns={"gene": "Species2"})

    pairs = s1.merge(s2, on="Ortholog_Group", how="inner")
    # keep only genes present in expression matrices
    pairs = pairs[pairs["Species1"].isin(s1_genes) &
                  pairs["Species2"].isin(s2_genes)]
    print(f"  {len(pairs):,} ortholog pairs "
          f"({pairs['Ortholog_Group'].nunique():,} orthogroups)", flush=True)
    return pairs.reset_index(drop=True)


# ── Mutual Rank ───────────────────────────────────────────────────────────────
def mutual_rank(expr_df, method="pearson"):
    """Return MR matrix (genes × genes) as float32 numpy array.

    Uses float32 throughout to halve peak memory vs np.corrcoef (float64).
    Ranking is done in chunks of 500 rows to bound the float64 intermediate.
    """
    N   = len(expr_df)
    mat = expr_df.values.astype(np.float32)

    if method == "spearman":
        mat = np.apply_along_axis(rankdata, 1, mat).astype(np.float32)

    t0 = time.time()
    # Manual float32 Pearson correlation — np.corrcoef uses float64 (2× memory)
    mat -= mat.mean(axis=1, keepdims=True)
    norms = np.sqrt((mat ** 2).sum(axis=1, keepdims=True))
    norms[norms < 1e-10] = 1.0
    mat /= norms
    C = np.dot(mat, mat.T).clip(-1.0, 1.0)   # float32, N×N
    del mat
    print(f"  correlation done ({time.time()-t0:.1f}s)", flush=True)

    t0 = time.time()
    # Row-wise ascending rank — chunked to keep float64 intermediate ≤ ~500 rows
    R = np.empty((N, N), dtype=np.float32)
    RCHUNK = 500
    for i in range(0, N, RCHUNK):
        j = min(i + RCHUNK, N)
        R[i:j] = rankdata(C[i:j], axis=1).astype(np.float32)
    del C
    MR = np.empty_like(R)
    np.multiply(R, R.T, out=MR)   # avoids N×N intermediate — peak is 2×N² not 3×N²
    np.sqrt(MR, out=MR)
    np.fill_diagonal(MR, 0)
    del R
    print(f"  MR normalisation done ({time.time()-t0:.1f}s)", flush=True)
    return MR


# ── Density threshold ─────────────────────────────────────────────────────────
def density_threshold(MR, density):
    upper = MR[np.triu_indices(len(MR), k=1)]
    upper_sorted = np.sort(upper)[::-1]
    idx = max(0, round(density * len(upper_sorted)) - 1)
    thr = upper_sorted[idx]
    print(f"  density threshold = {thr:.3f} "
          f"(top {density*100:.1f}% of {len(upper_sorted):,} edges)", flush=True)
    return thr


# The per-pair hypergeometric test is performed inline in the vectorised main loop
# below; a previously unused standalone hg_test() reference implementation was removed.


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    out  = Path(args.out_dir)
    out.joinpath("RData", "comparison_tables").mkdir(parents=True, exist_ok=True)
    out.joinpath("RData", "centrality").mkdir(parents=True, exist_ok=True)
    out.joinpath("correlationPlots").mkdir(parents=True, exist_ok=True)

    print(f"ComPlEx Python  |  {args.s1_name} vs {args.s2_name}", flush=True)
    print(f"density={args.density}  method={args.cor_method}", flush=True)

    # ── Load data ──────────────────────────────────────────────────────────
    print("\n[1] Loading expression data...", flush=True)
    s1 = load_expr(args.s1_expr, args.min_expr, args.min_samples)
    s2 = load_expr(args.s2_expr, args.min_expr, args.min_samples)

    print("\n[2] Building orthogroup pairs...", flush=True)
    ortho = build_ortho(args.orthologs, args.s1_name, args.s2_name,
                        set(s1.index), set(s2.index))

    # Filter expression to genes with orthologs
    s1 = s1.loc[s1.index.isin(ortho["Species1"])]
    s2 = s2.loc[s2.index.isin(ortho["Species2"])]
    print(f"  After ortholog filter: {len(s1):,} S1 genes, {len(s2):,} S2 genes",
          flush=True)

    # Index lookups
    idx1_to_s1 = dict(enumerate(s1.index))
    s1_to_idx1 = {g: i for i, g in idx1_to_s1.items()}
    idx2_to_s2 = dict(enumerate(s2.index))
    s2_to_idx2 = {g: i for i, g in idx2_to_s2.items()}

    # Ortholog maps
    s1_to_s2 = ortho.groupby("Species1")["Species2"].apply(list).to_dict()
    s2_to_s1 = ortho.groupby("Species2")["Species1"].apply(list).to_dict()

    N1, N2 = len(s1), len(s2)

    # ── MR network for S1 — compute, extract centrality + neighbourhood, free ─
    print(f"\n[3] Computing MR network for {args.s1_name}...", flush=True)
    MR1  = mutual_rank(s1, args.cor_method)
    thr1 = density_threshold(MR1, args.density)

    # Per-gene network degree (centrality): the threshold is the sorted top-1% value of
    # the MR upper triangle, matching the canonical rcomplex density-threshold pattern.
    # validate_complex.{py,R} confirm the Python and R per-gene degrees agree (Spearman
    # rho = 1.0). Degree is used only as a relative covariate (e.g. the pN/pS
    # degree/expression confound control).
    c1_thr = np.sort(MR1[np.triu_indices(N1, k=1)])[::-1][
                 round(0.01 * len(MR1[np.triu_indices(N1, k=1)])) - 1]
    deg1 = (MR1 >= c1_thr).sum(axis=1)
    pd.DataFrame({"Genes": s1.index, "Degree": deg1}).to_csv(
        out / "RData" / "centrality" / f"centrality_{args.s1_name}.tsv",
        sep="\t", index=False)

    neigh1 = (MR1 >= thr1).astype(np.uint8)
    np.fill_diagonal(neigh1, 0)
    m1_arr = neigh1.sum(axis=1).astype(np.int32)
    del MR1   # free ~3–5 GB before computing MR2

    # ── MR network for S2 — compute, extract centrality + neighbourhood, free ─
    print(f"\n[4] Computing MR network for {args.s2_name}...", flush=True)
    MR2  = mutual_rank(s2, args.cor_method)
    thr2 = density_threshold(MR2, args.density)

    c2_thr = np.sort(MR2[np.triu_indices(N2, k=1)])[::-1][
                 round(0.01 * len(MR2[np.triu_indices(N2, k=1)])) - 1]
    deg2 = (MR2 >= c2_thr).sum(axis=1)
    pd.DataFrame({"Genes": s2.index, "Degree": deg2}).to_csv(
        out / "RData" / "centrality" / f"centrality_{args.s2_name}.tsv",
        sep="\t", index=False)

    neigh2 = (MR2 >= thr2).astype(np.uint8)
    np.fill_diagonal(neigh2, 0)
    m2_arr = neigh2.sum(axis=1).astype(np.int32)
    del MR2   # free ~4–5 GB before mapped-matrix computation

    # ── Fully vectorised comparison ────────────────────────────────────────
    print(f"\n[5] Vectorised comparison of {len(ortho):,} ortholog pairs...",
          flush=True)
    t0 = time.time()

    from scipy.sparse import csr_matrix

    s1_idx_arr = np.array([s1_to_idx1[g] for g in ortho["Species1"]], dtype=np.int32)
    s2_idx_arr = np.array([s2_to_idx2[g] for g in ortho["Species2"]], dtype=np.int32)
    ones = np.ones(len(ortho), dtype=np.float32)
    ortho_fwd = csr_matrix((ones, (s2_idx_arr, s1_idx_arr)), shape=(N2, N1))
    ortho_rev = csr_matrix((ones, (s1_idx_arr, s2_idx_arr)), shape=(N1, N2))

    # Build mapped neighbourhood matrices in chunks to limit peak memory.
    # neigh2_in_s1[j, i] = 1 if any S2 neighbour of gene j has an S1 ortholog at i.
    # neigh1_in_s2[i, j] = 1 if any S1 neighbour of gene i has an S2 ortholog at j.
    print(f"  Building mapped neighbourhood matrices (chunked)...", flush=True)
    MCHUNK = 1000

    neigh2_in_s1 = np.zeros((N2, N1), dtype=np.uint8)
    for i in range(0, N2, MCHUNK):
        j = min(i + MCHUNK, N2)
        chunk = neigh2[i:j].astype(np.float32) @ ortho_fwd
        neigh2_in_s1[i:j] = (chunk > 0).astype(np.uint8)
        del chunk

    neigh1_in_s2 = np.zeros((N1, N2), dtype=np.uint8)
    for i in range(0, N1, MCHUNK):
        j = min(i + MCHUNK, N1)
        chunk = neigh1[i:j].astype(np.float32) @ ortho_rev
        neigh1_in_s2[i:j] = (chunk > 0).astype(np.uint8)
        del chunk

    print(f"  Mapped matrices ready ({(time.time()-t0)/60:.1f} min)", flush=True)

    CHUNK = 50_000
    n_pairs = len(ortho)
    og_arr = ortho["Ortholog_Group"].values
    s1_arr = ortho["Species1"].values
    s2_arr = ortho["Species2"].values
    res_rows = []

    for start in range(0, n_pairs, CHUNK):
        end = min(start + CHUNK, n_pairs)
        i1c = s1_idx_arr[start:end]
        i2c = s2_idx_arr[start:end]

        n1c = neigh1[i1c]; n2c = neigh2[i2c]
        n2m = neigh2_in_s1[i2c]; n1m = neigh1_in_s2[i1c]

        x1 = (n1c & n2m).sum(axis=1).astype(np.int32)
        k1 = n2m.astype(bool).sum(axis=1).astype(np.int32)
        x2 = (n2c & n1m).sum(axis=1).astype(np.int32)
        k2 = n1m.astype(bool).sum(axis=1).astype(np.int32)
        m1c = m1_arr[i1c]; m2c = m2_arr[i2c]

        # FDR-pool filter (see --min-overlap / --overlap-mode). Defines which
        # pairs count as tests in the BH correction below. p-values are still
        # only computable when x > 1 (x <= 1 -> p = 1), so this controls the FDR
        # denominator, not which pairs can be called co-expressologs.
        mo = args.min_overlap
        if args.overlap_mode == "and":
            keep = (x1 >= mo) & (x2 >= mo)
        else:
            keep = (x1 >= mo) | (x2 >= mo)
        for li in np.where(keep)[0]:
            gi = start + li
            xi1, xi2 = int(x1[li]), int(x2[li])
            ki1, ki2 = int(k1[li]), int(k2[li])
            mi1, mi2 = int(m1c[li]), int(m2c[li])
            p1 = hypergeom.sf(xi1-1, N1, mi1, ki1) if xi1 > 1 else 1.0
            p2 = hypergeom.sf(xi2-1, N2, mi2, ki2) if xi2 > 1 else 1.0
            res_rows.append((og_arr[gi], s1_arr[gi], s2_arr[gi],
                             xi1, p1, xi2, p2))

        elapsed = time.time() - t0
        print(f"  {end:,}/{n_pairs:,}  {elapsed/60:.1f} min  "
              f"~{(n_pairs-end)/max(end,1)*elapsed/60:.0f} min remaining  "
              f"{len(res_rows):,} candidates in FDR pool", flush=True)

    res = pd.DataFrame(res_rows, columns=[
        "OrthoGroup","Species1","Species2",
        "Species1.neigh.overlap","Species1.p.val",
        "Species2.neigh.overlap","Species2.p.val"])
    print(f"  Done in {(time.time()-t0)/60:.1f} min — "
          f"{len(res):,} candidate pairs in FDR pool "
          f"(min-overlap={args.min_overlap}, mode={args.overlap_mode})",
          flush=True)

    # ── FDR correction — filter to significant co-expressologs only ────────
    print("\n[6] FDR correction and saving...", flush=True)
    ct = pd.DataFrame(res)
    if len(ct) == 0:
        print("WARNING: no ortholog pairs with neighbourhood overlap found.")
        return

    _, p1_fdr, _, _ = multipletests(ct["Species1.p.val"], method="fdr_bh")
    _, p2_fdr, _, _ = multipletests(ct["Species2.p.val"], method="fdr_bh")
    ct["Species1.p.val"] = p1_fdr
    ct["Species2.p.val"] = p2_fdr
    ct["Max.p.val"]      = ct[["Species1.p.val","Species2.p.val"]].max(axis=1)
    ct = ct.sort_values("Max.p.val").reset_index(drop=True)

    n_candidates = len(ct)
    # Keep only significant co-expressologs — the full candidate table would contain
    # most pairs, since at 3% density a non-zero neighbourhood overlap arises by chance
    # for nearly all pairs.
    ct = ct[ct["Max.p.val"] < 0.05].reset_index(drop=True)

    out_tsv = (out / "RData" / "comparison_tables" /
               f"comparison_{args.s1_name}_{args.s2_name}.tsv")
    ct.to_csv(out_tsv, sep="\t", index=False)
    print(f"  Saved {out_tsv}", flush=True)
    print(f"  Co-expressologs at FDR<0.05: {len(ct):,} "
          f"(from {n_candidates:,} candidate pairs in the FDR pool)", flush=True)

    # Correlation plot of p-values
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(-np.log10(ct["Species1.p.val"] + 1e-300),
                   -np.log10(ct["Species2.p.val"] + 1e-300),
                   s=4, alpha=0.3)
        ax.set_xlabel(f"{args.s1_name} p-value (−log₁₀)")
        ax.set_ylabel(f"{args.s2_name} p-value (−log₁₀)")
        rho = np.corrcoef(-np.log10(ct["Species1.p.val"] + 1e-300),
                          -np.log10(ct["Species2.p.val"] + 1e-300))[0, 1]
        ax.set_title(f"ρ = {rho:.3f}")
        fig.savefig(out / "correlationPlots" /
                    f"ortholog_correlation_{args.s1_name}_{args.s2_name}.pdf",
                    bbox_inches="tight")
        plt.close()
    except Exception as e:
        print(f"  (plot skipped: {e})")

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
