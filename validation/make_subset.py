#!/usr/bin/env python3
"""
make_subset.py — build a reproducible subset for the R-vs-Python ComPlEx
validation (see validation/README.md).

Writes, into OUT_DIR (default ./complex_validation_output):
  s1_subset.tsv, s2_subset.tsv   expression subsets (gene x sample, "Genes" header)
                                 — read by BOTH validate_complex.R and complex_py.py
  ortho_subset.tsv               Ortholog_Group, Species1, Species2  (for validate_complex.R)
  ortho_for_python.tsv           gene, species, Ortholog_Group       (for complex_py.py)

The subset is the genes of a seeded random sample of cross-species ortholog
pairs for which both members pass the standard featureSelect filter
(>= MIN_EXPR in >= MIN_SAMPLES samples) in both expression matrices.

Usage:
  python3 make_subset.py \
      --s1-expr SCN_expression.txt --s2-expr PCN_expression.txt \
      --orthologs genes_ortholog_categories.tsv \
      [--n-pairs 3000] [--seed 42] [--out-dir complex_validation_output]
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

MIN_EXPR, MIN_SAMPLES = 1.0, 2

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--s1-expr", required=True, help="spruce expression matrix (e.g. SCN_expression.txt)")
    p.add_argument("--s2-expr", required=True, help="pine expression matrix (e.g. PCN_expression.txt)")
    p.add_argument("--orthologs", required=True, help="genes_ortholog_categories.tsv (gene, species, Ortholog_Group)")
    p.add_argument("--s1-species", default="Picea_abies")
    p.add_argument("--s2-species", default="Pinus_sylvestris")
    p.add_argument("--n-pairs", type=int, default=3000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out-dir", default="complex_validation_output")
    return p.parse_args()

def load_expr(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    sel = (df >= MIN_EXPR).sum(axis=1) >= MIN_SAMPLES
    df = df.loc[sel]
    return df.loc[df.std(axis=1) > 0]

def main():
    a = parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    print("Loading expression matrices...", flush=True)
    s1 = load_expr(a.s1_expr); s2 = load_expr(a.s2_expr)

    og = pd.read_csv(a.orthologs, sep="\t")
    # Exclude the unassigned-gene placeholder (matches the production parquet input,
    # which contains no "No orthogroup" rows). Without this, merging on the group
    # joins every unassigned gene to every other, producing millions of spurious
    # "pairs".
    og = og[~og.Ortholog_Group.isin(["No orthogroup", "No_orthogroup", ""])]
    p1 = og[og.species == a.s1_species][["gene", "Ortholog_Group"]].rename(columns={"gene": "Species1"})
    p2 = og[og.species == a.s2_species][["gene", "Ortholog_Group"]].rename(columns={"gene": "Species2"})
    p1 = p1[p1.Species1.isin(s1.index)]; p2 = p2[p2.Species2.isin(s2.index)]
    # Restrict the validation to 1:1 orthogroups (exactly one expressed gene per
    # species) so R (fixed pair list) and Python (re-merges on the group) operate
    # on identical pairs. The algorithm is identical for paralog-containing groups.
    one2one = (set(p1.groupby("Ortholog_Group").size()[lambda s: s == 1].index) &
               set(p2.groupby("Ortholog_Group").size()[lambda s: s == 1].index))
    p1 = p1[p1.Ortholog_Group.isin(one2one)]; p2 = p2[p2.Ortholog_Group.isin(one2one)]
    pairs = p1.merge(p2, on="Ortholog_Group").reset_index(drop=True)
    print(f"  {len(pairs):,} eligible 1:1 ortholog pairs (after excluding 'No orthogroup')", flush=True)

    rng = np.random.default_rng(a.seed)
    take = min(a.n_pairs, len(pairs))
    sub = pairs.iloc[rng.choice(len(pairs), size=take, replace=False)].reset_index(drop=True)
    s1_genes = sorted(sub.Species1.unique())
    s2_genes = sorted(sub.Species2.unique())

    s1.loc[s1_genes].to_csv(out / "s1_subset.tsv", sep="\t", index_label="Genes")
    s2.loc[s2_genes].to_csv(out / "s2_subset.tsv", sep="\t", index_label="Genes")
    sub[["Ortholog_Group", "Species1", "Species2"]].to_csv(out / "ortho_subset.tsv", sep="\t", index=False)

    # long format for complex_py.py (gene, species, Ortholog_Group)
    long = pd.concat([
        pd.DataFrame({"gene": s1_genes, "species": a.s1_species}).merge(
            sub[["Species1", "Ortholog_Group"]].drop_duplicates(), left_on="gene", right_on="Species1")[["gene", "species", "Ortholog_Group"]],
        pd.DataFrame({"gene": s2_genes, "species": a.s2_species}).merge(
            sub[["Species2", "Ortholog_Group"]].drop_duplicates(), left_on="gene", right_on="Species2")[["gene", "species", "Ortholog_Group"]],
    ], ignore_index=True)
    long.to_csv(out / "ortho_for_python.tsv", sep="\t", index=False)

    print(f"  wrote subset to {out}/ : {len(s1_genes)} spruce, {len(s2_genes)} pine genes, {take} pairs", flush=True)

if __name__ == "__main__":
    main()
