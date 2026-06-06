# ComPlEx Python implementation

Python reimplementation of the ComPlEx co-expressolog algorithm. It is intended
as a fast, faithful drop-in for the laboratory's R ComPlEx pipeline (the
RComPlEx `.Rmd` lineage, as captured by `validate_complex.R`), reproducing its
results while reducing runtime from hours to minutes for large co-expression
networks.

It is *not* a reproduction of the original publication (Netotea et al., 2014,
*BMC Genomics* 15:106, https://doi.org/10.1186/1471-2164-15-106). The R pipeline
it tracks has itself diverged from that paper in several deliberate ways; those
differences are documented in
[Relationship to the original publication](#relationship-to-the-original-publication)
below so that users understand what this implementation does and does not
reproduce.

## Algorithm

For each ortholog pair (gene_i in species 1, gene_j in species 2):

1. Compute Pearson correlation matrices for each species expression dataset.
2. Apply Mutual Rank (MR) normalisation: MR[i,j] = sqrt(rank_i(j) × rank_j(i)).
   MR reduces the influence of highly connected hub genes.
3. Retain the top `density` fraction of MR edges as gene neighbours (default 3%).
4. For each ortholog pair, compute the overlap between:
   - the S1 neighbourhood of gene_i, and
   - the S1-mapped neighbourhood of gene_j (S2 neighbours of j, mapped back to S1
     via the ortholog table).
5. Test overlap significance with a hypergeometric test (one-tailed); repeat in
   the S2 → S1 direction.
6. Apply Benjamini–Hochberg FDR correction (the denominator is selectable via
   `--fdr-denominator`; see below); report pairs with
   max(FDR_direction1, FDR_direction2) < 0.05 as co-expressologs.

## Relationship to the original publication

This implementation reproduces the laboratory's current R ComPlEx pipeline, not
the method exactly as published by Netotea et al. (2014). The core test — the
hypergeometric assessment of neighbourhood conservation — is implemented exactly
as described in the paper: for an ortholog pair, with `N` = genes in the focal
species, `n` = size of the focal gene's neighbourhood, `k` = the orthologs of the
partner gene's neighbours mapped back to the focal species, and `x` = their
overlap, the one-tailed probability `P(X ≥ x)` is computed from the
hypergeometric distribution and corrected with Benjamini–Hochberg FDR. The
points below list where the pipeline (and therefore this implementation) departs
from the publication. All are inherited from the R pipeline; none are introduced
here.

1. **Co-expression measure and normalisation.** The paper computed mutual
   information (B-spline estimator) and background-corrected it with the CLR
   (context likelihood of relatedness) method. This implementation, like the R
   pipeline, uses Pearson (default) or Spearman correlation normalised by Mutual
   Rank (MR), where `MR(i,j) = sqrt(rank_i(j) × rank_j(i))` (the standard
   element-wise definition). This matches `validate_complex.R`. Note that the
   Nextflow port `RComPlEx-NF` computes MR via a rank-matrix product rather than
   element-wise, so numerical agreement here is with the `.Rmd` /
   `validate_complex.R` reference, not necessarily with that port.

2. **Network gene universe.** The paper built genome-wide co-expression networks
   and only restricted the *reporting* universe to connected, ortholog-bearing
   genes. This implementation, like the R pipeline, restricts the correlation
   matrix to ortholog-bearing genes *before* computing correlations, so the
   neighbourhoods, the density threshold, and `N` are all defined within
   ortholog-only space. The rationale is cross-species comparability (equal gene
   universes, no species-specific genes inflating one side). The size of the
   effect relative to the paper depends on the fraction of each transcriptome
   that has orthologs.

3. **Directionality.** The paper reported neighbourhood conservation
   per direction. This implementation requires conservation to be *reciprocal*:
   a pair is reported as a co-expressolog only if it is significant in both
   directions, enforced by taking `max(FDR_S1→S2, FDR_S2→S1) < 0.05`. This is
   stricter than the paper and matches the R pipeline.

4. **FDR denominator.** The number of tests entering the BH correction is
   configurable; see [FDR denominator](#fdr-denominator---fdr-denominator) below.
   The default reproduces the original Python behaviour; `both-overlap`
   reproduces `validate_complex.R` exactly, and `all-pairs` reproduces the
   paper's convention of correcting across every ortholog pair.

5. **Scope — downstream stages not implemented.** This implementation produces
   the pairwise co-expressolog comparison table only. It does not perform the
   later RComPlEx pipeline stages — maximal-clique detection
   (`igraph::max_cliques`), CLR normalisation, or the signed/unsigned
   polarity-divergence analysis. Those are extensions in the R pipeline, not part
   of the 2014 method; add them downstream if required.

## Differences from the R implementation

### FDR denominator (`--fdr-denominator`)

The Benjamini–Hochberg correction is applied per direction, but the *number of
tests* entered into the correction (the denominator `m`) is configurable, since
this is the single choice that most affects how many co-expressologs are
called. A pair can only reach FDR < 0.05 if it has overlap > 1 in **both**
directions (the filter uses the max of the two adjusted p-values), and such
pairs are present under every mode; the modes differ only in `m`, which scales
the adjusted p-values. Larger `m` ⇒ more conservative.

| `--fdr-denominator` | tests entered into BH | matches |
|---|---|---|
| `candidates` (default) | pairs with overlap > 1 in at least one direction | original Python behaviour |
| `both-overlap` | pairs with overlap > 0 in both directions | `validate_complex.R` exactly |
| `all-pairs` | every tested ortholog pair (p = 1 where overlap ≤ 1) | the original ComPlEx paper (Netotea et al. 2014: "a p-value was computed for each ortholog pair … FDR … controlled at 0.05") |

`both-overlap` has been verified to reproduce the significant set and the
BH-adjusted p-values of the `validate_complex.R` algorithm to machine
precision. `all-pairs` reproduces the paper's convention of correcting across
every ortholog pair and is the most conservative. The non-significant pairs in
`both-overlap`/`all-pairs` all carry p = 1 and only enter via the denominator,
so they are accounted for by padding rather than by materialising a row each —
the result is identical to including them explicitly but uses no extra memory.

The default remains `candidates` for backward compatibility. For analyses
intended to reproduce the published method, use `--fdr-denominator all-pairs`;
to reproduce the R pipeline numerically, use `--fdr-denominator both-overlap`.

On an earlier 1,500-gene validation subset, the default `candidates` mode
entered 6,549 pairs into BH versus 7,130 for the R `x1 > 0 AND x2 > 0` filter
(the 581 extra pairs all had x = 1 in both directions and p = 1.0); all 28
R-significant pairs were a strict subset of the 43 `candidates` pairs, the 15
additional pairs having BH-adjusted p between 0.046 and 0.054. Selecting
`both-overlap` removes this discrepancy entirely.

### Memory and speed

| | R | Python (this implementation) |
|---|---|---|
| Correlation matrix | `cor()`, float64 | Manual dot-product, float32 |
| Ranking | `t(apply(net, 1, rank))` | `scipy.stats.rankdata`, chunked float32 |
| Neighbourhood overlap | Per-pair loop | Vectorised via sparse matrix × dense |
| Peak memory (N=30,000) | ~12 GB | ~5–6 GB |
| Runtime (N=30,000) | ~4 hours | ~10–15 minutes |

## Requirements

```
numpy
pandas
scipy
statsmodels
```

Install with:
```
pip install numpy pandas scipy statsmodels
```

## Usage

### Full analysis

```bash
python3 complex_py.py \
    --s1-expr  species1_expression.txt \
    --s2-expr  species2_expression.txt \
    --orthologs  genes_ortholog_categories.tsv \
    --s1-name  spruce \
    --s2-name  pine \
    --out-dir  results/ComPlEx/cold_needle \
    [--density 0.03] \
    [--cor-method pearson] \
    [--fdr-denominator candidates|both-overlap|all-pairs]
```

**Expression files**: tab-separated, first column = gene ID (row index),
remaining columns = VST or normalised count values.

**Orthologs file**: tab-separated with columns `gene`, `species`,
`Ortholog_Group`. The `species` column must match the values of the internal
species names (`Picea_abies` and `Pinus_sylvestris` by default; see
`build_ortho()` in `complex_py.py` to adapt).

**Output**: `<out-dir>/RData/comparison_tables/comparison_<s1>_<s2>.tsv`
with columns:
- `OrthoGroup` — orthogroup identifier
- `Species1`, `Species2` — gene identifiers
- `Species1.neigh.overlap`, `Species2.neigh.overlap` — overlap counts (x1, x2)
- `Species1.p.val`, `Species2.p.val` — BH-adjusted p-values per direction
- `Max.p.val` — max of the two adjusted p-values (filter criterion)

All rows in the output have `Max.p.val < 0.05`.

### Validation against R

1. Edit the data paths at the top of `validate_complex.py` (S1_EXPR_FILE,
   S2_EXPR_FILE, ORTHO_FILE, S1_SPECIES, S2_SPECIES, OUT_DIR).
2. Run the Python validation:
   ```bash
   python3 validate_complex.py
   ```
   This saves subset TSVs to OUT_DIR.

3. Set OUT_DIR in `validate_complex.R` to the same directory, then run:
   ```bash
   Rscript validate_complex.R
   ```

4. Re-run `validate_complex.py` to include the R comparison output.

**Expected result on a 1,500-gene subset**: optimised Python = loop-based
reference exactly (x1, x2, m1, m2, and raw p-values identical in all pairs;
same FDR<0.05 set). Float32 neighbourhood assignment matches float64 exactly
at 3% density. See the "Differences" section above for the R comparison.
