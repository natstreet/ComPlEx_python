# ComPlEx Python implementation

Python reimplementation of the ComPlEx co-expressolog algorithm. Produces
results numerically equivalent to the original R implementation while reducing
runtime from hours to minutes for large co-expression networks.

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
6. Apply Benjamini–Hochberg FDR correction across all pairs; report pairs with
   max(FDR_direction1, FDR_direction2) < 0.05 as co-expressologs.

## Differences from the R implementation

### FDR-pool definition (the default reproduces R)

By default this port reproduces the published R ComPlEx (Netotea *et al.* 2014)
for the FDR pool. With the shipped defaults `--min-overlap 1 --overlap-mode and`,
the pairs entering Benjamini–Hochberg correction are those with non-zero
neighbourhood overlap in both directions (`x1 > 0 AND x2 > 0`), exactly as in R.

An optional sensitivity mode, `--min-overlap 2 --overlap-mode or`, instead drops
pairs with overlap of at most one in both directions before BH correction (`x > 1`
in at least one direction). Because at large N and 3% density the expected random
overlap is approximately N × density² (about 22 genes for N > 20,000), `x > 0` is
satisfied for essentially all pairs; this mode therefore mainly shrinks the FDR
denominator and is slightly less conservative.

Neither setting changes which pairs can be *called* co-expressologs: a call requires
`Max.p < 0.05`, which already implies overlap of at least two in one direction
(overlap of at most one yields p = 1). The choice affects only the FDR denominator.

On a 1,500-gene validation subset the `x > 1` sensitivity mode enters 6,549 pairs
into BH against 7,130 under the default R behaviour (the 581 extra pairs all have
x = 1 in both directions and p = 1.0); the 28 pairs significant under R are a strict
subset of the 43 significant under that mode, the 15 additional pairs having
BH-adjusted p between 0.046 and 0.054. At full scale the practical effect is
negligible and both approaches are statistically valid.

### Network centrality (Degree)

The per-gene `Degree` column written to `RData/centrality/` uses the sorted top-1%
value of the mutual-rank upper triangle as its threshold, matching the canonical
rcomplex density-threshold pattern. The validation harness computes the degree under
both implementations on the same subset (`validate_complex.R` writes
`r_centrality_s1/s2.tsv`; `validate_complex.py` writes `py_centrality_s1/s2.tsv` and
reports the Spearman correlation between them): the Python and R per-gene degrees agree
exactly (Spearman rho = 1.0, identical median degree). `Degree` is used only as a
relative covariate (for example the pN/pS degree/expression confound control).

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
matplotlib
biopython
```

Install with:
```
pip install numpy pandas scipy statsmodels matplotlib biopython
```

See `requirements.txt` for version floors; record exact versions with `pip freeze` for the published run.

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
    [--min-overlap 1] [--overlap-mode and]
```

`--min-overlap` and `--overlap-mode` set the FDR pool; the defaults
(`1` and `and`) reproduce the published R behaviour (see "Differences" below).

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

The full validation — including a three-way cross-check against the canonical
published `rcomplex` — lives in [`validation/`](validation/); see
[`validation/VALIDATION_NOTE.md`](validation/VALIDATION_NOTE.md) for the write-up
and `validation/README.md` for how to reproduce it.

On a 5,000-gene-per-species cold-needle subset (seed 42; 92,537 ortholog pairs with
non-zero two-directional overlap), three implementations — this Python port
(`complex_py.py`, default `--min-overlap 1 --overlap-mode and`), an in-house R
transcription (`validation/validate_complex.R`), and the canonical published
`rcomplex` run verbatim (`validation/validate_against_canonical_rcomplex.R`) —
called an identical set of 11,254 co-expressologs (Jaccard = 1.000). Raw
hypergeometric p-values agreed to machine precision against the R transcription
(maximum absolute difference 2.1 × 10⁻¹⁵), and per-gene degree centrality was
identical (Spearman ρ = 1.000 in both species). Independent 1,500-gene cold-needle
and drought-root subsets (seed 123) gave the same agreement.

A lighter top-level quick-check (`validate_complex.py` + `validate_complex.R`,
recorded in `VALIDATION_RESULT.md`) reproduces the same conclusion on a 1,500-gene
subset, and separately confirms the float32-vectorised path matches the float64
loop reference exactly at 3% density.
