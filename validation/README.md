# Validation: Python ComPlEx reproduces the published R implementation

This package confirms that the parameterised Python re-implementation
(`complex_py.py`), run in **R-faithful mode**, reproduces the co-expressolog
calls of the original R ComPlEx (Netotea et al. 2014) on a shared subset.

R-faithful mode is the default and the headline setting for the manuscript:
`--min-overlap 1 --overlap-mode and`, i.e. the Benjamini–Hochberg pool is all
ortholog pairs with ≥1 shared neighbour in **both** directions (a p-value is
only computed when overlap > 1; pairs with overlap ≤ 1 carry p = 1 but remain in
the FDR denominator, exactly as in `ComPlEx.R`). The earlier optimisation
(`--min-overlap 2 --overlap-mode or`, which drops x = 1 pairs before FDR) is the
less-conservative variant used only for the sensitivity analysis.

`validate_complex.R` here is a faithful standalone transcription of the
`ComPlEx.R` core loop (Pearson → Mutual Rank → 3 % density threshold →
per-pair two-directional hypergeometric test → `x1>0 & x2>0` filter → BH FDR →
Max-p). It is line-equivalent to the production `ComPlEx.R` but reads the small
subset files so you do not have to reconfigure the full pipeline.

## Requirements
- Python: numpy, pandas, scipy, statsmodels
- R: dplyr, matrixStats

## Steps (run from this `validation/` directory)

```bash
# 1. Build a reproducible subset (one stress/tissue network pair, e.g. cold needles).
#    Point the paths at your copies of the expression matrices and ortholog table.
python3 make_subset.py \
    --s1-expr  /path/to/SCN_expression.txt \
    --s2-expr  /path/to/PCN_expression.txt \
    --orthologs /path/to/genes_ortholog_categories.tsv \
    --n-pairs 3000 --seed 42 \
    --out-dir complex_validation_output

# 2. Run the PRODUCTION Python implementation in R-faithful mode on the subset.
python3 ../complex_py.py \
    --s1-expr   complex_validation_output/s1_subset.tsv \
    --s2-expr   complex_validation_output/s2_subset.tsv \
    --orthologs complex_validation_output/ortho_for_python.tsv \
    --s1-name spruce --s2-name pine \
    --out-dir complex_validation_output \
    --min-overlap 1 --overlap-mode and

# 3. Run the R reference on the same subset (writes r_comparison.tsv).
Rscript validate_complex.R

# 4. Diff the two significant co-expressolog sets.
python3 compare.py \
    --python complex_validation_output/RData/comparison_tables/comparison_spruce_pine.tsv \
    --r      complex_validation_output/r_comparison.tsv
```

## Interpreting the result
- **PASS — exact reproduction**: identical significant sets (Jaccard = 1.0).
- **PASS (float-boundary)**: Jaccard ≥ 0.99; the few differing pairs sit exactly
  at the 3 % density cutoff and differ only because Python uses float32 and R
  uses float64. Inspect the `only-Python` / `only-R` pairs to confirm this.
- **CHECK**: materially different sets — do not proceed; report back.

For complete assurance you may additionally run the production `ComPlEx.R`
itself on the same subset and diff its output the same way; `validate_complex.R`
is provided because `ComPlEx.R` is not a stand-alone CLI.
