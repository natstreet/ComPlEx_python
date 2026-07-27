# validate_complex — Python vs R concordance (run 2026-07-06)

Subset: 1,500 genes/species, SEED=42; inputs spruce SC_expression.txt vs pine
PC_expression.txt, orthologs doc/genes_ortholog_categories.tsv (3% density).
Reproduce: `python3 validate_complex.py && Rscript validate_complex.R && python3 validate_complex.py`

Optimised (float32 vectorised) vs float64 loop reference:
  x1/x2 exact match 1.0000; m1/m2 1.0000; median rel p-error 0.00e+00
  FDR<0.05: opt 756 / ref 756 / intersection 756 (0 only-opt, 0 only-ref)

Python vs R (original ComPlEx algorithm):
  co-expressolog CALLS: py 756 / r 756 / intersection 756 (0 discordant)
  FDR pool: py 14,883 / r 14,878 (5 pairs differ at the x>0 boundary; denominator only)
  per-gene degree Spearman rho = 1.000 (both species)

Conclusion: identical co-expressolog calls; the manuscript's "reproduces the R
implementation exactly, identical calls" is substantiated on this subset.
