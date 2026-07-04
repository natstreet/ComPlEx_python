# Supplementary note: validation of the Python ComPlEx re-implementation

The co-expressolog analysis was performed with a Python/NumPy re-implementation
of ComPlEx (Netotea et al. 2014), developed for this study (`complex_py.py`).
To confirm that the re-implementation reproduces the published R implementation,
both were run on an identical validation subset and their co-expressolog calls
compared. Two independent R references were used: an in-house transcription of
the ComPlEx algorithm (`validate_complex.R`) and the canonical published
`rcomplex` code run verbatim (`validate_against_canonical_rcomplex.R`).

**Setup.** A reproducible subset of 5,000 genes per species (seed 42) was drawn
from the Norway spruce (SCN) and Scots pine (PCN) cold-needle expression matrices;
the cross-species ortholog pairs among them with non-zero neighbourhood overlap in
both directions (92,537 pairs) formed the validation set. All implementations used
the published settings: Pearson correlation, Mutual-Rank normalisation, a 3%
network density threshold, a two-directional hypergeometric test of neighbourhood
overlap, retention of pairs with overlap > 0 in both directions, Benjamini–Hochberg
FDR correction, and a Max-FDR < 0.05 significance threshold. In `complex_py.py`
this corresponds to the default, R-faithful mode (`--min-overlap 1 --overlap-mode
and`).

**Result.** All three implementations called an identical set of co-expressologs:

| | co-expressologs (Max FDR < 0.05) |
|---|---|
| R transcription (`validate_complex.R`) | 11,254 |
| Canonical published `rcomplex` | 11,254 |
| Python (`complex_py.py`, R-faithful mode) | 11,254 |
| shared (all three) | 11,254 (Jaccard = 1.000) |

Against the R transcription, raw hypergeometric p-values agreed to machine
precision across all 92,537 pairs (maximum absolute difference 2.1 × 10⁻¹⁵), and
degree centrality was identical (Spearman ρ = 1.000 in both species). Against the
canonical published `rcomplex`, the co-expressolog set was likewise identical
(only-canonical = 0, only-Python = 0, Jaccard = 1.000). Raw p-values agreed for
all but 3 of the 92,537 pairs (maximum absolute difference 1.2 × 10⁻²); those 3
pairs are all far from significance (Max FDR ≥ 0.82) and differ in only one of the
two test directions, reflecting a single borderline edge retained versus dropped at
the 3% density cutoff (a deterministic tie-break during network construction). No
pair's significance call is affected. The Python implementation therefore
reproduces the published R ComPlEx exactly at the level of the reported result.

**Note on the candidate / FDR pool.** A pair can only be called a co-expressolog
when its Max FDR is significant, which already requires overlap ≥ 2 in both
directions (overlap ≤ 1 yields p = 1). The published method (and the default
here) places all pairs with overlap > 0 in both directions in the FDR pool; an
optional, less conservative mode (`--min-overlap 2 --overlap-mode or`) excludes
overlap = 1 pairs before correction and is reported only as a sensitivity
analysis. The headline results use the conservative, R-faithful default.

A smaller 1,500-gene cold-needle subset and an independent 1,500-gene drought-root
subset (seed 123) gave the same agreement (identical co-expressolog sets).

Reproduce with the scripts in this directory (see `README.md`).
