# Supplementary note: validation of the Python ComPlEx re-implementation

The co-expressolog analysis was performed with a Python/NumPy re-implementation
of ComPlEx (Netotea et al. 2014), developed for this study (`complex_py.py`).
To confirm that the re-implementation reproduces the published R implementation,
both were run on an identical validation subset and their co-expressolog calls
compared.

**Setup.** A reproducible subset of 1,500 genes per species (seed 42) was drawn
from the Norway spruce (SCN) and Scots pine (PCN) cold-needle expression matrices;
the cross-species ortholog pairs among them with non-zero neighbourhood overlap in
both directions (7,130 pairs) formed the validation set. Both
implementations used the published settings: Pearson correlation, Mutual-Rank
normalisation, a 3% network density threshold, a two-directional hypergeometric
test of neighbourhood overlap, retention of pairs with overlap > 0 in both
directions, Benjamini–Hochberg FDR correction, and a Max-FDR < 0.05 significance
threshold. In `complex_py.py` this corresponds to the default, R-faithful mode
(`--min-overlap 1 --overlap-mode and`).

**Result.** The two implementations called an identical set of co-expressologs:

| | co-expressologs (Max FDR < 0.05) |
|---|---|
| R (`validate_complex.R`, transcription of `ComPlEx.R`) | 28 |
| Python (`complex_py.py`, R-faithful mode) | 28 |
| shared | 28 (Jaccard = 1.000) |

Across all 7,130 candidate pairs the raw hypergeometric p-values agreed to machine
precision (maximum absolute difference = 1.3 × 10⁻¹⁵; maximum |Δlog₁₀ p| = 3.6 ×
10⁻¹⁵), and the two implementations called an identical co-expressolog set (only-Python
= 0, only-R = 0). The
Python implementation therefore reproduces the R ComPlEx exactly.

**Note on the candidate / FDR pool.** A pair can only be called a co-expressolog
when its Max FDR is significant, which already requires overlap ≥ 2 in both
directions (overlap ≤ 1 yields p = 1). The published method (and the default
here) places all pairs with overlap > 0 in both directions in the FDR pool; an
optional, less conservative mode (`--min-overlap 2 --overlap-mode or`) excludes
overlap = 1 pairs before correction and is reported only as a sensitivity
analysis. The headline results use the conservative, R-faithful default.

Reproduce with the scripts in this directory (see `README.md`).
