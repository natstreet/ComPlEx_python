# Supplementary note: validation of the Python ComPlEx re-implementation

The co-expressolog analysis was performed with a Python/NumPy re-implementation
of ComPlEx (Netotea et al. 2014), developed for this study (`complex_py.py`).
To confirm that the re-implementation reproduces the published R implementation,
all were run on an identical validation subset and their co-expressolog calls
compared. Two R references were used: an in-house transcription of the ComPlEx
algorithm (`validate_complex.R`) and the canonical published `rcomplex` code run
verbatim (`validate_against_canonical_rcomplex.R`).

## Setup (reproducible)

A subset of 5,000 genes per species (seed 42) was drawn from the Norway spruce
(SCN) and Scots pine (PCN) cold-needle expression matrices; the cross-species
ortholog pairs among them with non-zero neighbourhood overlap in both directions
(92,537 pairs) formed the validation set. All implementations used the published
settings: Pearson correlation, Mutual-Rank normalisation (MR = sqrt(rank_i(j) x
rank_j(i))), a 3% network-density threshold, a two-directional hypergeometric
test of neighbourhood overlap, retention of pairs with overlap > 0 in both
directions, Benjamini-Hochberg FDR correction, and a Max-FDR < 0.05 significance
threshold. In `complex_py.py` this is the R-faithful default (`--min-overlap 1
--overlap-mode and`).

**Reproduce the exact subset and result** (from the repository root, with the
spruce-pine deposit checked out beside this repo):

    N_GENES=5000 SEED=42 \
    S1_EXPR=<deposit>/data/expression/SCN_expression.txt \
    S2_EXPR=<deposit>/data/expression/PCN_expression.txt \
    ORTHO_FILE=<deposit>/doc/genes_ortholog_categories.tsv \
    python3 validate_complex.py            # writes complex_validation_output/

    Rscript validation/validate_complex.R                 # in-house R transcription
    VALIDATION_DIR=complex_validation_output \
      Rscript validation/canonical_run.R                  # canonical published rcomplex

Output tables (subset, per-pair calls and raw p-values for each implementation)
are written to `complex_validation_output/` and committed here so the result is
checkable without re-running.

## Result

All implementations called an identical set of co-expressologs on the 92,537-pair
validation set:

| implementation | co-expressologs (Max FDR < 0.05) |
|---|---|
| Python port (`complex_py.py`, R-faithful, float32) | 11,254 |
| Python reference (explicit-loop, float64) | 11,254 |
| in-house R transcription (`validate_complex.R`) | 11,254 |
| canonical published `rcomplex` (`canonical_run.R`) | 11,254 |
| **shared (all four)** | **11,254 (Jaccard = 1.000, 0 discordant)** |

Raw hypergeometric p-values agreed to machine precision on 92,534 of the 92,537
pairs; the remaining 3 pairs differed by at most 1.2 x 10^-2, are all far from
significance (Max FDR >= 0.82), and reflect a single borderline edge retained
versus dropped at the 3% density cutoff (a deterministic float32/float64 tie-break
during network construction). No pair's significance call is affected.

As an additional, fully independent check, a from-scratch re-implementation of the
algorithm written without reference to `complex_py.py` (`independent_complex.py`)
was run on a 1:1-orthogroup subset, where the two hypergeometric test directions
are provably identical, so a single-direction implementation is exact. On that
subset it reproduced the port's co-expressolog calls exactly (2,588 of 2,588,
Jaccard = 1.000). Reproduce with:

    python3 validation/make_subset.py \
      --s1-expr <deposit>/data/expression/SCN_expression.txt \
      --s2-expr <deposit>/data/expression/PCN_expression.txt \
      --orthologs <deposit>/doc/genes_ortholog_categories.tsv \
      --n-pairs 5000 --seed 42 --out-dir repro_out
    python3 validation/independent_complex.py repro_out
    # -> 5,000 1:1 pairs, 4,434 candidates, 2,588 co-expressologs (== the port on
    #    the same subset). independent_complex.py assumes 1:1 alignment and is only
    #    valid on this make_subset output, not on the many-to-many headline set.

The Python implementation therefore reproduces the published R ComPlEx exactly at
the level of the reported result.

## Note on the reproducer scripts

`validate_complex.py` (repository root) builds the many-to-many validation set
used for the headline numbers above (all ortholog pairs among the sampled genes).
The scripts in `validation/` (`make_subset.py` + `compare.py`) build a smaller,
1:1-restricted subset for a fast, exactly-matched R-vs-Python cross-check; that
mode gives a correspondingly smaller co-expressolog count (it is a demonstration,
not the headline set). Both confirm identical calls.

Reproduce with the scripts in this directory (see `README.md`).
