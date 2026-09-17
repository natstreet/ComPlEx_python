#!/usr/bin/env python3
"""
Independent, from-scratch re-implementation of the ComPlEx co-expressolog
algorithm (Netotea et al. 2014), written to verify the complex_py.py port.
Not derived from complex_py.py — implements the published algorithm directly
with the stated settings: Pearson, Mutual-Rank normalisation, 3% density
threshold, two-directional hypergeometric neighbourhood-overlap test,
Benjamini-Hochberg FDR, Max-FDR < 0.05 significance.

Run on a 1:1 ortholog subset, where the two test directions are provably
identical (hypergeometric is symmetric in #successes vs #draws).
"""
import sys
import numpy as np, pandas as pd
from scipy.stats import rankdata, hypergeom
from statsmodels.stats.multitest import multipletests

OUT = sys.argv[1] if len(sys.argv) > 1 else "repro_out"
DENSITY = 0.03

def load(p):
    return pd.read_csv(p, sep="\t", index_col=0)

def mutual_rank(expr):
    C = np.corrcoef(expr.values.astype(np.float64))        # Pearson gene x gene
    R = rankdata(C, axis=1).astype(np.float64)             # ascending row rank
    MR = np.sqrt(R * R.T)
    np.fill_diagonal(MR, 0.0)
    return MR

def density_thr(MR, d):
    up = MR[np.triu_indices(MR.shape[0], k=1)]
    s = np.sort(up)[::-1]
    return s[max(0, round(d * len(s)) - 1)]

s1 = load(f"{OUT}/s1_subset.tsv"); s2 = load(f"{OUT}/s2_subset.tsv")
ortho = pd.read_csv(f"{OUT}/ortho_subset.tsv", sep="\t")

MR1 = mutual_rank(s1); MR2 = mutual_rank(s2)
n1 = (MR1 >= density_thr(MR1, DENSITY)).astype(np.uint8); np.fill_diagonal(n1, 0)
n2 = (MR2 >= density_thr(MR2, DENSITY)).astype(np.uint8); np.fill_diagonal(n2, 0)

# align neighbourhood matrices into 1:1 pair order
pos1 = {g: i for i, g in enumerate(s1.index)}
pos2 = {g: i for i, g in enumerate(s2.index)}
i1 = ortho.Species1.map(pos1).to_numpy()
i2 = ortho.Species2.map(pos2).to_numpy()
A1 = n1[np.ix_(i1, i1)]          # pair x pair : s1 neighbourhood
A2 = n2[np.ix_(i2, i2)]          # pair x pair : s2 neighbourhood
P = len(ortho); N = P

x = (A1 & A2).sum(axis=1).astype(int)      # overlap (both directions equal here)
m = A1.sum(axis=1).astype(int)
k = A2.sum(axis=1).astype(int)
pool = x >= 1                               # min-overlap 1, both directions
praw = np.ones(P)
idx = np.where(pool)[0]
praw[idx] = [hypergeom.sf(int(x[i]) - 1, N, int(m[i]), int(k[i])) for i in idx]
fdr = np.ones(P)
fdr[idx] = multipletests(praw[idx], method="fdr_bh")[1]

call = fdr < 0.05
res = pd.DataFrame({"OrthoGroup": ortho.Ortholog_Group, "Species1": ortho.Species1,
                    "Species2": ortho.Species2, "x": x, "p_raw": praw,
                    "max_fdr": fdr, "coexpressolog": call})
res.to_csv(f"{OUT}/independent_result.tsv", sep="\t", index=False)
print(f"pairs={P}  candidates(pool, x>=1)={pool.sum()}  co-expressologs(FDR<0.05)={int(call.sum())}")
