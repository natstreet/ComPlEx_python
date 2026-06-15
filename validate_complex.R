#!/usr/bin/env Rscript
# validate_complex.R
#
# Runs the original ComPlEx.R algorithm (exact copy of the loop logic) on the
# subset produced by validate_complex.py and writes r_comparison.tsv, plus the
# per-gene centrality (r_centrality_s1/s2.tsv) for the degree comparison.
# Uses x1>0 AND x2>0 post-loop filter as in the published R code.
#
# Usage:
#   1. Run validate_complex.py first to generate the subset TSVs in OUT_DIR.
#   2. Set OUT_DIR below to match the OUT_DIR in validate_complex.py.
#   3. Rscript validate_complex.R

suppressPackageStartupMessages({
  library(dplyr)
  library(matrixStats)
})

OUT_DIR <- "complex_validation_output"   # must match OUT_DIR in validate_complex.py
DENSITY <- 0.03

# ── Load subset produced by Python script ─────────────────────────────────────
cat("Loading subset...\n")
s1 <- read.table(file.path(OUT_DIR, "s1_subset.tsv"), sep="\t",
                 header=TRUE, row.names=1, check.names=FALSE)
s2 <- read.table(file.path(OUT_DIR, "s2_subset.tsv"), sep="\t",
                 header=TRUE, row.names=1, check.names=FALSE)
ortho <- read.table(file.path(OUT_DIR, "ortho_subset.tsv"), sep="\t",
                    header=TRUE, stringsAsFactors=FALSE)
# Column names: Ortholog_Group, Species1, Species2

cat(sprintf("  %d S1 genes  %d S2 genes  %d pairs\n",
            nrow(s1), nrow(s2), nrow(ortho)))

# ── Pearson correlation (R uses float64 by default) ───────────────────────────
cat("Computing Pearson correlations...\n")
s1_mat <- as.matrix(s1)
s2_mat <- as.matrix(s2)

# Remove zero-variance rows (should already be clean, but mirror R ComPlEx)
s1_mat <- s1_mat[rowSds(s1_mat) > 0, ]
s2_mat <- s2_mat[rowSds(s2_mat) > 0, ]

net1 <- cor(t(s1_mat), method="pearson")
dimnames(net1) <- list(rownames(s1_mat), rownames(s1_mat))
net2 <- cor(t(s2_mat), method="pearson")
dimnames(net2) <- list(rownames(s2_mat), rownames(s2_mat))

# ── Mutual Rank normalisation (exact R ComPlEx code) ─────────────────────────
cat("Computing Mutual Rank...\n")
R1 <- t(apply(net1, 1, rank))          # row-wise rank of correlations
net1_mr <- sqrt(R1 * t(R1))            # geometric mean of forward/reverse rank
diag(net1_mr) <- 0

R2 <- t(apply(net2, 1, rank))
net2_mr <- sqrt(R2 * t(R2))
diag(net2_mr) <- 0

# ── R centrality (per-gene degree), canonical rcomplex density-threshold ──────
# Matches the threshold pattern in Torgeir Hvidsten's canonical rcomplex
# (RComPlEx.Rmd): sort the mutual-rank upper triangle in decreasing order and take
# the value at the density-fraction position (here 1% for centrality). Degree is the
# number of MR edges >= that threshold per gene (diagonal is 0 and so excluded).
# Written for the Python-vs-R centrality comparison in validate_complex.py.
cat("Computing R centrality (degree, canonical sorted density threshold)...\n")
r_degree <- function(net_mr, genes, density = 0.01) {
  upper <- sort(net_mr[upper.tri(net_mr, diag = FALSE)], decreasing = TRUE)
  thr   <- upper[round(density * length(upper))]
  data.frame(Genes = genes, Degree = rowSums(net_mr >= thr))
}
write.table(r_degree(net1_mr, rownames(s1_mat)),
            file.path(OUT_DIR, "r_centrality_s1.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(r_degree(net2_mr, rownames(s2_mat)),
            file.path(OUT_DIR, "r_centrality_s2.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ── Density threshold ─────────────────────────────────────────────────────────
upper1 <- net1_mr[upper.tri(net1_mr)]
thr1   <- sort(upper1, decreasing=TRUE)[round(DENSITY * length(upper1))]

upper2 <- net2_mr[upper.tri(net2_mr)]
thr2   <- sort(upper2, decreasing=TRUE)[round(DENSITY * length(upper2))]

cat(sprintf("  S1 threshold: %.4f   S2 threshold: %.4f\n", thr1, thr2))

N1 <- nrow(net1_mr)
N2 <- nrow(net2_mr)

# ── Per-pair hypergeometric test (exact R ComPlEx loop) ───────────────────────
cat(sprintf("Running comparison loop over %d pairs...\n", nrow(ortho)))
n_pairs <- nrow(ortho)

x1_vec <- integer(n_pairs);  m1_vec <- integer(n_pairs);  k1_vec <- integer(n_pairs)
x2_vec <- integer(n_pairs);  m2_vec <- integer(n_pairs);  k2_vec <- integer(n_pairs)
p1_vec <- numeric(n_pairs);  p2_vec <- numeric(n_pairs)

for (i in seq_len(n_pairs)) {
  if (i %% 1000 == 0) cat(sprintf("  %d/%d\n", i, n_pairs))

  g1 <- ortho$Species1[i];  g2 <- ortho$Species2[i]

  # S1 → S2
  neigh1  <- names(net1_mr[g1, net1_mr[g1,] >= thr1])
  neigh2  <- names(net2_mr[g2, net2_mr[g2,] >= thr2])
  # Map S2 neighbourhood of g2 back to S1 gene names via ortho table
  ortho_n1 <- unique(ortho$Species1[ortho$Species2 %in% neigh2])

  m1 <- length(neigh1)
  k1 <- length(ortho_n1)
  x1 <- length(intersect(neigh1, ortho_n1))
  p1 <- ifelse(x1 > 1, phyper(x1-1, m1, N1-m1, k1, lower.tail=FALSE), 1.0)

  # S2 → S1
  # Map S1 neighbourhood of g1 back to S2 gene names via ortho table
  ortho_n2 <- unique(ortho$Species2[ortho$Species1 %in% neigh1])

  m2 <- length(neigh2)
  k2 <- length(ortho_n2)
  x2 <- length(intersect(neigh2, ortho_n2))
  p2 <- ifelse(x2 > 1, phyper(x2-1, m2, N2-m2, k2, lower.tail=FALSE), 1.0)

  x1_vec[i] <- x1;  m1_vec[i] <- m1;  k1_vec[i] <- k1;  p1_vec[i] <- p1
  x2_vec[i] <- x2;  m2_vec[i] <- m2;  k2_vec[i] <- k2;  p2_vec[i] <- p2
}

comparison <- data.frame(
  OrthoGroup = ortho$Ortholog_Group,
  Species1   = ortho$Species1,
  Species2   = ortho$Species2,
  x1 = x1_vec, m1 = m1_vec, k1 = k1_vec, p1_raw = p1_vec,
  x2 = x2_vec, m2 = m2_vec, k2 = k2_vec, p2_raw = p2_vec
)

# Apply the R ComPlEx post-loop filter: keep x1>0 AND x2>0
before <- nrow(comparison)
comparison <- comparison[comparison$x1 > 0 & comparison$x2 > 0, ]
cat(sprintf("  R filter (x1>0 AND x2>0): %d → %d pairs\n", before, nrow(comparison)))

# BH FDR correction
comparison$p1_fdr <- p.adjust(comparison$p1_raw, method="BH")
comparison$p2_fdr <- p.adjust(comparison$p2_raw, method="BH")
comparison$max_fdr <- pmax(comparison$p1_fdr, comparison$p2_fdr)
comparison <- comparison[order(comparison$max_fdr), ]

cat(sprintf("  FDR<0.05: %d pairs\n", sum(comparison$max_fdr < 0.05)))

write.table(comparison,
            file.path(OUT_DIR, "r_comparison.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("Saved to %s/r_comparison.tsv\n", OUT_DIR))
