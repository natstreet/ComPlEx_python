# validate_against_canonical_rcomplex.R
#
# UNBIASED cross-check: compares complex_py.py against the CANONICAL published
# ComPlEx (Torgeir Hvidsten's rcomplex, https://gitlab.com/hvidsten-lab/rcomplex,
# RComPlEx.Rmd), NOT the in-house transcription in validate_complex.R.
#
# The canonical algorithm below is copied verbatim from the "coexpression" and
# "comparison" chunks of RComPlEx.Rmd (Pearson -> MR sqrt(R*t(R)) -> 3%% density
# threshold -> phyper(x-1, m, N-m, k) for x>1 -> filter overlap>0 in both -> BH FDR
# -> Max.p < 0.05). It is run on the same validation subset that validate_complex.py
# writes, and its co-expressolog calls + raw p-values are compared to Python.
#
# Result (1500-gene cold-needle subset, seed 42): canonical = python = 28
# co-expressologs, Jaccard = 1.000, only-canonical = 0, only-python = 0; raw
# per-pair p-values agree to max abs diff 1.0e-15. The Python reimplementation
# therefore reproduces the canonical published ComPlEx exactly.
#
# Usage: run validate_complex.py first (writes the subset to its OUT_DIR), set V
# below to that directory, then: Rscript validate_against_canonical_rcomplex.R
#
suppressPackageStartupMessages(library(dplyr))
V <- Sys.getenv("VALIDATION_DIR", "complex_validation_output")  # dir written by validate_complex.py
s1 <- read.table(file.path(V,"s1_subset.tsv"), sep="\t", header=TRUE, row.names=1, check.names=FALSE)
s2 <- read.table(file.path(V,"s2_subset.tsv"), sep="\t", header=TRUE, row.names=1, check.names=FALSE)
ortho <- read.table(file.path(V,"ortho_subset.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
if ("Ortholog_Group" %in% names(ortho)) ortho$OrthoGroup <- ortho$Ortholog_Group
cat(sprintf("subset: %d S1 genes, %d S2 genes, %d pairs\n", nrow(s1), nrow(s2), nrow(ortho)))
# ===== CANONICAL rcomplex algorithm (verbatim from RComPlEx.Rmd) =====
density_thr <- 0.03
species1_net <- cor(t(s1), method="pearson"); dimnames(species1_net) <- list(rownames(s1), rownames(s1))
species2_net <- cor(t(s2), method="pearson"); dimnames(species2_net) <- list(rownames(s2), rownames(s2))
R <- t(apply(species1_net,1,rank)); species1_net <- sqrt(R * t(R))
R <- t(apply(species2_net,1,rank)); species2_net <- sqrt(R * t(R))
diag(species1_net) <- 0; diag(species2_net) <- 0
R <- sort(species1_net[upper.tri(species1_net, diag=FALSE)], decreasing=TRUE); species1_thr <- R[round(density_thr*length(R))]
R <- sort(species2_net[upper.tri(species2_net, diag=FALSE)], decreasing=TRUE); species2_thr <- R[round(density_thr*length(R))]
comparison <- ortho
comparison$s1o <- NA; comparison$s1p <- NA; comparison$s2o <- NA; comparison$s2p <- NA
for (i in 1:nrow(ortho)) {
  neigh <- names(which(species1_net[ortho$Species1[i],] >= species1_thr))
  on <- names(which(species2_net[ortho$Species2[i],] >= species2_thr))
  on <- unique(ortho$Species1[ortho$Species2 %in% on])
  N<-nrow(s1); m<-length(neigh); k<-length(on); x<-length(intersect(neigh,on))
  comparison$s1o[i]<-x; comparison$s1p[i]<- if(x>1) phyper(x-1,m,N-m,k,lower.tail=FALSE) else 1
  neigh <- names(which(species2_net[ortho$Species2[i],] >= species2_thr))
  on <- names(which(species1_net[ortho$Species1[i],] >= species1_thr))
  on <- unique(ortho$Species2[ortho$Species1 %in% on])
  N<-nrow(s2); m<-length(neigh); k<-length(on); x<-length(intersect(neigh,on))
  comparison$s2o[i]<-x; comparison$s2p[i]<- if(x>1) phyper(x-1,m,N-m,k,lower.tail=FALSE) else 1
}
comparison <- comparison %>% filter(s1o>0 & s2o>0)
comparison$s1fdr <- p.adjust(comparison$s1p, method="fdr")
comparison$s2fdr <- p.adjust(comparison$s2p, method="fdr")
comparison$maxfdr <- pmax(comparison$s1fdr, comparison$s2fdr)
cat(sprintf("CANONICAL: %d candidate pairs, %d co-expressologs (Max.p<0.05)\n", nrow(comparison), sum(comparison$maxfdr<0.05)))
# ===== compare to Python (R-faithful) =====
py <- read.table(file.path(V,"python_optimised.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat(sprintf("PYTHON:    %d candidate pairs, %d co-expressologs\n", nrow(py), sum(py$max_fdr<0.05)))
key <- function(d,a,b) paste(d[[a]], d[[b]])
canon_sig <- comparison$maxfdr<0.05; py_sig <- py$max_fdr<0.05
cs <- key(comparison[canon_sig,],"Species1","Species2"); ps <- key(py[py_sig,],"Species1","Species2")
cat(sprintf("\nco-expressolog set: canonical=%d python=%d shared=%d only-canon=%d only-py=%d Jaccard=%.4f\n",
    length(cs), length(ps), length(intersect(cs,ps)), length(setdiff(cs,ps)), length(setdiff(ps,cs)),
    length(intersect(cs,ps))/length(union(cs,ps))))
# raw p-value agreement on shared pairs
m <- merge(comparison[,c("Species1","Species2","s1p","s2p")], py[,c("Species1","Species2","p1_raw","p2_raw")], by=c("Species1","Species2"))
d1 <- abs(m$s1p - m$p1_raw); d2 <- abs(m$s2p - m$p2_raw)
cat(sprintf("raw p-value agreement (%d shared pairs): max abs diff = %.2e\n", nrow(m), max(c(d1,d2), na.rm=TRUE)))
