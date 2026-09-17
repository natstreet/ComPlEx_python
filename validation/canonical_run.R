# Standalone runner of the CANONICAL rcomplex algorithm (verbatim from
# validate_against_canonical_rcomplex.R / RComPlEx.Rmd), writing its own
# co-expressolog call set for independent comparison.
suppressPackageStartupMessages(library(dplyr))
V <- Sys.getenv("VALIDATION_DIR", "repro_out")
s1 <- read.table(file.path(V,"s1_subset.tsv"), sep="\t", header=TRUE, row.names=1, check.names=FALSE)
s2 <- read.table(file.path(V,"s2_subset.tsv"), sep="\t", header=TRUE, row.names=1, check.names=FALSE)
ortho <- read.table(file.path(V,"ortho_subset.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat(sprintf("subset: %d S1, %d S2, %d pairs\n", nrow(s1), nrow(s2), nrow(ortho)))
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
cat(sprintf("CANONICAL: %d candidate pairs, %d co-expressologs (Max.FDR<0.05)\n",
            nrow(comparison), sum(comparison$maxfdr<0.05)))
write.table(comparison[,c("Species1","Species2","s1p","s2p","maxfdr")],
            file.path(V,"canonical_result.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
