#!/usr/bin/env Rscript

options(warn=1)
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("bapred"))
source("config.R")
#
eset_tr_name <- "eset_gajewski_pittsburgh_sd0_wargo_zitvogel_sd0_mgs_ppm_cdd_mrg_tr"
eset_tr_file <- paste0("data/", eset_tr_name, ".Rda")
load(eset_tr_file)
Xtr <- t(exprs(get(eset_tr_name)))
ptr <- pData(get(eset_tr_name))
ytr <- as.factor(ptr$Class + 1)
btr <- ptr$Batch
butr <- sort(unique(btr))
for (j in 1:length(butr)) { if (j != butr[j]) { btr <- replace(btr, btr == butr[j], j) } }
btr <- as.factor(btr)
svg(file="results/pca_4_cohorts_ppm_cdd_v1.svg")
pcplot(Xtr, btr, col=as.numeric(unique(btr)), main="PCA 4 Cohorts PPM CDD")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
pc <- prcomp(Xtr, scale.=FALSE)
svg(file="results/pca_4_cohorts_ppm_cdd_v2.svg")
plot(pc$x[,1], pc$x[,2], col=as.numeric(btr), main="PCA 4 Cohorts PPM CDD", xlab="PC1", ylab="PC2")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
#
eset_tr_name <- "eset_gajewski_pittsburgh_sd0_wargo_zitvogel_sd0_mgs_ppm_ec_mrg_tr"
eset_tr_file <- paste0("data/", eset_tr_name, ".Rda")
load(eset_tr_file)
Xtr <- t(exprs(get(eset_tr_name)))
ptr <- pData(get(eset_tr_name))
ytr <- as.factor(ptr$Class + 1)
btr <- ptr$Batch
butr <- sort(unique(btr))
for (j in 1:length(butr)) { if (j != butr[j]) { btr <- replace(btr, btr == butr[j], j) } }
btr <- as.factor(btr)
svg(file="results/pca_4_cohorts_ppm_ec_v1.svg")
pcplot(Xtr, btr, col=as.numeric(unique(btr)), main="PCA 4 Cohorts PPM EC")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
pc <- prcomp(Xtr, scale.=FALSE)
svg(file="results/pca_4_cohorts_ppm_ec_v2.svg")
plot(pc$x[,1], pc$x[,2], col=as.numeric(btr), main="PCA 4 Cohorts PPM EC", xlab="PC1", ylab="PC2")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
#
eset_tr_name <- "eset_gajewski_pittsburgh_sd0_wargo_zitvogel_sd0_mgs_ppm_go_mrg_tr"
eset_tr_file <- paste0("data/", eset_tr_name, ".Rda")
load(eset_tr_file)
Xtr <- t(exprs(get(eset_tr_name)))
ptr <- pData(get(eset_tr_name))
ytr <- as.factor(ptr$Class + 1)
btr <- ptr$Batch
butr <- sort(unique(btr))
for (j in 1:length(butr)) { if (j != butr[j]) { btr <- replace(btr, btr == butr[j], j) } }
btr <- as.factor(btr)
svg(file="results/pca_4_cohorts_ppm_go_v1.svg")
pcplot(Xtr, btr, col=as.numeric(unique(btr)), main="PCA 4 Cohorts PPM GO")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
pc <- prcomp(Xtr, scale.=FALSE)
svg(file="results/pca_4_cohorts_ppm_go_v2.svg")
plot(pc$x[,1], pc$x[,2], col=as.numeric(btr), main="PCA 4 Cohorts PPM GO", xlab="PC1", ylab="PC2")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
#
eset_tr_name <- "eset_gajewski_pittsburgh_sd0_wargo_zitvogel_sd0_mgs_ppm_lkt_mrg_tr"
eset_tr_file <- paste0("data/", eset_tr_name, ".Rda")
load(eset_tr_file)
Xtr <- t(exprs(get(eset_tr_name)))
ptr <- pData(get(eset_tr_name))
ytr <- as.factor(ptr$Class + 1)
btr <- ptr$Batch
butr <- sort(unique(btr))
for (j in 1:length(butr)) { if (j != butr[j]) { btr <- replace(btr, btr == butr[j], j) } }
btr <- as.factor(btr)
svg(file="results/pca_4_cohorts_ppm_lkt_v1.svg")
pcplot(Xtr, btr, col=as.numeric(unique(btr)), main="PCA 4 Cohorts PPM LKT")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
pc <- prcomp(Xtr, scale.=FALSE)
svg(file="results/pca_4_cohorts_ppm_lkt_v2.svg")
plot(pc$x[,1], pc$x[,2], col=as.numeric(btr), main="PCA 4 Cohorts PPM LKT", xlab="PC1", ylab="PC2")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
#
eset_tr_name <- "eset_gajewski_pittsburgh_sd0_wargo_zitvogel_sd0_mgs_ppm_pfam_mrg_tr"
eset_tr_file <- paste0("data/", eset_tr_name, ".Rda")
load(eset_tr_file)
Xtr <- t(exprs(get(eset_tr_name)))
ptr <- pData(get(eset_tr_name))
ytr <- as.factor(ptr$Class + 1)
btr <- ptr$Batch
butr <- sort(unique(btr))
for (j in 1:length(butr)) { if (j != butr[j]) { btr <- replace(btr, btr == butr[j], j) } }
btr <- as.factor(btr)
svg(file="results/pca_4_cohorts_ppm_pfam_v1.svg")
pcplot(Xtr, btr, col=as.numeric(unique(btr)), main="PCA 4 Cohorts PPM PFAM")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
pc <- prcomp(Xtr, scale.=FALSE)
svg(file="results/pca_4_cohorts_ppm_pfam_v2.svg")
plot(pc$x[,1], pc$x[,2], col=as.numeric(btr), main="PCA 4 Cohorts PPM PFAM", xlab="PC1", ylab="PC2")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
#
eset_tr_name <- "eset_gajewski_pittsburgh_sd0_wargo_zitvogel_sd0_mgs_ppm_prints_mrg_tr"
eset_tr_file <- paste0("data/", eset_tr_name, ".Rda")
load(eset_tr_file)
Xtr <- t(exprs(get(eset_tr_name)))
ptr <- pData(get(eset_tr_name))
ytr <- as.factor(ptr$Class + 1)
btr <- ptr$Batch
butr <- sort(unique(btr))
for (j in 1:length(butr)) { if (j != butr[j]) { btr <- replace(btr, btr == butr[j], j) } }
btr <- as.factor(btr)
svg(file="results/pca_4_cohorts_ppm_prints_v1.svg")
pcplot(Xtr, btr, col=as.numeric(unique(btr)), main="PCA 4 Cohorts PPM PRINTS")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
pc <- prcomp(Xtr, scale.=FALSE)
svg(file="results/pca_4_cohorts_ppm_prints_v2.svg")
plot(pc$x[,1], pc$x[,2], col=as.numeric(btr), main="PCA 4 Cohorts PPM PRINTS", xlab="PC1", ylab="PC2")
legend("bottomleft", legend=c("gajewski","pittsburgh_sd0","wargo","zitvogel_sd0"), col=as.numeric(unique(btr)), pch=20, bty="n")
dev.off()
#
