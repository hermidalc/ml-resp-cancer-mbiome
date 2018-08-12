#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("tools"))
source("config.R")

data_file_prefix <- "YAMS_PD1ML"
pdata_file_basename <- paste0(c(data_file_prefix, "Metadata"), collapse="_")
cat("Loading:", pdata_file_basename, "\n")
pdata_all <- read.delim(paste0("data/", pdata_file_basename, ".txt"), row.names=1)
rownames(pdata_all) <- make.names(rownames(pdata_all))
for (feat_type in feat_types) {
    exprs_file_basename <- paste0(c("YAMS_PD1ML", toupper(feat_type), "PPM"), collapse="_")
    cat("Loading:", exprs_file_basename, "\n")
    eset_all <- ExpressionSet(
        assayData=as.matrix(read.delim(paste0("data/", exprs_file_basename, ".txt"), row.names=1)),
        phenoData=AnnotatedDataFrame(pdata_all),
        featureData=AnnotatedDataFrame(read.delim(paste0("data/YAMS_PD1ML_", toupper(feat_type), "_FeatTable.txt"), row.names=1))
    )
    for (dataset_name in dataset_names) {
        eset_name <- paste0(c("eset", dataset_name, feat_type), collapse="_")
        cat("Creating:", eset_name, "\n")
        eset <- eset_all[, eset_all$Study == toTitleCase(dataset_name)]
        assign(eset_name, eset)
        save(list=eset_name, file=paste0("data/", eset_name, ".Rda"))
    }
}
