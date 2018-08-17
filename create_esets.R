#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("tools"))
source("config.R")

parser <- ArgumentParser()
parser$add_argument("--datasets", type="character", nargs="+", help="datasets")
parser$add_argument("--feat-type", type="character", nargs="+", help="dataset feature type")
parser$add_argument("--norm-meth", type="character", nargs="+", help="preprocessing/normalization method")
args <- parser$parse_args()
if (!is.null(args$datasets)) {
    dataset_names <- intersect(dataset_names, args$datasets)
}
if (!is.null(args$feat_type)) {
    feat_types <- feat_types[feat_types %in% args$feat_type]
}
if (!is.null(args$norm_meth)) {
    norm_methods <- norm_methods[norm_methods %in% args$norm_meth]
}
pdata_file_basename <- paste0(c(toupper(data_file_prefix), "Metadata"), collapse="_")
cat("Loading:", pdata_file_basename, "\n")
pdata_all <- read.delim(paste0("data/", pdata_file_basename, ".txt"), row.names=1)
rownames(pdata_all) <- make.names(rownames(pdata_all))
for (feat_type in feat_types) {
    suffixes <- c(feat_type)
    for (norm_meth in norm_methods) {
        if (norm_meth != "none") suffixes <- c(suffixes, norm_meth)
        exprs_file_basename <- paste0(toupper(c(data_file_prefix, suffixes)), collapse="_")
        exprs_file <- paste0("data/", exprs_file_basename, ".txt")
        fdata_file_basename <- paste0(c(toupper(c(data_file_prefix, feat_type)), "FeatTable"), collapse="_")
        fdata_file <- paste0("data/", fdata_file_basename, ".txt")
        if (file.exists(exprs_file) && file.exists(fdata_file)) {
            cat("Loading:", exprs_file_basename, fdata_file_basename, "\n")
            eset_all <- ExpressionSet(
                assayData=as.matrix(read.delim(exprs_file, row.names=1)),
                phenoData=AnnotatedDataFrame(pdata_all),
                featureData=AnnotatedDataFrame(read.delim(fdata_file, row.names=1))
            )
            for (dataset_name in dataset_names) {
                eset_name <- paste0(c("eset", dataset_name, suffixes), collapse="_")
                cat("Creating:", eset_name, "\n")
                eset <- eset_all[, eset_all$Study == toTitleCase(dataset_name)]
                eset <- eset[rowSums(exprs(eset)) > 0, ]
                assign(eset_name, eset)
                save(list=eset_name, file=paste0("data/", eset_name, ".Rda"))
            }
        }
    }
}
