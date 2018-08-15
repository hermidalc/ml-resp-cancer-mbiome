#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Biobase"))
source("config.R")

parser <- ArgumentParser()
parser$add_argument("--num-tr-combo", type="integer", help="num datasets to combine")
parser$add_argument("--norm-meth", type="character", nargs="+", help="normalization method")
parser$add_argument("--feat-type", type="character", nargs="+", help="dataset feat type")
parser$add_argument("--load-only", action="store_true", default=FALSE, help="show search and eset load only")
args <- parser$parse_args()

num_tr_combo <- as.integer(args$num_tr_combo)
if (!is.null(args$norm_meth)) {
    norm_methods <- norm_methods[norm_methods %in% args$norm_meth]
}
if (!is.null(args$feat_type)) {
    feat_types <- feat_types[feat_types %in% args$feat_type]
}
for (norm_meth in norm_methods) {
    for (feat_type in feat_types) {
        suffixes <- c(norm_meth)
        if (feat_type != "none") suffixes <- c(suffixes, feat_type)
        for (dataset_name in dataset_names) {
            eset_name <- paste0(c("eset", dataset_name, suffixes), collapse="_")
            eset_file <- paste0("data/", eset_name, ".Rda")
            if (file.exists(eset_file)) {
                cat("Loading:", eset_name, "\n")
                load(eset_file)
            }
        }
        if (args$load_only) next
        dataset_name_combos <- combn(dataset_names, num_tr_combo)
        for (col in 1:ncol(dataset_name_combos)) {
            eset_merged_name <- paste0(c("eset", dataset_name_combos[,col], suffixes, "merged", "tr"), collapse="_")
            eset_1_name <- paste0(c("eset", dataset_name_combos[1,col], suffixes), collapse="_")
            eset_2_name <- paste0(c("eset", dataset_name_combos[2,col], suffixes), collapse="_")
            if (exists(eset_1_name) & exists(eset_2_name)) {
                cat("Creating:", eset_merged_name, "\n")
                eset_merged <- combine(get(eset_1_name), get(eset_2_name))
                if (nrow(dataset_name_combos) > 2) {
                    for (row in 3:nrow(dataset_name_combos)) {
                        eset_n_name <- paste0(c("eset", dataset_name_combos[row,col], suffixes), collapse="_")
                        eset_merged <- combine(eset_merged, get(eset_n_name))
                    }
                }
                assign(eset_merged_name, eset_merged)
                save(list=eset_merged_name, file=paste0("data/", eset_merged_name, ".Rda"))
                remove(list=c(eset_merged_name))
            }
        }
    }
}