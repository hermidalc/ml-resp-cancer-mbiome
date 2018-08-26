#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("tools"))
source("config.R")

parser <- ArgumentParser()
parser$add_argument("--base-datasets", type="character", nargs="+", help="base datasets")
parser$add_argument("--class-type", type="character", nargs="+", help="class type")
parser$add_argument("--data-type", type="character", nargs="+", help="data type")
parser$add_argument("--norm-meth", type="character", nargs="+", help="normalization method")
parser$add_argument("--feat-type", type="character", nargs="+", help="feature type")
args <- parser$parse_args()
if (!is.null(args$base_datasets)) {
    dataset_basenames <- intersect(dataset_basenames, args$base_datasets)
}
if (!is.null(args$class_type)) {
    class_types <- intersect(class_types, args$class_type)
}
if (!is.null(args$data_type)) {
    data_types <- intersect(data_types, args$data_type)
}
if (!is.null(args$norm_meth)) {
    norm_methods <- intersect(norm_methods, args$norm_meth)
}
if (!is.null(args$feat_type)) {
    feat_types <- intersect(feat_types, args$feat_type)
}
for (dataset_group in dataset_groups) {
    for (data_type in data_types) {
        pdata_file_basename <- paste0(c(toupper(dataset_group), "Metadata"), collapse="_")
        pdata_file <- paste0("data/", pdata_file_basename, ".txt")
        for (norm_meth in norm_methods) {
            if (norm_meth == 'none') next
            for (feat_type in feat_types) {
                if (feat_type == 'none') next
                suffixes <- c(data_type, norm_meth, feat_type)
                exprs_file_basename <- paste0(toupper(c(dataset_group, feat_type, norm_meth)), collapse="_")
                exprs_file <- paste0("data/", exprs_file_basename, ".txt")
                fdata_file_basename <- paste0(c(toupper(c(dataset_group, feat_type)), "FeatTable"), collapse="_")
                fdata_file <- paste0("data/", fdata_file_basename, ".txt")
                if (file.exists(pdata_file) && file.exists(exprs_file) && file.exists(fdata_file)) {
                    if (!exists("pdata")) {
                        cat(" Loading:", pdata_file_basename, "\n")
                        pdata <- read.delim(pdata_file, row.names=1)
                        rownames(pdata) <- make.names(rownames(pdata))
                    }
                    cat(" Loading:", exprs_file_basename, "\n")
                    cat(" Loading:", fdata_file_basename, "\n")
                    eset_all <- ExpressionSet(
                        assayData=as.matrix(read.delim(exprs_file, row.names=1)),
                        phenoData=AnnotatedDataFrame(pdata),
                        featureData=AnnotatedDataFrame(read.delim(fdata_file, row.names=1))
                    )
                    eset_all <- eset_all[rowSums(exprs(eset_all)) > 0,]
                    for (dataset_basename in dataset_basenames) {
                        for (class_type in class_types) {
                            eset_name <- paste0(c("eset", dataset_basename, class_type, suffixes), collapse="_")
                            cat("Creating:", eset_name, "\n")
                            eset <- eset_all[, eset_all$Study == toTitleCase(dataset_basename)]
                            for (label in varLabels(eset)) {
                                if (is.factor(eset[[label]])) eset[[label]] <- droplevels(eset[[label]])
                            }
                            if (class_type == "sd0") {
                                eset$Class <- ifelse(
                                    eset$Response %in% class_info$pos, 1, ifelse(
                                        eset$Response %in% c(class_info$neg, class_info$sd), 0, NA
                                    )
                                )
                            }
                            else if (class_type == "sd1") {
                                eset$Class <- ifelse(
                                    eset$Response %in% class_info$neg, 0, ifelse(
                                        eset$Response %in% c(class_info$pos, class_info$sd), 1, NA
                                    )
                                )
                            }
                            assign(eset_name, eset)
                            save(list=eset_name, file=paste0("data/", eset_name, ".Rda"))
                        }
                    }
                }
            }
        }
        if (exists("pdata")) remove(pdata)
    }
}
