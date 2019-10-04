#!/usr/bin/env Rscript

options(warn=1)
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Biobase"))
source("config.R")

parser <- ArgumentParser()
parser$add_argument(
    "--dataset-basename", type="character", nargs="+", help="dataset basename"
)
parser$add_argument(
    "--class-type", type="character", nargs="+", help="class type"
)
parser$add_argument(
    "--data-type", type="character", nargs="+", help="data type"
)
parser$add_argument(
    "--norm-meth", type="character", nargs="+", help="normalization method"
)
parser$add_argument(
    "--feat-type", type="character", nargs="+", help="feature type"
)
args <- parser$parse_args()
all_dataset_basenames <- dataset_basenames
if (!is.null(args$database_basename)) {
    dataset_basenames <- intersect(dataset_basenames, args$database_basename)
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
        pdata_file_basename <- paste0(
            c(dataset_group, "metadata"), collapse="_"
        )
        pdata_file <- paste0("data/", pdata_file_basename, ".tsv")
        for (norm_meth in norm_methods) {
            for (feat_type in feat_types) {
                suffixes <- c(data_type)
                for (suffix in c(norm_meth, feat_type)) {
                    if (suffix != "none") suffixes <- c(suffixes, suffix)
                }
                exprs_file_basename <- paste0(
                    c(dataset_group, feat_type, norm_meth), collapse="_"
                )
                exprs_file <- paste0("data/", exprs_file_basename, ".tsv")
                fdata_file_basename <- paste0(
                    c(dataset_group, feat_type, "features"), collapse="_"
                )
                fdata_file <- paste0("data/", fdata_file_basename, ".tsv")
                if (
                    file.exists(pdata_file) && file.exists(exprs_file) &&
                    file.exists(fdata_file)
                ) {
                    if (!exists("pdata")) {
                        cat(" Loading:", pdata_file_basename, "\n")
                        pdata <- read.delim(pdata_file, row.names=1)
                        rownames(pdata) <- make.names(rownames(pdata))
                    }
                    cat(" Loading:", exprs_file_basename, "\n")
                    counts <- as.matrix(read.delim(exprs_file, row.names=1))
                    counts <- ceiling(counts / max_read_length)
                    mode(counts) <- "integer"
                    cat(" Loading:", fdata_file_basename, "\n")
                    if (feat_type == "lkt") {
                        fdata <- read.delim(fdata_file, row.names=NULL)
                        row.names(fdata) <- as.vector(fdata[, ncol(fdata)])
                        fdata[, ncol(fdata)] <- NULL
                    } else {
                        fdata <- read.delim(fdata_file, row.names=1)
                    }
                    if (nrow(pdata) > ncol(counts)) {
                        cat(
                            "  Filter:", (nrow(pdata) - ncol(counts)),
                            "metadata samples\n"
                        )
                        pdata_feat_type <- pdata[colnames(counts), , drop=FALSE]
                    } else if (nrow(pdata) < ncol(counts)) {
                        cat(
                            "  Filter:", (ncol(counts) - nrow(pdata)),
                            "data matrix samples\n"
                        )
                        counts <- counts[, rownames(pdata), drop=FALSE]
                    } else {
                        pdata_feat_type <- pdata
                    }
                    eset_all <- ExpressionSet(
                        assayData=counts,
                        phenoData=AnnotatedDataFrame(pdata_feat_type),
                        featureData=AnnotatedDataFrame(fdata)
                    )
                    for (dataset_basename in dataset_basenames) {
                        if (dataset_basename %in% dataset_names) {
                            eset_name <- paste0(
                                c("eset", dataset_basename, suffixes),
                                collapse="_"
                            )
                            cat("Creating:", eset_name, "\n")
                            eset <- eset_all[,
                                tolower(eset_all$Study) == dataset_basename
                            ]
                            keep <- grep(
                                "_(none|_unclassified)$",
                                row.names(eset), perl=TRUE, invert=TRUE,
                                ignore.case=TRUE
                            )
                            if (dim(eset)[1] - length(keep) > 0) {
                                cat(
                                    "  Filter:",
                                    (dim(eset)[1] - length(keep)),
                                    "unmapped feature\n"
                                )
                                eset <- eset[keep, ]
                            }
                            # for (label in varLabels(eset)) {
                            #     if (is.factor(eset[[label]]))
                            #         eset[[label]] <- droplevels(eset[[label]])
                            # }
                            if (dim(eset)[2] > 0) {
                                eset$Class <- ifelse(
                                    eset$Clin_Response %in% class_info$pos, 1,
                                    ifelse(
                                        eset$Clin_Response %in% class_info$neg,
                                        0, NA
                                    )
                                )
                                eset <- eset[, !is.na(eset$Class)]
                                eset$Batch <- rep(which(
                                    dataset_basename == all_dataset_basenames
                                )[1])
                                assign(eset_name, eset)
                                save(
                                    list=eset_name,
                                    file=paste0("data/", eset_name, ".Rda")
                                )
                            } else {
                                cat("No samples\n")
                            }
                        }
                        else {
                            for (class_type in class_types) {
                                eset_name <- paste0(c(
                                    "eset", dataset_basename, class_type,
                                    suffixes
                                ), collapse="_")
                                cat("Creating:", eset_name, "\n")
                                eset <- eset_all[,
                                    tolower(eset_all$Study) == dataset_basename
                                ]
                                keep <- grep(
                                    "_(none|_unclassified)$",
                                    row.names(eset), perl=TRUE, invert=TRUE,
                                    ignore.case=TRUE
                                )
                                if (dim(eset)[1] - length(keep) > 0) {
                                    cat(
                                        "  Filter:",
                                        (dim(eset)[1] - length(keep)),
                                        "unmapped feature\n"
                                    )
                                    eset <- eset[keep, ]
                                }
                                # for (label in varLabels(eset)) {
                                #     if (is.factor(eset[[label]]))
                                #         eset[[label]] <-
                                #             droplevels(eset[[label]])
                                # }
                                if (class_type == "sd0") {
                                    eset$Class <- ifelse(
                                        eset$Clin_Response %in% class_info$pos,
                                        1, ifelse(
                                            eset$Clin_Response %in%
                                            c(class_info$neg, class_info$sd),
                                            0, NA
                                        )
                                    )
                                }
                                else if (class_type == "sd1") {
                                    eset$Class <- ifelse(
                                        eset$Clin_Response %in% class_info$neg,
                                        0, ifelse(
                                            eset$Clin_Response %in%
                                            c(class_info$pos, class_info$sd),
                                            1, NA
                                        )
                                    )
                                }
                                eset <- eset[, !is.na(eset$Class)]
                                eset$Batch <- rep(which(
                                    dataset_basename == all_dataset_basenames
                                )[1])
                                assign(eset_name, eset)
                                save(
                                    list=eset_name,
                                    file=paste0("data/", eset_name, ".Rda")
                                )
                            }
                        }
                    }
                }
            }
        }
        if (exists("pdata")) remove(pdata)
    }
}
