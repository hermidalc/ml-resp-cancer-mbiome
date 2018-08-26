# config
dataset_groups <- c(
    "yams_pd1ml"
)
dataset_basenames <- c(
    "gajewski",
    "pittsburgh",
    "wargo",
    "zitvogel"
)
class_types <- c(
    "sd0",
    "sd1"
)
dataset_names <- as.character(
    sapply(dataset_basenames, FUN=function(x) paste0(x, "-", class_types))
)
rna_seq_dataset_names <- c(
)
needs_log2_dataset_names <- c(
)
class_info <- list(
    pos = c(
        "Responder",
        "CR",
        "Complete response",
        "PR",
        "Partial response"
    ),
    neg = c(
        "Non-Responder",
        "PD",
        "Progression",
        "Dead"
    ),
    sd = c(
        "SD",
        "Stable"
    )
)
data_types <- c(
    "mgs"
)
norm_methods <- c(
    "ppm"
)
feat_types <- c(
    "cdd",
    "ec",
    "go",
    "lkt",
    "pfam",
    "prints"
)
prep_methods <- c(
    "none",
    "cff",
    "mrg"
)
bc_methods <- c(
    "none",
    "ctr",
    "std",
    "rta",
    "rtg",
    "qnorm",
    "cbt",
    "fab",
    "sva",
    "stica0",
    "stica025",
    "stica05",
    "stica1",
    "svd"
)
filter_types <- c(
    "none"
)
common_pheno_names <- c(
    "Study",
    "Cohort",
    "PatientID",
    "Batch",
    "Class",
    "Response",
    "Timepoint",
    "Clin_Response",
    "Host_disease_status",
    "GbNAHS",
    "PctAss"
)
matfact_k <- 20
