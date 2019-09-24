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
dataset_names <- c(
    "gajewski",
    "pittsburgh_sd0",
    "pittsburgh_sd1",
    "wargo",
    "zitvogel_sd0",
    "zitvogel_sd1"
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
    "none",
    "counts",
    "ppm",
    "tmm",
    "vst"
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
filt_types <- c(
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
pdata_cls_name <- "Class"
matfact_k <- 20
