# GSE50760 sample design preparation.
#
# The local metadata contains matched Stage IV colorectal cancer trios:
# primary tumor, adjacent/normal-looking colon epithelium, and liver metastasis.
# This script creates GSE50760_clinical_edit.csv with analysis_ columns that can
# be read directly by the downstream limma scripts.


# 0. Config -------------------------------------------------------------------

DATASET_ID <- "GSE50760"
DATA_TYPE <- "ngs"

SE_RDS_FILE <- "data/ngs/GSE50760/data_prepare/GSE50760_se_raw.rds"
CLINICAL_RAW_FILE <- "data/ngs/GSE50760/data_prepare/GSE50760_clinical_raw.csv"
CLINICAL_EDIT_FILE <- "data/ngs/GSE50760/data_prepare/GSE50760_clinical_edit.csv"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"

OUTPUT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID, "tables", "sample_design")

options(width = 200)


# 1. Packages and helpers ------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

source(REPORT_TABLE_FUNCTION_FILE)

extract_patient_id <- function(title) {
  title <- as.character(title)
  match_info <- regexec("AMC_([0-9]+)-[0-9]+", title)
  matched <- regmatches(title, match_info)

  vapply(matched, function(x) {
    if (length(x) >= 2) {
      return(paste0("AMC_", x[2]))
    }

    NA_character_
  }, character(1))
}

classify_tissue <- function(tissue, title) {
  text <- tolower(paste(tissue, title))

  if (grepl("metastatic|metastasized|liver", text)) {
    return("Liver_metastasis")
  }
  if (grepl("normal", text)) {
    return("Normal_colon")
  }
  if (grepl("primary", text)) {
    return("Primary_tumor")
  }

  "Unknown"
}

make_pair_summary <- function(clinical_data) {
  tissue_levels <- c("Primary_tumor", "Normal_colon", "Liver_metastasis")
  pair_table <- table(clinical_data$Patient_ID, clinical_data$tissue_class)
  missing_columns <- setdiff(tissue_levels, colnames(pair_table))
  if (length(missing_columns) > 0) {
    pair_table <- cbind(pair_table, matrix(
      0,
      nrow = nrow(pair_table),
      ncol = length(missing_columns),
      dimnames = list(rownames(pair_table), missing_columns)
    ))
  }
  pair_table <- pair_table[, tissue_levels, drop = FALSE]

  data.frame(
    Patient_ID = rownames(pair_table),
    Primary_tumor = as.integer(pair_table[, "Primary_tumor"]),
    Normal_colon = as.integer(pair_table[, "Normal_colon"]),
    Liver_metastasis = as.integer(pair_table[, "Liver_metastasis"]),
    Complete_Trio = rowSums(pair_table == 1) == length(tissue_levels),
    stringsAsFactors = FALSE
  )
}


# 2. Read sample information ---------------------------------------------------

dir.create(dirname(CLINICAL_EDIT_FILE), recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
stopifnot(inherits(se, "SummarizedExperiment"))

clinical_data <- read.csv(
  CLINICAL_RAW_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot("Title" %in% colnames(clinical_data))
stopifnot("tissue" %in% colnames(clinical_data))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))

missing_samples <- setdiff(colnames(se), clinical_data$Sample_ID)
stopifnot(length(missing_samples) == 0)

clinical_data <- clinical_data[
  match(colnames(se), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
stopifnot(all(clinical_data$Sample_ID == colnames(se)))


# 3. Derive sample annotations -------------------------------------------------

clinical_data$Patient_ID <- extract_patient_id(clinical_data$Title)
stopifnot(!any(is.na(clinical_data$Patient_ID)))

clinical_data$tissue_class <- mapply(
  classify_tissue,
  clinical_data$tissue,
  clinical_data$Title,
  USE.NAMES = FALSE
)

clinical_data$sample_system <- "matched_stage_iv_crc_primary_normal_liver_metastasis"
clinical_data$is_liver_metastasis_tumor_tissue <- clinical_data$tissue_class == "Liver_metastasis"
clinical_data$is_normal_tissue <- clinical_data$tissue_class == "Normal_colon"
clinical_data$is_primary_tumor_tissue <- clinical_data$tissue_class == "Primary_tumor"
clinical_data$is_tumor_tissue <- clinical_data$tissue_class %in% c("Primary_tumor", "Liver_metastasis")

pair_summary <- make_pair_summary(clinical_data)
complete_trio_count <- sum(pair_summary$Complete_Trio)


# 4. Define downstream analysis columns ---------------------------------------

clinical_data$analysis_Liver_metastasis_vs_normal <- ifelse(
  clinical_data$tissue_class == "Liver_metastasis",
  "Liver_metastasis_vs_normal",
  ifelse(clinical_data$tissue_class == "Normal_colon", "Normal_colon", "")
)

clinical_data$analysis_Liver_metastasis_vs_primary <- ifelse(
  clinical_data$tissue_class == "Liver_metastasis",
  "Liver_metastasis_vs_primary",
  ifelse(clinical_data$tissue_class == "Primary_tumor", "Primary_tumor", "")
)

clinical_data$analysis_Primary_tumor_vs_normal <- ifelse(
  clinical_data$tissue_class == "Primary_tumor",
  "Primary_tumor_vs_normal",
  ifelse(clinical_data$tissue_class == "Normal_colon", "Normal_colon", "")
)


# 5. Summaries -----------------------------------------------------------------

analysis_columns <- grep("^analysis_", colnames(clinical_data), value = TRUE)

design_notes <- c(
  Liver_metastasis_vs_normal = "Matched-pair contrast: liver metastasis tumor tissue versus normal-looking surrounding colonic epithelium.",
  Liver_metastasis_vs_primary = "Matched-pair contrast: liver metastasis tumor tissue versus primary colorectal cancer.",
  Primary_tumor_vs_normal = "Matched-pair contrast: primary colorectal cancer versus normal-looking surrounding colonic epithelium."
)

design_summary <- do.call(
  rbind,
  lapply(analysis_columns, function(column_name) {
    values <- trimws(as.character(clinical_data[[column_name]]))
    values[is.na(values)] <- ""
    used_values <- values[values != ""]
    group_counts <- table(used_values)
    analysis_name <- sub("^analysis_", "", column_name)
    experiment_group <- analysis_name
    control_group <- setdiff(names(group_counts), experiment_group)

    used_data <- clinical_data[values != "", , drop = FALSE]
    pair_table <- table(used_data$Patient_ID, used_values)
    complete_pairs <- sum(rowSums(pair_table > 0) == 2)

    data.frame(
      Dataset = DATASET_ID,
      Analysis_Column = column_name,
      Analysis_Name = analysis_name,
      Experiment_Group = experiment_group,
      Control_Group = paste(control_group, collapse = ";"),
      Samples_Used = length(used_values),
      Patients_Used = length(unique(used_data$Patient_ID)),
      Complete_Pairs = complete_pairs,
      Group_Counts = paste(
        paste(names(group_counts), as.integer(group_counts), sep = "="),
        collapse = ";"
      ),
      Ready_For_Two_Group_DE = length(group_counts) == 2 &&
        experiment_group %in% names(group_counts) &&
        length(control_group) == 1,
      Paired_Design_Recommended = complete_pairs > 0,
      Note = unname(design_notes[[analysis_name]]),
      stringsAsFactors = FALSE
    )
  })
)
rownames(design_summary) <- NULL

requested_comparison_audit <- data.frame(
  Dataset = DATASET_ID,
  Requested_Question = c(
    "ATF3 function in liver metastasis tumor tissue",
    "ATF3 difference between liver metastasis tumor tissue and normal tissue"
  ),
  Supported_By_Current_Sample_Info = c(TRUE, TRUE),
  Implemented_As = c(
    "Within-tissue ATF3 high/low DE, ATF3 correlation, and GSEA, including the liver metastasis tissue subset.",
    "Matched-pair limma contrast: Liver_metastasis_vs_normal."
  ),
  Liver_Metastasis_Tumor_Samples = sum(clinical_data$is_liver_metastasis_tumor_tissue),
  Normal_Tissue_Samples = sum(clinical_data$is_normal_tissue),
  Complete_Trio_Patients = complete_trio_count,
  stringsAsFactors = FALSE
)


# 6. Save outputs --------------------------------------------------------------

write.csv(clinical_data, CLINICAL_EDIT_FILE, row.names = FALSE)
write_csv_with_report_previews(
  clinical_data,
  file.path(OUTPUT_ROOT, "GSE50760_clinical_edit.csv")
)
write_csv_with_report_previews(
  pair_summary,
  file.path(OUTPUT_ROOT, "patient_pair_summary.csv")
)
write_csv_with_report_previews(
  design_summary,
  file.path(OUTPUT_ROOT, "analysis_design_summary.csv")
)
write_csv_with_report_previews(
  requested_comparison_audit,
  file.path(OUTPUT_ROOT, "requested_comparison_audit.csv")
)

cat("\nGSE50760 sample design prepared.\n")
cat("Clinical edit file: ", CLINICAL_EDIT_FILE, "\n", sep = "")
cat("Analysis designs: ", length(analysis_columns), "\n", sep = "")
cat("Complete matched trios: ", complete_trio_count, "\n\n", sep = "")
print(design_summary, row.names = FALSE)
