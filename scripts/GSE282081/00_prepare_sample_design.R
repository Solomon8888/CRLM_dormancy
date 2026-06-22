# GSE282081样本设计表整理
#
# 本脚本根据当前data_prepare中的临时样本信息，生成带analysis_列的
# GSE282081_clinical_edit.csv，供后续limma差异分析直接读取。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE282081"
DATA_TYPE <- "ngs"

SE_RDS_FILE <- "data/ngs/GSE282081/data_prepare/GSE282081_se_raw.rds"
CLINICAL_RAW_FILE <- "data/ngs/GSE282081/data_prepare/GSE282081_clinical_raw.csv"
CLINICAL_EDIT_FILE <- "data/ngs/GSE282081/data_prepare/GSE282081_clinical_edit.csv"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"

OUTPUT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID, "tables", "sample_design")

options(width = 200)


# 1. 加载包和函数 --------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

source(REPORT_TABLE_FUNCTION_FILE)


# 2. 读取样本信息 --------------------------------------------------------------

dir.create(dirname(CLINICAL_EDIT_FILE), recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
stopifnot(inherits(se, "SummarizedExperiment"))

if (file.exists(CLINICAL_RAW_FILE)) {
  clinical_data <- read.csv(
    CLINICAL_RAW_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else {
  clinical_data <- as.data.frame(colData(se), stringsAsFactors = FALSE)
}

stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot("Title" %in% colnames(clinical_data))
stopifnot("treatment" %in% colnames(clinical_data))
stopifnot(!any(duplicated(clinical_data$Sample_ID)))

missing_samples <- setdiff(colnames(se), clinical_data$Sample_ID)
stopifnot(length(missing_samples) == 0)

clinical_data <- clinical_data[
  match(colnames(se), clinical_data$Sample_ID),
  ,
  drop = FALSE
]
stopifnot(all(clinical_data$Sample_ID == colnames(se)))


# 3. 从Title/treatment中派生分析分组 ------------------------------------------

title_lower <- tolower(trimws(as.character(clinical_data$Title)))
treatment_lower <- tolower(trimws(as.character(clinical_data$treatment)))

is_3d <- grepl("^3d", title_lower)
is_2d <- grepl("^2d", title_lower)
has_j2 <- grepl("j2|3t3", title_lower) | grepl("j2|3t3", treatment_lower)
has_hepatocyte <- grepl("hepatocyte", title_lower) |
  grepl("hepatocyte", treatment_lower)

replicate_id <- sub(".*rep([0-9]+).*", "\\1", title_lower)
replicate_id[!grepl("rep[0-9]+", title_lower)] <- NA

clinical_data$culture_model <- ifelse(
  is_3d,
  "3D_organoid",
  ifelse(is_2d, "2D_culture", "Unknown")
)
clinical_data$J2_coculture <- ifelse(has_j2, "J2_coculture", "No_J2")
clinical_data$hepatocyte_coculture <- ifelse(
  has_hepatocyte,
  "Hepatocyte_coculture",
  "No_hepatocyte_coculture"
)
clinical_data$replicate <- replicate_id

clinical_data$treatment_class <- ifelse(
  has_j2 & has_hepatocyte,
  "J2_hepatocyte_coculture",
  ifelse(
    has_j2,
    "J2_coculture",
    ifelse(has_hepatocyte, "Hepatocyte_coculture", "None")
  )
)

clinical_data$liver_niche_model <- ifelse(
  clinical_data$culture_model == "3D_organoid" & has_hepatocyte,
  "3D_organoid_J2_hepatocyte",
  ifelse(
    clinical_data$culture_model == "3D_organoid",
    "3D_organoid_none",
    ifelse(
      clinical_data$culture_model == "2D_culture" & has_hepatocyte,
      "2D_J2_hepatocyte",
      ifelse(
        clinical_data$culture_model == "2D_culture" & has_j2,
        "2D_J2",
        "2D_none"
      )
    )
  )
)

# 当前本地样本表是体外SW480模型，不包含真实肝转移肿瘤组织或正常组织。
clinical_data$sample_system <- "in_vitro_SW480_model"
clinical_data$is_liver_metastasis_tumor_tissue <- FALSE
clinical_data$is_normal_tissue <- FALSE


# 4. 定义limma使用的analysis_列 ------------------------------------------------

clinical_data$analysis_Hepatocyte_coculture <- ifelse(
  has_hepatocyte,
  "Hepatocyte_coculture",
  "No_hepatocyte_coculture"
)

clinical_data$analysis_3D_organoid <- ifelse(
  clinical_data$culture_model == "3D_organoid",
  "3D_organoid",
  "2D_culture"
)

clinical_data$analysis_3D_hepatic_coculture <- ifelse(
  clinical_data$culture_model == "3D_organoid",
  ifelse(has_hepatocyte, "3D_hepatic_coculture", "3D_organoid_none"),
  ""
)

clinical_data$analysis_2D_J2_coculture <- ifelse(
  clinical_data$culture_model == "2D_culture" & !has_hepatocyte,
  ifelse(has_j2, "2D_J2_coculture", "2D_none"),
  ""
)

clinical_data$analysis_2D_J2_hepatocyte <- ifelse(
  clinical_data$culture_model == "2D_culture" & has_j2,
  ifelse(has_hepatocyte, "2D_J2_hepatocyte", "2D_J2_coculture"),
  ""
)

clinical_data$analysis_2D_J2_hepatocyte_vs_2D_none <- ifelse(
  clinical_data$culture_model == "2D_culture" & (!has_j2 | has_hepatocyte),
  ifelse(has_hepatocyte, "2D_J2_hepatocyte_vs_2D_none", "2D_none"),
  ""
)


# 5. 输出设计说明 --------------------------------------------------------------

analysis_columns <- grep("^analysis_", colnames(clinical_data), value = TRUE)

design_notes <- c(
  Hepatocyte_coculture = "All samples; approximates liver-cell niche exposure, but culture model and J2 status are partly confounded.",
  `3D_organoid` = "All samples; compares 3D organoid model against 2D culture.",
  `3D_hepatic_coculture` = "3D samples only; compares 3D organoids with J2/hepatocyte coculture against 3D organoids alone.",
  `2D_J2_coculture` = "2D samples without hepatocytes; compares J2 coculture against 2D culture alone.",
  `2D_J2_hepatocyte` = "2D J2 samples only; compares added hepatocytes against J2 coculture alone.",
  `2D_J2_hepatocyte_vs_2D_none` = "2D samples only; compares J2/hepatocyte coculture against 2D culture alone."
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

    data.frame(
      Dataset = DATASET_ID,
      Analysis_Column = column_name,
      Analysis_Name = analysis_name,
      Experiment_Group = experiment_group,
      Control_Group = paste(control_group, collapse = ";"),
      Samples_Used = length(used_values),
      Group_Counts = paste(
        paste(names(group_counts), as.integer(group_counts), sep = "="),
        collapse = ";"
      ),
      Ready_For_Two_Group_DE = length(group_counts) == 2 &&
        experiment_group %in% names(group_counts),
      Note = unname(design_notes[[analysis_name]]),
      stringsAsFactors = FALSE
    )
  })
)
rownames(design_summary) <- NULL

requested_comparison_audit <- data.frame(
  Dataset = DATASET_ID,
  Requested_Question = c(
    "ATF3 function in liver-metastasis-like tumor model",
    "ATF3 difference between liver metastasis tumor tissue and normal tissue"
  ),
  Supported_By_Current_Sample_Info = c(TRUE, FALSE),
  Implemented_As = c(
    "ATF3 expression, DE and GSEA across hepatocyte coculture and 3D organoid model comparisons.",
    "Not implemented as a formal DE contrast because current samples are in vitro SW480 cultures and contain no normal tissue samples."
  ),
  Tumor_Tissue_Samples = c(NA, sum(clinical_data$is_liver_metastasis_tumor_tissue)),
  Normal_Tissue_Samples = c(NA, sum(clinical_data$is_normal_tissue)),
  stringsAsFactors = FALSE
)


# 6. 保存结果 ------------------------------------------------------------------

write.csv(clinical_data, CLINICAL_EDIT_FILE, row.names = FALSE)
write_csv_with_report_previews(
  clinical_data,
  file.path(OUTPUT_ROOT, "GSE282081_clinical_edit.csv")
)
write_csv_with_report_previews(
  design_summary,
  file.path(OUTPUT_ROOT, "analysis_design_summary.csv")
)
write_csv_with_report_previews(
  requested_comparison_audit,
  file.path(OUTPUT_ROOT, "requested_comparison_audit.csv")
)

cat("\nGSE282081 sample design prepared.\n")
cat("Clinical edit file: ", CLINICAL_EDIT_FILE, "\n", sep = "")
cat("Analysis designs: ", length(analysis_columns), "\n\n", sep = "")
print(design_summary, row.names = FALSE)
