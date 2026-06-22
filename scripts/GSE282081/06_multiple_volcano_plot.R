# GSE282081多组差异表达火山图
#
# 按配置将多套DEG结果合并到同一张分面火山图中，绘图风格沿用GSE114012。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE282081"
DATA_TYPE <- "ngs"

DE_SCRIPT_FILE <- "scripts/GSE282081/02_limma_differential_expression.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
BATCH_VIS_FUNCTION_FILE <- "scripts/functions/gse_batch_visualization_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "multiple_volcano")
SUMMARY_ROOT <- file.path(RESULT_ROOT, "tables", "multiple_volcano_plot_summary")

SYNC_THRESHOLDS_FROM_DE_SCRIPT <- TRUE
P_VALUE_COLUMN <- "P.Value"
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5

MULTIPLE_VOLCANO_SCHEMES <- list(
  ALL_MAIN_DE = c(
    "Main_DE::Hepatocyte_coculture",
    "Main_DE::3D_organoid",
    "Main_DE::3D_hepatic_coculture",
    "Main_DE::2D_J2_coculture",
    "Main_DE::2D_J2_hepatocyte",
    "Main_DE::2D_J2_hepatocyte_vs_2D_none"
  ),
  COCULTURE_MODELS = c(
    "Main_DE::Hepatocyte_coculture",
    "Main_DE::2D_J2_coculture",
    "Main_DE::2D_J2_hepatocyte",
    "Main_DE::3D_hepatic_coculture"
  ),
  MODEL_STRUCTURE_CONTEXT = c(
    "Main_DE::3D_organoid",
    "Main_DE::3D_hepatic_coculture",
    "Main_DE::2D_J2_hepatocyte_vs_2D_none"
  )
)

CLEAN_MULTIPLE_VOLCANO_ROOT <- TRUE
OVERWRITE_SCHEME_OUTPUT <- TRUE
CUSTOM_LABEL_GENES <- character(0)

TOP_GENE_LABEL_NUDGE_Y <- 0.80

GROUP_LABEL_COLORS <- c(
  "Hepatocyte_coculture" = "#c70008",
  "3D_organoid" = "#0052a2",
  "3D_hepatic_coculture" = "#006b3f",
  "2D_J2_coculture" = "#eb7400",
  "2D_J2_hepatocyte" = "#0090c1",
  "2D_J2_hepatocyte_vs_2D_none" = "#805190"
)
GROUP_LABEL_ALPHA <- 0.24
GROUP_LABEL_WRAP_WIDTH <- 10
GROUP_LABEL_FONT_SIZE <- 5.3
GROUP_LABEL_LINE_HEIGHT <- 1.18
GROUP_LABEL_TEXT_DARKEN <- 0.78
GROUP_LABEL_BORDER_DARKEN <- 0.88
GROUP_LABEL_BOX_LOGFC_GAP <- 0.10
GROUP_LABEL_BOX_MIN_FRACTION <- 0.80
GROUP_LABEL_FONT_MIN_SIZE <- 3.9
GROUP_LABEL_FONT_LINE_SHRINK <- 0.55
GROUP_LABEL_LINE_HEIGHT_MIN <- 0.95
GROUP_LABEL_LINE_HEIGHT_SHRINK <- 0.08
GROUP_LABEL_BOX_X_MARGIN_FRACTION <- 0.05
GROUP_LABEL_BORDER_WIDTH <- 0.9

PANEL_SPACING_X_MM <- 4.6
LEGEND_TOP_MARGIN_PT <- 68
X_AXIS_PADDING_FRACTION <- 0.12
X_AXIS_PADDING_MIN <- 0.22
X_AXIS_LEFT_MIN_FRACTION <- 0.78
Y_AXIS_PADDING_FRACTION <- 0.07
LABEL_BOX_HEIGHT_PDF_SCALE <- 1.8

BASE_MULTIPLE_PDF_HEIGHT <- 6.0
GROUP_WIDTH_INCH <- 1.78
MULTIPLE_LEGEND_WIDTH_INCH <- 1.18
MIN_MULTIPLE_PDF_WIDTH <- 7.2
MAX_MULTIPLE_PDF_WIDTH <- 20.0
MAX_MULTIPLE_PDF_HEIGHT <- 8.2
LEGEND_TEXT_SIZE <- 13.5
LEGEND_POINT_SIZE_SCALE <- 1.45
LEGEND_KEY_SIZE_MM <- 6.6

options(width = 200)


# 1. 加载包和函数 --------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

source(REPORT_TABLE_FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)
source(BATCH_VIS_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()
USE_GG_REPEL <- requireNamespace("ggrepel", quietly = TRUE)


# 2. 同步阈值并检查输入 --------------------------------------------------------

SCHEMES_TO_RUN <- parse_env_vector(
  "GSE282081_MULTIPLE_VOLCANO_SCHEMES",
  names(MULTIPLE_VOLCANO_SCHEMES)
)

missing_schemes <- setdiff(SCHEMES_TO_RUN, names(MULTIPLE_VOLCANO_SCHEMES))
if (length(missing_schemes) > 0) {
  stop(
    "Undefined multiple volcano schemes: ",
    paste(missing_schemes, collapse = ", ")
  )
}

threshold_config <- sync_de_thresholds_from_script(
  script_file = DE_SCRIPT_FILE,
  p_value_column = P_VALUE_COLUMN,
  p_value_cutoff = P_VALUE_CUTOFF,
  logfc_cutoff = LOGFC_CUTOFF,
  sync = SYNC_THRESHOLDS_FROM_DE_SCRIPT
)
P_VALUE_COLUMN <- threshold_config$p_value_column
P_VALUE_CUTOFF <- threshold_config$p_value_cutoff
LOGFC_CUTOFF <- threshold_config$logfc_cutoff

file_info <- collect_batch_deg_file_info(TABLE_ROOT)
stopifnot(nrow(file_info) > 0)


# 3. 绘图并保存 ----------------------------------------------------------------

cat("\nRunning GSE282081 multiple volcano plot generation...\n")
cat("Schemes: ", paste(SCHEMES_TO_RUN, collapse = ", "), "\n", sep = "")
cat("P value column: ", P_VALUE_COLUMN, "\n", sep = "")
cat("P value cutoff: ", P_VALUE_CUTOFF, "\n", sep = "")
cat("logFC cutoff: ", LOGFC_CUTOFF, "\n", sep = "")

summary_table <- run_batch_multiple_volcano_plots(
  file_info = file_info,
  schemes = MULTIPLE_VOLCANO_SCHEMES,
  schemes_to_run = SCHEMES_TO_RUN,
  plot_root = PLOT_ROOT
)

dir.create(SUMMARY_ROOT, recursive = TRUE, showWarnings = FALSE)
summary_csv_file <- write_csv_with_report_previews(
  summary_table,
  file.path(SUMMARY_ROOT, "summary.csv"),
  n_rows = 21
)

cat("\nMultiple volcano plot summary:\n")
print(
  summary_table[
    ,
    c(
      "Plot_Name", "Analysis_Name", "Plot_Category", "Genes_Plotted",
      "Up", "Down", "NS", "Status", "PDF_Width", "PDF_Height"
    )
  ],
  row.names = FALSE
)

cat("\nMultiple volcano summary table: ", summary_csv_file, "\n", sep = "")
cat("Multiple volcano plot directory: ", PLOT_ROOT, "\n", sep = "")
cat("\nGSE282081 multiple volcano plot generation finished.\n")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
