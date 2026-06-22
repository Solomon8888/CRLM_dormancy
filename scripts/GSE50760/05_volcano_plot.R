# GSE50760 volcano plots.
#
# Reads main tissue DE and within-tissue ATF3 high/low DE all_genes.csv files,
# then writes PDF/PNG volcano plots with the same visual style as GSE114012.


# 0. Config -------------------------------------------------------------------

DATASET_ID <- "GSE50760"
DATA_TYPE <- "ngs"

DE_SCRIPT_FILE <- "scripts/GSE50760/02_limma_differential_expression.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
BATCH_VIS_FUNCTION_FILE <- "scripts/functions/gse_batch_visualization_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "volcano")
SUMMARY_ROOT <- file.path(RESULT_ROOT, "tables", "volcano_plot_summary")

SYNC_THRESHOLDS_FROM_DE_SCRIPT <- TRUE
P_VALUE_COLUMN <- "P.Value"
P_VALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5

CLEAN_VOLCANO_OUTPUT_DIR <- TRUE
NOT_SIGNIFICANT_COLOR <- "#B8B8B8"
CUSTOM_LABEL_GENES <- character(0)

THRESHOLD_LINE_COLOR <- "#333333"
THRESHOLD_LINE_WIDTH <- 0.45
THRESHOLD_LINE_TYPE <- "dashed"

BASE_PDF_HEIGHT <- 6.2
MAX_EXTRA_PDF_HEIGHT <- 0.8
LEGEND_WIDTH_INCH <- 0.95
RIGHT_LEGEND_GAP_INCH <- 0.18
MAX_PDF_WIDTH_HEIGHT_RATIO <- 1.28
PANEL_HEIGHT_WIDTH_RATIO <- 1.0

options(width = 200)


# 1. Packages and helpers ------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

source(REPORT_TABLE_FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)
source(BATCH_VIS_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()
USE_GG_REPEL <- requireNamespace("ggrepel", quietly = TRUE)


# 2. Sync thresholds and collect inputs ---------------------------------------

ANALYSES_TO_PLOT <- parse_env_vector("GSE50760_VOLCANO_ANALYSES", "all")
if (length(ANALYSES_TO_PLOT) == 1 && tolower(ANALYSES_TO_PLOT) == "all") {
  ANALYSES_TO_PLOT <- "all"
}
CATEGORIES_TO_PLOT <- parse_env_vector("GSE50760_VOLCANO_CATEGORIES", "all")

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
selected_file_info <- filter_deg_file_info(file_info, ANALYSES_TO_PLOT)
if (!identical(CATEGORIES_TO_PLOT, "all") &&
    !(length(CATEGORIES_TO_PLOT) == 1 && tolower(CATEGORIES_TO_PLOT) == "all")) {
  selected_file_info <- selected_file_info[
    selected_file_info$Plot_Category %in% CATEGORIES_TO_PLOT,
    ,
    drop = FALSE
  ]
}
stopifnot(nrow(selected_file_info) > 0)


# 3. Plot and save -------------------------------------------------------------

cat("\nRunning GSE50760 volcano plot generation...\n")
cat("Selected DEG tables: ", nrow(selected_file_info), "\n", sep = "")
cat("P value column: ", P_VALUE_COLUMN, "\n", sep = "")
cat("P value cutoff: ", P_VALUE_CUTOFF, "\n", sep = "")
cat("logFC cutoff: ", LOGFC_CUTOFF, "\n", sep = "")

summary_table <- run_batch_volcano_plots(
  file_info = selected_file_info,
  plot_root = PLOT_ROOT
)

dir.create(SUMMARY_ROOT, recursive = TRUE, showWarnings = FALSE)
summary_csv_file <- write_csv_with_report_previews(
  summary_table,
  file.path(SUMMARY_ROOT, "summary.csv"),
  n_rows = 21
)

cat("\nVolcano plot summary:\n")
print(
  summary_table[
    ,
    c(
      "Plot_Category", "Analysis_Name", "Genes_Plotted", "Up", "Down",
      "Not_Significant", "PDF_Width", "PDF_Height"
    )
  ],
  row.names = FALSE
)

cat("\nVolcano summary table: ", summary_csv_file, "\n", sep = "")
cat("Volcano plot directory: ", PLOT_ROOT, "\n", sep = "")
cat("\nGSE50760 volcano plot generation finished.\n")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
