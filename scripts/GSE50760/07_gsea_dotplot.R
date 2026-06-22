# GSE50760 batch GSEA dotplots.
#
# Scans main DE-GSEA, ATF3 correlation GSEA, and ATF3 high/low DE-GSEA tables,
# then saves dotplots by result category.


# 0. Config -------------------------------------------------------------------

DATASET_ID <- "GSE50760"
DATA_TYPE <- "ngs"

PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
GSEA_FUNCTION_FILE <- "scripts/functions/gsea_common_functions.R"
BATCH_VIS_FUNCTION_FILE <- "scripts/functions/gse_batch_visualization_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "GSEA")
SUMMARY_ROOT <- file.path(RESULT_ROOT, "tables", "GSEA_plot_summary")

SIMPLIFY_PATHWAY_PREFIX_IN_PLOT <- TRUE
REPLACE_UNDERSCORE_WITH_SPACE_IN_PLOT <- TRUE

GSEA_DOTPLOT_TOP_N <- 10
GSEA_DOTPLOT_P_COLUMN <- Sys.getenv("GSE50760_GSEA_DOTPLOT_P_COLUMN", unset = "pvalue")
GSEA_DOTPLOT_P_CUTOFF <- as.numeric(Sys.getenv("GSE50760_GSEA_DOTPLOT_P_CUTOFF", unset = "0.05"))
GSEA_DOTPLOT_LABEL_WIDTH <- 45
SHOW_NON_SIGNIFICANT_TOP_TERMS <- FALSE

GSEAVIS_POINT_SIZE_RANGE <- c(4.2, 9.8)
DOTPLOT_BODY_BASE_SIZE <- 4.8
DOTPLOT_LABEL_LINE_HEIGHT <- 0.34
DOTPLOT_TERM_GAP_HEIGHT <- 0.14
DOTPLOT_BODY_MIN_SIZE <- 5.2
DOTPLOT_BODY_MAX_SIZE <- 18.0
DOTPLOT_LABEL_BASE_WIDTH <- 1.8
DOTPLOT_LABEL_WIDTH_PER_CHARACTER <- 0.045
DOTPLOT_LABEL_MIN_WIDTH <- 2.4
DOTPLOT_LABEL_MAX_WIDTH <- 6.2
DOTPLOT_LEGEND_WIDTH <- 1.4
DOTPLOT_VERTICAL_PADDING <- 0.45

OVERWRITE_GSEA_DOTPLOT_OUTPUT <- TRUE

options(width = 200)
options(lifecycle_verbosity = "quiet")


# 1. Packages and helpers ------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

source(REPORT_TABLE_FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)
source(GSEA_FUNCTION_FILE)
source(BATCH_VIS_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()


# 2. Collect and filter inputs -------------------------------------------------

GSEA_CATEGORIES_TO_PLOT <- parse_env_vector("GSE50760_GSEA_DOTPLOT_CATEGORIES", "all")
GSEA_GENESETS_TO_PLOT <- parse_env_vector("GSE50760_GSEA_DOTPLOT_GENESETS", "all")
SHOW_NON_SIGNIFICANT_TOP_TERMS <- parse_env_logical(
  "GSE50760_GSEA_DOTPLOT_SHOW_NONSIG",
  SHOW_NON_SIGNIFICANT_TOP_TERMS
)
CLEAN_GSEA_DOTPLOT_ROOT <- parse_env_logical("GSE50760_GSEA_DOTPLOT_CLEAN", TRUE)

file_info <- collect_batch_gsea_file_info(TABLE_ROOT)
selected_file_info <- filter_gsea_file_info(
  file_info = file_info,
  categories_to_plot = GSEA_CATEGORIES_TO_PLOT,
  genesets_to_plot = GSEA_GENESETS_TO_PLOT
)
stopifnot(nrow(selected_file_info) > 0)

if (CLEAN_GSEA_DOTPLOT_ROOT && dir.exists(PLOT_ROOT)) {
  unlink(PLOT_ROOT, recursive = TRUE)
}
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)


# 3. Plot and save -------------------------------------------------------------

cat("\nRunning GSE50760 GSEA dotplot generation...\n")
cat("Selected GSEA tables: ", nrow(selected_file_info), "\n", sep = "")
cat("P value column: ", GSEA_DOTPLOT_P_COLUMN, "\n", sep = "")
cat("P value cutoff: ", GSEA_DOTPLOT_P_CUTOFF, "\n", sep = "")

summary_table <- run_batch_gsea_dotplots(
  file_info = selected_file_info,
  plot_root = PLOT_ROOT
)

dir.create(SUMMARY_ROOT, recursive = TRUE, showWarnings = FALSE)
summary_csv_file <- write_csv_with_report_previews(
  summary_table,
  file.path(SUMMARY_ROOT, "summary.csv"),
  n_rows = 21
)

cat("\nGSEA dotplot summary:\n")
print(
  summary_table[
    ,
    c(
      "Plot_Category", "Analysis_Name", "GeneSet_Name", "GSEA_Terms",
      "Terms_Plotted", "Positive_NES", "Negative_NES", "PDF_Width", "PDF_Height"
    )
  ],
  row.names = FALSE
)

cat("\nGSEA dotplot summary table: ", summary_csv_file, "\n", sep = "")
cat("GSEA dotplot directory: ", PLOT_ROOT, "\n", sep = "")
cat("\nGSE50760 GSEA dotplot generation finished.\n")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
