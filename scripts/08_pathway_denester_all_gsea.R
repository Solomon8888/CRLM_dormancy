# Project-wide PathwayDenester-style GSEA summary.
#
# This script scans GSEA result tables under results/ngs/<dataset>/tables,
# applies a PathwayDenester-style overlap dependency filter to each
# analysis x geneset result, and writes denested result tables plus summary
# plots following the project's PDF/PNG output convention.


# 0. Config -------------------------------------------------------------------

DATA_TYPE <- "ngs"

PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
GSEA_FUNCTION_FILE <- "scripts/functions/gsea_common_functions.R"
BATCH_VIS_FUNCTION_FILE <- "scripts/functions/gse_batch_visualization_functions.R"
DENESTER_FUNCTION_FILE <- "scripts/functions/pathway_denester_functions.R"

RESULT_BASE <- file.path("results", DATA_TYPE)
PROJECT_SUMMARY_ROOT <- file.path(RESULT_BASE, "GSEA_denester_project_summary")

SPECIES <- "human"
GENE_ID_TYPE <- "SYMBOL"
MSIGDB_REFERENCE_DIR <- file.path("data", "reference", "msigdb")
MSIGDB_REFERENCE_MAX_AGE_DAYS <- 3650

GSEA_DENESTER_P_COLUMN <- Sys.getenv("GSEA_DENESTER_P_COLUMN", unset = "pvalue")
GSEA_DENESTER_GSEA_P_CUTOFF <- as.numeric(Sys.getenv("GSEA_DENESTER_GSEA_P_CUTOFF", unset = "0.05"))
GSEA_DENESTER_P_CUTOFF <- as.numeric(Sys.getenv("GSEA_DENESTER_P_CUTOFF", unset = "0.05"))
GSEA_DENESTER_TO_TEST_THRESHOLD <- as.numeric(Sys.getenv("GSEA_DENESTER_TO_TEST_THRESHOLD", unset = "0"))
GSEA_DENESTER_MAX_UNEXPECTED_CORE_FRACTION <- as.numeric(
  Sys.getenv("GSEA_DENESTER_MAX_UNEXPECTED_CORE_FRACTION", unset = "0.10")
)

DENESTER_DOTPLOT_TOP_N <- as.integer(Sys.getenv("GSEA_DENESTER_DOTPLOT_TOP_N", unset = "10"))
DENESTER_OVERLAP_HEATMAP_TOP_N <- as.integer(Sys.getenv("GSEA_DENESTER_OVERLAP_TOP_N", unset = "30"))
SIMPLIFY_PATHWAY_PREFIX_IN_PLOT <- TRUE
REPLACE_UNDERSCORE_WITH_SPACE_IN_PLOT <- TRUE

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

options(width = 200)
options(lifecycle_verbosity = "quiet")


# 1. Packages and helpers ------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(qs2)
})

source(REPORT_TABLE_FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)
source(GSEA_FUNCTION_FILE)
source(BATCH_VIS_FUNCTION_FILE)
source(DENESTER_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

DATASETS_TO_RUN <- parse_env_vector("GSEA_DENESTER_DATASETS", "all")
if (length(DATASETS_TO_RUN) == 1 && tolower(DATASETS_TO_RUN) == "all") {
  DATASETS_TO_RUN <- basename(list.dirs(RESULT_BASE, recursive = FALSE, full.names = TRUE))
}
DATASETS_TO_RUN <- DATASETS_TO_RUN[dir.exists(file.path(RESULT_BASE, DATASETS_TO_RUN, "tables"))]
stopifnot(length(DATASETS_TO_RUN) > 0)

GSEA_DENESTER_CATEGORIES <- parse_env_vector("GSEA_DENESTER_CATEGORIES", "all")
GSEA_DENESTER_GENESETS <- parse_env_vector("GSEA_DENESTER_GENESETS", "all")
DRAW_OVERLAP_HEATMAP <- parse_env_logical("GSEA_DENESTER_DRAW_HEATMAP", TRUE)
CLEAN_DENESTER_OUTPUT <- parse_env_logical("GSEA_DENESTER_CLEAN", TRUE)


# 2. Load MSigDB references ----------------------------------------------------

cat("\nLoading MSigDB gene sets for denesting...\n")
MSIGDB_GENESET_CATALOG <- build_msigdb_geneset_catalog()
GSEA_GENESET_CONFIG <- select_msigdb_genesets(
  catalog = MSIGDB_GENESET_CATALOG,
  genesets_to_run = "all"
)

gene_set_cache <- lapply(names(GSEA_GENESET_CONFIG), function(geneset_name) {
  config <- GSEA_GENESET_CONFIG[[geneset_name]]
  terms <- load_msigdb_terms(geneset_name, config)
  make_gene_set_lookup(terms$term2gene)
})
names(gene_set_cache) <- names(GSEA_GENESET_CONFIG)


# 3. Collect GSEA result tables ------------------------------------------------

collect_dataset_gsea_info <- function(dataset_id) {
  result_root <- file.path(RESULT_BASE, dataset_id)
  table_root <- file.path(result_root, "tables")
  gsea_info <- collect_batch_gsea_file_info(table_root)

  if (nrow(gsea_info) == 0) {
    return(data.frame())
  }

  gsea_info <- filter_gsea_file_info(
    file_info = gsea_info,
    categories_to_plot = GSEA_DENESTER_CATEGORIES,
    genesets_to_plot = GSEA_DENESTER_GENESETS
  )

  if (nrow(gsea_info) == 0) {
    return(data.frame())
  }

  gsea_info$Dataset <- dataset_id
  gsea_info$Result_Root <- result_root
  gsea_info
}

gsea_info_list <- lapply(DATASETS_TO_RUN, collect_dataset_gsea_info)
gsea_info <- do.call(rbind, gsea_info_list)
rownames(gsea_info) <- NULL
stopifnot(nrow(gsea_info) > 0)

missing_reference <- setdiff(unique(gsea_info$GeneSet_Name), names(gene_set_cache))
if (length(missing_reference) > 0) {
  stop(
    "No MSigDB reference was loaded for geneset directories: ",
    paste(missing_reference, collapse = ", ")
  )
}

cat("\nPathway denester runtime configuration:\n")
cat("Datasets: ", paste(DATASETS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("GSEA result tables: ", nrow(gsea_info), "\n", sep = "")
cat("GSEA p column: ", GSEA_DENESTER_P_COLUMN, "\n", sep = "")
cat("GSEA p cutoff: ", GSEA_DENESTER_GSEA_P_CUTOFF, "\n", sep = "")
cat("Denester p cutoff: ", GSEA_DENESTER_P_CUTOFF, "\n", sep = "")
cat("To-test threshold: ", GSEA_DENESTER_TO_TEST_THRESHOLD, "\n", sep = "")
cat("Draw overlap heatmaps: ", DRAW_OVERLAP_HEATMAP, "\n", sep = "")

if (CLEAN_DENESTER_OUTPUT) {
  for (dataset_id in DATASETS_TO_RUN) {
    unlink(
      file.path(RESULT_BASE, dataset_id, "tables", "GSEA_denester"),
      recursive = TRUE,
      force = TRUE
    )
    unlink(
      file.path(RESULT_BASE, dataset_id, "tables", "GSEA_denester_summary"),
      recursive = TRUE,
      force = TRUE
    )
    unlink(
      file.path(RESULT_BASE, dataset_id, "plots", "GSEA_denester"),
      recursive = TRUE,
      force = TRUE
    )
  }
}


# 4. Run denesting -------------------------------------------------------------

cat("\nRunning PathwayDenester-style GSEA cleanup...\n")
summary_records <- vector("list", nrow(gsea_info))

for (i in seq_len(nrow(gsea_info))) {
  if (i == 1 || i %% 20 == 0 || i == nrow(gsea_info)) {
    cat(
      "Progress: ",
      i,
      "/",
      nrow(gsea_info),
      " ",
      gsea_info$Dataset[i],
      " / ",
      gsea_info$Plot_Category[i],
      " / ",
      gsea_info$Analysis_Name[i],
      " / ",
      gsea_info$GeneSet_Name[i],
      "\n",
      sep = ""
    )
  }

  summary_records[[i]] <- run_denester_for_gsea_file(
    gsea_info_row = gsea_info[i, , drop = FALSE],
    gene_sets = gene_set_cache[[gsea_info$GeneSet_Name[i]]],
    result_root = gsea_info$Result_Root[i],
    gsea_p_column = GSEA_DENESTER_P_COLUMN,
    gsea_p_cutoff = GSEA_DENESTER_GSEA_P_CUTOFF,
    denester_p_cutoff = GSEA_DENESTER_P_CUTOFF,
    to_test_threshold = GSEA_DENESTER_TO_TEST_THRESHOLD,
    max_unexpected_core_fraction = GSEA_DENESTER_MAX_UNEXPECTED_CORE_FRACTION,
    dotplot_top_n = DENESTER_DOTPLOT_TOP_N,
    overlap_heatmap_top_n = DENESTER_OVERLAP_HEATMAP_TOP_N,
    draw_overlap_heatmap = DRAW_OVERLAP_HEATMAP
  )
}

summary_table <- do.call(rbind, summary_records)
rownames(summary_table) <- NULL


# 5. Save dataset-level summaries ---------------------------------------------

for (dataset_id in unique(summary_table$Dataset)) {
  dataset_summary <- summary_table[summary_table$Dataset == dataset_id, , drop = FALSE]
  result_root <- file.path(RESULT_BASE, dataset_id)
  summary_table_dir <- file.path(result_root, "tables", "GSEA_denester_summary")
  summary_plot_dir <- file.path(result_root, "plots", "GSEA_denester", "summary")
  dir.create(summary_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summary_plot_dir, recursive = TRUE, showWarnings = FALSE)

  dataset_summary_file <- write_csv_with_report_previews(
    dataset_summary,
    file.path(summary_table_dir, "summary.csv"),
    n_rows = 21,
    na = "NA"
  )

  overview_plot <- make_denester_summary_barplot(dataset_summary)
  save_ggplot_pdf_png(
    plot = overview_plot,
    pdf_file = file.path(summary_plot_dir, "denester_excluded_percent_by_category.pdf"),
    width = max(7.2, length(unique(dataset_summary$Plot_Category)) * 2.0 + 3.2),
    height = 5.8
  )

  category_summary <- aggregate(
    cbind(Significant_Terms, Kept_Terms, Excluded_Terms) ~ Dataset + Plot_Category,
    data = dataset_summary,
    FUN = sum
  )
  category_summary$Excluded_Percent <- ifelse(
    category_summary$Significant_Terms > 0,
    category_summary$Excluded_Terms / category_summary$Significant_Terms * 100,
    0
  )
  write_csv_with_report_previews(
    category_summary,
    file.path(summary_table_dir, "summary_by_category.csv"),
    n_rows = 21,
    na = "NA"
  )

  cat("\nDataset denester summary saved: ", dataset_summary_file, "\n", sep = "")
}


# 6. Save project-level summary ------------------------------------------------

dir.create(PROJECT_SUMMARY_ROOT, recursive = TRUE, showWarnings = FALSE)
project_summary_file <- write_csv_with_report_previews(
  summary_table,
  file.path(PROJECT_SUMMARY_ROOT, "summary.csv"),
  n_rows = 21,
  na = "NA"
)

project_plot <- make_denester_summary_barplot(summary_table)
save_ggplot_pdf_png(
  plot = project_plot,
  pdf_file = file.path(PROJECT_SUMMARY_ROOT, "denester_excluded_percent_by_dataset_category.pdf"),
  width = max(8.2, length(unique(paste(summary_table$Dataset, summary_table$Plot_Category))) * 1.25 + 3.2),
  height = 6.2
)

project_category_summary <- aggregate(
  cbind(Significant_Terms, Kept_Terms, Excluded_Terms) ~ Dataset + Plot_Category,
  data = summary_table,
  FUN = sum
)
project_category_summary$Excluded_Percent <- ifelse(
  project_category_summary$Significant_Terms > 0,
  project_category_summary$Excluded_Terms / project_category_summary$Significant_Terms * 100,
  0
)
write_csv_with_report_previews(
  project_category_summary,
  file.path(PROJECT_SUMMARY_ROOT, "summary_by_dataset_category.csv"),
  n_rows = 21,
  na = "NA"
)

cat("\nProject denester summary:\n")
print(project_category_summary, row.names = FALSE)

cat("\nProject denester summary table: ", project_summary_file, "\n", sep = "")
cat("\nPathwayDenester-style GSEA cleanup finished.\n")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
