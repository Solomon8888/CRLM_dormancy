# GSE50760 ATF3 expression overview.
#
# Summarizes ATF3 expression across matched normal colon, primary tumor, and
# liver metastasis samples. Two-group tests use paired tests whenever complete
# patient pairs are available.


# 0. Config -------------------------------------------------------------------

DATASET_ID <- "GSE50760"
DATA_TYPE <- "ngs"
TARGET_GENE <- "ATF3"

SE_RDS_FILE <- "data/ngs/GSE50760/data_prepare/GSE50760_se_raw.rds"
CLINICAL_EDIT_FILE <- "data/ngs/GSE50760/data_prepare/GSE50760_clinical_edit.csv"
SAMPLE_DESIGN_SCRIPT <- "scripts/GSE50760/00_prepare_sample_design.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables", "ATF3_expression")
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "ATF3_expression")

options(width = 200)


# 1. Packages and helpers ------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggplot2)
})

source(REPORT_TABLE_FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)

get_target_index <- function(row_annotation, target_gene) {
  target_index <- which(toupper(trimws(as.character(row_annotation$Symbol))) == toupper(target_gene))
  if (length(target_index) == 0) {
    stop("Target gene was not found in rowData Symbol: ", target_gene)
  }

  target_index[1]
}

summarize_by_group <- function(dat, group_column) {
  groups <- unique(as.character(dat[[group_column]]))
  groups <- groups[!is.na(groups) & groups != ""]

  do.call(
    rbind,
    lapply(groups, function(group_value) {
      x <- dat$log2_tpm_plus_1[dat[[group_column]] == group_value]
      data.frame(
        Group_Variable = group_column,
        Group = group_value,
        N = length(x),
        Mean_Log2_TPM_Plus_1 = mean(x, na.rm = TRUE),
        Median_Log2_TPM_Plus_1 = median(x, na.rm = TRUE),
        SD_Log2_TPM_Plus_1 = sd(x, na.rm = TRUE),
        Mean_TPM = mean(dat$tpm[dat[[group_column]] == group_value], na.rm = TRUE),
        Median_TPM = median(dat$tpm[dat[[group_column]] == group_value], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
}

run_group_test <- function(analysis_column) {
  values <- trimws(as.character(clinical_data[[analysis_column]]))
  values[is.na(values)] <- ""
  keep <- values != ""
  values <- values[keep]
  patient_id <- clinical_data$Patient_ID[keep]
  y <- expression_table$log2_tpm_plus_1[keep]

  group_names <- unique(values)
  analysis_name <- sub("^analysis_", "", analysis_column)
  experiment_group <- analysis_name
  control_group <- setdiff(group_names, experiment_group)

  if (length(group_names) != 2 || !experiment_group %in% group_names || length(control_group) != 1) {
    return(data.frame(
      Analysis_Name = analysis_name,
      Contrast = NA_character_,
      N_Experiment = NA_integer_,
      N_Control = NA_integer_,
      Complete_Pairs = NA_integer_,
      Mean_Experiment = NA_real_,
      Mean_Control = NA_real_,
      Mean_Difference = NA_real_,
      Paired_T_Test_P = NA_real_,
      Paired_Wilcox_P = NA_real_,
      Unpaired_T_Test_P = NA_real_,
      Wilcox_P = NA_real_,
      Test_Mode = "Skipped",
      Status = "Skipped: expected exactly two groups with experiment group matching analysis name.",
      stringsAsFactors = FALSE
    ))
  }

  control_group <- control_group[1]
  x_exp <- y[values == experiment_group]
  x_ctrl <- y[values == control_group]

  pair_table <- data.frame(
    Patient_ID = patient_id,
    Group = values,
    Expression = y,
    stringsAsFactors = FALSE
  )
  pair_wide <- reshape(
    pair_table,
    idvar = "Patient_ID",
    timevar = "Group",
    direction = "wide"
  )
  exp_column <- paste0("Expression.", experiment_group)
  ctrl_column <- paste0("Expression.", control_group)
  complete_pair <- !is.na(pair_wide[[exp_column]]) & !is.na(pair_wide[[ctrl_column]])
  pair_wide <- pair_wide[complete_pair, , drop = FALSE]

  paired_t_p <- tryCatch(
    t.test(pair_wide[[exp_column]], pair_wide[[ctrl_column]], paired = TRUE)$p.value,
    error = function(e) NA_real_
  )
  paired_w_p <- tryCatch(
    wilcox.test(pair_wide[[exp_column]], pair_wide[[ctrl_column]], paired = TRUE, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  unpaired_t_p <- tryCatch(t.test(x_exp, x_ctrl)$p.value, error = function(e) NA_real_)
  wilcox_p <- tryCatch(
    wilcox.test(x_exp, x_ctrl, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )

  data.frame(
    Analysis_Name = analysis_name,
    Contrast = paste0(experiment_group, "_vs_", control_group),
    N_Experiment = length(x_exp),
    N_Control = length(x_ctrl),
    Complete_Pairs = nrow(pair_wide),
    Mean_Experiment = mean(x_exp, na.rm = TRUE),
    Mean_Control = mean(x_ctrl, na.rm = TRUE),
    Mean_Difference = mean(x_exp, na.rm = TRUE) - mean(x_ctrl, na.rm = TRUE),
    Paired_T_Test_P = paired_t_p,
    Paired_Wilcox_P = paired_w_p,
    Unpaired_T_Test_P = unpaired_t_p,
    Wilcox_P = wilcox_p,
    Test_Mode = ifelse(nrow(pair_wide) > 0, "paired", "unpaired"),
    Status = "OK",
    stringsAsFactors = FALSE
  )
}

get_significance_label <- function(p_value) {
  if (is.na(p_value)) return("NA")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  "ns"
}

run_pairwise_expression_tests <- function(dat, group_column, group_levels = NULL, paired_id_column = NULL) {
  groups <- as.character(dat[[group_column]])
  if (is.null(group_levels)) {
    group_levels <- unique(groups)
  }
  group_levels <- group_levels[group_levels %in% groups]

  if (length(group_levels) < 2) {
    return(data.frame())
  }

  pairs <- combn(group_levels, 2, simplify = FALSE)
  do.call(
    rbind,
    lapply(pairs, function(pair_value) {
      group_1 <- pair_value[1]
      group_2 <- pair_value[2]
      x1 <- dat$log2_tpm_plus_1[groups == group_1]
      x2 <- dat$log2_tpm_plus_1[groups == group_2]

      paired_t_p <- NA_real_
      paired_w_p <- NA_real_
      complete_pairs <- NA_integer_
      test_mode <- "unpaired"

      if (!is.null(paired_id_column) && paired_id_column %in% colnames(dat)) {
        pair_data <- data.frame(
          Pair_ID = dat[[paired_id_column]],
          Group = groups,
          Expression = dat$log2_tpm_plus_1,
          stringsAsFactors = FALSE
        )
        pair_data <- pair_data[pair_data$Group %in% c(group_1, group_2), , drop = FALSE]
        pair_wide <- reshape(
          pair_data,
          idvar = "Pair_ID",
          timevar = "Group",
          direction = "wide"
        )
        column_1 <- paste0("Expression.", group_1)
        column_2 <- paste0("Expression.", group_2)
        complete_index <- !is.na(pair_wide[[column_1]]) & !is.na(pair_wide[[column_2]])
        pair_wide <- pair_wide[complete_index, , drop = FALSE]
        complete_pairs <- nrow(pair_wide)

        if (complete_pairs > 0) {
          paired_t_p <- tryCatch(
            t.test(pair_wide[[column_2]], pair_wide[[column_1]], paired = TRUE)$p.value,
            error = function(e) NA_real_
          )
          paired_w_p <- tryCatch(
            wilcox.test(pair_wide[[column_2]], pair_wide[[column_1]], paired = TRUE, exact = FALSE)$p.value,
            error = function(e) NA_real_
          )
          test_mode <- "paired"
        }
      }

      unpaired_t_p <- tryCatch(t.test(x2, x1)$p.value, error = function(e) NA_real_)
      wilcox_p <- tryCatch(
        wilcox.test(x2, x1, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )

      label_p <- if (test_mode == "paired") paired_w_p else wilcox_p

      data.frame(
        Dataset = DATASET_ID,
        Target_Gene = TARGET_GENE,
        Group_Variable = group_column,
        Group_1 = group_1,
        Group_2 = group_2,
        Contrast = paste0(group_2, "_vs_", group_1),
        N_Group_1 = length(x1),
        N_Group_2 = length(x2),
        Complete_Pairs = complete_pairs,
        Mean_Group_1 = mean(x1, na.rm = TRUE),
        Mean_Group_2 = mean(x2, na.rm = TRUE),
        Mean_Difference_Group_2_minus_Group_1 = mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE),
        Paired_T_Test_P = paired_t_p,
        Paired_Wilcox_P = paired_w_p,
        Unpaired_T_Test_P = unpaired_t_p,
        Wilcox_P = wilcox_p,
        P_Value_For_Label = label_p,
        Significance_Label = get_significance_label(label_p),
        Test_Mode = test_mode,
        stringsAsFactors = FALSE
      )
    })
  )
}

make_significance_annotations <- function(pairwise_tests, group_levels, y_values) {
  if (nrow(pairwise_tests) == 0) {
    return(data.frame())
  }

  y_range <- diff(range(y_values, na.rm = TRUE))
  if (!is.finite(y_range) || y_range <= 0) {
    y_range <- 1
  }

  pairwise_tests$x_start <- match(pairwise_tests$Group_1, group_levels)
  pairwise_tests$x_end <- match(pairwise_tests$Group_2, group_levels)
  pairwise_tests$y_position <- max(y_values, na.rm = TRUE) +
    seq_len(nrow(pairwise_tests)) * y_range * 0.11
  pairwise_tests$label_y <- pairwise_tests$y_position + y_range * 0.025
  pairwise_tests$tip_y <- pairwise_tests$y_position - y_range * 0.025
  pairwise_tests
}

add_significance_annotations <- function(plot, annotation_table) {
  if (nrow(annotation_table) == 0) {
    return(plot)
  }

  plot +
    geom_segment(
      data = annotation_table,
      aes(x = x_start, xend = x_end, y = y_position, yend = y_position),
      inherit.aes = FALSE,
      linewidth = 0.32,
      color = TEXT_COLOR
    ) +
    geom_segment(
      data = annotation_table,
      aes(x = x_start, xend = x_start, y = tip_y, yend = y_position),
      inherit.aes = FALSE,
      linewidth = 0.32,
      color = TEXT_COLOR
    ) +
    geom_segment(
      data = annotation_table,
      aes(x = x_end, xend = x_end, y = tip_y, yend = y_position),
      inherit.aes = FALSE,
      linewidth = 0.32,
      color = TEXT_COLOR
    ) +
    geom_text(
      data = annotation_table,
      aes(x = (x_start + x_end) / 2, y = label_y, label = Significance_Label),
      inherit.aes = FALSE,
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      size = 3.4,
      color = TEXT_COLOR
    )
}


# 2. Read inputs ---------------------------------------------------------------

if (!file.exists(CLINICAL_EDIT_FILE)) {
  source(SAMPLE_DESIGN_SCRIPT)
}

dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)

se <- readRDS(SE_RDS_FILE)
clinical_data <- read.csv(
  CLINICAL_EDIT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(inherits(se, "SummarizedExperiment"))
stopifnot("Sample_ID" %in% colnames(clinical_data))
stopifnot(all(colnames(se) %in% clinical_data$Sample_ID))

clinical_data <- clinical_data[
  match(colnames(se), clinical_data$Sample_ID),
  ,
  drop = FALSE
]


# 3. Extract ATF3 expression ---------------------------------------------------

row_annotation <- as.data.frame(rowData(se), stringsAsFactors = FALSE)
target_index <- get_target_index(row_annotation, TARGET_GENE)

counts <- as.numeric(assay(se, "counts")[target_index, ])
tpm <- as.numeric(assay(se, "tpm")[target_index, ])
fpkm <- as.numeric(assay(se, "fpkm")[target_index, ])

expression_table <- data.frame(
  Dataset = DATASET_ID,
  Target_Gene = TARGET_GENE,
  Feature_ID = rownames(se)[target_index],
  GeneID = row_annotation$GeneID[target_index],
  Symbol = row_annotation$Symbol[target_index],
  Ensembl = row_annotation$Ensembl[target_index],
  Entrez = row_annotation$Entrez[target_index],
  Sample_ID = clinical_data$Sample_ID,
  Title = clinical_data$Title,
  Patient_ID = clinical_data$Patient_ID,
  tissue = clinical_data$tissue,
  tissue_class = clinical_data$tissue_class,
  ajcc_stage = clinical_data$ajcc_stage,
  counts = counts,
  fpkm = fpkm,
  tpm = tpm,
  log2_tpm_plus_1 = log2(tpm + 1),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# 4. Summaries and tests -------------------------------------------------------

group_summary <- summarize_by_group(expression_table, "tissue_class")
rownames(group_summary) <- NULL

analysis_columns <- grep("^analysis_", colnames(clinical_data), value = TRUE)
group_tests <- do.call(rbind, lapply(analysis_columns, run_group_test))
rownames(group_tests) <- NULL


# 5. Plot ----------------------------------------------------------------------

tissue_order <- c("Normal_colon", "Primary_tumor", "Liver_metastasis")
expression_table$tissue_class <- factor(
  expression_table$tissue_class,
  levels = tissue_order[tissue_order %in% expression_table$tissue_class]
)
tissue_levels <- levels(expression_table$tissue_class)
pairwise_tests <- run_pairwise_expression_tests(
  expression_table,
  group_column = "tissue_class",
  group_levels = tissue_levels,
  paired_id_column = "Patient_ID"
)
pairwise_plot_annotations <- make_significance_annotations(
  pairwise_tests,
  group_levels = tissue_levels,
  y_values = expression_table$log2_tpm_plus_1
)

fill_colors <- c(
  Normal_colon = "#4C78A8",
  Primary_tumor = "#8A7F2D",
  Liver_metastasis = "#C95F3F"
)

expression_plot <- ggplot(
  expression_table,
  aes(x = tissue_class, y = log2_tpm_plus_1, group = Patient_ID)
) +
  geom_line(color = "grey70", linewidth = 0.35, alpha = 0.75) +
  geom_boxplot(
    aes(fill = tissue_class, group = tissue_class),
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.72,
    color = TEXT_COLOR
  ) +
  geom_point(
    aes(fill = tissue_class),
    shape = 21,
    size = 2.8,
    stroke = 0.35,
    color = TEXT_COLOR,
    position = position_jitter(width = 0.08, height = 0, seed = 20260622),
    alpha = 0.95
  ) +
  scale_fill_manual(values = fill_colors) +
  scale_x_discrete(labels = function(x) wrap_label_by_underscore(x, width = 14)) +
  labs(
    x = NULL,
    y = paste0(TARGET_GENE, " log2(TPM + 1)"),
    fill = NULL
  ) +
  theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
    axis.text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
    axis.title = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
    legend.position = "none",
    text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
    plot.margin = margin(8, 12, 8, 8, unit = "pt")
  )

expression_plot <- add_significance_annotations(
  expression_plot,
  pairwise_plot_annotations
)

save_ggplot_pdf_png(
  plot = expression_plot,
  pdf_file = file.path(PLOT_ROOT, "ATF3_expression_by_tissue.pdf"),
  width = 7.4,
  height = 5.8
)


# 6. Save outputs --------------------------------------------------------------

write_csv_with_report_previews(
  expression_table,
  file.path(TABLE_ROOT, "ATF3_expression_by_sample.csv")
)
write_csv_with_report_previews(
  group_summary,
  file.path(TABLE_ROOT, "ATF3_group_summary.csv")
)
write_csv_with_report_previews(
  group_tests,
  file.path(TABLE_ROOT, "ATF3_group_tests.csv")
)
write_csv_with_report_previews(
  pairwise_tests,
  file.path(TABLE_ROOT, "ATF3_pairwise_tissue_class_tests.csv")
)

cat("\nGSE50760 ATF3 expression overview finished.\n")
cat("Output root: ", TABLE_ROOT, "\n\n", sep = "")
print(group_tests, row.names = FALSE)
