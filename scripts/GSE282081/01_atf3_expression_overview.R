# GSE282081 ATF3表达概览
#
# 读取00号脚本生成的样本设计，汇总ATF3在各培养/肝细胞共培养模型中的表达，
# 并输出样本级表格、分组统计、两组比较检验和表达图。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE282081"
DATA_TYPE <- "ngs"
TARGET_GENE <- "ATF3"

SE_RDS_FILE <- "data/ngs/GSE282081/data_prepare/GSE282081_se_raw.rds"
CLINICAL_EDIT_FILE <- "data/ngs/GSE282081/data_prepare/GSE282081_clinical_edit.csv"
SAMPLE_DESIGN_SCRIPT <- "scripts/GSE282081/00_prepare_sample_design.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables", "ATF3_expression")
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "ATF3_expression")

options(width = 200)


# 1. 加载包和函数 --------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggplot2)
})

source(REPORT_TABLE_FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)


# 2. 读取数据 ------------------------------------------------------------------

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


# 3. 提取ATF3表达 --------------------------------------------------------------

row_annotation <- as.data.frame(rowData(se), stringsAsFactors = FALSE)
target_index <- which(trimws(as.character(row_annotation$Symbol)) == TARGET_GENE)
stopifnot(length(target_index) >= 1)

if (length(target_index) > 1) {
  target_index <- target_index[1]
}

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
  culture_model = clinical_data$culture_model,
  treatment_class = clinical_data$treatment_class,
  liver_niche_model = clinical_data$liver_niche_model,
  hepatocyte_coculture = clinical_data$hepatocyte_coculture,
  J2_coculture = clinical_data$J2_coculture,
  counts = counts,
  fpkm = fpkm,
  tpm = tpm,
  log2_tpm_plus_1 = log2(tpm + 1),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# 4. 汇总和两组检验 ------------------------------------------------------------

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

group_summary <- do.call(
  rbind,
  lapply(
    c("culture_model", "treatment_class", "liver_niche_model", "hepatocyte_coculture", "J2_coculture"),
    function(group_column) summarize_by_group(expression_table, group_column)
  )
)
rownames(group_summary) <- NULL

analysis_columns <- grep("^analysis_", colnames(clinical_data), value = TRUE)

run_group_test <- function(analysis_column) {
  values <- trimws(as.character(clinical_data[[analysis_column]]))
  values[is.na(values)] <- ""
  keep <- values != ""
  values <- values[keep]
  y <- expression_table$log2_tpm_plus_1[keep]

  group_names <- unique(values)
  analysis_name <- sub("^analysis_", "", analysis_column)
  experiment_group <- analysis_name
  control_group <- setdiff(group_names, experiment_group)

  if (length(group_names) != 2 || !experiment_group %in% group_names || length(control_group) != 1) {
    return(data.frame(
      Analysis_Name = analysis_name,
      Contrast = NA,
      N_Experiment = NA,
      N_Control = NA,
      Mean_Experiment = NA,
      Mean_Control = NA,
      Mean_Difference = NA,
      T_Test_P = NA,
      Wilcox_P = NA,
      Status = "Skipped: expected exactly two groups with experiment group matching analysis name.",
      stringsAsFactors = FALSE
    ))
  }

  x_exp <- y[values == experiment_group]
  x_ctrl <- y[values == control_group]

  t_p <- tryCatch(t.test(x_exp, x_ctrl)$p.value, error = function(e) NA_real_)
  w_p <- tryCatch(
    wilcox.test(x_exp, x_ctrl, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )

  data.frame(
    Analysis_Name = analysis_name,
    Contrast = paste0(experiment_group, "_vs_", control_group),
    N_Experiment = length(x_exp),
    N_Control = length(x_ctrl),
    Mean_Experiment = mean(x_exp, na.rm = TRUE),
    Mean_Control = mean(x_ctrl, na.rm = TRUE),
    Mean_Difference = mean(x_exp, na.rm = TRUE) - mean(x_ctrl, na.rm = TRUE),
    T_Test_P = t_p,
    Wilcox_P = w_p,
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

run_pairwise_expression_tests <- function(dat, group_column, group_levels = NULL) {
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
      t_p <- tryCatch(t.test(x2, x1)$p.value, error = function(e) NA_real_)
      w_p <- tryCatch(
        wilcox.test(x2, x1, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )

      data.frame(
        Dataset = DATASET_ID,
        Target_Gene = TARGET_GENE,
        Group_Variable = group_column,
        Group_1 = group_1,
        Group_2 = group_2,
        Contrast = paste0(group_2, "_vs_", group_1),
        N_Group_1 = length(x1),
        N_Group_2 = length(x2),
        Mean_Group_1 = mean(x1, na.rm = TRUE),
        Mean_Group_2 = mean(x2, na.rm = TRUE),
        Mean_Difference_Group_2_minus_Group_1 = mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE),
        T_Test_P = t_p,
        Wilcox_P = w_p,
        P_Value_For_Label = w_p,
        Significance_Label = get_significance_label(w_p),
        Test_Mode = "unpaired",
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
    seq_len(nrow(pairwise_tests)) * y_range * 0.10
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
      size = 3.2,
      color = TEXT_COLOR
    )
}

group_tests <- do.call(rbind, lapply(analysis_columns, run_group_test))
rownames(group_tests) <- NULL

limitation_table <- data.frame(
  Dataset = DATASET_ID,
  Target_Gene = TARGET_GENE,
  Question = "ATF3 difference between liver metastasis tumor tissue and normal tissue",
  Status = "Not testable from current sample table",
  Reason = "Current GSE282081 local metadata contains in vitro SW480 2D/3D culture and coculture samples, but no normal tissue samples.",
  stringsAsFactors = FALSE
)


# 5. 绘图 ----------------------------------------------------------------------

model_order <- c(
  "2D_none",
  "2D_J2",
  "2D_J2_hepatocyte",
  "3D_organoid_none",
  "3D_organoid_J2_hepatocyte"
)
expression_table$liver_niche_model <- factor(
  expression_table$liver_niche_model,
  levels = model_order[model_order %in% expression_table$liver_niche_model]
)
model_levels <- levels(expression_table$liver_niche_model)
pairwise_tests <- run_pairwise_expression_tests(
  expression_table,
  group_column = "liver_niche_model",
  group_levels = model_levels
)
pairwise_plot_annotations <- make_significance_annotations(
  pairwise_tests,
  group_levels = model_levels,
  y_values = expression_table$log2_tpm_plus_1
)

fill_colors <- c(
  No_hepatocyte_coculture = "#7F8A6D",
  Hepatocyte_coculture = "#C95F3F"
)

expression_plot <- ggplot(
  expression_table,
  aes(
    x = liver_niche_model,
    y = log2_tpm_plus_1,
    fill = hepatocyte_coculture
  )
) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.70, color = TEXT_COLOR) +
  geom_point(
    aes(shape = culture_model),
    position = position_jitter(width = 0.10, height = 0, seed = 20260622),
    size = 3.0,
    alpha = 0.95,
    color = TEXT_COLOR
  ) +
  scale_fill_manual(values = fill_colors) +
  scale_x_discrete(labels = function(x) wrap_label_by_underscore(x, width = 16)) +
  labs(
    x = NULL,
    y = paste0(TARGET_GENE, " log2(TPM + 1)"),
    fill = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
    axis.text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
    axis.title = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
    legend.position = "right",
    legend.text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
    text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
    plot.margin = margin(8, 12, 8, 8, unit = "pt")
  )

expression_plot <- add_significance_annotations(
  expression_plot,
  pairwise_plot_annotations
)

save_ggplot_pdf_png(
  plot = expression_plot,
  pdf_file = file.path(PLOT_ROOT, "ATF3_expression_by_model.pdf"),
  width = 8.2,
  height = 6.8
)


# 6. 保存表格 ------------------------------------------------------------------

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
  file.path(TABLE_ROOT, "ATF3_pairwise_liver_niche_model_tests.csv")
)
write_csv_with_report_previews(
  limitation_table,
  file.path(TABLE_ROOT, "ATF3_unavailable_tissue_comparison.csv")
)

cat("\nGSE282081 ATF3 expression overview finished.\n")
cat("Output root: ", TABLE_ROOT, "\n\n", sep = "")
print(group_tests, row.names = FALSE)
