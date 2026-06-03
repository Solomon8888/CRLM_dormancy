# Common helper functions for limma differential expression analysis.

get_assay_matrix <- function(se, data_type) {
  assay_names <- names(SummarizedExperiment::assays(se))

  if (data_type == "ngs" && "counts" %in% assay_names) {
    return(as.matrix(SummarizedExperiment::assay(se, "counts")))
  }

  if ("expr" %in% assay_names) {
    return(as.matrix(SummarizedExperiment::assay(se, "expr")))
  }

  as.matrix(SummarizedExperiment::assay(se, 1))
}

sanitize_file_name <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- "analysis"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[x == ""] <- "analysis"
  x
}

get_analysis_designs <- function(clinical_data, prefix = "analysis_") {
  clinical_columns <- colnames(clinical_data)
  analysis_index <- which(startsWith(clinical_columns, prefix))

  if (length(analysis_index) == 0) {
    stop("No analysis design columns were found. Expected columns like analysis_OXLP.")
  }

  experiment_group <- sub(
    paste0("^", prefix),
    "",
    clinical_columns[analysis_index]
  )
  analysis_base_name <- sanitize_file_name(experiment_group)
  analysis_name <- analysis_base_name

  duplicated_name <- names(table(analysis_base_name))[table(analysis_base_name) > 1]

  for (name in duplicated_name) {
    duplicated_index <- which(analysis_base_name == name)
    analysis_name[duplicated_index] <- paste0(name, "_", seq_along(duplicated_index))
  }

  data.frame(
    Column_Index = analysis_index,
    Column_Name = clinical_columns[analysis_index],
    Analysis_Base_Name = analysis_base_name,
    Analysis_Name = analysis_name,
    Experiment_Group = experiment_group,
    stringsAsFactors = FALSE
  )
}

prepare_design_samples <- function(sample_info, group_column_index, experiment_group) {
  group_value <- trimws(as.character(sample_info[[group_column_index]]))
  group_value[is.na(group_value)] <- ""

  keep_sample <- group_value != ""
  sample_info_used <- sample_info[keep_sample, , drop = FALSE]
  group_value <- group_value[keep_sample]

  group_names <- unique(group_value)
  control_group <- setdiff(group_names, experiment_group)

  if (!experiment_group %in% group_names) {
    stop("Experiment group '", experiment_group, "' was not found in the design column.")
  }

  if (length(control_group) != 1) {
    stop(
      "Exactly one control group is required. Detected groups: ",
      paste(group_names, collapse = ", ")
    )
  }

  group_list <- factor(
    group_value,
    levels = c(control_group, experiment_group)
  )

  list(
    sample_info = sample_info_used,
    group_list = group_list,
    control_group = control_group
  )
}

needs_log2_transform <- function(exprSet) {
  q99 <- as.numeric(quantile(exprSet, 0.99, na.rm = TRUE))
  max_value <- max(exprSet, na.rm = TRUE)

  q99 > 50 || max_value > 100
}

get_distribution_diagnostics <- function(exprSet) {
  sample_median <- apply(exprSet, 2, median, na.rm = TRUE)
  sample_iqr <- apply(exprSet, 2, IQR, na.rm = TRUE)

  data.frame(
    Median_Spread = diff(range(sample_median)),
    IQR_Spread = diff(range(sample_iqr))
  )
}

prepare_microarray_data <- function(exprSet) {
  log2_transformed <- FALSE
  normalized_between_arrays <- FALSE

  if (needs_log2_transform(exprSet)) {
    if (min(exprSet, na.rm = TRUE) < 0) {
      stop("Expression values look unlogged, but contain negative values.")
    }

    exprSet <- log2(exprSet + 1)
    log2_transformed <- TRUE
  }

  diagnostics <- get_distribution_diagnostics(exprSet)

  if (diagnostics$Median_Spread > 0.5 || diagnostics$IQR_Spread > 0.5) {
    exprSet <- limma::normalizeBetweenArrays(exprSet, method = "quantile")
    normalized_between_arrays <- TRUE
    diagnostics <- get_distribution_diagnostics(exprSet)
  }

  list(
    data = exprSet,
    log2_transformed = log2_transformed,
    normalized_between_arrays = normalized_between_arrays,
    median_spread = diagnostics$Median_Spread,
    iqr_spread = diagnostics$IQR_Spread,
    filtered_genes = NA_integer_
  )
}

prepare_ngs_data <- function(counts, gene_annotation, group_list, design) {
  y <- edgeR::DGEList(
    counts = counts,
    group = group_list,
    genes = gene_annotation
  )

  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)

  v <- limma::voom(y, design = design, plot = FALSE)

  list(
    data = v,
    log2_transformed = NA,
    normalized_between_arrays = NA,
    median_spread = NA,
    iqr_spread = NA,
    filtered_genes = sum(!keep)
  )
}
