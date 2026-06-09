# IOBR分析公共函数
#
# 本文件服务于scripts/iobr目录下的IOBR快速分析脚本：
# - 统一项目路径、IOBR缓存、表格和图片保存逻辑；
# - 统一读取GSE114012 SummarizedExperiment与差异分析结果；
# - 统一调用parallel_runtime_functions.R的多进程进度条；
# - 每个IOBR模块独立记录success/skipped/failed，避免单个方法失败中断全流程。


# 0. 环境变量与路径 ------------------------------------------------------------

iobr_parse_env_vector <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  parts <- trimws(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
  parts[nzchar(parts)]
}

iobr_parse_env_logical <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(default)
  }

  value %in% c("1", "true", "t", "yes", "y")
}

iobr_parse_env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed)) {
    return(default)
  }

  parsed
}

iobr_find_project_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    marker <- file.path(current, "scripts", "functions", "iobr_common_functions.R")
    if (file.exists(marker)) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot locate project root from: ", getwd())
    }
    current <- parent
  }
}

iobr_source_project_functions <- function(project_root) {
  source(file.path(project_root, "scripts", "functions", "plotting_common_functions.R"))
  source(file.path(project_root, "scripts", "functions", "result_table_io_functions.R"))
  source(file.path(project_root, "scripts", "functions", "parallel_runtime_functions.R"))
  source(file.path(project_root, "scripts", "functions", "limma_de_functions.R"))
  invisible(TRUE)
}

iobr_sanitize <- function(x, default = "analysis") {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- default
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[x == ""] <- default
  x
}

iobr_format_p_value <- function(p_value, digits = 3L) {
  p_value <- suppressWarnings(as.numeric(p_value))
  vapply(p_value, function(p) {
    if (!is.finite(p)) {
      return("NA")
    }
    if (p < 0.001) {
      return("<0.001")
    }
    formatC(p, format = "f", digits = digits)
  }, character(1))
}

iobr_significance_stars <- function(p_value) {
  p_value <- suppressWarnings(as.numeric(p_value))
  ifelse(
    !is.finite(p_value), "",
    ifelse(
      p_value < 0.001, "***",
      ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns"))
    )
  )
}

iobr_make_p_label <- function(p_value, prefix = "p = ") {
  formatted <- iobr_format_p_value(p_value)
  formatted <- ifelse(grepl("^<", formatted), sub("^<", "< ", formatted), formatted)
  p_prefix <- ifelse(grepl("^<", formatted), "p ", prefix)
  stars <- iobr_significance_stars(p_value)
  label <- paste0(p_prefix, formatted)
  ifelse(nzchar(stars), paste(label, stars), label)
}

iobr_dataset_family <- function(dataset_id) {
  dataset_id <- as.character(dataset_id)[1]
  if (grepl("^TCGA_", dataset_id)) {
    return("TCGA")
  }
  if (grepl("^GTEx_", dataset_id)) {
    return("GTEX")
  }
  if (identical(dataset_id, "GSE114012") || grepl("^GSE114012", dataset_id)) {
    return("GSE114012")
  }
  "OTHER"
}

iobr_extract_design_from_stem <- function(file_stem) {
  stem <- iobr_sanitize(file_stem)
  if (grepl("_Target_Group($|_)", stem)) {
    return("Target_Group")
  }
  if (grepl("(^|_)analysis_", stem)) {
    design <- sub("^.*?(analysis_)", "\\1", stem)
    known_designs <- c(
      "analysis_DLD1_HCT15_SW48", "analysis_DLD1_HCT15",
      "analysis_SW948", "analysis_HCT15", "analysis_DLD1",
      "analysis_HT55", "analysis_SW48", "analysis_ALL", "analysis_RKO"
    )
    known_designs <- known_designs[order(nchar(known_designs), decreasing = TRUE)]
    matched_design <- known_designs[vapply(known_designs, function(candidate) {
      grepl(paste0("^", candidate, "($|_)"), design)
    }, logical(1))]
    if (length(matched_design) > 0L) {
      return(iobr_sanitize(matched_design[1], default = "analysis"))
    }
    design <- sub(
      "_(target_correlation|group_difference|group_tests|boxplot|pca|heatmap|all_features_heatmap|cell_barplot|scores|sig_gsea|batch_survival|target_group_roc|tme_cluster|cluster_composition|survival_forest|outlier_samples).*$",
      "",
      design
    )
    design <- sub("_(H|C[0-9]+|C[0-9]+_CP|C[0-9]+_CGP)$", "", design)
    return(iobr_sanitize(design, default = "analysis"))
  }
  "Target_Group"
}

iobr_extract_method_from_rest <- function(rest, module) {
  rest <- iobr_sanitize(rest)
  candidates <- c(
    "cibersort_abs", "mcpcounter", "quantiseq", "cibersort",
    "estimate", "timer", "xcell", "epic", "ips", "svr", "lsei",
    "integration", "ssgsea", "zscore", "pca",
    "C2_CGP", "C2_CP", "H", "C2", "C5"
  )
  candidates <- candidates[order(nchar(candidates), decreasing = TRUE)]
  matched <- candidates[vapply(candidates, function(candidate) {
    grepl(paste0("^", candidate, "($|_)"), rest)
  }, logical(1))]
  if (length(matched) > 0L) {
    return(iobr_sanitize(matched[1]))
  }
  parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
  if (length(parts) > 0L && nzchar(parts[1])) {
    return(iobr_sanitize(parts[1], default = module))
  }
  iobr_sanitize(module)
}

iobr_parse_output_context <- function(file_stem, module) {
  stem <- iobr_sanitize(file_stem)
  module <- iobr_sanitize(module)

  if (module %in% c("run_summary", "capability_catalog")) {
    return(list(
      is_global = TRUE,
      dataset_family = "00_metadata",
      dataset_id = "00_metadata",
      target_block = "global",
      design = "global",
      method = module
    ))
  }

  dataset_id <- "OTHER"
  target_gene <- NA_character_
  method <- module
  target_block <- "GENE_unknown"

  if (grepl("^TCGA_[A-Za-z0-9]+_01A_", stem)) {
    dataset_id <- sub("^(TCGA_[A-Za-z0-9]+_01A)_.*$", "\\1", stem)
    rest <- sub("^TCGA_[A-Za-z0-9]+_01A_", "", stem)
    parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
    target_gene <- parts[1]
    method <- iobr_extract_method_from_rest(sub(paste0("^", target_gene, "_?"), "", rest), module)
    target_block <- paste0("GENE_", iobr_sanitize(target_gene))
  } else if (grepl("^GTEx_[A-Za-z0-9]+_", stem)) {
    dataset_id <- sub("^(GTEx_[A-Za-z0-9]+)_.*$", "\\1", stem)
    rest <- sub("^GTEx_[A-Za-z0-9]+_", "", stem)
    parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
    target_gene <- parts[1]
    method <- iobr_extract_method_from_rest(sub(paste0("^", target_gene, "_?"), "", rest), module)
    target_block <- paste0("GENE_", iobr_sanitize(target_gene))
  } else if (grepl("^GSE114012_", stem)) {
    dataset_id <- "GSE114012"
    rest <- sub("^GSE114012_", "", stem)
    parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
    if (length(parts) >= 2L && !identical(parts[1], "analysis")) {
      target_gene <- parts[1]
      method <- iobr_extract_method_from_rest(sub(paste0("^", target_gene, "_?"), "", rest), module)
      target_block <- paste0("GENE_", iobr_sanitize(target_gene))
    } else {
      target_block <- "DEG"
      method <- sub("^.*_([A-Z][A-Za-z0-9]*)_sig_gsea.*$", "\\1", stem)
      if (identical(method, stem)) {
        method <- module
      }
    }
  }

  list(
    is_global = FALSE,
    dataset_family = iobr_dataset_family(dataset_id),
    dataset_id = iobr_sanitize(dataset_id),
    target_block = iobr_sanitize(target_block),
    design = iobr_extract_design_from_stem(stem),
    method = iobr_sanitize(method, default = module)
  )
}


# 1. 包安装、能力目录与输出 ----------------------------------------------------

iobr_install_if_needed <- function(auto_install = TRUE) {
  if (requireNamespace("IOBR", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  if (!auto_install) {
    stop(
      "IOBR is not installed. Install with: ",
      "remotes::install_github('IOBR/IOBR', dependencies = TRUE)"
    )
  }

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("IOBR/IOBR", dependencies = TRUE, upgrade = "never")

  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("IOBR installation was attempted but the package is still unavailable.")
  }

  invisible(TRUE)
}

iobr_setup_runtime <- function(
    project_root,
    iobr_cache_dir,
    auto_install = TRUE,
    parallel_backend = "auto",
    quiet_strategy = TRUE) {
  dir.create(iobr_cache_dir, recursive = TRUE, showWarnings = FALSE)
  options(IOBR.cache_dir = iobr_cache_dir)
  options(parallel_runtime_backend = parallel_backend)
  options(parallel_runtime_force_single_line_progress = TRUE)
  options(parallel_runtime_quiet_strategy = quiet_strategy)

  iobr_install_if_needed(auto_install = auto_install)
  suppressPackageStartupMessages(library(IOBR))

  if (exists("set_iobr_cache_dir", envir = asNamespace("IOBR"), mode = "function")) {
    try(IOBR::set_iobr_cache_dir(iobr_cache_dir), silent = TRUE)
  }

  invisible(TRUE)
}

iobr_prepare_output_tree <- function(output_root, temporary_root, clear_previous = TRUE) {
  if (clear_previous) {
    unlink(output_root, recursive = TRUE, force = TRUE)
    unlink(temporary_root, recursive = TRUE, force = TRUE)
  }

  dirs <- c(
    output_root,
    file.path(output_root, "tables"),
    file.path(output_root, "plots"),
    file.path(output_root, "logs"),
    temporary_root
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

iobr_write_csv <- function(dat, file, n_rows = 21, na = "NA") {
  write_csv_with_report_previews(dat, file, n_rows = n_rows, na = na)
}

iobr_clamp <- function(x, lower, upper, default) {
  x <- suppressWarnings(as.numeric(x)[1])
  if (!is.finite(x)) {
    x <- default
  }
  min(max(x, lower), upper)
}

iobr_dynamic_plot_size <- function(
    plot_type = "generic",
    n_items = NA_integer_,
    n_samples = NA_integer_,
    n_groups = NA_integer_,
    max_label_chars = NA_integer_,
    title = "",
    width = NULL,
    height = NULL) {
  plot_type <- iobr_sanitize(plot_type, default = "generic")
  n_items <- iobr_clamp(n_items, 1, 5000, 8)
  n_samples <- iobr_clamp(n_samples, 1, 5000, 20)
  n_groups <- iobr_clamp(n_groups, 1, 200, 2)
  max_label_chars <- iobr_clamp(max_label_chars, 1, 200, 24)

  title <- paste(as.character(title), collapse = " ")
  title_lines <- max(1L, length(strwrap(title, width = 64)))
  title_extra <- 0.25 + 0.30 * max(0L, title_lines - 1L)

  size <- switch(
    plot_type,
    bar = list(
      width = iobr_clamp(7.0 + 0.080 * max_label_chars, 8.8, 22.0, 10.0),
      height = iobr_clamp(3.4 + 0.42 * n_items + title_extra, 5.8, 24.0, 7.0)
    ),
    box = list(
      width = iobr_clamp(5.4 + 0.86 * n_groups + 0.028 * max_label_chars, 6.4, 14.0, 7.2),
      height = iobr_clamp(5.2 + title_extra, 6.0, 9.5, 6.4)
    ),
    pca = list(
      width = iobr_clamp(6.5 + 0.22 * n_groups + 0.012 * max_label_chars, 7.2, 12.0, 7.8),
      height = iobr_clamp(5.8 + title_extra, 6.4, 9.5, 6.8)
    ),
    heatmap = list(
      width = if (n_samples > 140) {
        iobr_clamp(8.0 + 0.018 * n_samples, 10.0, 22.0, 14.0)
      } else {
        iobr_clamp(5.8 + 0.075 * n_samples + 0.010 * max_label_chars, 7.0, 18.0, 10.0)
      },
      height = iobr_clamp(3.6 + 0.190 * n_items + title_extra, 5.8, 46.0, 8.0)
    ),
    all_signature_heatmap = list(
      width = if (n_samples > 140) {
        iobr_clamp(8.6 + 0.016 * n_samples, 11.0, 22.0, 15.0)
      } else {
        iobr_clamp(6.4 + 0.070 * n_samples + 0.010 * max_label_chars, 8.0, 18.0, 12.0)
      },
      height = iobr_clamp(4.4 + 0.095 * n_items + title_extra, 8.0, 52.0, 18.0)
    ),
    heatmap_chunk = list(
      width = if (n_samples > 140) {
        iobr_clamp(8.4 + 0.016 * n_samples, 10.5, 21.0, 14.0)
      } else {
        iobr_clamp(6.2 + 0.070 * n_samples + 0.010 * max_label_chars, 8.0, 18.0, 11.0)
      },
      height = iobr_clamp(3.8 + 0.260 * n_items + title_extra, 7.2, 18.0, 10.0)
    ),
    cell_bar = list(
      width = iobr_clamp(7.4 + 0.040 * max_label_chars, 8.6, 15.5, 9.5),
      height = iobr_clamp(4.8 + 0.38 * n_items + title_extra, 6.4, 18.0, 7.2)
    ),
    forest = list(
      width = iobr_clamp(7.4 + 0.045 * max_label_chars, 8.2, 15.5, 9.0),
      height = iobr_clamp(3.8 + 0.44 * n_items + title_extra, 6.2, 24.0, 8.0)
    ),
    composition = list(
      width = iobr_clamp(5.6 + 0.70 * n_groups + 0.018 * max_label_chars, 6.8, 12.5, 7.2),
      height = iobr_clamp(5.0 + title_extra, 5.8, 8.5, 6.0)
    ),
    list(
      width = iobr_clamp(7.2 + 0.015 * max_label_chars, 6.0, 14.0, 8.0),
      height = iobr_clamp(5.0 + title_extra, 4.8, 10.0, 5.8)
    )
  )

  if (!is.null(width)) {
    size$width <- iobr_clamp(width, 3.5, 30.0, size$width)
  }
  if (!is.null(height)) {
    size$height <- iobr_clamp(height, 3.5, 60.0, size$height)
  }

  size
}

iobr_common_ggplot_theme <- function(
    base_size = BASE_FONT_SIZE,
    axis_text_size = NULL,
    axis_title_size = NULL,
    title_size = NULL,
    legend_text_size = NULL,
    legend_title_size = NULL,
    rotate_x = FALSE,
    hide_x_text = FALSE,
    hide_y_text = FALSE,
    margin = ggplot2::margin(10, 12, 10, 10, unit = "pt")) {
  if (is.null(axis_text_size)) {
    axis_text_size <- base_size
  }
  if (is.null(axis_title_size)) {
    axis_title_size <- base_size + 1
  }
  if (is.null(title_size)) {
    title_size <- base_size + 1
  }
  if (is.null(legend_text_size)) {
    legend_text_size <- max(base_size - 1, 7)
  }
  if (is.null(legend_title_size)) {
    legend_title_size <- base_size
  }

  x_text <- if (hide_x_text) {
    ggplot2::element_blank()
  } else {
    ggplot2::element_text(
      color = TEXT_COLOR,
      face = TEXT_FONT_FACE,
      size = axis_text_size,
      angle = if (rotate_x) 90 else 0,
      hjust = if (rotate_x) 1 else 0.5,
      vjust = if (rotate_x) 0.5 else 0.5
    )
  }

  y_text <- if (hide_y_text) {
    ggplot2::element_blank()
  } else {
    ggplot2::element_text(
      color = TEXT_COLOR,
      face = TEXT_FONT_FACE,
      size = axis_text_size
    )
  }

  ggplot2::theme_bw(base_size = base_size, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      text = ggplot2::element_text(
        color = TEXT_COLOR,
        family = TEXT_FONT_FAMILY,
        face = TEXT_FONT_FACE
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = TEXT_FONT_FACE,
        size = title_size,
        margin = ggplot2::margin(b = 10, unit = "pt")
      ),
      axis.text.x = x_text,
      axis.text.y = y_text,
      axis.title = ggplot2::element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        size = axis_title_size
      ),
      axis.line = ggplot2::element_line(color = TEXT_COLOR, linewidth = AXIS_LINE_WIDTH),
      axis.ticks = ggplot2::element_line(color = TEXT_COLOR, linewidth = AXIS_LINE_WIDTH * 0.7),
      panel.border = ggplot2::element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
      panel.grid.major = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        size = legend_title_size
      ),
      legend.text = ggplot2::element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        size = legend_text_size
      ),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      strip.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      plot.margin = margin
    )
}

iobr_evenly_spaced_values <- function(values, max_labels) {
  values <- unique(as.character(values))
  n <- length(values)
  max_labels <- suppressWarnings(as.integer(max_labels)[1])
  if (!is.finite(max_labels) || max_labels <= 0L) {
    return(character(0))
  }

  if (n <= max_labels) {
    return(values)
  }

  unique(values[round(seq(1, n, length.out = max_labels))])
}

iobr_make_sparse_axis_breaks <- function(values, max_labels = 25L) {
  values <- as.character(values)
  n <- length(values)
  max_labels <- suppressWarnings(as.integer(max_labels)[1])
  if (!is.finite(max_labels) || max_labels <= 0L) {
    return(character(0))
  }
  if (n <= max_labels) {
    return(unique(values))
  }

  unique(values[round(seq(1, n, length.out = max_labels))])
}

iobr_max_sample_axis_labels <- function(n_samples, max_label_chars = 20L) {
  n_samples <- suppressWarnings(as.integer(n_samples)[1])
  max_label_chars <- suppressWarnings(as.integer(max_label_chars)[1])
  if (!is.finite(n_samples)) {
    n_samples <- 1L
  }
  if (!is.finite(max_label_chars)) {
    max_label_chars <- 20L
  }
  n_samples <- max(n_samples, 1L)
  max_label_chars <- max(max_label_chars, 1L)
  base <- if (n_samples <= 60L) {
    60L
  } else if (n_samples <= 90L) {
    32L
  } else if (n_samples <= 140L) {
    18L
  } else {
    0L
  }

  if (max_label_chars > 24L) {
    base <- if (base <= 0L) 0L else max(6L, floor(base * 0.65))
  }
  base
}

iobr_max_feature_axis_labels <- function(n_features, max_label_chars = 24L) {
  n_features <- suppressWarnings(as.integer(n_features)[1])
  max_label_chars <- suppressWarnings(as.integer(max_label_chars)[1])
  if (!is.finite(n_features)) {
    n_features <- 1L
  }
  if (!is.finite(max_label_chars)) {
    max_label_chars <- 24L
  }
  n_features <- max(n_features, 1L)
  max_label_chars <- max(max_label_chars, 1L)
  base <- if (n_features <= 45L) {
    n_features
  } else if (n_features <= 90L) {
    42L
  } else if (n_features <= 180L) {
    34L
  } else {
    28L
  }

  if (max_label_chars > 34L) {
    base <- max(18L, floor(base * 0.75))
  }
  base
}

iobr_heatmap_feature_chunk_size <- function(n_samples) {
  n_samples <- suppressWarnings(as.integer(n_samples)[1])
  if (!is.finite(n_samples)) {
    n_samples <- 1L
  }
  n_samples <- max(n_samples, 1L)
  if (n_samples > 140L) {
    return(35L)
  }
  42L
}

iobr_save_ggplot <- function(
    plot,
    plot_root,
    module,
    file_stem,
    width = NULL,
    height = NULL,
    plot_type = "generic",
    n_items = NA_integer_,
    n_samples = NA_integer_,
    n_groups = NA_integer_,
    max_label_chars = NA_integer_,
    title = file_stem) {
  pdf_dir <- file.path(plot_root, module, "pdf")
  png_dir <- file.path(plot_root, module, "png")
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

  file_stem <- iobr_sanitize(file_stem)
  pdf_file <- file.path(pdf_dir, paste0(file_stem, ".pdf"))
  size <- iobr_dynamic_plot_size(
    plot_type = plot_type,
    n_items = n_items,
    n_samples = n_samples,
    n_groups = n_groups,
    max_label_chars = max_label_chars,
    title = title,
    width = width,
    height = height
  )
  save_ggplot_pdf_png(plot, pdf_file = pdf_file, width = size$width, height = size$height)
}

iobr_empty_result <- function(module, task, status, message = "", output_file = "") {
  data.frame(
    Module = module,
    Task = task,
    Status = status,
    Message = as.character(message),
    Output_File = as.character(output_file),
    stringsAsFactors = FALSE
  )
}

iobr_safe_task <- function(expr, module, task) {
  tryCatch(
    expr,
    error = function(error) {
      iobr_empty_result(module, task, "failed", conditionMessage(error))
    }
  )
}

iobr_get_exported_functions <- function() {
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    return(character(0))
  }
  sort(getNamespaceExports("IOBR"))
}

iobr_get_tme_methods <- function() {
  if (requireNamespace("IOBR", quietly = TRUE) &&
      exists("tme_deconvolution_methods", envir = asNamespace("IOBR"))) {
    methods <- get("tme_deconvolution_methods", envir = asNamespace("IOBR"))
    return(unique(unname(as.character(methods))))
  }

  c(
    "mcpcounter", "epic", "xcell", "cibersort", "cibersort_abs",
    "ips", "estimate", "timer", "quantiseq", "svr", "lsei"
  )
}

iobr_get_signature_methods <- function() {
  if (requireNamespace("IOBR", quietly = TRUE) &&
      exists("signature_score_calculation_methods", envir = asNamespace("IOBR"))) {
    methods <- get("signature_score_calculation_methods", envir = asNamespace("IOBR"))
    return(unique(unname(as.character(methods))))
  }
  c("pca", "zscore", "ssgsea", "integration")
}

iobr_write_capability_catalog <- function(output_root, scope = "core") {
  exported <- iobr_get_exported_functions()
  core_functions <- c(
    "deconvo_tme", "deconvo_mcpcounter", "deconvo_epic", "deconvo_xcell",
    "deconvo_cibersort", "deconvo_estimate", "deconvo_ips",
    "deconvo_timer", "deconvo_quantiseq", "calculate_sig_score",
    "sigScore", "batch_wilcoxon", "batch_kruskal", "batch_cor",
    "sig_box", "sig_box_batch", "sig_gsea", "iobr_pca", "tme_cluster"
  )
  optional_functions <- c(
    "LR_cal", "sig_roc", "batch_surv", "sig_surv_plot", "PrognosticModel",
    "PrognosticResult", "RegressionResult", "lasso_select", "add_riskscore",
    "find_mutations", "make_mut_matrix", "RemoveBatchEffect", "remove_batcheffect",
    "iobr_deg", "limma.dif", "feature_select", "find_outlier_samples"
  )

  catalog <- data.frame(
    Function = sort(unique(c(core_functions, optional_functions, exported))),
    Exported_By_IOBR = sort(unique(c(core_functions, optional_functions, exported))) %in% exported,
    Planned_Scope = ifelse(
      sort(unique(c(core_functions, optional_functions, exported))) %in% core_functions,
      "core_directly_relevant",
      ifelse(
        sort(unique(c(core_functions, optional_functions, exported))) %in% optional_functions,
        "optional_or_input_dependent",
        "available_not_prioritized"
      )
    ),
    Current_Script = scope,
    stringsAsFactors = FALSE
  )

  iobr_write_csv(
    catalog,
    file.path(output_root, "tables", "capability_catalog", "000_iobr_function_catalog.csv")
  )
}


# 2. 数据读取与整理 ------------------------------------------------------------

iobr_load_gse114012_inputs <- function(
    se_file,
    clinical_file,
    assay_name = "tpm",
    symbol_column = "Symbol",
    biotype_filter = "protein_coding") {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package SummarizedExperiment is required.")
  }
  se <- readRDS(se_file)
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(assay_name %in% names(SummarizedExperiment::assays(se)))
  stopifnot("Sample_ID" %in% colnames(clinical))

  expr <- as.matrix(SummarizedExperiment::assay(se, assay_name))
  row_data <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE)

  missing_samples <- setdiff(colnames(expr), clinical$Sample_ID)
  if (length(missing_samples) > 0L) {
    stop("Samples missing from clinical table: ", paste(missing_samples, collapse = ", "))
  }
  clinical <- clinical[match(colnames(expr), clinical$Sample_ID), , drop = FALSE]

  expr_symbol <- iobr_collapse_expression_by_symbol(
    expr = expr,
    row_data = row_data,
    symbol_column = symbol_column,
    biotype_filter = biotype_filter
  )

  list(
    se = se,
    clinical = clinical,
    expression = expr_symbol,
    row_data = row_data,
    se_file = se_file,
    clinical_file = clinical_file,
    assay_name = assay_name
  )
}

iobr_collapse_expression_by_symbol <- function(
    expr,
    row_data,
    symbol_column = "Symbol",
    biotype_filter = "protein_coding") {
  stopifnot(nrow(expr) == nrow(row_data))
  stopifnot(symbol_column %in% colnames(row_data))

  keep <- !is.na(row_data[[symbol_column]]) & trimws(row_data[[symbol_column]]) != ""
  if (!is.null(biotype_filter) &&
      "Biotype" %in% colnames(row_data) &&
      length(biotype_filter) > 0L) {
    keep <- keep & row_data$Biotype %in% biotype_filter
  }

  expr <- expr[keep, , drop = FALSE]
  symbols <- trimws(as.character(row_data[[symbol_column]][keep]))
  symbols <- make.unique(symbols, sep = "__dup__")
  base_symbols <- sub("__dup__[0-9]+$", "", symbols)

  storage.mode(expr) <- "double"
  collapsed <- rowsum(expr, group = base_symbols, reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(factor(base_symbols, levels = rownames(collapsed))))
  collapsed <- sweep(collapsed, 1, counts, "/")
  collapsed <- collapsed[order(rownames(collapsed)), , drop = FALSE]
  collapsed
}

iobr_detect_sample_id_column <- function(col_data) {
  candidates <- c("Sample_ID", "sample_id", "ID", "Run", "Title")
  detected <- candidates[candidates %in% colnames(col_data)][1]
  if (is.na(detected)) {
    return(NULL)
  }
  detected
}

iobr_filter_samples <- function(
    col_data,
    sample_id,
    detail_values = NULL,
    detail_column = "group_detail",
    barcode_pattern = NULL) {
  keep <- rep(TRUE, nrow(col_data))
  detail_filter_applied <- FALSE

  if (!is.null(detail_values) &&
      length(detail_values) > 0L &&
      detail_column %in% colnames(col_data)) {
    detail <- tolower(trimws(as.character(col_data[[detail_column]])))
    keep <- keep & detail %in% tolower(detail_values)
    detail_filter_applied <- TRUE
  }

  if (!is.null(barcode_pattern) &&
      nzchar(barcode_pattern) &&
      (!detail_filter_applied || !any(keep)) &&
      length(sample_id) == nrow(col_data)) {
    keep <- grepl(barcode_pattern, sample_id, ignore.case = TRUE)
  }

  keep
}

iobr_prepare_local_se_input <- function(
    dataset_id,
    se_file,
    target_gene,
    project_id = dataset_id,
    assay_name = "tpm",
    deconvolution_assay = assay_name,
    score_assay = assay_name,
    score_log2 = TRUE,
    sample_detail_values = NULL,
    sample_detail_column = "group_detail",
    sample_barcode_pattern = NULL,
    symbol_column = "Symbol",
    biotype_filter = "protein_coding",
    high_group_label = "High",
    low_group_label = "Low",
    tumor = TRUE,
    timer_indication = NULL,
    clinical_file = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package SummarizedExperiment is required.")
  }
  if (!file.exists(se_file)) {
    stop("SE file not found: ", se_file)
  }

  se <- readRDS(se_file)
  stopifnot(inherits(se, "SummarizedExperiment"))
  assay_names <- names(SummarizedExperiment::assays(se))
  if (!deconvolution_assay %in% assay_names) {
    stop("Assay not found for deconvolution: ", deconvolution_assay)
  }
  if (!score_assay %in% assay_names) {
    stop("Assay not found for signature scoring: ", score_assay)
  }

  col_data <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  sample_id_column <- iobr_detect_sample_id_column(col_data)
  sample_id <- if (is.null(sample_id_column)) {
    colnames(se)
  } else {
    as.character(col_data[[sample_id_column]])
  }
  sample_id[is.na(sample_id) | sample_id == ""] <- colnames(se)[is.na(sample_id) | sample_id == ""]
  stopifnot(!any(duplicated(sample_id)))

  keep_sample <- iobr_filter_samples(
    col_data = col_data,
    sample_id = sample_id,
    detail_values = sample_detail_values,
    detail_column = sample_detail_column,
    barcode_pattern = sample_barcode_pattern
  )
  if (!any(keep_sample)) {
    stop("No samples remain after filtering for dataset: ", dataset_id)
  }

  deconvolution_expr <- as.matrix(SummarizedExperiment::assay(se, deconvolution_assay))
  score_expr <- as.matrix(SummarizedExperiment::assay(se, score_assay))
  row_data <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE)

  deconvolution_expr <- deconvolution_expr[, keep_sample, drop = FALSE]
  score_expr <- score_expr[, keep_sample, drop = FALSE]
  col_data <- col_data[keep_sample, , drop = FALSE]
  sample_id <- sample_id[keep_sample]
  colnames(deconvolution_expr) <- sample_id
  colnames(score_expr) <- sample_id

  deconvolution_expr <- iobr_collapse_expression_by_symbol(
    expr = deconvolution_expr,
    row_data = row_data,
    symbol_column = symbol_column,
    biotype_filter = biotype_filter
  )
  score_expr <- iobr_collapse_expression_by_symbol(
    expr = score_expr,
    row_data = row_data,
    symbol_column = symbol_column,
    biotype_filter = biotype_filter
  )

  if (score_log2) {
    score_expr <- log2(score_expr + 1)
  }

  target_gene <- trimws(as.character(target_gene)[1])
  if (!target_gene %in% rownames(score_expr)) {
    stop("Target gene not found in ", dataset_id, ": ", target_gene)
  }

  target_expression <- as.numeric(score_expr[target_gene, ])
  target_median <- median(target_expression, na.rm = TRUE)
  target_group <- ifelse(
    target_expression >= target_median,
    high_group_label,
    low_group_label
  )
  target_group <- factor(target_group, levels = c(low_group_label, high_group_label))

  pdata <- cbind(
    data.frame(
      ID = sample_id,
      Sample_ID = sample_id,
      Dataset_ID = dataset_id,
      ProjectID = project_id,
      Target_Gene = target_gene,
      Target_Expression = target_expression,
      Target_Median = target_median,
      Target_Group = as.character(target_group),
      stringsAsFactors = FALSE
    ),
    col_data,
    stringsAsFactors = FALSE
  )
  pdata <- pdata[, !duplicated(colnames(pdata)), drop = FALSE]

  if (!is.null(clinical_file) && file.exists(clinical_file)) {
    clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
    if ("Sample_ID" %in% colnames(clinical)) {
      pdata <- merge(pdata, clinical, by = "Sample_ID", all.x = TRUE, suffixes = c("", "_clinical"))
      pdata <- pdata[match(sample_id, pdata$Sample_ID), , drop = FALSE]
    }
  }

  list(
    dataset_id = dataset_id,
    project_id = project_id,
    se_file = se_file,
    assay_name = assay_name,
    deconvolution_assay = deconvolution_assay,
    score_assay = score_assay,
    score_log2 = score_log2,
    deconvolution_expr = deconvolution_expr,
    score_expr = score_expr,
    pdata = pdata,
    target_gene = target_gene,
    tumor = tumor,
    timer_indication = timer_indication,
    sample_count = ncol(score_expr),
    gene_count = nrow(score_expr)
  )
}

iobr_add_survival_columns <- function(pdata) {
  pdata <- as.data.frame(pdata, stringsAsFactors = FALSE, check.names = FALSE)
  if ("days_to_death" %in% colnames(pdata) || "days_to_last_follow_up" %in% colnames(pdata)) {
    death <- if ("days_to_death" %in% colnames(pdata)) {
      suppressWarnings(as.numeric(pdata$days_to_death))
    } else {
      rep(NA_real_, nrow(pdata))
    }
    follow <- if ("days_to_last_follow_up" %in% colnames(pdata)) {
      suppressWarnings(as.numeric(pdata$days_to_last_follow_up))
    } else {
      rep(NA_real_, nrow(pdata))
    }
    pdata$time <- ifelse(is.na(death), follow, death)
  }

  if ("OS_status" %in% colnames(pdata)) {
    pdata$status <- suppressWarnings(as.numeric(pdata$OS_status))
  } else if ("vital_status" %in% colnames(pdata)) {
    vital <- tolower(trimws(as.character(pdata$vital_status)))
    pdata$status <- ifelse(vital %in% c("dead", "deceased", "1"), 1, 0)
  }

  pdata
}

iobr_bind_rows <- function(items) {
  items <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, items)
  if (length(items) == 0L) {
    return(data.frame())
  }
  dplyr::bind_rows(items)
}

iobr_feature_columns <- function(dat, exclude = character(0)) {
  default_exclude <- c(
    "ID", "Sample_ID", "Dataset_ID", "ProjectID", "Target_Gene",
    "Target_Expression", "Target_Median", "Target_Group",
    "time", "status"
  )
  iobr_numeric_feature_columns(dat, exclude = unique(c(default_exclude, exclude)))
}

iobr_status_row <- function(
    dataset_id,
    module,
    task,
    status,
    message = "",
    output_file = "",
    n_rows = NA_integer_,
    n_cols = NA_integer_) {
  data.frame(
    Dataset_ID = dataset_id,
    Module = module,
    Task = task,
    Status = status,
    Message = as.character(message),
    Rows = n_rows,
    Columns = n_cols,
    Output_File = as.character(output_file),
    stringsAsFactors = FALSE
  )
}

iobr_prepare_signature_collection <- function(selected_signatures = NULL) {
  data(signature_collection, package = "IOBR", envir = environment())
  signatures <- signature_collection
  if (!is.null(selected_signatures) && length(selected_signatures) > 0L) {
    selected_signatures <- intersect(selected_signatures, names(signatures))
    signatures <- signatures[selected_signatures]
  }
  signatures
}

iobr_load_lm22_reference <- function() {
  ref <- IOBR::load_data("lm22")
  as.data.frame(ref, check.names = FALSE)
}

iobr_make_sample_pdata <- function(clinical) {
  pdata <- clinical
  pdata$ID <- pdata$Sample_ID
  pdata
}

iobr_get_analysis_design_table <- function(clinical) {
  get_analysis_designs(clinical)
}

iobr_make_design_group_table <- function(clinical, design_row) {
  design_samples <- prepare_design_samples(
    sample_info = clinical,
    group_column_index = design_row$Column_Index,
    experiment_group = design_row$Experiment_Group
  )
  data.frame(
    ID = design_samples$sample_info$Sample_ID,
    Analysis_Name = design_row$Analysis_Name,
    Group = as.character(design_samples$group_list),
    Control_Group = design_samples$control_group,
    Experiment_Group = design_row$Experiment_Group,
    stringsAsFactors = FALSE
  )
}

iobr_get_deg_file_info_safe <- function(table_root) {
  get_deg_file_info(table_root)
}

iobr_read_deg_for_sig_gsea <- function(file) {
  dat <- read_report_csv(file)
  symbol_col <- intersect(c("Symbol", "symbol", "Gene", "gene"), colnames(dat))[1]
  logfc_col <- intersect(c("logFC", "log2FoldChange", "avg_log2FC"), colnames(dat))[1]
  if (is.na(symbol_col) || is.na(logfc_col)) {
    stop("DEG file lacks Symbol/logFC columns: ", file)
  }
  dat <- dat[!is.na(dat[[symbol_col]]) & dat[[symbol_col]] != "", , drop = FALSE]
  data.frame(
    symbol = as.character(dat[[symbol_col]]),
    log2FoldChange = as.numeric(dat[[logfc_col]]),
    stringsAsFactors = FALSE
  )
}


# 3. 结果矩阵、统计和绘图 ------------------------------------------------------

iobr_normalize_id_column <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  id_candidates <- intersect(c("ID", "Sample", "sample", "Sample_ID", "Mixture"), colnames(dat))
  if (length(id_candidates) == 0L) {
    dat$ID <- rownames(dat)
  } else if (!identical(id_candidates[1], "ID")) {
    colnames(dat)[match(id_candidates[1], colnames(dat))] <- "ID"
  }
  dat$ID <- as.character(dat$ID)
  dat
}

iobr_numeric_feature_columns <- function(dat, exclude = character(0)) {
  dat <- as.data.frame(dat)
  numeric_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  setdiff(numeric_cols, exclude)
}

iobr_merge_feature_tables <- function(tables) {
  tables <- Filter(function(x) is.data.frame(x) && "ID" %in% colnames(x), tables)
  if (length(tables) == 0L) {
    return(data.frame())
  }
  Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), tables)
}

iobr_add_atf3_expression <- function(feature_table, expr, target_gene) {
  if (!target_gene %in% rownames(expr)) {
    stop("Target gene not found in expression matrix: ", target_gene)
  }
  atf3 <- data.frame(
    ID = colnames(expr),
    Target_Gene = target_gene,
    Target_Expression = as.numeric(expr[target_gene, ]),
    stringsAsFactors = FALSE
  )
  merge(atf3, feature_table, by = "ID", all.x = TRUE)
}

iobr_batch_cor <- function(dat, target, features, method = "spearman") {
  if (requireNamespace("IOBR", quietly = TRUE)) {
    res <- try(IOBR::batch_cor(data = dat, target = target, feature = features, method = method), silent = TRUE)
    if (!inherits(res, "try-error")) {
      res <- as.data.frame(res, stringsAsFactors = FALSE)
      res$Method_Source <- "IOBR::batch_cor"
      return(res)
    }
  }

  rows <- lapply(features, function(feature) {
    ok <- is.finite(dat[[target]]) & is.finite(dat[[feature]])
    if (sum(ok) < 3L) {
      return(NULL)
    }
    test <- suppressWarnings(cor.test(dat[[target]][ok], dat[[feature]][ok], method = method))
    data.frame(
      sig_names = feature,
      p.value = test$p.value,
      statistic = unname(test$estimate),
      Method_Source = "stats::cor.test",
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, rows)
  if (is.null(res)) {
    return(data.frame())
  }
  res$p.adj <- p.adjust(res$p.value, method = "BH")
  res[order(res$p.value), , drop = FALSE]
}

iobr_group_test <- function(dat, group_col, features) {
  prepared <- iobr_prepare_group_data(
    dat = dat,
    group_col = group_col,
    require_two_groups = iobr_is_gse_analysis_group(group_col)
  )
  dat <- prepared$data
  test_group_col <- prepared$group_col
  groups <- prepared$levels
  if (length(groups) < 2L) {
    return(data.frame())
  }

  if (requireNamespace("IOBR", quietly = TRUE)) {
    res <- tryCatch({
      if (length(groups) == 2L) {
        IOBR::batch_wilcoxon(data = dat, target = test_group_col, feature = features)
      } else {
        IOBR::batch_kruskal(data = dat, group = test_group_col, feature = features)
      }
    }, error = function(error) error)
    if (!inherits(res, "error")) {
      res <- as.data.frame(res, stringsAsFactors = FALSE)
      res$Method_Source <- if (length(groups) == 2L) "IOBR::batch_wilcoxon" else "IOBR::batch_kruskal"
      res$Comparison <- paste(groups, collapse = " vs ")
      return(res)
    }
  }

  rows <- lapply(features, function(feature) {
    formula <- stats::as.formula(paste0("`", feature, "` ~ `", test_group_col, "`"))
    if (length(groups) == 2L) {
      test <- suppressWarnings(wilcox.test(formula, data = dat))
      stat <- unname(test$statistic)
      source <- "stats::wilcox.test"
    } else {
      test <- suppressWarnings(kruskal.test(formula, data = dat))
      stat <- unname(test$statistic)
      source <- "stats::kruskal.test"
    }
    data.frame(
      sig_names = feature,
      p.value = test$p.value,
      statistic = stat,
      Comparison = paste(groups, collapse = " vs "),
      Method_Source = source,
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, rows)
  res$p.adj <- p.adjust(res$p.value, method = "BH")
  res[order(res$p.value), , drop = FALSE]
}

iobr_make_top_barplot <- function(
    dat,
    feature_col,
    value_col,
    title,
    xlab = NULL,
    top_n = 30,
    label_col = NULL,
    ylab = value_col) {
  dat <- dat[is.finite(dat[[value_col]]), , drop = FALSE]
  dat <- dat[order(abs(dat[[value_col]]), decreasing = TRUE), , drop = FALSE]
  dat <- head(dat, top_n)
  dat[[feature_col]] <- factor(dat[[feature_col]], levels = rev(dat[[feature_col]]))
  if (!is.null(label_col) && label_col %in% colnames(dat)) {
    dat$IOBR_Bar_Label <- as.character(dat[[label_col]])
  } else {
    dat$IOBR_Bar_Label <- ""
  }

  ggplot2::ggplot(dat, ggplot2::aes(x = .data[[feature_col]], y = .data[[value_col]])) +
    ggplot2::geom_col(fill = "#2166AC", alpha = 0.86) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = IOBR_Bar_Label,
        hjust = ifelse(.data[[value_col]] >= 0, -0.03, 1.03)
      ),
      size = 4.2,
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      color = TEXT_COLOR,
      na.rm = TRUE
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_x_discrete(labels = function(x) wrap_label_by_underscore(x, width = 30)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.14, 0.48))) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    iobr_common_ggplot_theme(
      base_size = BASE_FONT_SIZE,
      axis_text_size = 12,
      axis_title_size = 14,
      title_size = 15,
      margin = ggplot2::margin(12, 70, 12, 18, unit = "pt")
    )
}

iobr_make_feature_boxplot <- function(dat, feature, group_col, title) {
  prepared <- iobr_prepare_group_data(
    dat = dat,
    group_col = group_col,
    require_two_groups = iobr_is_gse_analysis_group(group_col)
  )
  dat <- prepared$data
  plot_group_col <- prepared$group_col
  groups <- prepared$levels
  if (nrow(dat) == 0L || length(groups) < 2L) {
    stop("No valid groups for boxplot: ", group_col)
  }

  dat[[feature]] <- suppressWarnings(as.numeric(dat[[feature]]))
  dat <- dat[is.finite(dat[[feature]]) & !is.na(dat[[plot_group_col]]), , drop = FALSE]
  if (nrow(dat) == 0L) {
    stop("No finite values for boxplot feature: ", feature)
  }

  formula <- stats::as.formula(paste0("`", feature, "` ~ `", plot_group_col, "`"))
  p_value <- if (length(groups) == 2L) {
    suppressWarnings(stats::wilcox.test(formula, data = dat)$p.value)
  } else {
    suppressWarnings(stats::kruskal.test(formula, data = dat)$p.value)
  }
  p_label <- iobr_make_p_label(p_value)
  y_range <- range(dat[[feature]], na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- max(abs(y_range), 1, na.rm = TRUE)
  }
  y_pos <- y_range[2] + y_span * 0.12

  ggplot2::ggplot(dat, ggplot2::aes(x = .data[[plot_group_col]], y = .data[[feature]], fill = .data[[plot_group_col]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.82) +
    ggplot2::geom_jitter(width = 0.12, size = 1.8, alpha = 0.78) +
    ggplot2::annotate(
      "text",
      x = mean(seq_along(groups)),
      y = y_pos,
      label = p_label,
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      size = 4.8,
      color = TEXT_COLOR
    ) +
    ggplot2::expand_limits(y = y_pos + y_span * 0.10) +
    ggplot2::labs(x = iobr_group_axis_label(group_col), y = feature, title = title) +
    ggplot2::scale_x_discrete(labels = function(x) wrap_label(x, width = 18)) +
    iobr_common_ggplot_theme(
      base_size = BASE_FONT_SIZE + 1,
      axis_text_size = 12,
      axis_title_size = 14,
      title_size = 15,
      margin = ggplot2::margin(16, 18, 14, 18, unit = "pt")
    ) +
    ggplot2::theme(
      legend.position = "none"
    )
}

iobr_make_pca_plot <- function(dat, features, color_col, title) {
  features <- features[vapply(features, function(feature) {
    values <- suppressWarnings(as.numeric(dat[[feature]]))
    sum(is.finite(values)) >= 3L && stats::sd(values, na.rm = TRUE) > 0
  }, logical(1))]
  if (length(features) < 2L) {
    stop("Too few non-constant features for PCA.")
  }

  complete <- stats::complete.cases(dat[, features, drop = FALSE])
  pca_data <- dat[complete, c("ID", color_col, features), drop = FALSE]
  if (nrow(pca_data) < 3L || length(features) < 2L) {
    stop("Too few complete rows/features for PCA.")
  }
  matrix_data <- scale(as.matrix(pca_data[, features, drop = FALSE]))
  pca <- prcomp(matrix_data, center = FALSE, scale. = FALSE)
  plot_data <- data.frame(
    ID = pca_data$ID,
    Group = pca_data[[color_col]],
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

  ggplot2::ggplot(plot_data, ggplot2::aes(PC1, PC2, color = Group)) +
    ggplot2::geom_point(size = 3.2, alpha = 0.86) +
    ggplot2::labs(
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)"),
      title = title
    ) +
    iobr_common_ggplot_theme(
      base_size = BASE_FONT_SIZE,
      axis_text_size = 10.5,
      axis_title_size = 13,
      title_size = 13
    )
}

iobr_make_heatmap_plot <- function(dat, features, sample_order, title) {
  features <- features[features %in% colnames(dat)]
  dat <- dat[match(sample_order, dat$ID), c("ID", features), drop = FALSE]
  sample_label_count <- iobr_max_sample_axis_labels(
    length(sample_order),
    max(nchar(as.character(sample_order)), na.rm = TRUE)
  )
  sample_breaks <- iobr_make_sparse_axis_breaks(
    dat$ID,
    max_labels = sample_label_count
  )
  feature_label_count <- iobr_max_feature_axis_labels(
    length(features),
    max(nchar(as.character(features)), na.rm = TRUE)
  )
  feature_breaks <- iobr_make_sparse_axis_breaks(
    features,
    max_labels = feature_label_count
  )
  show_sample_labels <- length(sample_breaks) > 0L
  show_feature_labels <- length(feature_breaks) > 0L
  row_text_size <- iobr_clamp(
    9.2 - 0.030 * length(feature_breaks),
    4.8,
    8.4,
    6.8
  )
  col_text_size <- iobr_clamp(
    8.0 - 0.030 * length(sample_breaks),
    5.2,
    7.4,
    6.2
  )
  x_axis_note <- if (!show_sample_labels) {
    "; sample IDs hidden to avoid overlap"
  } else if (length(sample_breaks) < length(sample_order)) {
    "; sample IDs shown sparsely"
  } else {
    ""
  }
  y_axis_note <- if (length(feature_breaks) < length(features)) {
    "; feature labels shown sparsely"
  } else {
    ""
  }
  bottom_margin <- if (show_sample_labels) 42 else 18
  left_margin <- if (show_feature_labels) 24 else 16
  long <- reshape(
    dat,
    varying = features,
    v.names = "Value",
    timevar = "Feature",
    times = features,
    direction = "long"
  )
  long$ID <- factor(long$ID, levels = sample_order)
  long$Feature <- factor(long$Feature, levels = rev(features))

  ggplot2::ggplot(long, ggplot2::aes(x = ID, y = Feature, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(
      breaks = sample_breaks,
      labels = function(x) wrap_label_by_underscore(x, width = 16)
    ) +
    ggplot2::scale_y_discrete(
      breaks = feature_breaks,
      labels = function(x) wrap_label_by_underscore(x, width = 34)
    ) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D73027") +
    ggplot2::labs(
      x = paste0("Samples ordered by group and target expression (n = ", length(sample_order), x_axis_note, ")"),
      y = paste0("IOBR features / signatures (n = ", length(features), y_axis_note, ")"),
      title = title,
      fill = "Score"
    ) +
    iobr_common_ggplot_theme(
      base_size = 10,
      axis_text_size = 7,
      axis_title_size = 11.5,
      title_size = 12,
      rotate_x = TRUE,
      hide_x_text = !show_sample_labels,
      hide_y_text = !show_feature_labels,
      margin = ggplot2::margin(14, 16, bottom_margin, left_margin, unit = "pt")
    ) +
    ggplot2::theme(
      axis.text.x = if (show_sample_labels) {
        ggplot2::element_text(
          color = TEXT_COLOR,
          face = TEXT_FONT_FACE,
          size = col_text_size,
          angle = 90,
          hjust = 1,
          vjust = 0.5
        )
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x = if (show_sample_labels) {
        ggplot2::element_line(color = TEXT_COLOR, linewidth = AXIS_LINE_WIDTH * 0.7)
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = if (show_feature_labels) {
        ggplot2::element_text(size = row_text_size, face = TEXT_FONT_FACE)
      } else {
        ggplot2::element_blank()
      },
      panel.grid = ggplot2::element_blank()
    )
}

iobr_make_cell_fraction_summary_plot <- function(
    feature_table,
    features,
    group_col = "Target_Group",
    title,
    top_n = 20L) {
  features <- intersect(features, colnames(feature_table))
  features <- head(features, min(top_n, length(features)))
  if (length(features) == 0L) {
    stop("No features available for cell fraction summary plot.")
  }

  if (!group_col %in% colnames(feature_table)) {
    feature_table$All_Samples <- "All samples"
    group_col <- "All_Samples"
  }

  dat <- feature_table[, c("ID", group_col, features), drop = FALSE]
  long <- reshape(
    dat,
    varying = features,
    v.names = "Value",
    timevar = "Feature",
    times = features,
    direction = "long"
  )
  long$Value <- suppressWarnings(as.numeric(long$Value))
  long[[group_col]] <- as.character(long[[group_col]])
  long <- long[is.finite(long$Value) & !is.na(long[[group_col]]), , drop = FALSE]
  if (nrow(long) == 0L) {
    stop("No finite feature values available for cell fraction summary plot.")
  }

  summary_dat <- stats::aggregate(
    Value ~ Feature + Group,
    data = data.frame(
      Feature = long$Feature,
      Group = long[[group_col]],
      Value = long$Value,
      stringsAsFactors = FALSE
    ),
    FUN = stats::median,
    na.rm = TRUE
  )
  summary_dat$Abs_Value <- abs(summary_dat$Value)
  feature_order <- stats::aggregate(
    Abs_Value ~ Feature,
    data = summary_dat,
    FUN = max,
    na.rm = TRUE
  )
  feature_order <- feature_order$Feature[order(feature_order[[2]], decreasing = TRUE)]
  summary_dat$Feature <- factor(summary_dat$Feature, levels = rev(feature_order))

  ggplot2::ggplot(summary_dat, ggplot2::aes(x = Feature, y = Value, fill = Group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.64, alpha = 0.88) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(labels = function(x) wrap_label_by_underscore(x, width = 34)) +
    ggplot2::labs(
      x = "TME features / scores",
      y = "Median score / fraction",
      fill = group_col,
      title = title
    ) +
    iobr_common_ggplot_theme(
      base_size = BASE_FONT_SIZE,
      axis_text_size = 10,
      axis_title_size = 12.5,
      title_size = 13,
      legend_text_size = 9.5,
      legend_title_size = 10.5,
      margin = ggplot2::margin(12, 18, 12, 18, unit = "pt")
    )
}

iobr_feature_table_from_score <- function(score_table, prefix = NULL) {
  score_table <- iobr_normalize_id_column(score_table)
  numeric_cols <- iobr_numeric_feature_columns(score_table)
  if (!is.null(prefix) && nzchar(prefix)) {
    rename_cols <- setdiff(colnames(score_table), "ID")
    colnames(score_table)[match(rename_cols, colnames(score_table))] <- paste0(prefix, "__", rename_cols)
  }
  score_table
}


# 4. IOBR快速分析缓存与通用任务函数 ------------------------------------------

iobr_cache_write <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qs_save(object, file)
  } else {
    saveRDS(object, file)
  }
  invisible(file)
}

iobr_cache_read <- function(file) {
  if (requireNamespace("qs2", quietly = TRUE)) {
    return(qs2::qs_read(file))
  }
  readRDS(file)
}

iobr_make_cache_file <- function(cache_root, file_stem) {
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  extension <- if (requireNamespace("qs2", quietly = TRUE)) ".qs2" else ".rds"
  file.path(cache_root, paste0(iobr_sanitize(file_stem), extension))
}

iobr_table_dir <- function(output_root, module, file_stem = NULL) {
  module <- iobr_sanitize(module)
  if (is.null(file_stem) || module %in% c("run_summary", "capability_catalog")) {
    path <- file.path(output_root, "tables", module)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    return(path)
  }

  context <- iobr_parse_output_context(file_stem, module)
  path <- file.path(
    output_root,
    context$dataset_family,
    context$dataset_id,
    context$target_block,
    context$design,
    module,
    context$method,
    "tables"
  )
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

iobr_plot_dir <- function(output_root, module, format = c("pdf", "png"), file_stem = NULL) {
  format <- match.arg(format)
  module <- iobr_sanitize(module)
  if (is.null(file_stem) || module %in% c("run_summary", "capability_catalog")) {
    path <- file.path(output_root, "plots", module, format)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    return(path)
  }

  context <- iobr_parse_output_context(file_stem, module)
  path <- file.path(
    output_root,
    context$dataset_family,
    context$dataset_id,
    context$target_block,
    context$design,
    module,
    context$method,
    "plots",
    format
  )
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

iobr_write_module_csv <- function(dat, output_root, module, file_stem, n_rows = 21) {
  iobr_write_csv(
    dat,
    file.path(iobr_table_dir(output_root, module, file_stem), paste0(iobr_sanitize(file_stem), ".csv")),
    n_rows = n_rows,
    na = "NA"
  )
}

iobr_save_module_plot <- function(
    plot,
    output_root,
    module,
    file_stem,
    width = NULL,
    height = NULL,
    plot_type = "generic",
    n_items = NA_integer_,
    n_samples = NA_integer_,
    n_groups = NA_integer_,
    max_label_chars = NA_integer_,
    title = file_stem) {
  pdf_file <- file.path(
    iobr_plot_dir(output_root, module, "pdf", file_stem),
    paste0(iobr_sanitize(file_stem), ".pdf")
  )
  size <- iobr_dynamic_plot_size(
    plot_type = plot_type,
    n_items = n_items,
    n_samples = n_samples,
    n_groups = n_groups,
    max_label_chars = max_label_chars,
    title = title,
    width = width,
    height = height
  )
  save_ggplot_pdf_png(plot, pdf_file = pdf_file, width = size$width, height = size$height)
}

iobr_append_title <- function(plot, title) {
  if (inherits(plot, "ggplot")) {
    return(
      plot +
        ggplot2::ggtitle(title) +
        iobr_common_ggplot_theme(
          base_size = BASE_FONT_SIZE + 1,
          axis_text_size = 12,
          axis_title_size = 14,
          title_size = 15,
          legend_text_size = 11,
          legend_title_size = 12,
          margin = ggplot2::margin(14, 18, 14, 18, unit = "pt")
        )
    )
  }
  plot
}

iobr_empty_feature_result <- function() {
  data.frame(
    sig_names = character(0),
    p.value = numeric(0),
    statistic = numeric(0),
    p.adj = numeric(0),
    Method_Source = character(0),
    stringsAsFactors = FALSE
  )
}

iobr_is_gse_analysis_group <- function(group_col) {
  grepl("^analysis_", as.character(group_col)[1])
}

iobr_group_design_name <- function(group_col) {
  group_col <- as.character(group_col)[1]
  if (iobr_is_gse_analysis_group(group_col)) {
    return(iobr_sanitize(group_col))
  }
  iobr_sanitize(group_col)
}

iobr_group_axis_label <- function(group_col) {
  if (iobr_is_gse_analysis_group(group_col)) {
    return("BULK vs LRC")
  }
  if (identical(as.character(group_col)[1], "Target_Group")) {
    return("Target expression group")
  }
  iobr_sanitize(group_col)
}

iobr_prepare_group_data <- function(dat, group_col, require_two_groups = FALSE) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  if (!group_col %in% colnames(dat)) {
    return(list(data = dat[0, , drop = FALSE], group_col = group_col, levels = character(0)))
  }

  raw_values <- trimws(as.character(dat[[group_col]]))
  keep <- !is.na(raw_values) & raw_values != ""
  dat <- dat[keep, , drop = FALSE]
  raw_values <- raw_values[keep]

  display_col <- paste0(group_col, "__IOBR_Display")
  if (iobr_is_gse_analysis_group(group_col)) {
    display_values <- ifelse(toupper(raw_values) == "BULK", "BULK", "LRC")
    levels <- intersect(c("BULK", "LRC"), unique(display_values))
  } else if (identical(group_col, "Target_Group")) {
    display_values <- raw_values
    levels <- intersect(c("Low", "High"), unique(display_values))
    if (length(levels) == 0L) {
      levels <- sort(unique(display_values))
    }
  } else {
    display_values <- raw_values
    levels <- sort(unique(display_values))
  }

  dat[[display_col]] <- factor(display_values, levels = levels)
  dat <- dat[!is.na(dat[[display_col]]), , drop = FALSE]
  levels <- levels(dat[[display_col]])[levels(dat[[display_col]]) %in% as.character(dat[[display_col]])]

  if (require_two_groups && length(levels) != 2L) {
    return(list(data = dat[0, , drop = FALSE], group_col = display_col, levels = levels))
  }

  list(data = dat, group_col = display_col, levels = levels)
}

iobr_standardize_test_result <- function(dat, feature_col = NULL) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(dat) == 0L) {
    return(iobr_empty_feature_result())
  }

  if (is.null(feature_col)) {
    feature_col <- intersect(
      c("sig_names", "feature", "Feature", "signature", "Signature", "variables", "Variable"),
      colnames(dat)
    )[1]
  }
  if (is.na(feature_col) || is.null(feature_col)) {
    feature_col <- colnames(dat)[1]
  }
  if (!identical(feature_col, "sig_names")) {
    colnames(dat)[match(feature_col, colnames(dat))] <- "sig_names"
  }

  p_col <- intersect(c("p.value", "p.value.Wilcox", "pvalue", "p", "P.Value"), colnames(dat))[1]
  if (is.na(p_col)) {
    dat$p.value <- NA_real_
  } else if (!identical(p_col, "p.value")) {
    colnames(dat)[match(p_col, colnames(dat))] <- "p.value"
  }

  stat_col <- intersect(c("statistic", "estimate", "cor", "Correlation", "rho", "R"), colnames(dat))[1]
  if (is.na(stat_col)) {
    dat$statistic <- NA_real_
  } else if (!identical(stat_col, "statistic")) {
    colnames(dat)[match(stat_col, colnames(dat))] <- "statistic"
  }

  if (!"p.adj" %in% colnames(dat)) {
    padj_col <- intersect(c("p.adj", "p.adjust", "p.adjusted", "FDR", "adj.P.Val"), colnames(dat))[1]
    if (!is.na(padj_col)) {
      colnames(dat)[match(padj_col, colnames(dat))] <- "p.adj"
    } else {
      dat$p.adj <- p.adjust(suppressWarnings(as.numeric(dat$p.value)), method = "BH")
    }
  }

  dat$sig_names <- as.character(dat$sig_names)
  dat$p.value <- suppressWarnings(as.numeric(dat$p.value))
  dat$statistic <- suppressWarnings(as.numeric(dat$statistic))
  dat$p.adj <- suppressWarnings(as.numeric(dat$p.adj))
  dat$p.value.label <- iobr_format_p_value(dat$p.value)
  dat$p.adj.label <- iobr_format_p_value(dat$p.adj)
  dat$Significance <- iobr_significance_stars(dat$p.value)
  dat$Adjusted_Significance <- iobr_significance_stars(dat$p.adj)
  dat[order(dat$p.value, na.last = TRUE), , drop = FALSE]
}

iobr_get_signature_groups <- function(selected_groups = NULL) {
  e <- new.env(parent = emptyenv())
  data(sig_group, package = "IOBR", envir = e)
  groups <- e$sig_group
  if (!is.null(selected_groups) && length(selected_groups) > 0L) {
    selected_groups <- intersect(selected_groups, names(groups))
    groups <- groups[selected_groups]
  }
  groups
}

iobr_get_signature_collection_list <- function(selected_signatures = NULL) {
  e <- new.env(parent = emptyenv())
  data(signature_collection, package = "IOBR", envir = e)
  signatures <- e$signature_collection
  if (!is.null(selected_signatures) && length(selected_signatures) > 0L) {
    selected_signatures <- intersect(selected_signatures, names(signatures))
    signatures <- signatures[selected_signatures]
  }
  signatures
}

iobr_select_signatures <- function(
    selected_groups = NULL,
    selected_signatures = NULL,
    max_signatures = Inf) {
  signatures <- iobr_get_signature_collection_list(selected_signatures)

  if (!is.null(selected_groups) && length(selected_groups) > 0L) {
    groups <- iobr_get_signature_groups(selected_groups)
    wanted <- unique(unlist(groups, use.names = FALSE))
    signatures <- signatures[intersect(wanted, names(signatures))]
  }

  if (is.finite(max_signatures) && length(signatures) > max_signatures) {
    signatures <- signatures[seq_len(max_signatures)]
  }

  signatures
}

iobr_filter_signatures_for_expression <- function(
    signatures,
    expression_matrix,
    mini_gene_count = 3L,
    method = "pca") {
  stopifnot(is.list(signatures))
  if (length(signatures) == 0L) {
    return(signatures)
  }

  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "double"
  gene_sd <- apply(expression_matrix, 1, stats::sd, na.rm = TRUE)
  usable_genes <- rownames(expression_matrix)[is.finite(gene_sd) & gene_sd > 0]
  mini_gene_count <- suppressWarnings(as.integer(mini_gene_count)[1])
  if (!is.finite(mini_gene_count)) {
    mini_gene_count <- 3L
  }
  min_genes <- if (identical(tolower(method), "pca")) {
    max(mini_gene_count, 2L)
  } else {
    max(mini_gene_count, 1L)
  }

  filtered <- lapply(signatures, function(genes) {
    unique(intersect(as.character(genes), usable_genes))
  })
  keep <- vapply(filtered, length, integer(1)) >= min_genes
  filtered[keep]
}

iobr_write_signature_catalog <- function(output_root) {
  signatures <- iobr_get_signature_collection_list()
  groups <- iobr_get_signature_groups()

  signature_table <- data.frame(
    Signature = names(signatures),
    Gene_Count = vapply(signatures, length, integer(1)),
    stringsAsFactors = FALSE
  )
  group_table <- do.call(rbind, lapply(names(groups), function(group_name) {
    data.frame(
      Signature_Group = group_name,
      Signature = groups[[group_name]],
      stringsAsFactors = FALSE
    )
  }))

  iobr_write_module_csv(signature_table, output_root, "capability_catalog", "001_iobr_signature_collection")
  iobr_write_module_csv(group_table, output_root, "capability_catalog", "002_iobr_signature_groups")
}

iobr_prepare_and_cache_local_input <- function(
    project_root,
    cache_root,
    dataset_id,
    se_file,
    target_gene,
    project_id = dataset_id,
    assay_name = "tpm",
    deconvolution_assay = assay_name,
    score_assay = assay_name,
    score_log2 = TRUE,
    sample_detail_values = NULL,
    sample_detail_column = "group_detail",
    sample_barcode_pattern = NULL,
    clinical_file = NULL,
    tumor = TRUE,
    timer_indication = NULL) {
  input <- iobr_prepare_local_se_input(
    dataset_id = dataset_id,
    se_file = se_file,
    target_gene = target_gene,
    project_id = project_id,
    assay_name = assay_name,
    deconvolution_assay = deconvolution_assay,
    score_assay = score_assay,
    score_log2 = score_log2,
    sample_detail_values = sample_detail_values,
    sample_detail_column = sample_detail_column,
    sample_barcode_pattern = sample_barcode_pattern,
    clinical_file = clinical_file,
    tumor = tumor,
    timer_indication = timer_indication
  )
  input$pdata <- iobr_add_survival_columns(input$pdata)

  cache_file <- iobr_make_cache_file(
    cache_root,
    paste(dataset_id, target_gene, "iobr_input", sep = "_")
  )
  iobr_cache_write(input, cache_file)

  data.frame(
    Dataset_ID = dataset_id,
    ProjectID = project_id,
    Target_Gene = target_gene,
    Sample_Count = input$sample_count,
    Gene_Count = input$gene_count,
    SE_File = normalizePath(se_file, winslash = "/", mustWork = TRUE),
    Cache_File = cache_file,
    Tumor = tumor,
    TIMER_Indication = ifelse(is.null(timer_indication), NA_character_, timer_indication),
    stringsAsFactors = FALSE
  )
}

iobr_get_group_columns <- function(input, include_target_group = TRUE, gse_designs = TRUE) {
  pdata <- input$pdata
  group_columns <- character(0)

  if (include_target_group && "Target_Group" %in% colnames(pdata)) {
    group_columns <- c(group_columns, "Target_Group")
  }
  if (gse_designs) {
    design_cols <- grep("^analysis_", colnames(pdata), value = TRUE)
    group_columns <- c(group_columns, design_cols)
  }

  unique(group_columns[vapply(group_columns, function(column) {
    prepared <- iobr_prepare_group_data(
      dat = pdata,
      group_col = column,
      require_two_groups = iobr_is_gse_analysis_group(column)
    )
    length(prepared$levels) >= 2L
  }, logical(1))])
}

iobr_make_group_test_table <- function(feature_table, group_columns, features) {
  if (length(group_columns) == 0L || length(features) == 0L) {
    return(data.frame())
  }

  results <- lapply(group_columns, function(group_col) {
    res <- iobr_standardize_test_result(
      iobr_group_test(feature_table, group_col = group_col, features = features)
    )
    if (nrow(res) == 0L) {
      return(NULL)
    }
    res$Group_Column <- group_col
    res$Group_Display <- iobr_group_axis_label(group_col)
    res$Analysis_Design <- iobr_group_design_name(group_col)
    res
  })

  iobr_bind_rows(results)
}

iobr_make_correlation_table <- function(feature_table, target = "Target_Expression", features, method = "spearman") {
  if (!target %in% colnames(feature_table) || length(features) == 0L) {
    return(data.frame())
  }
  iobr_standardize_test_result(
    iobr_batch_cor(feature_table, target = target, features = features, method = method)
  )
}

iobr_save_feature_summary_plots <- function(
    feature_table,
    output_root,
    module,
    file_prefix,
    group_columns = character(0),
    correlation_table = data.frame(),
    group_test_table = data.frame(),
    top_n = 15L,
    draw_all_feature_heatmap = TRUE,
    all_feature_heatmap_max_features = Inf) {
  features <- unique(c(
    if ("sig_names" %in% colnames(correlation_table)) correlation_table$sig_names else character(0),
    if ("sig_names" %in% colnames(group_test_table)) group_test_table$sig_names else character(0)
  ))
  features <- intersect(features, colnames(feature_table))
  if (length(features) == 0L) {
    features <- iobr_feature_columns(feature_table)
  }
  status <- list()

  try({
  if (nrow(correlation_table) > 0L && "statistic" %in% colnames(correlation_table)) {
    plot_data <- head(
      correlation_table[order(abs(correlation_table$statistic), decreasing = TRUE), , drop = FALSE],
      top_n
    )
    if (nrow(plot_data) > 0L) {
      plot_data$P_Label <- iobr_make_p_label(plot_data$p.value)
      p <- iobr_make_top_barplot(
        plot_data,
        feature_col = "sig_names",
        value_col = "statistic",
        title = paste0(file_prefix, " correlation with target gene"),
        xlab = NULL,
        top_n = top_n,
        label_col = "P_Label",
        ylab = "Correlation statistic"
      )
      files <- iobr_save_module_plot(
        p,
        output_root = output_root,
        module = module,
        file_stem = paste0(file_prefix, "_target_correlation_top", top_n),
        plot_type = "bar",
        n_items = nrow(plot_data),
        max_label_chars = max(nchar(as.character(plot_data$sig_names)), na.rm = TRUE),
        title = paste0(file_prefix, " correlation with target gene")
      )
      status[[length(status) + 1L]] <- files
    }
  }
  }, silent = TRUE)

  try({
	  if (nrow(group_test_table) > 0L) {
	    for (group_col in unique(group_test_table$Group_Column)) {
	      current <- group_test_table[group_test_table$Group_Column == group_col, , drop = FALSE]
	      try(
	        iobr_write_module_csv(
	          current,
	          output_root = output_root,
	          module = module,
	          file_stem = paste0(file_prefix, "_", group_col, "_group_tests")
	        ),
	        silent = TRUE
	      )
	      current$minus_log10_p <- -log10(pmax(current$p.value, .Machine$double.xmin))
	      current$P_Label <- iobr_make_p_label(current$p.value)
      current <- head(current[order(current$p.value, na.last = TRUE), , drop = FALSE], top_n)
      if (nrow(current) == 0L) {
        next
      }

      p <- iobr_make_top_barplot(
        current,
        feature_col = "sig_names",
        value_col = "minus_log10_p",
        title = paste0(file_prefix, " group differences: ", iobr_group_axis_label(group_col)),
        xlab = NULL,
        top_n = top_n,
        label_col = "P_Label",
        ylab = "-log10(p value)"
      )
      iobr_save_module_plot(
        p,
        output_root = output_root,
        module = module,
        file_stem = paste0(file_prefix, "_", group_col, "_group_difference_top", top_n),
        plot_type = "bar",
        n_items = nrow(current),
        max_label_chars = max(nchar(as.character(current$sig_names)), na.rm = TRUE),
        title = paste0(file_prefix, " group differences: ", iobr_group_axis_label(group_col))
      )

      box_features <- intersect(current$sig_names, features)
      if (length(box_features) > 0L) {
        for (feature in head(box_features, min(6L, length(box_features)))) {
          p_box <- try(
            iobr_make_feature_boxplot(
              feature_table,
              feature = feature,
              group_col = group_col,
              title = paste(file_prefix, feature, "by", iobr_group_axis_label(group_col))
            ),
            silent = TRUE
          )
          if (inherits(p_box, "try-error") || !inherits(p_box, "ggplot")) {
            next
          }
          prepared_box <- iobr_prepare_group_data(
            feature_table,
            group_col = group_col,
            require_two_groups = iobr_is_gse_analysis_group(group_col)
          )
          iobr_save_module_plot(
            p_box,
            output_root = output_root,
            module = module,
            file_stem = paste0(file_prefix, "_", group_col, "_", feature, "_boxplot"),
            plot_type = "box",
            n_groups = length(prepared_box$levels),
            max_label_chars = max(nchar(c(prepared_box$levels, feature)), na.rm = TRUE),
            title = paste(file_prefix, feature, "by", iobr_group_axis_label(group_col))
          )
        }
      }
    }
  }
  }, silent = TRUE)

  pca_features <- character(0)
  if (nrow(correlation_table) > 0L) {
    pca_features <- correlation_table$sig_names[order(abs(correlation_table$statistic), decreasing = TRUE)]
  }
  if (length(pca_features) < 2L && nrow(group_test_table) > 0L) {
    pca_features <- c(pca_features, group_test_table$sig_names[order(group_test_table$p.value, na.last = TRUE)])
  }
  pca_features <- unique(intersect(pca_features, features))
  pca_features <- head(pca_features, min(top_n, length(pca_features)))
  pca_features <- pca_features[vapply(pca_features, function(feature) {
    values <- suppressWarnings(as.numeric(feature_table[[feature]]))
    sum(is.finite(values)) >= 3L && stats::sd(values, na.rm = TRUE) > 0
  }, logical(1))]

  try({
  if (length(pca_features) >= 2L && length(group_columns) > 0L) {
    group_col <- group_columns[1]
    prepared_group <- iobr_prepare_group_data(
      feature_table,
      group_col = group_col,
      require_two_groups = iobr_is_gse_analysis_group(group_col)
    )
    feature_table_for_group <- prepared_group$data
    plot_group_col <- prepared_group$group_col
    if (nrow(feature_table_for_group) == 0L || length(prepared_group$levels) < 2L) {
      stop("No valid group data for PCA/heatmap: ", group_col)
    }
    p_pca <- try(
      IOBR::iobr_pca(
        data = t(as.matrix(feature_table_for_group[, pca_features, drop = FALSE])),
        is.matrix = TRUE,
        scale = TRUE,
        is.log = FALSE,
        pdata = feature_table_for_group[, c("ID", plot_group_col), drop = FALSE],
        id_pdata = "ID",
        group = plot_group_col,
        repel = FALSE,
        addEllipses = FALSE
      ),
      silent = TRUE
    )
    if (inherits(p_pca, "try-error") || !inherits(p_pca, "ggplot")) {
      p_pca <- iobr_make_pca_plot(
        feature_table_for_group,
        features = pca_features,
        color_col = plot_group_col,
        title = paste0(file_prefix, " PCA")
      )
    } else {
      p_pca <- iobr_append_title(p_pca, paste0(file_prefix, " PCA"))
    }
    iobr_save_module_plot(
      p_pca,
      output_root = output_root,
      module = module,
      file_stem = paste0(file_prefix, "_", group_col, "_pca"),
      plot_type = "pca",
      n_groups = length(prepared_group$levels),
      max_label_chars = max(nchar(c(as.character(feature_table[[group_col]]), pca_features)), na.rm = TRUE),
      title = paste0(file_prefix, " PCA")
    )

    sample_order <- feature_table_for_group$ID[order(feature_table_for_group[[plot_group_col]], feature_table_for_group$Target_Expression)]
    p_heat <- iobr_make_heatmap_plot(
      feature_table_for_group,
      features = pca_features,
      sample_order = sample_order,
      title = paste0(file_prefix, " top feature heatmap")
    )
    iobr_save_module_plot(
      p_heat,
      output_root = output_root,
      module = module,
      file_stem = paste0(file_prefix, "_", group_col, "_heatmap"),
      plot_type = "heatmap",
      n_items = length(pca_features),
      n_samples = length(sample_order),
      max_label_chars = max(nchar(c(pca_features, sample_order)), na.rm = TRUE),
      title = paste0(file_prefix, " top feature heatmap")
    )
  }
  }, silent = TRUE)

  try({
  if (isTRUE(draw_all_feature_heatmap) && length(features) >= 2L && length(group_columns) > 0L) {
    heat_features <- features
    if (nrow(correlation_table) > 0L && "statistic" %in% colnames(correlation_table)) {
      ordered_features <- correlation_table$sig_names[order(abs(correlation_table$statistic), decreasing = TRUE)]
      heat_features <- unique(c(intersect(ordered_features, features), features))
    }
    if (is.finite(all_feature_heatmap_max_features)) {
      heat_features <- head(heat_features, all_feature_heatmap_max_features)
    }
    heat_features <- heat_features[vapply(heat_features, function(feature) {
      values <- suppressWarnings(as.numeric(feature_table[[feature]]))
      sum(is.finite(values)) >= 3L && stats::sd(values, na.rm = TRUE) > 0
    }, logical(1))]

    if (length(heat_features) >= 2L) {
      group_col <- group_columns[1]
      prepared_group <- iobr_prepare_group_data(
        feature_table,
        group_col = group_col,
        require_two_groups = iobr_is_gse_analysis_group(group_col)
      )
      feature_table_for_group <- prepared_group$data
      plot_group_col <- prepared_group$group_col
      if (nrow(feature_table_for_group) == 0L || length(prepared_group$levels) < 2L) {
        stop("No valid group data for all-feature heatmap: ", group_col)
      }
      sample_order <- feature_table_for_group$ID[order(feature_table_for_group[[plot_group_col]], feature_table_for_group$Target_Expression)]
      p_all_heat <- iobr_make_heatmap_plot(
        feature_table_for_group,
        features = heat_features,
        sample_order = sample_order,
        title = paste0(file_prefix, " all IOBR features heatmap")
      )
      iobr_save_module_plot(
        p_all_heat,
        output_root = output_root,
        module = module,
        file_stem = paste0(file_prefix, "_", group_col, "_all_features_heatmap"),
        plot_type = "all_signature_heatmap",
        n_items = length(heat_features),
        n_samples = length(sample_order),
        max_label_chars = max(nchar(c(heat_features, sample_order)), na.rm = TRUE),
        title = paste0(file_prefix, " all IOBR features heatmap")
      )

      chunk_size <- iobr_heatmap_feature_chunk_size(length(sample_order))
      if (length(heat_features) > chunk_size) {
        chunk_index <- ceiling(seq_along(heat_features) / chunk_size)
        heat_feature_chunks <- split(heat_features, chunk_index)
        for (chunk_id in seq_along(heat_feature_chunks)) {
          chunk_features <- heat_feature_chunks[[chunk_id]]
          chunk_title <- paste0(
            file_prefix,
            " all IOBR features heatmap part ",
            sprintf("%03d", chunk_id),
            "/",
            sprintf("%03d", length(heat_feature_chunks))
          )
          p_chunk_heat <- iobr_make_heatmap_plot(
            feature_table_for_group,
            features = chunk_features,
            sample_order = sample_order,
            title = chunk_title
          )
          iobr_save_module_plot(
            p_chunk_heat,
            output_root = output_root,
            module = module,
            file_stem = paste0(
              file_prefix,
              "_",
              group_col,
              "_all_features_heatmap_part",
              sprintf("%03d", chunk_id)
            ),
            plot_type = "heatmap_chunk",
            n_items = length(chunk_features),
            n_samples = length(sample_order),
            max_label_chars = max(nchar(c(chunk_features, sample_order)), na.rm = TRUE),
            title = chunk_title
          )
        }
      }
    }
  }
  }, silent = TRUE)

  invisible(status)
}

iobr_run_deconvolution_task <- function(
    task,
    output_root,
    correlation_method = "spearman",
    top_n = 15L,
    draw_all_feature_heatmap = TRUE,
    all_feature_heatmap_max_features = Inf) {
  input <- iobr_cache_read(task$Input_File)
  method <- task$Method
  module <- "deconvolution"
  file_prefix <- iobr_sanitize(paste(task$Dataset_ID, task$Target_Gene, method, sep = "_"))

  result <- iobr_safe_task({
    reference <- NULL
    if (method %in% c("svr", "lsei")) {
      reference <- iobr_load_lm22_reference()
    }
    group_list <- NULL
    if (identical(method, "timer") && !is.null(input$timer_indication)) {
      group_list <- rep(input$timer_indication, ncol(input$deconvolution_expr))
    }

    deconv <- suppressWarnings(suppressMessages(IOBR::deconvo_tme(
      eset = input$deconvolution_expr,
      project = input$project_id,
      method = method,
      tumor = input$tumor,
      perm = task$Perm,
      reference = reference,
      group_list = group_list,
      arrays = FALSE,
      plot = FALSE
    )))
    deconv <- iobr_normalize_id_column(deconv)
    feature_table <- merge(input$pdata, deconv, by = "ID", all.x = TRUE)
    features <- iobr_feature_columns(feature_table, exclude = setdiff(colnames(input$pdata), "ID"))
    group_columns <- iobr_get_group_columns(input)
    cor_table <- iobr_make_correlation_table(
      feature_table,
      target = "Target_Expression",
      features = features,
      method = correlation_method
    )
    group_table <- iobr_make_group_test_table(
      feature_table,
      group_columns = group_columns,
      features = features
    )

    feature_file <- iobr_write_module_csv(feature_table, output_root, module, paste0(file_prefix, "_scores"))
    cor_file <- iobr_write_module_csv(cor_table, output_root, module, paste0(file_prefix, "_target_correlation"))
    group_file <- iobr_write_module_csv(group_table, output_root, module, paste0(file_prefix, "_group_tests"))

    iobr_save_feature_summary_plots(
      feature_table = feature_table,
      output_root = output_root,
      module = module,
      file_prefix = file_prefix,
      group_columns = group_columns,
      correlation_table = cor_table,
      group_test_table = group_table,
      top_n = top_n,
      draw_all_feature_heatmap = draw_all_feature_heatmap,
      all_feature_heatmap_max_features = all_feature_heatmap_max_features
    )

    if (length(features) > 1L && nrow(feature_table) > 1L) {
      group_col_for_cell <- if (length(group_columns) > 0L) group_columns[1] else "Target_Group"
      cell_plot <- try(
        iobr_make_cell_fraction_summary_plot(
          feature_table = feature_table,
          features = features,
          group_col = group_col_for_cell,
          title = paste0(file_prefix, " cell fractions by group"),
          top_n = min(length(features), 20L)
        ),
        silent = TRUE
      )
      if (!inherits(cell_plot, "try-error") && inherits(cell_plot, "ggplot")) {
        iobr_save_module_plot(
          cell_plot,
          output_root = output_root,
          module = module,
          file_stem = paste0(file_prefix, "_cell_barplot"),
          plot_type = "cell_bar",
          n_items = min(length(features), 20L),
          n_samples = nrow(feature_table),
          max_label_chars = max(nchar(c(head(features, min(length(features), 20L)), feature_table$ID)), na.rm = TRUE),
          title = paste0(file_prefix, " cell fractions by group")
        )
      }
    }

    iobr_status_row(
      dataset_id = task$Dataset_ID,
      module = module,
      task = paste("deconvo_tme", method, sep = "::"),
      status = "success",
      output_file = feature_file,
      n_rows = nrow(feature_table),
      n_cols = ncol(feature_table)
    )
  }, module = module, task = paste(task$Dataset_ID, method, sep = "_"))

  result
}

iobr_run_signature_score_task <- function(
    task,
    output_root,
    signatures,
    correlation_method = "spearman",
    top_n = 15L,
    draw_all_feature_heatmap = TRUE,
    all_feature_heatmap_max_features = Inf) {
  input <- iobr_cache_read(task$Input_File)
  method <- task$Method
  module <- "signature_score"
  file_prefix <- iobr_sanitize(paste(task$Dataset_ID, task$Target_Gene, method, sep = "_"))

  result <- iobr_safe_task({
    task_signatures <- iobr_filter_signatures_for_expression(
      signatures = signatures,
      expression_matrix = input$score_expr,
      mini_gene_count = task$Mini_Gene_Count,
      method = method
    )
    if (length(task_signatures) == 0L) {
      return(iobr_status_row(
        dataset_id = task$Dataset_ID,
        module = module,
        task = paste("calculate_sig_score", method, sep = "::"),
        status = "skipped",
        message = "No signatures retained after zero-variance gene filtering."
      ))
    }

    score <- suppressWarnings(suppressMessages(IOBR::calculate_sig_score(
      pdata = input$pdata[, c("ID"), drop = FALSE],
      eset = input$score_expr,
      signature = task_signatures,
      method = method,
      mini_gene_count = task$Mini_Gene_Count,
      column_of_sample = "ID",
      print_gene_proportion = FALSE,
      print_filtered_signatures = FALSE,
      adjust_eset = FALSE,
      parallel.size = task$Inner_Workers
    )))
    score <- iobr_normalize_id_column(score)
    feature_table <- merge(input$pdata, score, by = "ID", all.x = TRUE)
    features <- iobr_feature_columns(feature_table, exclude = setdiff(colnames(input$pdata), "ID"))
    group_columns <- iobr_get_group_columns(input)
    cor_table <- iobr_make_correlation_table(
      feature_table,
      target = "Target_Expression",
      features = features,
      method = correlation_method
    )
    group_table <- iobr_make_group_test_table(
      feature_table,
      group_columns = group_columns,
      features = features
    )

    feature_file <- iobr_write_module_csv(feature_table, output_root, module, paste0(file_prefix, "_scores"))
    cor_file <- iobr_write_module_csv(cor_table, output_root, module, paste0(file_prefix, "_target_correlation"))
    group_file <- iobr_write_module_csv(group_table, output_root, module, paste0(file_prefix, "_group_tests"))

    iobr_save_feature_summary_plots(
      feature_table = feature_table,
      output_root = output_root,
      module = module,
      file_prefix = file_prefix,
      group_columns = group_columns,
      correlation_table = cor_table,
      group_test_table = group_table,
      top_n = top_n,
      draw_all_feature_heatmap = draw_all_feature_heatmap,
      all_feature_heatmap_max_features = all_feature_heatmap_max_features
    )

    iobr_status_row(
      dataset_id = task$Dataset_ID,
      module = module,
      task = paste("calculate_sig_score", method, sep = "::"),
      status = "success",
      output_file = feature_file,
      n_rows = nrow(feature_table),
      n_cols = ncol(feature_table)
    )
  }, module = module, task = paste(task$Dataset_ID, method, sep = "_"))

  result
}
