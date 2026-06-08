# 本地TCGA/GTEx quickanalysis共用函数
#
# 服务于quickanalysis 03/04/05/06号脚本：
# - 读取本项目已经准备好的TCGA或GTEx SummarizedExperiment对象；
# - 按脚本配置统一筛选样本；
# - 统一筛选protein-coding基因；
# - 复用NGS流程中的GSEA计算、绘图、表格保存、并行进度条和耗时统计函数。


# 0. 命令行覆盖配置的小工具 ---------------------------------------------------

parse_env_vector <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  parts <- unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

parse_env_logical <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(default)
  }

  value %in% c("1", "true", "t", "yes", "y")
}

parse_env_integer <- function(name, default) {
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

qa_sanitize_file_name <- function(x, default = "analysis") {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- default
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[x == ""] <- default
  x
}


# 1. 项目路径与依赖 -----------------------------------------------------------

qa_get_current_script_file <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(normalizePath(frames[[i]]$ofile, winslash = "/", mustWork = TRUE))
    }
  }

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }

  NA_character_
}

qa_find_project_root <- function(script_file = NA_character_) {
  start_points <- unique(c(
    getwd(),
    if (!is.na(script_file)) dirname(script_file) else character(0)
  ))

  for (start_point in start_points) {
    current <- normalizePath(start_point, winslash = "/", mustWork = TRUE)
    repeat {
      marker <- file.path(current, "scripts", "functions", "limma_de_functions.R")
      if (file.exists(marker) && dir.exists(file.path(current, "data"))) {
        return(current)
      }

      parent <- dirname(current)
      if (identical(parent, current)) {
        break
      }
      current <- parent
    }
  }

  stop("Cannot locate project root. Please run this script from the project directory.")
}

qa_require_packages <- function(packages) {
  is_available <- vapply(packages, function(package_name) {
    suppressWarnings(
      suppressPackageStartupMessages(
        requireNamespace(package_name, quietly = TRUE)
      )
    )
  }, logical(1))

  missing_packages <- packages[!is_available]
  if (length(missing_packages) > 0L) {
    stop(
      "Please install required R packages before running this script: ",
      paste(missing_packages, collapse = ", ")
    )
  }

  invisible(TRUE)
}

qa_source_project_functions <- function(project_root) {
  source(file.path(project_root, "scripts", "functions", "limma_de_functions.R"))
  source(file.path(project_root, "scripts", "functions", "plotting_common_functions.R"))
  source(file.path(project_root, "scripts", "functions", "report_table_functions.R"))
  source(file.path(project_root, "scripts", "functions", "result_table_io_functions.R"))
  source(file.path(project_root, "scripts", "functions", "parallel_runtime_functions.R"))
  source(file.path(project_root, "scripts", "functions", "gsea_common_functions.R"))
  invisible(TRUE)
}


# 2. 本地SE读取与筛选 ---------------------------------------------------------

qa_apply_se_file_template <- function(
    project_root,
    se_file_template,
    data_source,
    dataset_id) {
  se_file <- se_file_template
  replacements <- c(
    data_source = data_source,
    dataset = dataset_id,
    cancer = dataset_id,
    tissue = dataset_id
  )

  for (key in names(replacements)) {
    se_file <- gsub(
      paste0("{", key, "}"),
      replacements[[key]],
      se_file,
      fixed = TRUE
    )
  }

  if (!grepl("^/", se_file)) {
    se_file <- file.path(project_root, se_file)
  }

  se_file
}

qa_get_local_se_file <- function(
    project_root,
    data_source,
    dataset_id,
    se_file_template = NULL) {
  data_source <- trimws(as.character(data_source))
  dataset_id <- trimws(as.character(dataset_id))

  if (is.null(se_file_template)) {
    return(file.path(
      project_root,
      "data",
      data_source,
      dataset_id,
      "data_prepare",
      paste0(dataset_id, "_se_raw.rds")
    ))
  }

  qa_apply_se_file_template(
    project_root = project_root,
    se_file_template = se_file_template,
    data_source = data_source,
    dataset_id = dataset_id
  )
}

qa_get_tcga_se_file <- function(project_root, cancer, se_file_template = NULL) {
  qa_get_local_se_file(
    project_root = project_root,
    data_source = "TCGA",
    dataset_id = cancer,
    se_file_template = se_file_template
  )
}

qa_get_assay_matrix <- function(se, assay_name) {
  assay_names <- names(SummarizedExperiment::assays(se))
  if (!assay_name %in% assay_names) {
    stop(
      "Assay '",
      assay_name,
      "' was not found. Available assays: ",
      paste(assay_names, collapse = ", ")
    )
  }

  as.matrix(SummarizedExperiment::assay(se, assay_name))
}

qa_plain_data_frame <- function(x) {
  # S4Vectors::DataFrame的as.data.frame方法会忽略stringsAsFactors/check.names等参数；
  # 这里先做基础转换，再统一把factor列转成字符，避免批量运行时产生无意义warning。
  dat <- as.data.frame(x)
  dat <- data.frame(dat, check.names = FALSE, stringsAsFactors = FALSE)

  for (column_name in colnames(dat)) {
    if (is.factor(dat[[column_name]])) {
      dat[[column_name]] <- as.character(dat[[column_name]])
    }
  }

  dat
}

qa_prepare_expression_matrix <- function(se, assay_name, log2_transform = TRUE) {
  expr <- qa_get_assay_matrix(se, assay_name)
  if (log2_transform) {
    if (min(expr, na.rm = TRUE) < 0) {
      stop("Cannot apply log2(x + 1): assay contains negative values.")
    }
    expr <- log2(expr + 1)
  }

  expr
}

qa_make_sample_id <- function(se, sample_info) {
  if (!"Sample_ID" %in% colnames(sample_info)) {
    sample_info$Sample_ID <- colnames(se)
  }

  sample_info$Sample_ID <- trimws(as.character(sample_info$Sample_ID))
  empty_sample_id <- is.na(sample_info$Sample_ID) | sample_info$Sample_ID == ""
  sample_info$Sample_ID[empty_sample_id] <- colnames(se)[empty_sample_id]
  sample_info
}

qa_select_local_samples <- function(
    se,
    sample_filter_column = NULL,
    sample_filter_values = NULL,
    sample_filter_regex = NULL,
    sample_filter_label = "all samples") {
  sample_info <- qa_plain_data_frame(SummarizedExperiment::colData(se))
  sample_info <- qa_make_sample_id(se, sample_info)

  if (!any(duplicated(sample_info$Sample_ID))) {
    rownames(sample_info) <- sample_info$Sample_ID
  }

  keep_sample <- rep(TRUE, nrow(sample_info))
  filter_source <- sample_filter_label

  has_column_filter <- !is.null(sample_filter_column) &&
    length(sample_filter_column) == 1L &&
    nzchar(sample_filter_column) &&
    sample_filter_column %in% colnames(sample_info) &&
    !is.null(sample_filter_values) &&
    length(sample_filter_values) > 0L

  if (has_column_filter) {
    keep_sample <- trimws(as.character(sample_info[[sample_filter_column]])) %in%
      as.character(sample_filter_values)
    filter_source <- paste0(
      sample_filter_column,
      " in ",
      paste(sample_filter_values, collapse = ", ")
    )
  } else if (!is.null(sample_filter_regex) &&
             length(sample_filter_regex) == 1L &&
             nzchar(sample_filter_regex)) {
    keep_sample <- grepl(sample_filter_regex, sample_info$Sample_ID)
    filter_source <- paste0("Sample_ID pattern ", sample_filter_regex)
  }

  if (!any(keep_sample)) {
    stop("No samples were selected by ", filter_source, ".")
  }

  sample_info <- sample_info[keep_sample, , drop = FALSE]
  se <- se[, sample_info$Sample_ID, drop = FALSE]

  stopifnot(all(sample_info$Sample_ID == colnames(se)))

  list(
    se = se,
    sample_info = sample_info,
    filter_source = filter_source,
    selected_sample_count = ncol(se)
  )
}

qa_select_01a_samples <- function(
    se,
    sample_detail_filter = "tumor_primary_01a",
    sample_barcode_pattern = "-01A$") {
  qa_select_local_samples(
    se = se,
    sample_filter_column = "group_detail",
    sample_filter_values = sample_detail_filter,
    sample_filter_regex = sample_barcode_pattern,
    sample_filter_label = "TCGA 01A primary tumor samples"
  )
}

qa_load_local_se_inputs <- function(
    project_root,
    data_source,
    dataset_id,
    se_file_template = NULL,
    count_assay_name = "counts",
    expression_assay_name = "tpm",
    expression_log2_transform = TRUE,
    gene_biotype_filter = "coding",
    sample_filter_column = NULL,
    sample_filter_values = NULL,
    sample_filter_regex = NULL,
    sample_filter_label = "all samples") {
  se_file <- qa_get_local_se_file(
    project_root = project_root,
    data_source = data_source,
    dataset_id = dataset_id,
    se_file_template = se_file_template
  )
  if (!file.exists(se_file)) {
    stop(data_source, " SE file does not exist: ", se_file)
  }

  se <- readRDS(se_file)
  stopifnot(inherits(se, "SummarizedExperiment"))

  sample_selection <- qa_select_local_samples(
    se = se,
    sample_filter_column = sample_filter_column,
    sample_filter_values = sample_filter_values,
    sample_filter_regex = sample_filter_regex,
    sample_filter_label = sample_filter_label
  )
  se <- sample_selection$se
  sample_info <- sample_selection$sample_info

  counts_all <- qa_get_assay_matrix(se, count_assay_name)
  expression_all <- qa_prepare_expression_matrix(
    se = se,
    assay_name = expression_assay_name,
    log2_transform = expression_log2_transform
  )

  feature_id <- rownames(counts_all)
  if (is.null(feature_id)) {
    feature_id <- paste0("Feature_", seq_len(nrow(counts_all)))
    rownames(counts_all) <- feature_id
    rownames(expression_all) <- feature_id
  }

  gene_annotation <- data.frame(
    Feature_ID = feature_id,
    qa_plain_data_frame(SummarizedExperiment::rowData(se)),
    check.names = FALSE
  )
  rownames(gene_annotation) <- rownames(counts_all)

  gene_filter <- filter_genes_by_biotype(
    exprSet = counts_all,
    gene_annotation = gene_annotation,
    biotype_filter = gene_biotype_filter
  )

  counts_coding <- gene_filter$exprSet
  gene_annotation_coding <- gene_filter$gene_annotation
  expression_coding <- expression_all[rownames(counts_coding), , drop = FALSE]

  list(
    data_source = data_source,
    dataset_id = dataset_id,
    cancer = dataset_id,
    se_file = se_file,
    se = se,
    sample_info = sample_info,
    counts = counts_coding,
    expression = expression_coding,
    gene_annotation = gene_annotation_coding,
    sample_filter_source = sample_selection$filter_source,
    selected_sample_count = sample_selection$selected_sample_count,
    original_gene_count = nrow(counts_all),
    selected_gene_count = gene_filter$selected_gene_count,
    gene_biotype_filter = gene_filter$filter,
    biotype_column = gene_filter$biotype_column
  )
}

qa_load_tcga_01a_inputs <- function(
    project_root,
    cancer,
    se_file_template = NULL,
    count_assay_name = "counts",
    expression_assay_name = "tpm",
    expression_log2_transform = TRUE,
    gene_biotype_filter = "coding",
    sample_detail_filter = "tumor_primary_01a",
    sample_barcode_pattern = "-01A$") {
  qa_load_local_se_inputs(
    project_root = project_root,
    data_source = "TCGA",
    dataset_id = cancer,
    se_file_template = se_file_template,
    count_assay_name = count_assay_name,
    expression_assay_name = expression_assay_name,
    expression_log2_transform = expression_log2_transform,
    gene_biotype_filter = gene_biotype_filter,
    sample_filter_column = "group_detail",
    sample_filter_values = sample_detail_filter,
    sample_filter_regex = sample_barcode_pattern,
    sample_filter_label = "TCGA 01A primary tumor samples"
  )
}

qa_make_gene_result_slug <- function(target_genes, suffix) {
  gene_slug <- qa_sanitize_file_name(paste(unique(target_genes), collapse = "_"))
  paste0(gene_slug, "-", suffix)
}


# 3. 基因定位、分组和输出整理 --------------------------------------------------

qa_strip_ensembl_version <- function(x) {
  sub("[.][0-9]+(_PAR_Y)?$", "\\1", as.character(x))
}

qa_find_target_feature <- function(gene_annotation, expression_matrix, target_gene) {
  target_gene <- trimws(as.character(target_gene))
  if (!nzchar(target_gene)) {
    stop("TARGET_GENE cannot be empty.")
  }

  candidates <- rep(FALSE, nrow(gene_annotation))

  if ("Symbol" %in% colnames(gene_annotation)) {
    candidates <- candidates |
      toupper(trimws(as.character(gene_annotation$Symbol))) == toupper(target_gene)
  }

  if ("GeneID" %in% colnames(gene_annotation)) {
    candidates <- candidates |
      trimws(as.character(gene_annotation$GeneID)) == target_gene |
      qa_strip_ensembl_version(gene_annotation$GeneID) == qa_strip_ensembl_version(target_gene)
  }

  if ("Ensembl" %in% colnames(gene_annotation)) {
    candidates <- candidates |
      trimws(as.character(gene_annotation$Ensembl)) == target_gene |
      qa_strip_ensembl_version(gene_annotation$Ensembl) == qa_strip_ensembl_version(target_gene)
  }

  if ("Entrez" %in% colnames(gene_annotation)) {
    candidates <- candidates | trimws(as.character(gene_annotation$Entrez)) == target_gene
  }

  candidate_index <- which(candidates)
  if (length(candidate_index) == 0L) {
    stop("Target gene was not found in protein-coding gene annotation: ", target_gene)
  }

  candidate_expression <- expression_matrix[candidate_index, , drop = FALSE]
  mean_expression <- rowMeans(candidate_expression, na.rm = TRUE)
  candidate_symbol <- if ("Symbol" %in% colnames(gene_annotation)) {
    as.character(gene_annotation$Symbol[candidate_index])
  } else {
    rownames(gene_annotation)[candidate_index]
  }

  selected_index <- candidate_index[order(-mean_expression, candidate_symbol)[1]]
  selected_row <- gene_annotation[selected_index, , drop = FALSE]

  list(
    feature_index = selected_index,
    feature_id = rownames(gene_annotation)[selected_index],
    annotation = selected_row,
    candidate_count = length(candidate_index),
    mean_expression = mean_expression[match(selected_index, candidate_index)]
  )
}

qa_make_analysis_design <- function(
    sample_info,
    analysis_name,
    group_values,
    experiment_group,
    group_column_name = NULL) {
  if (is.null(group_column_name)) {
    group_column_name <- paste0("analysis_", sanitize_file_name(analysis_name))
  }

  sample_info[[group_column_name]] <- as.character(group_values)
  column_index <- match(group_column_name, colnames(sample_info))

  analysis_designs <- data.frame(
    Analysis_Order = 1L,
    Column_Index = column_index,
    Column_Name = group_column_name,
    Analysis_Base_Name = sanitize_file_name(analysis_name),
    Duplicate_Order = 1L,
    Analysis_Name = sanitize_file_name(analysis_name),
    Experiment_Group = experiment_group,
    stringsAsFactors = FALSE
  )

  list(
    sample_info = sample_info,
    analysis_designs = analysis_designs
  )
}

qa_prepare_ranked_output_table <- function(dat, drop_columns = c("Feature_ID", "Biotype", "Length")) {
  keep_columns <- setdiff(colnames(dat), drop_columns)
  dat <- dat[, keep_columns, drop = FALSE]

  if ("Entrez" %in% colnames(dat)) {
    if (is.numeric(dat$Entrez)) {
      dat$Entrez <- ifelse(
        is.na(dat$Entrez),
        "",
        format(dat$Entrez, scientific = FALSE, trim = TRUE)
      )
    } else {
      dat$Entrez <- as.character(dat$Entrez)
    }
  }

  dat
}

qa_count_regulation <- function(dat, p_value_column, p_value_cutoff, effect_column, effect_cutoff) {
  stopifnot(p_value_column %in% colnames(dat))
  stopifnot(effect_column %in% colnames(dat))

  effect <- as.numeric(dat[[effect_column]])
  p_value <- as.numeric(dat[[p_value_column]])

  significant <- is.finite(effect) &
    is.finite(p_value) &
    abs(effect) > effect_cutoff &
    p_value < p_value_cutoff

  up <- significant & effect > effect_cutoff
  down <- significant & effect < -effect_cutoff

  list(
    significant = significant,
    up = up,
    down = down
  )
}


# 4. GSEA计算和绘图 -----------------------------------------------------------

qa_make_gsea_geneset_cache <- function() {
  msigdb_catalog <- build_msigdb_geneset_catalog()
  runtime_genesets <- get_runtime_genesets_to_run()
  geneset_config <- select_msigdb_genesets(
    catalog = msigdb_catalog,
    genesets_to_run = runtime_genesets
  )

  cat("\nLoading MSigDB gene sets...\n")
  geneset_cache <- lapply(names(geneset_config), function(geneset_name) {
    config <- geneset_config[[geneset_name]]
    terms <- load_msigdb_terms(geneset_name, config)

    list(
      config = config,
      term2gene = terms$term2gene,
      term2name = terms$term2name,
      cache_source = terms$Cache_Source
    )
  })
  names(geneset_cache) <- names(geneset_config)

  geneset_summary <- do.call(
    rbind,
    lapply(names(geneset_cache), function(geneset_name) {
      cache <- geneset_cache[[geneset_name]]
      data.frame(
        GeneSet_Name = geneset_name,
        Output_Name = cache$config$output_name,
        Terms = length(unique(cache$term2gene$term)),
        Term_Gene_Links = nrow(cache$term2gene),
        Source = cache$cache_source,
        stringsAsFactors = FALSE
      )
    })
  )
  print(geneset_summary, row.names = FALSE)

  geneset_cache
}

qa_prepare_gsea_analysis_cache <- function(
    ranked_file_info,
    expression_input_list) {
  analysis_cache <- lapply(seq_len(nrow(ranked_file_info)), function(i) {
    analysis_name <- ranked_file_info$Analysis_Name[i]
    ranked_file <- ranked_file_info$All_Genes_File[i]
    ranked_table <- read.csv(
      ranked_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    gene_list <- prepare_gene_list(
      deg_table = ranked_table,
      gene_id_type = GENE_ID_TYPE,
      rank_metric_column = RANK_METRIC_COLUMN
    )

    expression_data <- NULL
    if (isTRUE(DRAW_SINGLE_PATHWAY_GSEA)) {
      expression_input <- expression_input_list[[analysis_name]]
      if (is.null(expression_input)) {
        stop("Missing expression input for GSEA plotting: ", analysis_name)
      }
      expression_data <- prepare_single_pathway_expression_table(
        analysis_name = analysis_name,
        expression_matrix_all = expression_input$expression_matrix,
        gene_annotation_all = expression_input$gene_annotation,
        analysis_designs = expression_input$analysis_designs,
        sample_info_all = expression_input$sample_info
      )
    }

    list(
      deg_file = ranked_file,
      gene_list = gene_list,
      expression_data = expression_data
    )
  })

  names(analysis_cache) <- ranked_file_info$Analysis_Name
  analysis_cache
}

qa_clean_gsea_outputs <- function(table_output_root, plot_root, analysis_names) {
  unlink(file.path(table_output_root, analysis_names, "GSEA"), recursive = TRUE, force = TRUE)
  unlink(file.path(plot_root, sanitize_file_name(analysis_names)), recursive = TRUE, force = TRUE)
  unlink(file.path(table_output_root, "GSEA_summary"), recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

qa_run_gsea_compute_and_plot <- function(
    ranked_file_info,
    expression_input_list,
    table_output_root,
    plot_root,
    parallel_workers,
    clean_outputs = TRUE,
    run_plots = TRUE) {
  if (nrow(ranked_file_info) == 0L) {
    stop("No ranked gene files were supplied for GSEA.")
  }

  stopifnot(all(c("Analysis_Name", "All_Genes_File") %in% colnames(ranked_file_info)))
  stopifnot(all(file.exists(ranked_file_info$All_Genes_File)))

  if (clean_outputs) {
    qa_clean_gsea_outputs(
      table_output_root = table_output_root,
      plot_root = plot_root,
      analysis_names = ranked_file_info$Analysis_Name
    )
  }

  geneset_cache <- qa_make_gsea_geneset_cache()
  analysis_cache <- qa_prepare_gsea_analysis_cache(
    ranked_file_info = ranked_file_info,
    expression_input_list = expression_input_list
  )

  task_table <- expand.grid(
    Analysis_Name = ranked_file_info$Analysis_Name,
    GeneSet_Name = names(geneset_cache),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  total_tasks <- nrow(task_table)

  parallel_strategy <- setup_parallel_strategy(
    total_tasks = total_tasks,
    max_workers = parallel_workers,
    inner_label = "GSEA nproc per task",
    nested_label = "Single-pathway workers"
  )
  gsea_task_workers <- parallel_strategy$task_workers
  assign("GSEA_INNER_NPROC", parallel_strategy$inner_workers, envir = .GlobalEnv)
  assign("SINGLE_PATHWAY_PLOT_WORKERS", parallel_strategy$nested_workers, envir = .GlobalEnv)

  run_one_gsea_task <- function(task_id) {
    analysis_name <- task_table$Analysis_Name[task_id]
    geneset_name <- task_table$GeneSet_Name[task_id]
    analysis_input <- analysis_cache[[analysis_name]]
    cache <- geneset_cache[[geneset_name]]
    output_name <- cache$config$output_name
    output_dir_name <- sanitize_file_name(output_name)

    table_output_dir <- file.path(
      table_output_root,
      analysis_name,
      "GSEA",
      output_dir_name
    )
    dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

    gsea_run <- load_or_run_gsea(
      analysis_name = analysis_name,
      deg_file = analysis_input$deg_file,
      geneset_name = geneset_name,
      config = cache$config,
      gene_list = analysis_input$gene_list,
      term2gene = cache$term2gene,
      term2name = cache$term2name
    )

    gsea_result <- gsea_run$result
    csv_file <- file.path(table_output_dir, "gsea_result.csv")
    result_table <- write_gsea_result_tables(gsea_result, csv_file)
    csv_file <- resolve_report_csv_file(csv_file)

    plot_files <- list(pdf_file = "", png_file = "")
    single_pathway_count <- 0L

    if (run_plots) {
      plot_output_dir <- file.path(
        plot_root,
        sanitize_file_name(analysis_name),
        output_dir_name
      )
      dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)

      if (is_gsea_result_object(gsea_result)) {
        plot_gsea_result <- prepare_gsea_result_for_plot(gsea_result)
        plot_result_table <- as.data.frame(plot_gsea_result)
        dotplot_result <- make_gsea_dotplot(
          gsea_result = plot_gsea_result,
          result_table = plot_result_table,
          analysis_name = analysis_name,
          geneset_name = geneset_name
        )
      } else {
        dotplot_result <- list(
          plot = make_empty_gsea_plot(),
          shown_terms = 0L,
          plot_labels = character(0)
        )
      }

      plot_size <- get_gsea_dotplot_size(
        shown_terms = dotplot_result$shown_terms,
        plot_labels = dotplot_result$plot_labels
      )
      plot_files <- with_gsea_warnings_suppressed(
        save_ggplot_pdf_png(
          plot = dotplot_result$plot,
          pdf_file = file.path(plot_output_dir, "dotplot.pdf"),
          width = plot_size$width,
          height = plot_size$height
        )
      )

      if (DRAW_SINGLE_PATHWAY_GSEA && is_gsea_result_object(gsea_result)) {
        single_pathway_count <- with_gsea_warnings_suppressed(
          save_single_pathway_gsea_plots(
            gsea_result = gsea_result,
            result_table = result_table,
            analysis_name = analysis_name,
            geneset_name = geneset_name,
            plot_output_dir = plot_output_dir,
            expression_data = analysis_input$expression_data
          )
        )
      }
    }

    data.frame(
      Analysis_Name = analysis_name,
      GeneSet_Name = geneset_name,
      Source = gsea_run$source,
      Ranked_Genes = length(analysis_input$gene_list),
      GSEA_Terms = nrow(result_table),
      Positive_NES = count_nes_direction(result_table, "positive"),
      Negative_NES = count_nes_direction(result_table, "negative"),
      Single_Pathway_Plots = single_pathway_count,
      CSV_File = csv_file,
      PDF_File = plot_files$pdf_file,
      PNG_File = plot_files$png_file,
      stringsAsFactors = FALSE
    )
  }

  cat("\nRunning GSEA compute and plotting tasks...\n")
  task_ids <- seq_len(total_tasks)
  summary_records <- run_parallel_tasks_with_progress(
    task_ids = task_ids,
    task_function = run_one_gsea_task,
    workers = gsea_task_workers,
    progress_label = "GSEA"
  )
  stop_on_parallel_errors(summary_records, task_ids = task_ids, label = "GSEA tasks")

  summary_table <- do.call(rbind, summary_records)
  rownames(summary_table) <- NULL

  summary_output_dir <- file.path(table_output_root, "GSEA_summary")
  dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_csv_file <- write_csv_with_report_previews(
    summary_table,
    file.path(summary_output_dir, "summary.csv"),
    n_rows = 21
  )

  cat("\nGSEA summary:\n")
  print(
    summary_table[
      ,
      c(
        "Analysis_Name", "GeneSet_Name", "Source", "Ranked_Genes",
        "GSEA_Terms", "Positive_NES", "Negative_NES",
        "Single_Pathway_Plots"
      )
    ],
    row.names = FALSE
  )

  list(
    summary = summary_table,
    summary_csv_file = summary_csv_file
  )
}
