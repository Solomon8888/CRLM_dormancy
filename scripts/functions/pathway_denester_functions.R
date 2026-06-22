# Pathway denesting helpers for GSEA result tables.
#
# This is a project-local R implementation of the PathwayDenester idea:
# enriched terms are checked against more significant overlapping terms to flag
# pathways whose leading-edge genes are mostly explained by a stronger pathway.


# 0. Small helpers -------------------------------------------------------------

split_core_genes <- function(x) {
  x <- as.character(x)
  if (length(x) == 0 || is.na(x[1]) || x[1] == "") {
    return(character(0))
  }

  genes <- unlist(strsplit(x[1], "/", fixed = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  unique(genes[!is.na(genes) & genes != ""])
}

make_gene_set_lookup <- function(term2gene) {
  stopifnot(all(c("term", "gene") %in% colnames(term2gene)))
  term2gene$term <- as.character(term2gene$term)
  term2gene$gene <- as.character(term2gene$gene)
  term2gene <- term2gene[
    !is.na(term2gene$term) & term2gene$term != "" &
      !is.na(term2gene$gene) & term2gene$gene != "",
    ,
    drop = FALSE
  ]

  split(term2gene$gene, term2gene$term)
}

safe_numeric <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(x))
  value[!is.finite(value)] <- default
  value
}

get_first_or_default <- function(x, default = NA_character_) {
  if (length(x) == 0 || is.na(x[1])) {
    return(default)
  }

  as.character(x[1])
}

denester_hypergeom_p <- function(deg_count, deg_in_intersection, intersection_size, term_size) {
  deg_count <- as.integer(deg_count)
  deg_in_intersection <- as.integer(deg_in_intersection)
  intersection_size <- as.integer(intersection_size)
  term_size <- as.integer(term_size)

  if (!is.finite(deg_count) || !is.finite(deg_in_intersection) ||
      !is.finite(intersection_size) || !is.finite(term_size)) {
    return(NA_real_)
  }
  if (term_size <= 0 || deg_count <= 0 || intersection_size <= 0 ||
      deg_in_intersection <= 0) {
    return(1)
  }
  if (deg_count > term_size || intersection_size > term_size) {
    return(NA_real_)
  }

  phyper(
    q = deg_in_intersection - 1,
    m = deg_count,
    n = term_size - deg_count,
    k = intersection_size,
    lower.tail = FALSE
  )
}

normalize_denester_status <- function(x) {
  x <- as.character(x)
  x[x == "keep"] <- "Keep"
  x[x == "exclude"] <- "Exclude"
  x
}


# 1. PathwayDenester-style calculation ----------------------------------------

prepare_denester_input_table <- function(
    result_table,
    gene_sets,
    gsea_p_column = "pvalue",
    gsea_p_cutoff = 0.05,
    max_unexpected_core_fraction = 0.10) {
  required_columns <- c("ID", "Description", "NES", gsea_p_column, "core_enrichment")
  missing_columns <- setdiff(required_columns, colnames(result_table))
  if (length(missing_columns) > 0) {
    stop(
      "Missing columns in GSEA result table: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  result_table$ID <- as.character(result_table$ID)
  result_table$Description <- as.character(result_table$Description)
  result_table$NES <- safe_numeric(result_table$NES)
  result_table[[gsea_p_column]] <- safe_numeric(result_table[[gsea_p_column]])

  keep_index <- !is.na(result_table$ID) &
    result_table$ID != "" &
    !is.na(result_table$Description) &
    result_table$Description != "" &
    is.finite(result_table$NES) &
    is.finite(result_table[[gsea_p_column]]) &
    result_table[[gsea_p_column]] <= gsea_p_cutoff &
    result_table$ID %in% names(gene_sets)

  result_table <- result_table[keep_index, , drop = FALSE]
  if (nrow(result_table) == 0) {
    return(list(pathways = list(), skipped = data.frame()))
  }

  pathway_list <- vector("list", nrow(result_table))
  skipped_list <- list()
  kept_count <- 0L

  for (i in seq_len(nrow(result_table))) {
    term_id <- result_table$ID[i]
    all_genes <- unique(as.character(gene_sets[[term_id]]))
    all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
    core_genes_raw <- split_core_genes(result_table$core_enrichment[i])
    core_genes <- intersect(core_genes_raw, all_genes)

    unexpected_fraction <- if (length(core_genes_raw) > 0) {
      length(setdiff(core_genes_raw, all_genes)) / length(core_genes_raw)
    } else {
      NA_real_
    }

    if (length(all_genes) == 0 || length(core_genes) == 0 ||
        (!is.na(unexpected_fraction) && unexpected_fraction > max_unexpected_core_fraction)) {
      skipped_list[[length(skipped_list) + 1L]] <- data.frame(
        Pathway_ID = term_id,
        Description = result_table$Description[i],
        Skip_Reason = ifelse(
          length(core_genes) == 0,
          "No core enrichment genes matched the reference gene set.",
          "Too many core enrichment genes were absent from the reference gene set."
        ),
        Core_Gene_Count_Raw = length(core_genes_raw),
        Core_Gene_Count_Matched = length(core_genes),
        Term_Size = length(all_genes),
        Unexpected_Core_Gene_Fraction = unexpected_fraction,
        stringsAsFactors = FALSE
      )
      next
    }

    kept_count <- kept_count + 1L
    p_adjust <- if ("p.adjust" %in% colnames(result_table)) {
      safe_numeric(result_table$`p.adjust`[i])
    } else {
      NA_real_
    }
    q_value <- if ("qvalue" %in% colnames(result_table)) {
      safe_numeric(result_table$qvalue[i])
    } else {
      NA_real_
    }

    pathway_list[[kept_count]] <- list(
      id = term_id,
      name = result_table$Description[i],
      p_value = result_table[[gsea_p_column]][i],
      p_adjust = p_adjust,
      q_value = q_value,
      nes = result_table$NES[i],
      all_genes = all_genes,
      core_genes = core_genes,
      raw_core_genes = core_genes_raw,
      density = length(core_genes) / length(all_genes),
      unexpected_fraction = unexpected_fraction,
      result = 1,
      reciprocal = 0,
      filter = "keep",
      vs = "itself",
      vs_name = "",
      intersection_size = 0L,
      core_genes_in_intersection = 0L,
      jaccard_with_vs = 0,
      fraction_core_in_intersection = 0,
      fraction_term_overlap_with_vs = 0
    )
  }

  pathway_list <- pathway_list[seq_len(kept_count)]
  if (length(pathway_list) > 0) {
    p_values <- vapply(pathway_list, `[[`, numeric(1), "p_value")
    density <- vapply(pathway_list, `[[`, numeric(1), "density")
    core_count <- vapply(pathway_list, function(x) length(x$core_genes), integer(1))
    order_index <- order(p_values, -density, -core_count)
    pathway_list <- pathway_list[order_index]
  }

  skipped_table <- if (length(skipped_list) > 0) {
    do.call(rbind, skipped_list)
  } else {
    data.frame()
  }

  list(
    pathways = pathway_list,
    skipped = skipped_table
  )
}

run_pathway_denester_on_table <- function(
    result_table,
    gene_sets,
    gsea_p_column = "pvalue",
    gsea_p_cutoff = 0.05,
    denester_p_cutoff = 0.05,
    to_test_threshold = 0,
    max_unexpected_core_fraction = 0.10) {
  prepared <- prepare_denester_input_table(
    result_table = result_table,
    gene_sets = gene_sets,
    gsea_p_column = gsea_p_column,
    gsea_p_cutoff = gsea_p_cutoff,
    max_unexpected_core_fraction = max_unexpected_core_fraction
  )
  pathways <- prepared$pathways

  if (length(pathways) == 0) {
    return(list(
      result_table = data.frame(),
      skipped_table = prepared$skipped
    ))
  }

  if (length(pathways) >= 2) {
    for (current_line in seq(2, length(pathways))) {
      degs_in_current <- length(pathways[[current_line]]$core_genes)
      size_current <- length(pathways[[current_line]]$all_genes)

      for (test_line in seq_len(current_line - 1L)) {
        if (pathways[[test_line]]$filter != "keep") {
          next
        }

        degs_in_test <- length(pathways[[test_line]]$core_genes)
        intersection_genes <- intersect(
          pathways[[current_line]]$all_genes,
          pathways[[test_line]]$all_genes
        )
        degs_in_intersection_current <- length(intersect(
          pathways[[current_line]]$core_genes,
          intersection_genes
        ))

        if (degs_in_intersection_current <=
            to_test_threshold * min(degs_in_current, degs_in_test)) {
          next
        }

        intersection_size <- length(intersection_genes)
        size_test <- length(pathways[[test_line]]$all_genes)

        current_result <- denester_hypergeom_p(
          deg_count = degs_in_current,
          deg_in_intersection = degs_in_intersection_current,
          intersection_size = intersection_size,
          term_size = size_current
        )
        reverse_result <- denester_hypergeom_p(
          deg_count = degs_in_test,
          deg_in_intersection = degs_in_intersection_current,
          intersection_size = intersection_size,
          term_size = size_test
        )

        union_size <- length(union(
          pathways[[current_line]]$all_genes,
          pathways[[test_line]]$all_genes
        ))
        jaccard <- ifelse(union_size > 0, intersection_size / union_size, 0)
        fraction_core <- degs_in_intersection_current / degs_in_current
        fraction_term <- intersection_size / size_current

        if (is.finite(current_result) &&
            current_result < denester_p_cutoff &&
            (!is.finite(reverse_result) || reverse_result > denester_p_cutoff)) {
          pathways[[current_line]]$result <- current_result
          pathways[[current_line]]$reciprocal <- reverse_result
          pathways[[current_line]]$filter <- "exclude"
          pathways[[current_line]]$vs <- pathways[[test_line]]$id
          pathways[[current_line]]$vs_name <- pathways[[test_line]]$name
          pathways[[current_line]]$intersection_size <- intersection_size
          pathways[[current_line]]$core_genes_in_intersection <- degs_in_intersection_current
          pathways[[current_line]]$jaccard_with_vs <- jaccard
          pathways[[current_line]]$fraction_core_in_intersection <- fraction_core
          pathways[[current_line]]$fraction_term_overlap_with_vs <- fraction_term
          break
        }

        if (is.finite(current_result) &&
            current_result < pathways[[current_line]]$result &&
            pathways[[current_line]]$filter == "keep") {
          pathways[[current_line]]$result <- current_result
          pathways[[current_line]]$reciprocal <- reverse_result
          pathways[[current_line]]$vs <- pathways[[test_line]]$id
          pathways[[current_line]]$vs_name <- pathways[[test_line]]$name
          pathways[[current_line]]$intersection_size <- intersection_size
          pathways[[current_line]]$core_genes_in_intersection <- degs_in_intersection_current
          pathways[[current_line]]$jaccard_with_vs <- jaccard
          pathways[[current_line]]$fraction_core_in_intersection <- fraction_core
          pathways[[current_line]]$fraction_term_overlap_with_vs <- fraction_term
        }
      }
    }
  }

  result_table <- do.call(rbind, lapply(seq_along(pathways), function(i) {
    pathway <- pathways[[i]]
    data.frame(
      Denester_Rank = i,
      Pathway_ID = pathway$id,
      Description = pathway$name,
      NES = pathway$nes,
      GSEA_P_Value = pathway$p_value,
      GSEA_P_Adjust = pathway$p_adjust,
      GSEA_Q_Value = pathway$q_value,
      Core_Gene_Count = length(pathway$core_genes),
      Term_Size = length(pathway$all_genes),
      Core_Density = pathway$density,
      Unexpected_Core_Gene_Fraction = pathway$unexpected_fraction,
      Intersection_Size_With_MSP = pathway$intersection_size,
      Core_Genes_In_MSP_Intersection = pathway$core_genes_in_intersection,
      Fraction_Core_In_MSP_Intersection = pathway$fraction_core_in_intersection,
      Fraction_Term_Overlap_With_MSP = pathway$fraction_term_overlap_with_vs,
      Jaccard_With_MSP = pathway$jaccard_with_vs,
      Denester_P_Value = pathway$result,
      Reciprocal_P_Value = pathway$reciprocal,
      Denester_Status = normalize_denester_status(pathway$filter),
      Versus_Pathway_ID = pathway$vs,
      Versus_Pathway_Name = pathway$vs_name,
      Is_Reciprocal = pathway$reciprocal < denester_p_cutoff &&
        pathway$p_value < denester_p_cutoff,
      Core_Enrichment_Matched = paste(pathway$core_genes, collapse = "/"),
      Core_Enrichment_Raw = paste(pathway$raw_core_genes, collapse = "/"),
      stringsAsFactors = FALSE
    )
  }))
  rownames(result_table) <- NULL

  list(
    result_table = result_table,
    skipped_table = prepared$skipped
  )
}


# 2. Output path helpers -------------------------------------------------------

get_denester_table_dir <- function(result_root, category, analysis_name, geneset_name) {
  file.path(
    result_root,
    "tables",
    "GSEA_denester",
    sanitize_file_name(category),
    sanitize_file_name(analysis_name),
    sanitize_file_name(geneset_name)
  )
}

get_denester_plot_dir <- function(result_root, category, analysis_name, geneset_name) {
  file.path(
    result_root,
    "plots",
    "GSEA_denester",
    sanitize_file_name(category),
    sanitize_file_name(analysis_name),
    sanitize_file_name(geneset_name)
  )
}


# 3. Plotting ------------------------------------------------------------------

select_denester_terms_for_dotplot <- function(denester_table, top_n = 10) {
  if (nrow(denester_table) == 0) {
    return(denester_table)
  }

  keep_table <- denester_table[denester_table$Denester_Status == "Keep", , drop = FALSE]
  exclude_table <- denester_table[denester_table$Denester_Status == "Exclude", , drop = FALSE]

  if (nrow(keep_table) > 0) {
    keep_table <- keep_table[order(keep_table$GSEA_P_Value, -abs(keep_table$NES)), , drop = FALSE]
    keep_table <- keep_table[seq_len(min(top_n, nrow(keep_table))), , drop = FALSE]
  }
  if (nrow(exclude_table) > 0) {
    exclude_table <- exclude_table[order(exclude_table$GSEA_P_Value, -abs(exclude_table$NES)), , drop = FALSE]
    exclude_table <- exclude_table[seq_len(min(top_n, nrow(exclude_table))), , drop = FALSE]
  }

  plot_table <- rbind(keep_table, exclude_table)
  plot_table <- plot_table[order(plot_table$Denester_Status, plot_table$NES), , drop = FALSE]
  rownames(plot_table) <- NULL
  plot_table
}

make_denester_empty_plot <- function(label = "No denester terms to display") {
  ggplot2::ggplot() +
    ggplot2::geom_text(
      ggplot2::aes(x = 0, y = 0, label = label),
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      size = 4.2
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void(base_family = TEXT_FONT_FAMILY)
}

make_denester_dotplot <- function(denester_table, top_n = 10) {
  plot_table <- select_denester_terms_for_dotplot(denester_table, top_n = top_n)
  if (nrow(plot_table) == 0) {
    return(list(
      plot = make_denester_empty_plot(),
      shown_terms = 0L,
      plot_labels = character(0)
    ))
  }

  plot_table$Description_For_Plot <- format_gsea_description_for_plot(plot_table$Description)
  plot_table$Description_For_Plot <- wrap_label(plot_table$Description_For_Plot, width = 45)
  plot_table$Minus_Log10_GSEA_P <- -log10(plot_table$GSEA_P_Value)
  plot_table$Denester_Status <- factor(
    plot_table$Denester_Status,
    levels = c("Keep", "Exclude")
  )
  plot_table$Description_For_Plot <- factor(
    plot_table$Description_For_Plot,
    levels = plot_table$Description_For_Plot
  )

  status_colors <- c(
    Keep = UP_COLOR,
    Exclude = "#8F8F8F"
  )

  plot <- ggplot2::ggplot(
    plot_table,
    ggplot2::aes(x = NES, y = Description_For_Plot)
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey72",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = Core_Gene_Count, fill = Denester_Status),
      shape = 21,
      color = TEXT_COLOR,
      stroke = 0.35,
      alpha = 0.90
    ) +
    ggplot2::scale_fill_manual(values = status_colors, name = NULL) +
    ggplot2::scale_size(range = c(3.8, 9.2), name = "Core genes") +
    ggplot2::labs(
      x = "Normalized enrichment score",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.22),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
      axis.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.title = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      text = ggplot2::element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        family = TEXT_FONT_FAMILY
      ),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )

  list(
    plot = plot,
    shown_terms = nrow(plot_table),
    plot_labels = as.character(plot_table$Description_For_Plot)
  )
}

get_denester_dotplot_size <- function(shown_terms, plot_labels) {
  if (exists("get_gsea_dotplot_size", mode = "function")) {
    return(get_gsea_dotplot_size(shown_terms = shown_terms, plot_labels = plot_labels))
  }

  list(
    width = max(7.2, 5.2 + max(nchar(plot_labels), 12) * 0.045),
    height = max(5.2, shown_terms * 0.35 + 2)
  )
}

make_denester_overlap_heatmap <- function(denester_table, gene_sets, top_n = 30) {
  if (nrow(denester_table) < 2) {
    return(list(
      plot = make_denester_empty_plot("Not enough terms for overlap heatmap"),
      width = 6.2,
      height = 5.2
    ))
  }

  plot_table <- denester_table[
    order(denester_table$GSEA_P_Value, -abs(denester_table$NES)),
    ,
    drop = FALSE
  ]
  plot_table <- plot_table[seq_len(min(top_n, nrow(plot_table))), , drop = FALSE]
  term_ids <- plot_table$Pathway_ID

  overlap_matrix <- matrix(0, nrow = length(term_ids), ncol = length(term_ids))
  rownames(overlap_matrix) <- term_ids
  colnames(overlap_matrix) <- term_ids

  for (i in seq_along(term_ids)) {
    genes_i <- unique(gene_sets[[term_ids[i]]])
    for (j in seq_along(term_ids)) {
      genes_j <- unique(gene_sets[[term_ids[j]]])
      overlap_matrix[i, j] <- ifelse(
        length(genes_i) > 0,
        length(intersect(genes_i, genes_j)) / length(genes_i),
        0
      )
    }
  }

  if (length(term_ids) > 2) {
    distance_matrix <- as.dist(1 - pmax(overlap_matrix, t(overlap_matrix)))
    order_index <- tryCatch(
      hclust(distance_matrix, method = "average")$order,
      error = function(e) seq_along(term_ids)
    )
  } else {
    order_index <- seq_along(term_ids)
  }

  plot_table <- plot_table[order_index, , drop = FALSE]
  term_ids <- plot_table$Pathway_ID
  overlap_matrix <- overlap_matrix[term_ids, term_ids, drop = FALSE]

  display_labels <- format_gsea_description_for_plot(plot_table$Description)
  display_labels <- wrap_label(display_labels, width = 34)
  names(display_labels) <- term_ids

  heatmap_data <- expand.grid(
    Row_ID = term_ids,
    Column_ID = term_ids,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  heatmap_data$Overlap <- mapply(
    function(row_id, col_id) overlap_matrix[row_id, col_id],
    heatmap_data$Row_ID,
    heatmap_data$Column_ID
  )
  heatmap_data$Row_Label <- factor(
    display_labels[heatmap_data$Row_ID],
    levels = rev(display_labels[term_ids])
  )
  heatmap_data$Column_Index <- match(heatmap_data$Column_ID, term_ids)

  status_data <- data.frame(
    Row_Label = factor(display_labels[term_ids], levels = rev(display_labels[term_ids])),
    Status_X = 0,
    Denester_Status = factor(plot_table$Denester_Status, levels = c("Keep", "Exclude")),
    stringsAsFactors = FALSE
  )

  status_colors <- c(
    Keep = UP_COLOR,
    Exclude = "#8F8F8F"
  )

  heatmap_plot <- ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(x = Column_Index, y = Row_Label, fill = Overlap)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.12) +
    ggplot2::geom_point(
      data = status_data,
      ggplot2::aes(x = Status_X, y = Row_Label, color = Denester_Status),
      inherit.aes = FALSE,
      shape = 15,
      size = 3.8
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = DOWN_COLOR,
      limits = c(0, 1),
      name = "Overlap"
    ) +
    ggplot2::scale_color_manual(values = status_colors, name = NULL) +
    ggplot2::scale_x_continuous(
      breaks = c(0, seq_along(term_ids)),
      labels = c("Status", seq_along(term_ids)),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
      axis.text.x = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE, angle = 0),
      axis.text.y = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE, size = 8.5),
      legend.title = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      text = ggplot2::element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        family = TEXT_FONT_FAMILY
      ),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )

  list(
    plot = heatmap_plot,
    width = max(7.2, length(term_ids) * 0.32 + 4.6),
    height = max(5.4, length(term_ids) * 0.30 + 2.2)
  )
}

make_denester_summary_barplot <- function(summary_table) {
  if (nrow(summary_table) == 0) {
    return(make_denester_empty_plot("No denester summary to display"))
  }

  plot_table <- summary_table
  plot_table$Dataset_Group <- paste(plot_table$Dataset, plot_table$Plot_Category, sep = " / ")
  plot_table$Excluded_Percent <- as.numeric(plot_table$Excluded_Percent)

  ggplot2::ggplot(
    plot_table,
    ggplot2::aes(x = Dataset_Group, y = Excluded_Percent, fill = Plot_Category)
  ) +
    ggplot2::geom_boxplot(
      width = 0.62,
      outlier.shape = NA,
      alpha = 0.72,
      color = TEXT_COLOR
    ) +
    ggplot2::geom_jitter(
      width = 0.16,
      height = 0,
      size = 1.8,
      alpha = 0.62,
      color = TEXT_COLOR
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Excluded pathways (%)",
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
      axis.text.x = ggplot2::element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        angle = 35,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.position = "right",
      legend.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE)
    )
}


# 4. Batch runner --------------------------------------------------------------

run_denester_for_gsea_file <- function(
    gsea_info_row,
    gene_sets,
    result_root,
    gsea_p_column,
    gsea_p_cutoff,
    denester_p_cutoff,
    to_test_threshold,
    max_unexpected_core_fraction,
    dotplot_top_n,
    overlap_heatmap_top_n,
    draw_overlap_heatmap = TRUE) {
  result_table <- read_result_csv(gsea_info_row$GSEA_Result_File)
  denester_result <- run_pathway_denester_on_table(
    result_table = result_table,
    gene_sets = gene_sets,
    gsea_p_column = gsea_p_column,
    gsea_p_cutoff = gsea_p_cutoff,
    denester_p_cutoff = denester_p_cutoff,
    to_test_threshold = to_test_threshold,
    max_unexpected_core_fraction = max_unexpected_core_fraction
  )

  denester_table <- denester_result$result_table
  skipped_table <- denester_result$skipped_table

  if (nrow(denester_table) > 0) {
    denester_table$Dataset <- gsea_info_row$Dataset
    denester_table$Plot_Category <- gsea_info_row$Plot_Category
    denester_table$Analysis_Name <- gsea_info_row$Analysis_Name
    denester_table$GeneSet_Name <- gsea_info_row$GeneSet_Name
    leading_columns <- c(
      "Dataset", "Plot_Category", "Analysis_Name", "GeneSet_Name"
    )
    denester_table <- denester_table[
      ,
      c(leading_columns, setdiff(colnames(denester_table), leading_columns)),
      drop = FALSE
    ]
  }

  table_dir <- get_denester_table_dir(
    result_root = result_root,
    category = gsea_info_row$Plot_Category,
    analysis_name = gsea_info_row$Analysis_Name,
    geneset_name = gsea_info_row$GeneSet_Name
  )
  plot_dir <- get_denester_plot_dir(
    result_root = result_root,
    category = gsea_info_row$Plot_Category,
    analysis_name = gsea_info_row$Analysis_Name,
    geneset_name = gsea_info_row$GeneSet_Name
  )
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  denester_csv <- write_csv_with_report_previews(
    denester_table,
    file.path(table_dir, "denester_result.csv"),
    n_rows = 21,
    na = "NA"
  )
  kept_csv <- write_csv_with_report_previews(
    denester_table[denester_table$Denester_Status == "Keep", , drop = FALSE],
    file.path(table_dir, "denester_kept_terms.csv"),
    n_rows = 21,
    na = "NA"
  )
  excluded_csv <- write_csv_with_report_previews(
    denester_table[denester_table$Denester_Status == "Exclude", , drop = FALSE],
    file.path(table_dir, "denester_excluded_terms.csv"),
    n_rows = 21,
    na = "NA"
  )

  skipped_csv <- NA_character_
  if (nrow(skipped_table) > 0) {
    skipped_csv <- write_csv_with_report_previews(
      skipped_table,
      file.path(table_dir, "denester_skipped_terms.csv"),
      n_rows = 21,
      na = "NA"
    )
  }

  dotplot_result <- make_denester_dotplot(
    denester_table = denester_table,
    top_n = dotplot_top_n
  )
  dotplot_size <- get_denester_dotplot_size(
    shown_terms = dotplot_result$shown_terms,
    plot_labels = dotplot_result$plot_labels
  )
  dotplot_files <- save_ggplot_pdf_png(
    plot = dotplot_result$plot,
    pdf_file = file.path(plot_dir, "denester_dotplot.pdf"),
    width = dotplot_size$width,
    height = dotplot_size$height
  )

  heatmap_pdf <- NA_character_
  heatmap_png <- NA_character_
  if (isTRUE(draw_overlap_heatmap)) {
    heatmap_result <- make_denester_overlap_heatmap(
      denester_table = denester_table,
      gene_sets = gene_sets,
      top_n = overlap_heatmap_top_n
    )
    heatmap_files <- save_ggplot_pdf_png(
      plot = heatmap_result$plot,
      pdf_file = file.path(plot_dir, "denester_overlap_heatmap.pdf"),
      width = heatmap_result$width,
      height = heatmap_result$height
    )
    heatmap_pdf <- heatmap_files$pdf_file
    heatmap_png <- heatmap_files$png_file
  }

  significant_terms <- nrow(denester_table)
  kept_terms <- sum(denester_table$Denester_Status == "Keep", na.rm = TRUE)
  excluded_terms <- sum(denester_table$Denester_Status == "Exclude", na.rm = TRUE)

  data.frame(
    Dataset = gsea_info_row$Dataset,
    Plot_Category = gsea_info_row$Plot_Category,
    Analysis_Name = gsea_info_row$Analysis_Name,
    GeneSet_Name = gsea_info_row$GeneSet_Name,
    Significant_Terms = significant_terms,
    Kept_Terms = kept_terms,
    Excluded_Terms = excluded_terms,
    Excluded_Percent = ifelse(significant_terms > 0, excluded_terms / significant_terms * 100, 0),
    Reciprocal_Terms = sum(denester_table$Is_Reciprocal, na.rm = TRUE),
    Skipped_Terms = nrow(skipped_table),
    GSEA_Result_File = resolve_table_file(gsea_info_row$GSEA_Result_File),
    Denester_Result_File = denester_csv,
    Kept_Terms_File = kept_csv,
    Excluded_Terms_File = excluded_csv,
    Skipped_Terms_File = skipped_csv,
    Dotplot_PDF_File = dotplot_files$pdf_file,
    Dotplot_PNG_File = dotplot_files$png_file,
    Heatmap_PDF_File = heatmap_pdf,
    Heatmap_PNG_File = heatmap_png,
    stringsAsFactors = FALSE
  )
}
