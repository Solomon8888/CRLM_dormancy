# Batch visualization helpers for GEO NGS analysis scripts.
#
# These helpers keep the plotting style aligned with the GSE114012 volcano and
# GSEA scripts, while adding support for nested ATF3_function result folders.


# 0. Runtime config helpers ----------------------------------------------------

parse_env_vector <- function(variable_name, default_value) {
  value <- trimws(Sys.getenv(variable_name, unset = ""))
  if (value == "") {
    return(default_value)
  }

  unique(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

parse_env_logical <- function(variable_name, default_value) {
  value <- tolower(trimws(Sys.getenv(variable_name, unset = "")))
  if (value == "") {
    return(default_value)
  }

  value %in% c("1", "true", "t", "yes", "y")
}

resolve_table_file <- function(file) {
  if (exists("resolve_report_csv_file", mode = "function")) {
    return(resolve_report_csv_file(file))
  }

  file
}

read_result_csv <- function(file) {
  read.csv(
    resolve_table_file(file),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

get_relative_path_parts <- function(file, root) {
  root_norm <- normalizePath(root, mustWork = FALSE)
  file_norm <- normalizePath(file, mustWork = FALSE)
  root_prefix <- paste0(root_norm, .Platform$file.sep)

  relative_path <- if (startsWith(file_norm, root_prefix)) {
    substring(file_norm, nchar(root_prefix) + 1L)
  } else {
    file
  }

  strsplit(relative_path, .Platform$file.sep, fixed = TRUE)[[1]]
}


# 1. DEG file discovery --------------------------------------------------------

is_deg_all_genes_file <- function(file) {
  if (basename(file) != "all_genes.csv") {
    return(FALSE)
  }

  parent_dir <- basename(dirname(file))
  grandparent_dir <- basename(dirname(dirname(file)))

  parent_dir == "DEG" || (parent_dir == "csv" && grandparent_dir == "DEG")
}

get_deg_plot_metadata <- function(file, table_root) {
  parts <- get_relative_path_parts(file, table_root)

  if (length(parts) >= 3 &&
      parts[1] != "ATF3_function" &&
      parts[2] == "DEG") {
    return(list(
      plot_category = "Main_DE",
      analysis_name = parts[1],
      display_name = parts[1]
    ))
  }

  if (length(parts) >= 5 &&
      parts[1] == "ATF3_function" &&
      parts[3] == "ATF3_high_low_DE" &&
      parts[4] == "DEG") {
    return(list(
      plot_category = "ATF3_high_low_DE",
      analysis_name = parts[2],
      display_name = parts[2]
    ))
  }

  NULL
}

collect_batch_deg_file_info <- function(table_root) {
  all_gene_files <- list.files(
    table_root,
    pattern = "^all_genes[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  all_gene_files <- all_gene_files[
    vapply(all_gene_files, is_deg_all_genes_file, logical(1))
  ]

  if (length(all_gene_files) == 0) {
    return(data.frame())
  }

  get_file_key <- function(file) {
    metadata <- get_deg_plot_metadata(file, table_root)
    if (is.null(metadata)) {
      return(NA_character_)
    }

    paste(metadata$plot_category, metadata$analysis_name, sep = "::")
  }

  file_keys <- vapply(all_gene_files, get_file_key, character(1))
  all_gene_files <- all_gene_files[!is.na(file_keys)]

  if (length(all_gene_files) == 0) {
    return(data.frame())
  }

  if (exists("prefer_report_csv_files", mode = "function")) {
    all_gene_files <- prefer_report_csv_files(all_gene_files, get_file_key)
  }

  info_list <- lapply(all_gene_files, function(file) {
    metadata <- get_deg_plot_metadata(file, table_root)
    data.frame(
      Plot_Category = metadata$plot_category,
      Analysis_Name = metadata$analysis_name,
      Display_Name = metadata$display_name,
      Plot_Key = paste(metadata$plot_category, metadata$analysis_name, sep = "::"),
      All_Genes_File = file,
      stringsAsFactors = FALSE
    )
  })

  file_info <- do.call(rbind, info_list)
  file_info <- file_info[order(file_info$Plot_Category, file_info$Analysis_Name), , drop = FALSE]
  rownames(file_info) <- NULL

  duplicated_keys <- unique(file_info$Plot_Key[duplicated(file_info$Plot_Key)])
  if (length(duplicated_keys) > 0) {
    stop(
      "More than one all_genes.csv file was found for: ",
      paste(duplicated_keys, collapse = ", ")
    )
  }

  file_info
}

filter_deg_file_info <- function(file_info, analyses_to_plot) {
  if (nrow(file_info) == 0) {
    return(file_info)
  }

  if (identical(analyses_to_plot, "all") ||
      (length(analyses_to_plot) == 1 && tolower(analyses_to_plot) == "all")) {
    return(file_info)
  }

  wanted <- trimws(as.character(analyses_to_plot))
  matched <- file_info$Plot_Key %in% wanted |
    file_info$Analysis_Name %in% wanted |
    file_info$Display_Name %in% wanted

  missing <- setdiff(
    wanted,
    unique(c(
      file_info$Plot_Key[matched],
      file_info$Analysis_Name[matched],
      file_info$Display_Name[matched]
    ))
  )
  missing <- missing[!missing %in% wanted[matched]]

  if (!any(matched)) {
    stop("No DEG all_genes.csv file was selected.")
  }

  if (length(setdiff(wanted, unique(c(
    file_info$Plot_Key,
    file_info$Analysis_Name,
    file_info$Display_Name
  )))) > 0) {
    stop(
      "No DEG all_genes.csv file was found for: ",
      paste(setdiff(wanted, unique(c(
        file_info$Plot_Key,
        file_info$Analysis_Name,
        file_info$Display_Name
      ))), collapse = ", ")
    )
  }

  file_info[matched, , drop = FALSE]
}


# 2. Single volcano plots ------------------------------------------------------

get_batch_volcano_axis_limits <- function(plot_data) {
  x_abs <- max(abs(plot_data$logFC), LOGFC_CUTOFF, na.rm = TRUE)
  x_limit <- ceiling(x_abs * 1.08 * 2) / 2
  x_limit <- max(x_limit, 1)

  y_threshold <- -log10(P_VALUE_CUTOFF)
  y_limit <- ceiling(max(plot_data$Neg_Log10_P, y_threshold, na.rm = TRUE) * 1.08)
  y_limit <- max(y_limit, 2)

  list(
    x = c(-x_limit, x_limit),
    y = c(0, y_limit)
  )
}

get_batch_volcano_pdf_size <- function(axis_limits) {
  y_span <- diff(axis_limits$y)
  extra_height <- min(max((y_span - 6) * 0.04, 0), MAX_EXTRA_PDF_HEIGHT)
  pdf_height <- BASE_PDF_HEIGHT + extra_height

  pdf_width <- pdf_height + LEGEND_WIDTH_INCH + RIGHT_LEGEND_GAP_INCH
  pdf_width <- min(pdf_width, pdf_height * MAX_PDF_WIDTH_HEIGHT_RATIO)
  pdf_width <- max(pdf_width, pdf_height * 1.08)

  list(
    width = pdf_width,
    height = pdf_height
  )
}

make_batch_volcano_plot <- function(plot_data, axis_limits) {
  label_data <- get_volcano_label_data(
    plot_data = plot_data,
    custom_label_genes = CUSTOM_LABEL_GENES,
    symbol_column = TOP_GENE_SYMBOL_COLUMN,
    match_columns = CUSTOM_LABEL_MATCH_COLUMNS,
    p_value_column = P_VALUE_COLUMN,
    top_up_n = TOP_UP_LABEL_N,
    top_down_n = TOP_DOWN_LABEL_N
  )
  label_colors <- get_regulation_label_colors(
    label_data = label_data,
    up_color = UP_COLOR,
    down_color = DOWN_COLOR,
    darken_fraction = TOP_GENE_LABEL_COLOR_DARKEN
  )

  use_gg_repel <- exists("USE_GG_REPEL", mode = "logical") && isTRUE(USE_GG_REPEL)

  volcano_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = logFC, y = Neg_Log10_P, color = Regulation)
  ) +
    ggplot2::geom_point(
      size = POINT_SIZE,
      alpha = POINT_ALPHA,
      shape = 16,
      stroke = 0
    ) +
    ggplot2::geom_vline(
      xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF),
      linewidth = THRESHOLD_LINE_WIDTH,
      linetype = THRESHOLD_LINE_TYPE,
      color = THRESHOLD_LINE_COLOR
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(P_VALUE_CUTOFF),
      linewidth = THRESHOLD_LINE_WIDTH,
      linetype = THRESHOLD_LINE_TYPE,
      color = THRESHOLD_LINE_COLOR
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Not significant" = NOT_SIGNIFICANT_COLOR,
        "Down" = DOWN_COLOR,
        "Up" = UP_COLOR
      ),
      breaks = c("Up", "Down", "Not significant"),
      labels = c("Sig_Up", "Sig_Down", "Not_Sig")
    ) +
    ggplot2::scale_x_continuous(
      limits = axis_limits$x,
      breaks = pretty(axis_limits$x, n = 7),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = axis_limits$y,
      breaks = pretty(axis_limits$y, n = 6),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      x = "log2 fold change",
      y = paste0("-log10(", P_VALUE_COLUMN, ")"),
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
      axis.line = ggplot2::element_line(color = TEXT_COLOR, linewidth = AXIS_LINE_WIDTH),
      axis.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      aspect.ratio = PANEL_HEIGHT_WIDTH_RATIO,
      legend.position = "right",
      legend.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.key = ggplot2::element_blank(),
      legend.key.height = grid::unit(5.5, "mm"),
      legend.key.width = grid::unit(5.5, "mm"),
      legend.box.spacing = grid::unit(8, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 4, unit = "pt"),
      strip.text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      text = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      plot.margin = ggplot2::margin(10, 12, 10, 10, unit = "pt")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = POINT_SIZE * 1.15, alpha = 0.85)
      )
    )

  add_volcano_gene_label_layer(
    plot = volcano_plot,
    label_data = label_data,
    label_colors = label_colors,
    use_gg_repel = use_gg_repel,
    text_family = TEXT_FONT_FAMILY,
    fontface = TOP_GENE_LABEL_FONT_FACE,
    font_size = TOP_GENE_LABEL_FONT_SIZE,
    box_padding = TOP_GENE_LABEL_BOX_PADDING,
    point_padding = TOP_GENE_LABEL_POINT_PADDING,
    segment_width = TOP_GENE_LABEL_SEGMENT_WIDTH,
    force = TOP_GENE_LABEL_FORCE,
    force_pull = TOP_GENE_LABEL_FORCE_PULL,
    max_overlaps = TOP_GENE_LABEL_MAX_OVERLAPS,
    fallback_vjust = -0.8
  )
}

run_batch_volcano_plots <- function(file_info, plot_root) {
  if (nrow(file_info) == 0) {
    return(data.frame())
  }

  summary_list <- lapply(seq_len(nrow(file_info)), function(i) {
    current_info <- file_info[i, , drop = FALSE]
    dat <- read_result_csv(current_info$All_Genes_File)

    plot_data <- prepare_volcano_data(
      dat = dat,
      analysis_name = current_info$Display_Name,
      p_value_column = P_VALUE_COLUMN,
      p_value_cutoff = P_VALUE_CUTOFF,
      logfc_cutoff = LOGFC_CUTOFF,
      ns_label = "Not significant",
      regulation_levels = c("Not significant", "Down", "Up")
    )
    axis_limits <- get_batch_volcano_axis_limits(plot_data)
    pdf_size <- get_batch_volcano_pdf_size(axis_limits)
    volcano_plot <- make_batch_volcano_plot(plot_data, axis_limits)

    output_dir <- file.path(
      plot_root,
      sanitize_file_name(current_info$Plot_Category),
      sanitize_file_name(current_info$Analysis_Name)
    )
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    if (isTRUE(CLEAN_VOLCANO_OUTPUT_DIR)) {
      unlink(list.files(output_dir, pattern = "[.](pdf|png)$", full.names = TRUE))
    }

    output_files <- save_ggplot_pdf_png(
      plot = volcano_plot,
      pdf_file = file.path(output_dir, "volcano_plot.pdf"),
      width = pdf_size$width,
      height = pdf_size$height
    )

    status_counts <- table(plot_data$Regulation)
    data.frame(
      Plot_Category = current_info$Plot_Category,
      Analysis_Name = current_info$Analysis_Name,
      Display_Name = current_info$Display_Name,
      Genes_Plotted = nrow(plot_data),
      Up = count_status(status_counts, "Up"),
      Down = count_status(status_counts, "Down"),
      Not_Significant = count_status(status_counts, "Not significant"),
      X_Min = axis_limits$x[1],
      X_Max = axis_limits$x[2],
      Y_Max = axis_limits$y[2],
      PDF_Width = round(pdf_size$width, 2),
      PDF_Height = round(pdf_size$height, 2),
      PDF_File = output_files$pdf_file,
      PNG_File = output_files$png_file,
      stringsAsFactors = FALSE
    )
  })

  summary_table <- do.call(rbind, summary_list)
  rownames(summary_table) <- NULL
  summary_table
}


# 3. Multiple volcano plots ----------------------------------------------------

get_batch_group_label_colors <- function(selected_analyses) {
  group_colors <- GROUP_LABEL_COLORS[selected_analyses]
  missing_index <- is.na(group_colors)
  duplicated_index <- duplicated(group_colors) & !is.na(group_colors)

  if (any(missing_index | duplicated_index)) {
    replacement_colors <- grDevices::hcl.colors(
      sum(missing_index | duplicated_index),
      palette = "Dark 3"
    )
    group_colors[missing_index | duplicated_index] <- replacement_colors
  }

  if (any(duplicated(group_colors))) {
    group_colors <- grDevices::hcl.colors(
      length(selected_analyses),
      palette = "Dark 3"
    )
  }

  names(group_colors) <- selected_analyses
  group_colors
}

get_batch_group_label_box_y_limits <- function() {
  box_half_height <- max(
    LOGFC_CUTOFF - GROUP_LABEL_BOX_LOGFC_GAP,
    LOGFC_CUTOFF * GROUP_LABEL_BOX_MIN_FRACTION
  )
  box_half_height <- min(box_half_height, LOGFC_CUTOFF)

  c(-box_half_height, box_half_height)
}

get_batch_group_label_text_style <- function(max_line_count) {
  line_extra <- max(max_line_count - 1, 0)

  list(
    font_size = max(
      GROUP_LABEL_FONT_SIZE - line_extra * GROUP_LABEL_FONT_LINE_SHRINK,
      GROUP_LABEL_FONT_MIN_SIZE
    ),
    line_height = max(
      GROUP_LABEL_LINE_HEIGHT - line_extra * GROUP_LABEL_LINE_HEIGHT_SHRINK,
      GROUP_LABEL_LINE_HEIGHT_MIN
    )
  )
}

get_batch_group_layout_data <- function(plot_data, selected_analyses) {
  group_labels <- wrap_label_by_underscore(
    selected_analyses,
    width = GROUP_LABEL_WRAP_WIDTH
  )
  names(group_labels) <- selected_analyses

  threshold_x <- -log10(P_VALUE_CUTOFF)
  max_line_count <- max(vapply(
    group_labels,
    get_label_line_count,
    integer(1)
  ))

  label_style <- get_batch_group_label_text_style(max_line_count)
  label_box_y_limits <- get_batch_group_label_box_y_limits()
  box_height <- diff(label_box_y_limits)

  do.call(rbind, lapply(selected_analyses, function(analysis_name) {
    dat <- plot_data[
      plot_data$Analysis_Name == analysis_name &
        plot_data$Regulation %in% c("Up", "Down"),
      ,
      drop = FALSE
    ]

    if (nrow(dat) == 0) {
      dat <- plot_data[plot_data$Analysis_Name == analysis_name, , drop = FALSE]
    }

    x_max <- max(
      dat$Neg_Log10_P,
      threshold_x,
      na.rm = TRUE
    )

    x_range <- max(x_max - threshold_x, 0.5)
    x_padding <- max(
      x_range * X_AXIS_PADDING_FRACTION,
      X_AXIS_PADDING_MIN
    )
    x_limit_min <- max(
      threshold_x * X_AXIS_LEFT_MIN_FRACTION,
      threshold_x - x_padding
    )
    x_limit_max <- x_max + x_padding
    box_x_margin <- (x_limit_max - x_limit_min) *
      GROUP_LABEL_BOX_X_MARGIN_FRACTION

    data.frame(
      Analysis_Name = analysis_name,
      Label = group_labels[analysis_name],
      Label_X = (x_limit_min + box_x_margin + x_limit_max - box_x_margin) / 2,
      Label_Y = 0,
      Box_X_Min = x_limit_min + box_x_margin,
      Box_X_Max = x_limit_max - box_x_margin,
      Box_Y_Min = label_box_y_limits[1],
      Box_Y_Max = label_box_y_limits[2],
      X_Min = x_limit_min,
      X_Max = x_limit_max,
      X_Data_Min = threshold_x,
      X_Data_Max = x_max,
      X_Padding = x_padding,
      Box_X_Margin = box_x_margin,
      Label_Box_Height = box_height,
      Label_Font_Size = label_style$font_size,
      Label_Line_Height = label_style$line_height,
      stringsAsFactors = FALSE
    )
  }))
}

get_batch_multiple_axis_info <- function(plot_data) {
  plotted_data <- plot_data[
    plot_data$Regulation %in% c("Up", "Down"),
    ,
    drop = FALSE
  ]

  true_y_min <- min(plotted_data$logFC, -LOGFC_CUTOFF, na.rm = TRUE)
  true_y_max <- max(plotted_data$logFC, LOGFC_CUTOFF, na.rm = TRUE)
  true_y_span <- max(true_y_max - true_y_min, 1)
  true_y_padding <- true_y_span * Y_AXIS_PADDING_FRACTION
  true_limits <- c(true_y_min - true_y_padding, true_y_max + true_y_padding)

  true_breaks <- pretty(true_limits, n = 7)
  true_breaks <- true_breaks[
    true_breaks >= true_limits[1] &
      true_breaks <= true_limits[2]
  ]
  true_breaks <- unique(c(true_breaks, -LOGFC_CUTOFF, 0, LOGFC_CUTOFF))
  true_breaks <- true_breaks[
    true_breaks >= true_limits[1] &
      true_breaks <= true_limits[2]
  ]
  true_breaks <- sort(unique(round(true_breaks, 6)))

  list(
    true_limits = true_limits,
    display_limits = true_limits,
    breaks = true_breaks,
    labels = format(true_breaks, trim = TRUE, scientific = FALSE)
  )
}

get_batch_multiple_pdf_size <- function(group_count, display_y_limits, group_layout) {
  y_span <- diff(display_y_limits)
  label_box_height <- max(group_layout$Label_Box_Height, na.rm = TRUE)
  label_height_extra <- max(
    label_box_height - LOGFC_CUTOFF,
    0
  ) * LABEL_BOX_HEIGHT_PDF_SCALE

  pdf_height <- BASE_MULTIPLE_PDF_HEIGHT +
    min(max((y_span - 5) * 0.08, 0), 1.2) +
    label_height_extra
  pdf_height <- min(pdf_height, MAX_MULTIPLE_PDF_HEIGHT)

  pdf_width <- GROUP_WIDTH_INCH * group_count + MULTIPLE_LEGEND_WIDTH_INCH + 1.2
  pdf_width <- max(pdf_width, MIN_MULTIPLE_PDF_WIDTH)
  pdf_width <- max(pdf_width, pdf_height * 1.12)
  pdf_width <- min(pdf_width, MAX_MULTIPLE_PDF_WIDTH)

  list(
    width = pdf_width,
    height = pdf_height
  )
}

make_batch_multiple_volcano_plot <- function(plot_data, group_layout, selected_analyses, axis_info) {
  plot_data$Analysis_Name <- factor(plot_data$Analysis_Name, levels = selected_analyses)
  group_layout$Analysis_Name <- factor(group_layout$Analysis_Name, levels = selected_analyses)

  point_data <- plot_data[
    plot_data$Regulation %in% c("Up", "Down"),
    ,
    drop = FALSE
  ]
  stopifnot(nrow(point_data) > 0)

  top_label_data <- get_volcano_label_data(
    plot_data = plot_data,
    custom_label_genes = CUSTOM_LABEL_GENES,
    symbol_column = TOP_GENE_SYMBOL_COLUMN,
    match_columns = CUSTOM_LABEL_MATCH_COLUMNS,
    p_value_column = P_VALUE_COLUMN,
    top_up_n = TOP_UP_LABEL_N,
    top_down_n = TOP_DOWN_LABEL_N
  )
  if (nrow(top_label_data) > 0) {
    top_label_data$Analysis_Name <- factor(top_label_data$Analysis_Name, levels = selected_analyses)
    top_label_data$Label_Nudge_Y <- ifelse(
      top_label_data$Regulation == "Up",
      TOP_GENE_LABEL_NUDGE_Y,
      -TOP_GENE_LABEL_NUDGE_Y
    )
  }
  top_label_colors <- get_regulation_label_colors(
    label_data = top_label_data,
    up_color = UP_COLOR,
    down_color = DOWN_COLOR,
    darken_fraction = TOP_GENE_LABEL_COLOR_DARKEN
  )

  x_axis_anchor_data <- data.frame(
    Analysis_Name = factor(
      rep(as.character(group_layout$Analysis_Name), each = 2),
      levels = selected_analyses
    ),
    Neg_Log10_P = c(as.vector(rbind(group_layout$X_Min, group_layout$X_Max))),
    Plot_LogFC = 0,
    stringsAsFactors = FALSE
  )

  group_label_colors <- get_batch_group_label_colors(selected_analyses)
  group_label_text_colors <- darken_color(
    group_label_colors,
    fraction = GROUP_LABEL_TEXT_DARKEN
  )
  group_label_border_colors <- darken_color(
    group_label_colors,
    fraction = GROUP_LABEL_BORDER_DARKEN
  )

  point_colors <- c(
    "Down" = DOWN_COLOR,
    "Up" = UP_COLOR
  )

  use_gg_repel <- exists("USE_GG_REPEL", mode = "logical") && isTRUE(USE_GG_REPEL)

  volcano_plot <- ggplot2::ggplot(
    point_data,
    ggplot2::aes(x = Neg_Log10_P, y = Plot_LogFC, color = Regulation)
  ) +
    ggplot2::geom_point(
      size = POINT_SIZE,
      alpha = POINT_ALPHA,
      shape = 16,
      stroke = 0
    ) +
    ggplot2::geom_blank(
      data = x_axis_anchor_data,
      ggplot2::aes(x = Neg_Log10_P, y = Plot_LogFC),
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(Analysis_Name),
      scales = "free_x",
      space = "fixed"
    ) +
    ggplot2::scale_color_manual(
      values = point_colors,
      breaks = c("Up", "Down"),
      labels = c("Sig_Up", "Sig_Down")
    ) +
    ggplot2::scale_x_continuous(
      breaks = function(x) {
        x_breaks <- pretty(x, n = 4)
        x_breaks[x_breaks >= 0]
      },
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = axis_info$breaks,
      labels = axis_info$labels,
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::coord_cartesian(
      ylim = axis_info$display_limits,
      clip = "off"
    ) +
    ggplot2::labs(
      x = paste0("-log10(", P_VALUE_COLUMN, ")"),
      y = "log2 fold change",
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.22),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
      panel.spacing.x = grid::unit(PANEL_SPACING_X_MM, "mm"),
      axis.line = ggplot2::element_line(color = TEXT_COLOR, linewidth = AXIS_LINE_WIDTH),
      axis.text.x = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.text.y = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = ggplot2::element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.position = "right",
      legend.justification = "top",
      legend.box.just = "top",
      legend.text = ggplot2::element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        size = LEGEND_TEXT_SIZE
      ),
      legend.key = ggplot2::element_blank(),
      legend.key.height = grid::unit(LEGEND_KEY_SIZE_MM, "mm"),
      legend.key.width = grid::unit(LEGEND_KEY_SIZE_MM, "mm"),
      legend.box.spacing = grid::unit(8, "pt"),
      legend.margin = ggplot2::margin(LEGEND_TOP_MARGIN_PT, 0, 0, 4, unit = "pt"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      text = ggplot2::element_text(
        color = TEXT_COLOR,
        face = TEXT_FONT_FACE,
        family = TEXT_FONT_FAMILY
      ),
      plot.margin = ggplot2::margin(10, 12, 10, 10, unit = "pt")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          size = POINT_SIZE * LEGEND_POINT_SIZE_SCALE,
          alpha = POINT_ALPHA
        )
      )
    )

  for (analysis_name in selected_analyses) {
    current_label <- group_layout[
      group_layout$Analysis_Name == analysis_name,
      ,
      drop = FALSE
    ]

    volcano_plot <- volcano_plot +
      ggplot2::geom_rect(
        data = current_label,
        ggplot2::aes(
          xmin = Box_X_Min,
          xmax = Box_X_Max,
          ymin = Box_Y_Min,
          ymax = Box_Y_Max
        ),
        inherit.aes = FALSE,
        color = group_label_border_colors[analysis_name],
        fill = grDevices::adjustcolor(
          group_label_colors[analysis_name],
          alpha.f = GROUP_LABEL_ALPHA
        ),
        linewidth = GROUP_LABEL_BORDER_WIDTH
      ) +
      ggplot2::geom_text(
        data = current_label,
        ggplot2::aes(x = Label_X, y = Label_Y, label = Label),
        inherit.aes = FALSE,
        family = TEXT_FONT_FAMILY,
        fontface = TEXT_FONT_FACE,
        size = current_label$Label_Font_Size[1],
        lineheight = current_label$Label_Line_Height[1],
        color = group_label_text_colors[analysis_name]
      )
  }

  add_volcano_gene_label_layer(
    plot = volcano_plot,
    label_data = top_label_data,
    label_colors = top_label_colors,
    use_gg_repel = use_gg_repel,
    text_family = TEXT_FONT_FAMILY,
    fontface = TOP_GENE_LABEL_FONT_FACE,
    font_size = TOP_GENE_LABEL_FONT_SIZE,
    box_padding = TOP_GENE_LABEL_BOX_PADDING,
    point_padding = TOP_GENE_LABEL_POINT_PADDING,
    segment_width = TOP_GENE_LABEL_SEGMENT_WIDTH,
    force = TOP_GENE_LABEL_FORCE,
    force_pull = TOP_GENE_LABEL_FORCE_PULL,
    max_overlaps = TOP_GENE_LABEL_MAX_OVERLAPS,
    nudge_y = top_label_data$Label_Nudge_Y,
    fallback_vjust = ifelse(top_label_data$Regulation == "Up", -0.7, 1.2)
  )
}

run_batch_multiple_volcano_scheme <- function(scheme_name, selected_keys, file_info, plot_root) {
  missing_keys <- setdiff(selected_keys, file_info$Plot_Key)
  if (length(missing_keys) > 0) {
    stop(
      "No DEG all_genes.csv file was found for scheme ",
      scheme_name,
      ": ",
      paste(missing_keys, collapse = ", ")
    )
  }

  selected_info <- file_info[match(selected_keys, file_info$Plot_Key), , drop = FALSE]
  selected_analyses <- make.unique(selected_info$Display_Name, sep = "_")

  plot_data_list <- vector("list", length(selected_keys))
  names(plot_data_list) <- selected_analyses

  for (analysis_index in seq_along(selected_keys)) {
    dat <- read_result_csv(selected_info$All_Genes_File[analysis_index])

    current_plot_data <- prepare_volcano_data(
      dat = dat,
      analysis_name = selected_analyses[analysis_index],
      p_value_column = P_VALUE_COLUMN,
      p_value_cutoff = P_VALUE_CUTOFF,
      logfc_cutoff = LOGFC_CUTOFF,
      ns_label = "NS",
      regulation_levels = c("NS", "Down", "Up")
    )
    current_plot_data$Plot_LogFC <- current_plot_data$logFC
    plot_data_list[[analysis_index]] <- current_plot_data
  }

  plot_data <- do.call(rbind, plot_data_list)
  plot_data <- plot_data[order(plot_data$Regulation), , drop = FALSE]

  plot_data_for_axis <- plot_data[
    plot_data$Regulation %in% c("Up", "Down"),
    ,
    drop = FALSE
  ]
  if (nrow(plot_data_for_axis) == 0) {
    return(data.frame(
      Plot_Name = scheme_name,
      Analysis_Name = paste(selected_analyses, collapse = ";"),
      Total_Genes = nrow(plot_data),
      Genes_Plotted = 0L,
      Up = 0L,
      Down = 0L,
      NS = nrow(plot_data),
      X_Min = NA_real_,
      X_Max = NA_real_,
      Y_Min = NA_real_,
      Y_Max = NA_real_,
      PDF_Width = NA_real_,
      PDF_Height = NA_real_,
      PDF_File = NA_character_,
      PNG_File = NA_character_,
      Status = "Skipped: no significant Up/Down genes.",
      stringsAsFactors = FALSE
    ))
  }

  group_layout <- get_batch_group_layout_data(
    plot_data = plot_data,
    selected_analyses = selected_analyses
  )
  axis_info <- get_batch_multiple_axis_info(plot_data = plot_data)
  pdf_size <- get_batch_multiple_pdf_size(
    group_count = length(selected_analyses),
    display_y_limits = axis_info$display_limits,
    group_layout = group_layout
  )

  multiple_volcano_plot <- make_batch_multiple_volcano_plot(
    plot_data = plot_data,
    group_layout = group_layout,
    selected_analyses = selected_analyses,
    axis_info = axis_info
  )

  output_dir <- file.path(plot_root, sanitize_file_name(scheme_name))
  if (isTRUE(OVERWRITE_SCHEME_OUTPUT) && dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output_files <- save_ggplot_pdf_png(
    plot = multiple_volcano_plot,
    pdf_file = file.path(output_dir, "multiple_volcano_plot.pdf"),
    width = pdf_size$width,
    height = pdf_size$height
  )

  summary_table <- do.call(rbind, lapply(seq_along(selected_analyses), function(i) {
    analysis_name <- selected_analyses[i]
    dat <- plot_data[plot_data$Analysis_Name == analysis_name, , drop = FALSE]
    plotted_dat <- dat[dat$Regulation %in% c("Up", "Down"), , drop = FALSE]
    status_counts <- table(dat$Regulation)
    layout_row <- group_layout[group_layout$Analysis_Name == analysis_name, , drop = FALSE]

    data.frame(
      Plot_Name = scheme_name,
      Analysis_Name = analysis_name,
      Plot_Category = selected_info$Plot_Category[i],
      Source_Analysis = selected_info$Analysis_Name[i],
      Total_Genes = nrow(dat),
      Genes_Plotted = nrow(plotted_dat),
      Up = count_status(status_counts, "Up"),
      Down = count_status(status_counts, "Down"),
      NS = count_status(status_counts, "NS"),
      X_Min = round(layout_row$X_Min, 2),
      X_Max = round(layout_row$X_Max, 2),
      Y_Min = round(axis_info$true_limits[1], 2),
      Y_Max = round(axis_info$true_limits[2], 2),
      PDF_Width = round(pdf_size$width, 2),
      PDF_Height = round(pdf_size$height, 2),
      PDF_File = output_files$pdf_file,
      PNG_File = output_files$png_file,
      Status = "OK",
      stringsAsFactors = FALSE
    )
  }))

  rownames(summary_table) <- NULL
  summary_table
}

run_batch_multiple_volcano_plots <- function(file_info, schemes, schemes_to_run, plot_root) {
  if (isTRUE(CLEAN_MULTIPLE_VOLCANO_ROOT) && dir.exists(plot_root)) {
    unlink(plot_root, recursive = TRUE)
  }
  dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)

  summary_list <- lapply(schemes_to_run, function(scheme_name) {
    run_batch_multiple_volcano_scheme(
      scheme_name = scheme_name,
      selected_keys = schemes[[scheme_name]],
      file_info = file_info,
      plot_root = plot_root
    )
  })
  names(summary_list) <- schemes_to_run

  summary_table <- do.call(rbind, summary_list)
  rownames(summary_table) <- NULL
  summary_table
}


# 4. GSEA result discovery and dotplots ---------------------------------------

get_gsea_plot_metadata <- function(file, table_root) {
  parts <- get_relative_path_parts(file, table_root)

  if (length(parts) >= 4 &&
      parts[1] != "ATF3_function" &&
      parts[2] == "GSEA") {
    main_deg_file <- file.path(table_root, parts[1], "DEG", "all_genes.csv")
    main_deg_report_file <- file.path(table_root, parts[1], "DEG", "csv", "all_genes.csv")
    if (!file.exists(main_deg_file) && !file.exists(main_deg_report_file)) {
      return(NULL)
    }

    return(list(
      plot_category = "Main_DE",
      analysis_name = parts[1],
      geneset_name = parts[3],
      display_name = parts[1]
    ))
  }

  if (length(parts) >= 5 &&
      parts[1] == "ATF3_function" &&
      parts[3] == "GSEA") {
    return(list(
      plot_category = "ATF3_correlation",
      analysis_name = parts[2],
      geneset_name = parts[4],
      display_name = parts[2]
    ))
  }

  if (length(parts) >= 6 &&
      parts[1] == "ATF3_function" &&
      parts[3] == "ATF3_high_low_DE" &&
      parts[4] == "GSEA") {
    return(list(
      plot_category = "ATF3_high_low_DE",
      analysis_name = parts[2],
      geneset_name = parts[5],
      display_name = parts[2]
    ))
  }

  NULL
}

collect_batch_gsea_file_info <- function(table_root) {
  gsea_files <- list.files(
    table_root,
    pattern = "^gsea_result[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(gsea_files) == 0) {
    return(data.frame())
  }

  info_list <- lapply(gsea_files, function(file) {
    metadata <- get_gsea_plot_metadata(file, table_root)
    if (is.null(metadata)) {
      return(NULL)
    }

    data.frame(
      Plot_Category = metadata$plot_category,
      Analysis_Name = metadata$analysis_name,
      Display_Name = metadata$display_name,
      GeneSet_Name = metadata$geneset_name,
      Plot_Key = paste(
        metadata$plot_category,
        metadata$analysis_name,
        metadata$geneset_name,
        sep = "::"
      ),
      GSEA_Result_File = file,
      stringsAsFactors = FALSE
    )
  })

  info_list <- info_list[!vapply(info_list, is.null, logical(1))]
  if (length(info_list) == 0) {
    return(data.frame())
  }

  file_info <- do.call(rbind, info_list)
  file_info <- file_info[
    order(file_info$Plot_Category, file_info$Analysis_Name, file_info$GeneSet_Name),
    ,
    drop = FALSE
  ]
  rownames(file_info) <- NULL
  file_info
}

filter_gsea_file_info <- function(file_info, categories_to_plot, genesets_to_plot) {
  if (nrow(file_info) == 0) {
    return(file_info)
  }

  keep <- rep(TRUE, nrow(file_info))

  if (!identical(categories_to_plot, "all") &&
      !(length(categories_to_plot) == 1 && tolower(categories_to_plot) == "all")) {
    keep <- keep & file_info$Plot_Category %in% categories_to_plot
  }

  if (!identical(genesets_to_plot, "all") &&
      !(length(genesets_to_plot) == 1 && tolower(genesets_to_plot) == "all")) {
    keep <- keep & file_info$GeneSet_Name %in% genesets_to_plot
  }

  file_info[keep, , drop = FALSE]
}

parse_gene_ratio_value <- function(x) {
  x <- as.character(x)

  vapply(x, function(value) {
    if (is.na(value) || value == "") {
      return(NA_real_)
    }

    if (grepl("/", value, fixed = TRUE)) {
      parts <- strsplit(value, "/", fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        numerator <- as.numeric(parts[1])
        denominator <- as.numeric(parts[2])
        if (is.finite(numerator) && is.finite(denominator) && denominator > 0) {
          return(numerator / denominator)
        }
      }
    }

    as.numeric(value)
  }, numeric(1))
}

get_gsea_point_size_value <- function(result_table) {
  if ("GeneRatio" %in% colnames(result_table)) {
    ratio <- parse_gene_ratio_value(result_table$GeneRatio)
    if (any(is.finite(ratio))) {
      return(ratio)
    }
  }

  if (all(c("Count", "setSize") %in% colnames(result_table))) {
    count <- as.numeric(result_table$Count)
    set_size <- as.numeric(result_table$setSize)
    ratio <- count / set_size
    if (any(is.finite(ratio))) {
      return(ratio)
    }
  }

  if ("setSize" %in% colnames(result_table)) {
    return(as.numeric(result_table$setSize))
  }

  rep(1, nrow(result_table))
}

select_gsea_terms_for_dotplot <- function(result_table) {
  required_columns <- c("Description", "NES", "pvalue", "p.adjust")
  if (nrow(result_table) == 0 || any(!required_columns %in% colnames(result_table))) {
    return(result_table[0, , drop = FALSE])
  }

  result_table$NES <- as.numeric(result_table$NES)
  result_table$pvalue <- as.numeric(result_table$pvalue)
  result_table$`p.adjust` <- as.numeric(result_table$`p.adjust`)

  valid_index <- !is.na(result_table$Description) &
    result_table$Description != "" &
    is.finite(result_table$NES) &
    is.finite(result_table$pvalue) &
    is.finite(result_table$`p.adjust`)
  result_table <- result_table[valid_index, , drop = FALSE]

  if (nrow(result_table) == 0) {
    return(result_table)
  }

  p_column <- GSEA_DOTPLOT_P_COLUMN
  if (!p_column %in% colnames(result_table)) {
    p_column <- "pvalue"
  }

  display_table <- result_table[
    !is.na(result_table[[p_column]]) &
      result_table[[p_column]] <= GSEA_DOTPLOT_P_CUTOFF,
    ,
    drop = FALSE
  ]

  if (nrow(display_table) == 0 && isTRUE(SHOW_NON_SIGNIFICANT_TOP_TERMS)) {
    display_table <- result_table
  }

  if (nrow(display_table) == 0) {
    return(display_table)
  }

  positive_table <- display_table[display_table$NES > 0, , drop = FALSE]
  negative_table <- display_table[display_table$NES < 0, , drop = FALSE]

  if (nrow(positive_table) > 0) {
    positive_table <- positive_table[
      order(positive_table[[p_column]], -abs(positive_table$NES)),
      ,
      drop = FALSE
    ]
    positive_table <- positive_table[seq_len(min(GSEA_DOTPLOT_TOP_N, nrow(positive_table))), , drop = FALSE]
  }

  if (nrow(negative_table) > 0) {
    negative_table <- negative_table[
      order(negative_table[[p_column]], -abs(negative_table$NES)),
      ,
      drop = FALSE
    ]
    negative_table <- negative_table[seq_len(min(GSEA_DOTPLOT_TOP_N, nrow(negative_table))), , drop = FALSE]
  }

  selected_table <- rbind(positive_table, negative_table)
  selected_table <- selected_table[order(selected_table$NES), , drop = FALSE]
  rownames(selected_table) <- NULL
  selected_table
}

make_gsea_table_dotplot <- function(result_table) {
  plot_table <- select_gsea_terms_for_dotplot(result_table)

  if (nrow(plot_table) == 0) {
    return(list(
      plot = make_empty_gsea_plot(),
      shown_terms = 0L,
      plot_labels = character(0)
    ))
  }

  plot_table$Description_For_Plot <- format_gsea_description_for_plot(plot_table$Description)
  plot_table$Description_For_Plot <- wrap_label(plot_table$Description_For_Plot, width = GSEA_DOTPLOT_LABEL_WIDTH)
  plot_table$Point_Size_Value <- get_gsea_point_size_value(plot_table)
  plot_table$Minus_Log10_P <- -log10(plot_table[[GSEA_DOTPLOT_P_COLUMN]])
  plot_table$Direction <- ifelse(plot_table$NES >= 0, "Positive NES", "Negative NES")
  plot_table$Description_For_Plot <- factor(
    plot_table$Description_For_Plot,
    levels = plot_table$Description_For_Plot
  )

  dotplot <- ggplot2::ggplot(
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
      ggplot2::aes(size = Point_Size_Value, color = Minus_Log10_P),
      alpha = 0.88
    ) +
    ggplot2::scale_color_gradient(
      low = DOWN_COLOR,
      high = UP_COLOR,
      name = paste0("-log10(", GSEA_DOTPLOT_P_COLUMN, ")")
    ) +
    ggplot2::scale_size(
      range = GSEAVIS_POINT_SIZE_RANGE,
      name = "GeneRatio"
    ) +
    ggplot2::labs(
      x = "Normalized enrichment score",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.22),
      panel.grid.major.x = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.22),
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
    plot = dotplot,
    shown_terms = nrow(plot_table),
    plot_labels = as.character(plot_table$Description_For_Plot)
  )
}

run_batch_gsea_dotplots <- function(file_info, plot_root) {
  if (nrow(file_info) == 0) {
    return(data.frame())
  }

  summary_list <- lapply(seq_len(nrow(file_info)), function(i) {
    current_info <- file_info[i, , drop = FALSE]
    result_table <- read_result_csv(current_info$GSEA_Result_File)
    dotplot_result <- make_gsea_table_dotplot(result_table)
    plot_size <- get_gsea_dotplot_size(
      shown_terms = dotplot_result$shown_terms,
      plot_labels = dotplot_result$plot_labels
    )

    output_dir <- file.path(
      plot_root,
      sanitize_file_name(current_info$Plot_Category),
      sanitize_file_name(current_info$Analysis_Name),
      sanitize_file_name(current_info$GeneSet_Name)
    )
    if (isTRUE(OVERWRITE_GSEA_DOTPLOT_OUTPUT) && dir.exists(output_dir)) {
      unlink(output_dir, recursive = TRUE)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    output_files <- save_ggplot_pdf_png(
      plot = dotplot_result$plot,
      pdf_file = file.path(output_dir, "dotplot.pdf"),
      width = plot_size$width,
      height = plot_size$height
    )

    data.frame(
      Plot_Category = current_info$Plot_Category,
      Analysis_Name = current_info$Analysis_Name,
      Display_Name = current_info$Display_Name,
      GeneSet_Name = current_info$GeneSet_Name,
      GSEA_Terms = nrow(result_table),
      Terms_Plotted = dotplot_result$shown_terms,
      Positive_NES = count_nes_direction(result_table, "positive"),
      Negative_NES = count_nes_direction(result_table, "negative"),
      PDF_Width = round(plot_size$width, 2),
      PDF_Height = round(plot_size$height, 2),
      CSV_File = resolve_table_file(current_info$GSEA_Result_File),
      PDF_File = output_files$pdf_file,
      PNG_File = output_files$png_file,
      stringsAsFactors = FALSE
    )
  })

  summary_table <- do.call(rbind, summary_list)
  rownames(summary_table) <- NULL
  summary_table
}
