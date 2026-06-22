# Plot PathwayDenester summaries.
#
# Reads normalized PathwayDenester result tables created by
# 01_run_pathway_denester.py and saves project-style PDF/PNG plots.


# 0. Config -------------------------------------------------------------------

RESULTS_ROOT <- "results"
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"

CLEAN_PATHWAY_DENESTER_PLOTS <- TRUE
TOP_HITCHHIKER_N <- 20L

GENESET_ORDER <- c(
  "hallmark",
  "CP_BIOCARTA",
  "CP_KEGG_MEDICUS",
  "CP_KEGG_LEGACY",
  "CP_REACTOME",
  "CP_WIKIPATHWAYS",
  "TFT_TFT_LEGACY",
  "TFT_GTRD",
  "GO_BP",
  "GO_CC",
  "GO_MF",
  "HPO",
  "C6",
  "IMMUNESIGDB"
)

options(width = 200)


# 1. Packages and helpers ------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

source(REPORT_TABLE_FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)

get_numeric_column <- function(dat, column_name) {
  if (!column_name %in% colnames(dat)) {
    return(rep(NA_real_, nrow(dat)))
  }

  suppressWarnings(as.numeric(dat[[column_name]]))
}

make_plot_label <- function(dat) {
  ifelse(
    dat$Plot_Category == "Main_DE",
    dat$Analysis_Name,
    paste(dat$Plot_Category, dat$Analysis_Name, sep = "::")
  )
}

wrap_pathway_label <- function(x, width = 42) {
  x <- as.character(x)
  vapply(x, function(label) {
    label <- gsub("_", " ", label)
    wrap_label(label, width = width)
  }, character(1))
}

make_empty_pathway_denester_plot <- function(label = "No excluded hitchhiker terms") {
  ggplot() +
    geom_text(
      aes(x = 0, y = 0, label = label),
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      size = 4.2,
      color = TEXT_COLOR
    ) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void(base_family = TEXT_FONT_FAMILY) +
    theme(
      text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      plot.margin = margin(8, 8, 8, 8)
    )
}

apply_pathway_denester_theme <- function(plot) {
  plot +
    theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    theme(
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.22),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
      axis.text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.title = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.key = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "grey90", color = TEXT_COLOR, linewidth = AXIS_LINE_WIDTH),
      strip.text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      plot.title = element_blank(),
      plot.margin = margin(8, 10, 8, 8)
    )
}

get_safe_neglog10 <- function(p_value) {
  p_value <- as.numeric(p_value)
  positive_p <- p_value[is.finite(p_value) & p_value > 0]
  if (length(positive_p) == 0) {
    positive_p <- 1e-300
  }
  safe_p <- p_value
  safe_p[!is.finite(safe_p) | safe_p <= 0] <- min(positive_p, na.rm = TRUE) * 0.1
  -log10(safe_p)
}

rbind_fill <- function(...) {
  dat_list <- list(...)
  dat_list <- dat_list[!vapply(dat_list, is.null, logical(1))]
  if (length(dat_list) == 0) {
    return(data.frame())
  }

  all_columns <- unique(unlist(lapply(dat_list, colnames), use.names = FALSE))
  dat_list <- lapply(dat_list, function(dat) {
    missing_columns <- setdiff(all_columns, colnames(dat))
    for (column_name in missing_columns) {
      dat[[column_name]] <- NA
    }
    dat[, all_columns, drop = FALSE]
  })

  do.call(rbind, dat_list)
}


# 2. Summary plots -------------------------------------------------------------

make_dataset_heatmap <- function(summary_dat) {
  plot_dat <- summary_dat[summary_dat$Status %in% c("OK", "Skipped: existing result"), , drop = FALSE]
  plot_dat$Plot_Label <- make_plot_label(plot_dat)
  plot_dat$Plot_Label <- factor(
    plot_dat$Plot_Label,
    levels = rev(unique(plot_dat$Plot_Label))
  )
  plot_dat$GeneSet_Name <- factor(
    plot_dat$GeneSet_Name,
    levels = GENESET_ORDER[GENESET_ORDER %in% plot_dat$GeneSet_Name]
  )

  plot <- ggplot(
    plot_dat,
    aes(x = GeneSet_Name, y = Plot_Label, fill = Excluded_Percent)
  ) +
    geom_tile(color = "white", linewidth = 0.45) +
    scale_fill_gradient(
      low = "#F7FBFF",
      high = UP_COLOR,
      limits = c(0, 100),
      name = "Excluded %"
    ) +
    scale_x_discrete(labels = function(x) wrap_label_by_underscore(x, width = 11)) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off")

  apply_pathway_denester_theme(plot) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank()
    )
}

make_dataset_count_barplot <- function(summary_dat) {
  plot_dat <- summary_dat[summary_dat$Status %in% c("OK", "Skipped: existing result"), , drop = FALSE]
  plot_dat$Plot_Label <- make_plot_label(plot_dat)
  plot_dat$Plot_Label <- factor(
    plot_dat$Plot_Label,
    levels = unique(plot_dat$Plot_Label)
  )

  count_dat <- rbind(
    data.frame(
      Plot_Label = plot_dat$Plot_Label,
      Plot_Category = plot_dat$Plot_Category,
      Count_Type = "Kept",
      Count = as.numeric(plot_dat$Kept_Terms),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Plot_Label = plot_dat$Plot_Label,
      Plot_Category = plot_dat$Plot_Category,
      Count_Type = "Excluded",
      Count = as.numeric(plot_dat$Excluded_Terms),
      stringsAsFactors = FALSE
    )
  )

  plot <- ggplot(
    count_dat,
    aes(x = Plot_Label, y = Count, fill = Count_Type)
  ) +
    geom_col(width = 0.72, color = TEXT_COLOR, linewidth = 0.18) +
    facet_wrap(~ Plot_Category, scales = "free_x", nrow = 1) +
    scale_fill_manual(
      values = c(Kept = DOWN_COLOR, Excluded = UP_COLOR),
      breaks = c("Kept", "Excluded"),
      name = NULL
    ) +
    labs(x = NULL, y = "Pathway count") +
    scale_x_discrete(labels = function(x) wrap_label_by_underscore(x, width = 12))

  apply_pathway_denester_theme(plot) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right"
    )
}

save_dataset_summary_plots <- function(summary_dat, result_root) {
  output_dir <- file.path(
    RESULTS_ROOT,
    summary_dat$Data_Type[1],
    summary_dat$Dataset_ID[1],
    "plots",
    "PathwayDenester",
    "summary"
  )

  if (CLEAN_PATHWAY_DENESTER_PLOTS && dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  plot_label_count <- length(unique(make_plot_label(summary_dat)))
  heatmap_width <- max(9.5, length(unique(summary_dat$GeneSet_Name)) * 0.52 + 3.2)
  heatmap_height <- max(5.4, plot_label_count * 0.30 + 2.1)

  heatmap_files <- save_ggplot_pdf_png(
    plot = make_dataset_heatmap(summary_dat),
    pdf_file = file.path(output_dir, "pathway_denester_exclusion_heatmap.pdf"),
    width = heatmap_width,
    height = heatmap_height
  )

  bar_width <- max(9.2, plot_label_count * 0.30 + 4.0)
  bar_height <- 6.4
  bar_files <- save_ggplot_pdf_png(
    plot = make_dataset_count_barplot(summary_dat),
    pdf_file = file.path(output_dir, "pathway_denester_kept_excluded_barplot.pdf"),
    width = bar_width,
    height = bar_height
  )

  data.frame(
    Data_Type = summary_dat$Data_Type[1],
    Dataset_ID = summary_dat$Dataset_ID[1],
    Plot_Type = c("exclusion_heatmap", "kept_excluded_barplot"),
    PDF_File = c(heatmap_files$pdf_file, bar_files$pdf_file),
    PNG_File = c(heatmap_files$png_file, bar_files$png_file),
    stringsAsFactors = FALSE
  )
}


# 3. Per-result plots ----------------------------------------------------------

make_top_hitchhiker_plot <- function(result_dat) {
  if (!"Filtered" %in% colnames(result_dat)) {
    return(make_empty_pathway_denester_plot("PathwayDenester result is missing Filtered column"))
  }

  excluded_dat <- result_dat[
    tolower(as.character(result_dat$Filtered)) == "exclude",
    ,
    drop = FALSE
  ]

  if (nrow(excluded_dat) == 0) {
    return(make_empty_pathway_denester_plot())
  }

  excluded_dat$Result_Numeric <- get_numeric_column(excluded_dat, "Result")
  excluded_dat$Reciprocal_Numeric <- get_numeric_column(excluded_dat, "Reciprocal pvalue")
  excluded_dat$DEGs_in_Intersection <- get_numeric_column(excluded_dat, "DEGs in Intersection")
  excluded_dat$Neg_Log10_Result <- get_safe_neglog10(excluded_dat$Result_Numeric)
  excluded_dat <- excluded_dat[
    order(excluded_dat$Result_Numeric, -excluded_dat$DEGs_in_Intersection),
    ,
    drop = FALSE
  ]
  excluded_dat <- excluded_dat[seq_len(min(TOP_HITCHHIKER_N, nrow(excluded_dat))), , drop = FALSE]

  excluded_dat$Pathway_Label <- wrap_pathway_label(excluded_dat$Name, width = 45)
  excluded_dat$Driver_Label <- wrap_pathway_label(excluded_dat$`Versus Name`, width = 32)
  excluded_dat$Pathway_Label <- factor(
    excluded_dat$Pathway_Label,
    levels = rev(excluded_dat$Pathway_Label)
  )

  plot <- ggplot(
    excluded_dat,
    aes(x = Neg_Log10_Result, y = Pathway_Label)
  ) +
    geom_col(
      aes(fill = DEGs_in_Intersection),
      width = 0.72,
      color = TEXT_COLOR,
      linewidth = 0.18
    ) +
    scale_fill_gradient(
      low = "#FEE8C8",
      high = UP_COLOR,
      name = "Leading-edge\noverlap"
    ) +
    labs(
      x = "-log10(PathwayDenester p value)",
      y = NULL
    )

  apply_pathway_denester_theme(plot)
}

save_top_hitchhiker_plot <- function(summary_row, result_dat, output_dir) {
  plot <- make_top_hitchhiker_plot(result_dat)

  shown_terms <- min(TOP_HITCHHIKER_N, as.numeric(summary_row$Excluded_Terms))
  plot_height <- max(5.2, shown_terms * 0.36 + 1.8)
  plot_width <- 8.8

  files <- save_ggplot_pdf_png(
    plot = plot,
    pdf_file = file.path(output_dir, "top_hitchhiker_terms.pdf"),
    width = plot_width,
    height = plot_height
  )

  data.frame(
    Data_Type = summary_row$Data_Type,
    Dataset_ID = summary_row$Dataset_ID,
    Plot_Category = summary_row$Plot_Category,
    Analysis_Name = summary_row$Analysis_Name,
    GeneSet_Name = summary_row$GeneSet_Name,
    Plot_Type = "top_hitchhiker_terms",
    Terms_Plotted = shown_terms,
    PDF_File = files$pdf_file,
    PNG_File = files$png_file,
    stringsAsFactors = FALSE
  )
}

plot_one_pathway_denester_result <- function(summary_row) {
  result_file <- summary_row$PathwayDenester_Result_File
  result_dat <- read.csv(result_file, stringsAsFactors = FALSE, check.names = FALSE)

  output_dir <- file.path(
    RESULTS_ROOT,
    summary_row$Data_Type,
    summary_row$Dataset_ID,
    "plots",
    "PathwayDenester",
    sanitize_file_name(summary_row$Plot_Category),
    sanitize_file_name(summary_row$Analysis_Name),
    sanitize_file_name(summary_row$GeneSet_Name)
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Pathway-overlap heatmaps are drawn as PathwayDenester-style combo plots
  # by 03_plot_pathway_overlap_combo_heatmap.py.
  save_top_hitchhiker_plot(
    summary_row = summary_row,
    result_dat = result_dat,
    output_dir = output_dir
  )
}


# 4. Run -----------------------------------------------------------------------

summary_file <- file.path(RESULTS_ROOT, "PathwayDenester_summary", "summary.csv")
if (!file.exists(summary_file)) {
  stop("PathwayDenester summary was not found. Run 01_run_pathway_denester.py first: ", summary_file)
}

summary_dat <- read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)
summary_dat <- summary_dat[summary_dat$Status %in% c("OK", "Skipped: existing result"), , drop = FALSE]
stopifnot(nrow(summary_dat) > 0)

cat("\nRunning PathwayDenester plotting...\n")
cat("Result tables: ", nrow(summary_dat), "\n", sep = "")

dataset_keys <- unique(paste(summary_dat$Data_Type, summary_dat$Dataset_ID, sep = "\t"))
summary_plot_records <- do.call(
  rbind,
  lapply(dataset_keys, function(dataset_key) {
    parts <- strsplit(dataset_key, "\t", fixed = TRUE)[[1]]
    dat <- summary_dat[
      summary_dat$Data_Type == parts[1] &
        summary_dat$Dataset_ID == parts[2],
      ,
      drop = FALSE
    ]
    save_dataset_summary_plots(dat, RESULTS_ROOT)
  })
)

per_result_records <- do.call(
  rbind,
  lapply(seq_len(nrow(summary_dat)), function(i) {
    plot_one_pathway_denester_result(summary_dat[i, , drop = FALSE])
  })
)

plot_summary <- rbind_fill(summary_plot_records, per_result_records)
rownames(plot_summary) <- NULL

global_plot_summary_dir <- file.path(RESULTS_ROOT, "PathwayDenester_plot_summary")
dir.create(global_plot_summary_dir, recursive = TRUE, showWarnings = FALSE)
write_csv_with_report_previews(
  plot_summary,
  file.path(global_plot_summary_dir, "summary.csv"),
  n_rows = 21
)

for (dataset_key in dataset_keys) {
  parts <- strsplit(dataset_key, "\t", fixed = TRUE)[[1]]
  dat <- plot_summary[
    plot_summary$Data_Type == parts[1] &
      plot_summary$Dataset_ID == parts[2],
    ,
    drop = FALSE
  ]
  out_dir <- file.path(RESULTS_ROOT, parts[1], parts[2], "tables", "PathwayDenester_plot_summary")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv_with_report_previews(
    dat,
    file.path(out_dir, "summary.csv"),
    n_rows = 21
  )
}

cat("\nPathwayDenester plotting finished.\n")
cat("Global plot summary: ", file.path(global_plot_summary_dir, "summary.csv"), "\n", sep = "")
