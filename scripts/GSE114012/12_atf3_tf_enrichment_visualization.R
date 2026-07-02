# GSE114012 ATF3-focused TF enrichment visualization
#
# Draw bubble plots for DoRothEA, ChEA3 and CollecTRI method-final results,
# then draw 3-set Venn diagrams for their TF overlaps.


# 0. Config -------------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"
TARGET_TF <- "ATF3"

PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TF_DEG_SUMMARY_ROOT <- file.path(RESULT_ROOT, "TF_summary", "deg")

PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "TF_ATF3_visualization")
BUBBLE_PLOT_ROOT <- file.path(PLOT_ROOT, "bubble")
VENN_PLOT_ROOT <- file.path(PLOT_ROOT, "venn")

TABLE_ROOT <- file.path(RESULT_ROOT, "tables", "TF_ATF3_visualization")

# The three methods selected from the ATF3 ranking strategy search.
SELECTED_METHODS <- c("dorothea", "chea3", "collectri")
TF_METHOD_LABELS <- c(
  dorothea = "DoRothEA",
  chea3 = "ChEA3",
  collectri = "CollecTRI"
)

# Set to "all" for all DEG analysis schemes, or provide selected names.
ANALYSES_TO_PLOT <- "all"
# ANALYSES_TO_PLOT <- c("ALL", "DLD1_HCT15", "DLD1_HCT15_SW48", "HCT15", "SW948")

# Bubble plot settings.
TOP_N_TF_TO_PLOT <- 20
LABEL_TOP_N_TF <- 3
INCLUDE_TARGET_TF_IF_OUTSIDE_TOP_N <- TRUE

# User-customizable colors.
BUBBLE_COLORS <- c(
  "ATF3" = "#D73027",
  "Top candidate" = "#2C7FB8",
  "Other" = "#BDBDBD"
)

VENN_COLORS <- c(
  "DoRothEA" = "#00A087",
  "ChEA3" = "#3C5488",
  "CollecTRI" = "#E64B35"
)
VENN_ALPHA <- 0.34
VENN_BORDER_WIDTH <- 1.1

TARGET_LABEL_FILL <- "#FFF4F0"
TARGET_LABEL_COLOR <- "#B2182B"

P_VALUE_REFERENCE <- 0.05
MIN_P_VALUE_FOR_LOG <- 1e-300
MIN_SCORE_FOR_LOG <- 1e-300

BUBBLE_SIZE_RANGE <- c(3.2, 11.0)
BUBBLE_POINT_ALPHA <- 0.78
BUBBLE_STROKE_WIDTH <- 0.35

BUBBLE_PDF_WIDTH <- 6.8
BUBBLE_PDF_HEIGHT <- 6.2
VENN_PDF_WIDTH <- 6.4
VENN_PDF_HEIGHT <- 5.7

CLEAN_TF_ATF3_VISUALIZATION_OUTPUT <- TRUE

options(width = 200)


# 1. Packages and helpers ------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

source(PLOTTING_FUNCTION_FILE)
source(REPORT_TABLE_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()
USE_GG_REPEL <- requireNamespace("ggrepel", quietly = TRUE)

safe_read_csv <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
  }

  tryCatch(
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

to_number <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

standardize_tf <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == "" | is.na(x) | x == "NA"] <- NA_character_
  x
}

first_existing_column <- function(dat, column_names) {
  hit <- column_names[column_names %in% colnames(dat)]
  if (length(hit) == 0L) {
    return(NA_character_)
  }
  hit[1]
}

get_numeric_column <- function(dat, column_names, default = NA_real_) {
  column <- first_existing_column(dat, column_names)
  if (is.na(column)) {
    return(rep(default, nrow(dat)))
  }

  to_number(dat[[column]])
}

get_character_column <- function(dat, column_names, default = "") {
  column <- first_existing_column(dat, column_names)
  if (is.na(column)) {
    return(rep(default, nrow(dat)))
  }

  as.character(dat[[column]])
}

safe_neg_log10 <- function(x, min_value) {
  x <- to_number(x)
  x[is.na(x) | x <= 0] <- min_value
  -log10(x)
}

get_method_result_file <- function(analysis_name, method_name) {
  method_dir <- file.path(TF_DEG_SUMMARY_ROOT, analysis_name, "method_final", method_name)
  if (!dir.exists(method_dir)) {
    return(NA_character_)
  }

  files <- list.files(method_dir, pattern = "[.]csv$", full.names = TRUE)
  files <- files[!grepl("/summary/", files)]
  if (length(files) == 0L) {
    return(NA_character_)
  }

  files[1]
}

read_method_final_table <- function(analysis_name, method_name) {
  result_file <- get_method_result_file(analysis_name, method_name)
  if (is.na(result_file)) {
    return(NULL)
  }

  dat <- safe_read_csv(result_file)
  if (is.null(dat) || nrow(dat) == 0L || !"TF" %in% colnames(dat)) {
    return(NULL)
  }

  dat$TF <- standardize_tf(dat$TF)
  dat <- dat[!is.na(dat$TF) & dat$TF != "", , drop = FALSE]
  if (nrow(dat) == 0L) {
    return(NULL)
  }

  dat$Analysis_Name <- analysis_name
  dat$Method <- method_name
  dat$Method_Label <- TF_METHOD_LABELS[[method_name]]
  dat$Method_Rank <- get_numeric_column(dat, c("Method_Rank", "Rank", "Best_Rank"))
  dat$Method_Score <- get_numeric_column(dat, c("Method_Score", "Score", "score", "Best_Combined_Score", "Activity_Score_Mean"))
  dat$Method_P_Value <- get_numeric_column(dat, c("Method_P_Value", "p_value", "P.Value", "p.value", "Best_P_Value"))
  dat$Method_Adjusted_P_Value <- get_numeric_column(dat, c("Method_Adjusted_P_Value", "Best_Adjusted_P_Value", "adj.P.Val", "p.adjust", "FDR"))
  dat$Method_Direction <- get_character_column(dat, c("Method_Direction", "Direction", "direction"))
  dat$CheA3_Library_Count <- get_numeric_column(dat, "CheA3_Library_Count", default = 0)
  dat$CheA3_Integrated_TopRank <- get_numeric_column(dat, "CheA3_Integrated_TopRank")
  dat$Source_File <- result_file

  dat$Method_Rank[is.na(dat$Method_Rank)] <- seq_len(nrow(dat))[is.na(dat$Method_Rank)]
  dat$Method_Direction[is.na(dat$Method_Direction)] <- ""
  dat$CheA3_Library_Count[is.na(dat$CheA3_Library_Count)] <- 0

  dat <- dat[order(dat$Method_Rank, dat$TF, na.last = TRUE), , drop = FALSE]
  dat[!duplicated(dat$TF), , drop = FALSE]
}

prepare_bubble_data <- function(method_table) {
  stopifnot(nrow(method_table) > 0L)
  method_name <- method_table$Method[1]

  method_table$Bubble_Group <- ifelse(
    method_table$TF == TARGET_TF,
    TARGET_TF,
    ifelse(method_table$Method_Rank <= LABEL_TOP_N_TF, "Top candidate", "Other")
  )
  method_table$Bubble_Group <- factor(
    method_table$Bubble_Group,
    levels = c(TARGET_TF, "Top candidate", "Other")
  )

  if (method_name == "chea3") {
    method_table$Evidence_Strength <- safe_neg_log10(method_table$Method_Score, MIN_SCORE_FOR_LOG)
    method_table$Bubble_Size_Value <- method_table$Evidence_Strength
    method_table$Evidence_Axis_Label <- "-log10(ChEA3 integrated score)"
    method_table$Size_Legend_Label <- "-log10(score)"
    method_table$Reference_X <- NA_real_
  } else {
    method_table$Evidence_Strength <- safe_neg_log10(method_table$Method_P_Value, MIN_P_VALUE_FOR_LOG)
    method_table$Bubble_Size_Value <- abs(method_table$Method_Score)
    method_table$Evidence_Axis_Label <- "-log10(P value)"
    method_table$Size_Legend_Label <- "TF score"
    method_table$Reference_X <- -log10(P_VALUE_REFERENCE)
  }

  top_index <- method_table$Method_Rank <= TOP_N_TF_TO_PLOT
  if (INCLUDE_TARGET_TF_IF_OUTSIDE_TOP_N) {
    top_index <- top_index | method_table$TF == TARGET_TF
  }

  plot_data <- method_table[top_index, , drop = FALSE]
  plot_data <- plot_data[order(plot_data$Method_Rank, plot_data$TF), , drop = FALSE]
  plot_data$Label <- ifelse(
    plot_data$TF == TARGET_TF | plot_data$Method_Rank <= LABEL_TOP_N_TF,
    plot_data$TF,
    ""
  )
  plot_data$Target_TF_In_Top_N <- any(
    plot_data$TF == TARGET_TF &
      plot_data$Method_Rank <= TOP_N_TF_TO_PLOT
  )

  plot_data
}

make_tf_bubble_plot <- function(plot_data, analysis_name, method_name) {
  method_label <- TF_METHOD_LABELS[[method_name]]
  x_label <- plot_data$Evidence_Axis_Label[1]
  size_label <- plot_data$Size_Legend_Label[1]
  reference_x <- plot_data$Reference_X[1]

  y_max <- max(plot_data$Method_Rank, na.rm = TRUE)
  y_breaks <- sort(unique(c(
    1,
    pretty(c(1, y_max), n = 6),
    plot_data$Method_Rank[plot_data$TF == TARGET_TF]
  )))
  y_breaks <- y_breaks[y_breaks >= 1 & y_breaks <= y_max]

  p <- ggplot(
    plot_data,
    aes(
      x = Evidence_Strength,
      y = Method_Rank,
      size = Bubble_Size_Value,
      fill = Bubble_Group
    )
  ) +
    geom_point(
      shape = 21,
      color = TEXT_COLOR,
      alpha = BUBBLE_POINT_ALPHA,
      stroke = BUBBLE_STROKE_WIDTH
    ) +
    scale_fill_manual(values = BUBBLE_COLORS, drop = FALSE) +
    scale_size_continuous(range = BUBBLE_SIZE_RANGE) +
    scale_y_reverse(
      breaks = y_breaks,
      limits = c(y_max + 0.6, 0.4),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_x_continuous(
      breaks = pretty(plot_data$Evidence_Strength, n = 5),
      expand = expansion(mult = c(0.04, 0.12))
    ) +
    labs(
      x = x_label,
      y = "TF rank (1 = top)",
      fill = NULL,
      size = size_label,
      title = paste0(analysis_name, " / ", method_label)
    ) +
    theme_bw(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    theme(
      panel.grid.major = element_line(color = "#E7E7E7", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = TEXT_COLOR, fill = NA, linewidth = AXIS_LINE_WIDTH),
      axis.line = element_line(color = TEXT_COLOR, linewidth = AXIS_LINE_WIDTH),
      axis.text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      axis.title = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      plot.title = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE, hjust = 0.5, size = BASE_FONT_SIZE + 1),
      legend.position = "right",
      legend.text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.title = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      legend.key = element_blank(),
      text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      plot.margin = margin(10, 14, 10, 10, unit = "pt")
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = 5.2, alpha = 0.85)),
      size = guide_legend(override.aes = list(fill = "#A6CEE3", alpha = 0.85))
    )

  if (!is.na(reference_x)) {
    p <- p +
      geom_vline(
        xintercept = reference_x,
        linewidth = 0.45,
        linetype = "dashed",
        color = "#555555"
      ) +
      annotate(
        "text",
        x = reference_x,
        y = y_max + 0.35,
        label = "P = 0.05",
        hjust = -0.04,
        vjust = 1,
        size = 3.1,
        family = TEXT_FONT_FAMILY,
        fontface = TEXT_FONT_FACE,
        color = "#555555"
      )
  }

  label_data <- plot_data[plot_data$Label != "", , drop = FALSE]
  if (nrow(label_data) > 0L) {
    if (USE_GG_REPEL) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_data,
          aes(label = Label),
          size = 3.4,
          family = TEXT_FONT_FAMILY,
          fontface = TEXT_FONT_FACE,
          color = ifelse(label_data$TF == TARGET_TF, TARGET_LABEL_COLOR, TEXT_COLOR),
          box.padding = 0.30,
          point.padding = 0.18,
          min.segment.length = 0,
          segment.size = 0.25,
          max.overlaps = Inf,
          seed = 114012,
          show.legend = FALSE
        )
    } else {
      p <- p +
        geom_text(
          data = label_data,
          aes(label = Label),
          nudge_y = -0.35,
          size = 3.2,
          family = TEXT_FONT_FAMILY,
          fontface = TEXT_FONT_FACE,
          color = ifelse(label_data$TF == TARGET_TF, TARGET_LABEL_COLOR, TEXT_COLOR),
          show.legend = FALSE
        )
    }
  }

  p
}

get_tf_set <- function(method_table) {
  if (is.null(method_table) || nrow(method_table) == 0L) {
    return(character(0))
  }

  sort(unique(method_table$TF[!is.na(method_table$TF) & method_table$TF != ""]))
}

get_venn_counts <- function(sets) {
  a <- sets[[1]]
  b <- sets[[2]]
  c <- sets[[3]]

  data.frame(
    Region = c("A_only", "B_only", "C_only", "AB_only", "AC_only", "BC_only", "ABC"),
    Count = c(
      length(setdiff(a, union(b, c))),
      length(setdiff(b, union(a, c))),
      length(setdiff(c, union(a, b))),
      length(setdiff(intersect(a, b), c)),
      length(setdiff(intersect(a, c), b)),
      length(setdiff(intersect(b, c), a)),
      length(Reduce(intersect, sets))
    ),
    stringsAsFactors = FALSE
  )
}

make_circle_polygon <- function(center_x, center_y, radius, n = 240) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = center_x + radius * cos(theta),
    y = center_y + radius * sin(theta)
  )
}

make_venn_plot <- function(sets, analysis_name) {
  method_labels <- unname(TF_METHOD_LABELS[SELECTED_METHODS])
  names(sets) <- method_labels

  circle_params <- data.frame(
    Set = method_labels,
    Center_X = c(-0.78, 0.78, 0),
    Center_Y = c(0.35, 0.35, -0.55),
    Radius = c(1.22, 1.22, 1.22),
    stringsAsFactors = FALSE
  )

  circle_data <- do.call(rbind, lapply(seq_len(nrow(circle_params)), function(i) {
    dat <- make_circle_polygon(
      center_x = circle_params$Center_X[i],
      center_y = circle_params$Center_Y[i],
      radius = circle_params$Radius[i]
    )
    dat$Set <- circle_params$Set[i]
    dat
  }))

  counts <- get_venn_counts(sets)
  count_positions <- data.frame(
    Region = c("A_only", "B_only", "C_only", "AB_only", "AC_only", "BC_only", "ABC"),
    x = c(-1.34, 1.34, 0, 0, -0.56, 0.56, 0),
    y = c(0.40, 0.40, -1.24, 0.76, -0.38, -0.38, 0.08),
    stringsAsFactors = FALSE
  )
  counts <- merge(counts, count_positions, by = "Region", all.x = TRUE, sort = FALSE)

  set_labels <- data.frame(
    Set = method_labels,
    Label = paste0(method_labels, "\n", "n = ", vapply(sets, length, integer(1))),
    x = c(-1.45, 1.45, 0),
    y = c(1.70, 1.70, -1.96),
    stringsAsFactors = FALSE
  )

  in_three_way <- TARGET_TF %in% Reduce(intersect, sets)
  target_label <- data.frame(
    x = 0,
    y = -0.18,
    Label = ifelse(in_three_way, paste0(TARGET_TF, "\n3-way overlap"), paste0(TARGET_TF, "\nnot in 3-way overlap")),
    stringsAsFactors = FALSE
  )

  ggplot() +
    geom_polygon(
      data = circle_data,
      aes(x = x, y = y, group = Set, fill = Set, color = Set),
      alpha = VENN_ALPHA,
      linewidth = VENN_BORDER_WIDTH
    ) +
    geom_text(
      data = counts,
      aes(x = x, y = y, label = Count),
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      color = TEXT_COLOR,
      size = 4.4
    ) +
    geom_label(
      data = target_label,
      aes(x = x, y = y, label = Label),
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      color = TARGET_LABEL_COLOR,
      fill = TARGET_LABEL_FILL,
      linewidth = 0.35,
      size = 3.4,
      lineheight = 0.95
    ) +
    geom_text(
      data = set_labels,
      aes(x = x, y = y, label = Label, color = Set),
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      size = 4.2,
      lineheight = 1.02
    ) +
    scale_fill_manual(values = VENN_COLORS[method_labels], drop = FALSE) +
    scale_color_manual(values = darken_color(VENN_COLORS[method_labels], 0.72), drop = FALSE) +
    coord_equal(xlim = c(-2.25, 2.25), ylim = c(-2.15, 2.0), expand = FALSE) +
    labs(title = paste0(analysis_name, " / DoRothEA + ChEA3 + CollecTRI")) +
    theme_void(base_size = BASE_FONT_SIZE, base_family = TEXT_FONT_FAMILY) +
    theme(
      plot.title = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE, hjust = 0.5, size = BASE_FONT_SIZE + 1),
      legend.position = "none",
      text = element_text(color = TEXT_COLOR, face = TEXT_FONT_FACE),
      plot.margin = margin(10, 10, 10, 10, unit = "pt")
    )
}

make_intersection_candidate_table <- function(method_tables) {
  sets <- lapply(method_tables, get_tf_set)
  intersected_tfs <- Reduce(intersect, sets)
  if (length(intersected_tfs) == 0L) {
    return(data.frame())
  }

  do.call(rbind, lapply(sort(intersected_tfs), function(tf) {
    ranks <- vapply(SELECTED_METHODS, function(method_name) {
      tab <- method_tables[[method_name]]
      tab$Method_Rank[match(tf, tab$TF)]
    }, numeric(1))

    p_values <- vapply(SELECTED_METHODS, function(method_name) {
      tab <- method_tables[[method_name]]
      tab$Method_P_Value[match(tf, tab$TF)]
    }, numeric(1))

    scores <- vapply(SELECTED_METHODS, function(method_name) {
      tab <- method_tables[[method_name]]
      tab$Method_Score[match(tf, tab$TF)]
    }, numeric(1))

    first_tab <- method_tables[[SELECTED_METHODS[1]]]
    first_idx <- match(tf, first_tab$TF)

    data.frame(
      TF = tf,
      Mean_Selected_Rank = mean(ranks, na.rm = TRUE),
      Best_Selected_Rank = min(ranks, na.rm = TRUE),
      Worst_Selected_Rank = max(ranks, na.rm = TRUE),
      DoRothEA_Rank = ranks[["dorothea"]],
      ChEA3_Rank = ranks[["chea3"]],
      CollecTRI_Rank = ranks[["collectri"]],
      DoRothEA_P_Value = p_values[["dorothea"]],
      CollecTRI_P_Value = p_values[["collectri"]],
      DoRothEA_Score = scores[["dorothea"]],
      ChEA3_Score = scores[["chea3"]],
      CollecTRI_Score = scores[["collectri"]],
      CheA3_Library_Count = if (!is.na(first_idx)) first_tab$CheA3_Library_Count[first_idx] else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

rank_intersection_candidates <- function(candidate_table) {
  if (is.null(candidate_table) || nrow(candidate_table) == 0L) {
    return(candidate_table)
  }

  candidate_table <- candidate_table[order(
    candidate_table$Mean_Selected_Rank,
    candidate_table$Best_Selected_Rank,
    -candidate_table$CheA3_Library_Count,
    candidate_table$TF,
    na.last = TRUE
  ), , drop = FALSE]
  candidate_table$Consensus_Rank <- seq_len(nrow(candidate_table))
  candidate_table
}


# 2. Discover inputs -----------------------------------------------------------

all_analysis_names <- basename(list.dirs(TF_DEG_SUMMARY_ROOT, recursive = FALSE, full.names = TRUE))
all_analysis_names <- all_analysis_names[all_analysis_names != ""]
stopifnot(length(all_analysis_names) > 0L)

selected_analyses <- if (identical(ANALYSES_TO_PLOT, "all")) {
  all_analysis_names
} else {
  intersect(ANALYSES_TO_PLOT, all_analysis_names)
}
stopifnot(length(selected_analyses) > 0L)

if (CLEAN_TF_ATF3_VISUALIZATION_OUTPUT) {
  unlink(PLOT_ROOT, recursive = TRUE, force = TRUE)
  unlink(TABLE_ROOT, recursive = TRUE, force = TRUE)
}
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)

method_tables_by_analysis <- setNames(vector("list", length(selected_analyses)), selected_analyses)
for (analysis_name in selected_analyses) {
  method_tables_by_analysis[[analysis_name]] <- setNames(
    lapply(SELECTED_METHODS, function(method_name) {
      read_method_final_table(analysis_name, method_name)
    }),
    SELECTED_METHODS
  )
}


# 3. Plot bubbles --------------------------------------------------------------

cat("\nRunning ATF3-focused TF bubble plots...\n")
cat("Methods: ", paste(unname(TF_METHOD_LABELS[SELECTED_METHODS]), collapse = ", "), "\n", sep = "")
cat("Analyses: ", paste(selected_analyses, collapse = ", "), "\n", sep = "")

bubble_tasks <- expand.grid(
  Analysis_Name = selected_analyses,
  Method = SELECTED_METHODS,
  stringsAsFactors = FALSE
)

run_one_bubble_plot <- function(i) {
  analysis_name <- bubble_tasks$Analysis_Name[i]
  method_name <- bubble_tasks$Method[i]
  method_table <- method_tables_by_analysis[[analysis_name]][[method_name]]
  if (is.null(method_table) || nrow(method_table) == 0L) {
    return(data.frame(
      Analysis_Name = analysis_name,
      Method = method_name,
      Method_Label = TF_METHOD_LABELS[[method_name]],
      Status = "missing_method_table",
      stringsAsFactors = FALSE
    ))
  }

  plot_data <- prepare_bubble_data(method_table)
  bubble_plot <- make_tf_bubble_plot(plot_data, analysis_name, method_name)

  output_dir <- file.path(
    BUBBLE_PLOT_ROOT,
    sanitize_file_name(method_name),
    sanitize_file_name(analysis_name)
  )
  pdf_file <- file.path(output_dir, "tf_bubble_plot.pdf")
  output_files <- save_ggplot_pdf_png(
    plot = bubble_plot,
    pdf_file = pdf_file,
    width = BUBBLE_PDF_WIDTH,
    height = BUBBLE_PDF_HEIGHT
  )

  atf3_row <- method_table[method_table$TF == TARGET_TF, , drop = FALSE]
  if (nrow(atf3_row) == 0L) {
    atf3_row <- method_table[0, , drop = FALSE]
  }

  data.frame(
    Analysis_Name = analysis_name,
    Method = method_name,
    Method_Label = TF_METHOD_LABELS[[method_name]],
    Status = "ok",
    TFs_Plotted = nrow(plot_data),
    ATF3_Found = nrow(atf3_row) > 0L,
    ATF3_Rank = if (nrow(atf3_row) > 0L) atf3_row$Method_Rank[1] else NA_real_,
    ATF3_P_Value = if (nrow(atf3_row) > 0L) atf3_row$Method_P_Value[1] else NA_real_,
    ATF3_Adjusted_P_Value = if (nrow(atf3_row) > 0L) atf3_row$Method_Adjusted_P_Value[1] else NA_real_,
    ATF3_Score = if (nrow(atf3_row) > 0L) atf3_row$Method_Score[1] else NA_real_,
    ATF3_Evidence_Strength = plot_data$Evidence_Strength[match(TARGET_TF, plot_data$TF)],
    ATF3_In_Top_N = any(plot_data$TF == TARGET_TF & plot_data$Method_Rank <= TOP_N_TF_TO_PLOT),
    PDF_File = output_files$pdf_file,
    PNG_File = output_files$png_file,
    stringsAsFactors = FALSE
  )
}

bubble_summary_list <- run_indexed_tasks_with_progress(
  total_tasks = nrow(bubble_tasks),
  workers = 1L,
  progress_label = "TF bubble plots",
  task_function = run_one_bubble_plot
)
stop_on_parallel_errors(bubble_summary_list, task_ids = seq_len(nrow(bubble_tasks)), label = "TF bubble plots")

bubble_summary <- do.call(rbind, bubble_summary_list)
bubble_summary_file <- write_csv_with_report_previews(
  bubble_summary,
  file.path(TABLE_ROOT, "tf_bubble_plot_summary.csv"),
  n_rows = 21
)


# 4. Plot Venn diagrams --------------------------------------------------------

cat("\nRunning DoRothEA + ChEA3 + CollecTRI Venn plots...\n")

run_one_venn_plot <- function(i) {
  analysis_name <- selected_analyses[i]
  method_tables <- method_tables_by_analysis[[analysis_name]]
  if (any(vapply(method_tables, is.null, logical(1)))) {
    return(data.frame(
      Analysis_Name = analysis_name,
      Status = "missing_method_table",
      stringsAsFactors = FALSE
    ))
  }

  sets <- lapply(method_tables, get_tf_set)
  names(sets) <- unname(TF_METHOD_LABELS[SELECTED_METHODS])
  counts <- get_venn_counts(sets)

  venn_plot <- make_venn_plot(sets, analysis_name)
  output_dir <- file.path(VENN_PLOT_ROOT, sanitize_file_name(analysis_name))
  pdf_file <- file.path(output_dir, "dorothea_chea3_collectri_venn.pdf")
  output_files <- save_ggplot_pdf_png(
    plot = venn_plot,
    pdf_file = pdf_file,
    width = VENN_PDF_WIDTH,
    height = VENN_PDF_HEIGHT
  )

  candidate_table <- rank_intersection_candidates(
    make_intersection_candidate_table(method_tables)
  )
  candidate_output_dir <- file.path(TABLE_ROOT, "intersection_candidates", sanitize_file_name(analysis_name))
  candidate_file <- write_csv_with_report_previews(
    candidate_table,
    file.path(candidate_output_dir, "dorothea_chea3_collectri_candidates.csv"),
    n_rows = 21
  )
  atf3_row <- candidate_table[candidate_table$TF == TARGET_TF, , drop = FALSE]

  data.frame(
    Analysis_Name = analysis_name,
    Status = "ok",
    DoRothEA_TF_Count = length(sets[["DoRothEA"]]),
    ChEA3_TF_Count = length(sets[["ChEA3"]]),
    CollecTRI_TF_Count = length(sets[["CollecTRI"]]),
    Three_Way_Intersection_Count = counts$Count[counts$Region == "ABC"],
    ATF3_In_Three_Way_Intersection = nrow(atf3_row) > 0L,
    ATF3_Consensus_Rank = if (nrow(atf3_row) > 0L) atf3_row$Consensus_Rank[1] else NA_real_,
    ATF3_Mean_Selected_Rank = if (nrow(atf3_row) > 0L) atf3_row$Mean_Selected_Rank[1] else NA_real_,
    Candidate_Table = candidate_file,
    PDF_File = output_files$pdf_file,
    PNG_File = output_files$png_file,
    stringsAsFactors = FALSE
  )
}

venn_summary_list <- run_indexed_tasks_with_progress(
  total_tasks = length(selected_analyses),
  workers = 1L,
  progress_label = "TF Venn plots",
  task_function = run_one_venn_plot
)
stop_on_parallel_errors(venn_summary_list, task_ids = selected_analyses, label = "TF Venn plots")

venn_summary <- do.call(rbind, venn_summary_list)
venn_summary_file <- write_csv_with_report_previews(
  venn_summary,
  file.path(TABLE_ROOT, "dorothea_chea3_collectri_venn_summary.csv"),
  n_rows = 21
)


# 5. Console summary -----------------------------------------------------------

cat("\nBubble plot summary table: ", bubble_summary_file, "\n", sep = "")
cat("Venn plot summary table:   ", venn_summary_file, "\n", sep = "")
cat("Bubble plot directory:     ", BUBBLE_PLOT_ROOT, "\n", sep = "")
cat("Venn plot directory:       ", VENN_PLOT_ROOT, "\n", sep = "")

cat("\nATF3 consensus ranks in DoRothEA + ChEA3 + CollecTRI intersection:\n")
print(
  venn_summary[
    ,
    c(
      "Analysis_Name",
      "Three_Way_Intersection_Count",
      "ATF3_In_Three_Way_Intersection",
      "ATF3_Consensus_Rank",
      "ATF3_Mean_Selected_Rank"
    )
  ],
  row.names = FALSE
)

cat("\nATF3-focused TF visualization finished.\n")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
