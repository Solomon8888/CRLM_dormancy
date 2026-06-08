# TCGAplot快速单基因与泛癌分析脚本
#
# 设计目的：
# 1. 先以TCGAplot包为主，快速完成ATF3的TCGA单癌种和泛癌扫描；
# 2. 所有正式结果统一保存到results/quickanalysis/tcga，不让TCGAplot把文件散落到项目根目录；
# 3. 下载或从包内提取的大型数据统一保存到data/tcgaplot；
# 4. TCGAplot运行时产生的中间文件统一放到temporary/tcgaplot；
# 5. 额外审计TCGAplot内置COAD样本和本项目本地SE对象是否一致。
#
# 当前默认设计：
# - 主分析基因：ATF3
# - 主分析癌种：COAD，不包含READ
# - 泛癌分析：使用TCGAplot内置33癌种pan-cancer数据，不受TARGET_CANCERS限制


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


# 1. 分析设计区：日常主要修改这里 ---------------------------------------------

# 单基因分析的目标基因。默认只分析ATF3。
# 临时换基因可在终端使用：
# TCGAPLOT_TARGET_GENES=ATF3,MYC Rscript scripts/quickanalysis/01_tcgaplot_quick_analysis.R
TARGET_GENES <- parse_env_vector("TCGAPLOT_TARGET_GENES", c("ATF3"))

# 单癌种分析的TCGA癌种。当前项目主线只分析COAD，不纳入READ。
# 注意：pan_boxplot、pan_forest、gene_TMB_radar等泛癌函数仍然自动使用TCGAplot内置全部癌种。
TARGET_CANCERS <- parse_env_vector("TCGAPLOT_TARGET_CANCERS", c("COAD"))

# 本地SE文件只用于样本一致性审计，不参与TCGAplot绘图运算。
LOCAL_TCGA_SE_FILES <- c(
  COAD = "data/TCGA/COAD/data_prepare/COAD_se_raw.rds",
  READ = "data/TCGA/READ/data_prepare/READ_se_raw.rds"
)

# gene_gene_scatter需要两个基因。默认先不运行该类图，避免把ATF3和任意基因硬配对。
# 如果后续需要看ATF3与免疫检查点或候选通路基因的相关性，在这里增加行即可。
GENE_PAIR_DESIGNS <- data.frame(
  Cancer = character(0),
  Gene1 = character(0),
  Gene2 = character(0),
  Density = character(0),
  stringsAsFactors = FALSE
)
# 示例：
# GENE_PAIR_DESIGNS <- data.frame(
#   Cancer = c("COAD", "COAD"),
#   Gene1 = c("ATF3", "ATF3"),
#   Gene2 = c("PDCD1", "CD274"),
#   Density = c("F", "F"),
#   stringsAsFactors = FALSE
# )

# GO/KEGG网络富集可以输入一个基因或一组基因。默认用TARGET_GENES，也就是ATF3。
NETWORK_GENES <- TARGET_GENES

# gene_network_go/gene_network_kegg本质是富集网络图，单个基因通常没有可用富集条目。
# 默认至少2个基因才运行；ATF3单基因分析中这两个任务会记录为skipped而不是failed。
NETWORK_MIN_GENES <- parse_env_integer("TCGAPLOT_NETWORK_MIN_GENES", 2L)

# 基因集分析配置。TCGAplot的gs_*函数需要一个MSigDB基因集名或用户自定义基因向量，
# 同时需要geneset_alias作为图中展示名称。默认不跑基因集分析。
TARGET_GENESETS <- data.frame(
  GeneSet = character(0),
  Alias = character(0),
  stringsAsFactors = FALSE
)
# Example:
# TARGET_GENESETS <- data.frame(
#   GeneSet = c("ALONSO_METASTASIS_EMT_DN"),
#   Alias = c("METASTASIS EMT"),
#   stringsAsFactors = FALSE
# )

# 只有这里为TRUE且TARGET_GENESETS不为空时，才会运行TCGAplot的gs_*基因集分析。
RUN_GENESET_ANALYSES <- parse_env_logical("TCGAPLOT_RUN_GENESETS", FALSE)

# 分析开关：
# - "literature_quick_scan"：默认值，按近年单基因TCGA文章常用模块快速跑一版；
# - "all"：运行本脚本封装的全部TCGAplot绘图/分析任务；
# - 也可以在终端指定逗号分隔的任务名，例如：
#   TCGAPLOT_ANALYSES=sample_audit,data_summary,pan_boxplot,tcga_boxplot,tcga_roc Rscript ...
ANALYSES_TO_RUN <- parse_env_vector(
  "TCGAPLOT_ANALYSES",
  "literature_quick_scan"
)

# TCGAplot通用参数。尽量集中在头部，方便快速调整。
CORRELATION_METHOD <- "pearson"
BOXPLOT_PALETTE <- "jco"
SURVIVAL_PALETTE <- "jco"
LEGEND_POSITION <- "none"
PAN_LEGEND_POSITION <- "right"
BOXPLOT_ADD_LAYER <- "jitter"
P_VALUE_LABEL <- "p.signif"
GROUP_TEST_METHOD <- "wilcox.test"
AGE_CUTOFF <- 60
AGE_YOUNG_CUTOFF <- 40
AGE_OLD_CUTOFF <- 60
TOP_N_GENES <- 20
HEATMAP_LOW_COLOR <- "#2166AC"
HEATMAP_HIGH_COLOR <- "#D73027"
HEATMAP_CLUSTER_ROW <- TRUE
HEATMAP_CLUSTER_COL <- TRUE
HEATMAP_LEGEND <- TRUE

# pan_forest的原包adjust=TRUE只校正age。
# 这里按TCGAplot内置meta的全部可用临床协变量逐一扩展：age、gender、stage。
PAN_FOREST_ADJUST_VARIABLES <- parse_env_vector(
  "TCGAPLOT_PAN_FOREST_ADJUST_VARIABLES",
  c("age", "gender", "stage")
)

# 输出目录。这里统一使用绝对路径，避免TCGAplot在任务临时目录运行时写错位置。
PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
SCRIPT_FILE <- file.path(PROJECT_ROOT, "scripts", "quickanalysis", "01_tcgaplot_quick_analysis.R")
DATASET_ID <- "tcga"
RESULT_ROOT <- file.path(PROJECT_ROOT, "results", "quickanalysis", DATASET_ID)
PLOT_ROOT <- file.path(RESULT_ROOT, "plots", "TCGAplot")
TABLE_ROOT <- file.path(RESULT_ROOT, "tables", "TCGAplot")
PLOT_PDF_DIR <- file.path(PLOT_ROOT, "pdf")
PLOT_PNG_DIR <- file.path(PLOT_ROOT, "png")
SUMMARY_ROOT <- TABLE_ROOT
DATA_ROOT <- file.path(PROJECT_ROOT, "data", "tcgaplot")
TEMP_ROOT <- file.path(PROJECT_ROOT, "temporary", "tcgaplot")
TCGAPLOT_REFERENCE_CACHE_ROOT <- file.path(DATA_ROOT, "reference_cache")
TCGAPLOT_TASK_CACHE_ROOT <- file.path(TCGAPLOT_REFERENCE_CACHE_ROOT, "task_manifest")
OMNIPATHR_CACHE_DIR <- file.path(TCGAPLOT_REFERENCE_CACHE_ROOT, "omnipathr_cache")
OMNIPATHR_LOG_DIR <- file.path(TEMP_ROOT, "omnipathr-log")

# 是否把TCGAplot内置大矩阵额外提取保存到data/tcgaplot。
# 默认FALSE，因为get_all_tpm/get_all_promoter_methy等对象较大；
# 如果需要保存用于后续复查，可在终端设置：
# TCGAPLOT_SAVE_BUILTIN_DATA=1 TCGAPLOT_ANALYSES=data_summary Rscript ...
SAVE_BUILTIN_DATA_EXTRACTS <- parse_env_logical(
  "TCGAPLOT_SAVE_BUILTIN_DATA",
  FALSE
)

# 默认在运行前清空本脚本前次运行产生的全部结果与中间文件。
# 这会删除results/quickanalysis/tcga、temporary/tcgaplot，
# 以及本脚本的任务manifest缓存；不会删除data/tcgaplot中保存的TCGAplot数据提取结果。
CLEAR_PREVIOUS_RUN_OUTPUTS <- parse_env_logical("TCGAPLOT_CLEAR_PREVIOUS_OUTPUTS", TRUE)

# 单个任务运行前清理自己的临时工作目录，避免TCGAplot内部写出的文件相互混入。
CLEAN_TASK_OUTPUT_DIR <- parse_env_logical("TCGAPLOT_CLEAN_OUTPUT", TRUE)

# 并行配置。默认使用parallel_runtime_functions.R自动识别的可用核心数；
# 如需限制CPU占用，可在终端设置：
# TCGAPLOT_PARALLEL_WORKERS=4 Rscript scripts/quickanalysis/01_tcgaplot_quick_analysis.R
MAX_PARALLEL_WORKERS <- parse_env_integer("TCGAPLOT_PARALLEL_WORKERS", NA_integer_)
QUICKANALYSIS_VERBOSE <- parse_env_logical("TCGAPLOT_VERBOSE", FALSE)
DISABLE_FORK_PARALLEL <- parse_env_logical("TCGAPLOT_DISABLE_FORK", interactive())
PARALLEL_BACKEND <- Sys.getenv("TCGAPLOT_PARALLEL_BACKEND", unset = "auto")

# 对网络/富集类TCGAplot任务启用任务输出缓存。
# 这不会改变正式结果目录，只是在data/tcgaplot/reference_cache下保存输出manifest，
# 同一参数重复运行且结果文件仍存在时可直接复用，避免重复跑耗时网络/富集图。
USE_TCGAPLOT_TASK_CACHE <- parse_env_logical("TCGAPLOT_USE_TASK_CACHE", TRUE)
TCGAPLOT_TASK_CACHE_MAX_AGE_DAYS <- parse_env_integer("TCGAPLOT_TASK_CACHE_MAX_AGE_DAYS", 30L)

# gene_gsea_kegg会经由clusterProfiler访问KEGG REST。
# 在macOS fork并行子进程中，KEGG下载偶尔会触发R进程级segfault；
# 这类任务默认放入独立Rscript进程运行，即使底层包崩溃也不会带崩主会话。
PROCESS_ISOLATED_ANALYSES <- parse_env_vector(
  "TCGAPLOT_PROCESS_ISOLATED_ANALYSES",
  c("gene_gsea_kegg")
)
DISABLE_PROCESS_ISOLATION <- parse_env_logical("TCGAPLOT_DISABLE_PROCESS_ISOLATION", FALSE)

options(width = 200)
options(bitmapType = "cairo")
options(quickanalysis_verbose = QUICKANALYSIS_VERBOSE)
options(parallel_runtime_force_single_line_progress = TRUE)
options(parallel_runtime_quiet_strategy = !QUICKANALYSIS_VERBOSE)
options(parallel_runtime_disable_fork = DISABLE_FORK_PARALLEL)
options(parallel_runtime_backend = PARALLEL_BACKEND)


# 2. 加载R包和项目共用函数 ----------------------------------------------------

# TCGAplot负责主要分析；SummarizedExperiment只用于读取本项目本地SE对象做样本审计。
# ggplot2/grid/Cairo用于把TCGAplot返回的图形对象统一保存成PDF和PNG。

required_packages <- c(
  "TCGAplot",
  "SummarizedExperiment",
  "ggplot2",
  "grid",
  "Cairo"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    function(package) {
      suppressPackageStartupMessages(
        suppressWarnings(requireNamespace(package, quietly = TRUE))
      )
    },
    logical(1)
  )
]
if (length(missing_packages) > 0) {
  stop(
    "Please install required R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  suppressWarnings(library(TCGAplot))
  library(SummarizedExperiment)
})

PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
TABLE_IO_FUNCTION_FILE <- "scripts/functions/result_table_io_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"
NETWORK_CACHE_FUNCTION_FILE <- "scripts/functions/network_cache_functions.R"

source(PLOTTING_FUNCTION_FILE)
source(TABLE_IO_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)
source(NETWORK_CACHE_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

PARALLEL_WORKERS <- if (is.na(MAX_PARALLEL_WORKERS)) {
  get_available_worker_count()
} else {
  max(1L, MAX_PARALLEL_WORKERS)
}


# 3. 通用工具函数 --------------------------------------------------------------

safe_name <- function(x, default = "analysis") {
  sanitize_file_name(x, default = default)
}

paste_compact <- function(x, collapse = "_", default = "all") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(default)
  }
  paste(x, collapse = collapse)
}

write_table <- function(dat, output_dir, file_stem) {
  # 调用项目通用表格输出函数，同时生成报告预览文件，保持和其他脚本一致。
  write_csv_with_report_previews(
    dat = dat,
    csv_file = file.path(output_dir, paste0(file_stem, ".csv")),
    na = "NA"
  )
}

format_task_timestamp <- function(time) {
  if (is.null(time) || is.na(time)) {
    return("")
  }

  format(time, "%Y-%m-%d %H:%M:%S")
}

get_tcgaplot_function <- function(function_name) {
  if (exists(function_name, envir = .GlobalEnv, mode = "function")) {
    return(get(function_name, envir = .GlobalEnv, mode = "function"))
  }

  if (!exists(function_name, envir = asNamespace("TCGAplot"), mode = "function")) {
    stop("TCGAplot function is unavailable: ", function_name)
  }

  get(function_name, envir = asNamespace("TCGAplot"), mode = "function")
}

add_title_to_ggplot_if_needed <- function(plot, title) {
  if (is.null(title) || !nzchar(title)) {
    return(plot)
  }

  if (inherits(plot, "patchwork") && requireNamespace("patchwork", quietly = TRUE)) {
    return(plot + patchwork::plot_annotation(title = title))
  }

  if (!inherits(plot, "ggplot")) {
    return(plot)
  }

  plot + ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

draw_outer_title <- function(title) {
  if (is.null(title) || !nzchar(title)) {
    return(invisible(FALSE))
  }

  grid::grid.text(
    label = title,
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.985, "npc"),
    gp = grid::gpar(fontface = "bold", fontsize = 13)
  )
  invisible(TRUE)
}

draw_in_title_viewport <- function(title, draw_fun) {
  if (is.null(title) || !nzchar(title)) {
    draw_fun()
    return(invisible(NULL))
  }

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.0, "npc"),
    width = grid::unit(1, "npc"),
    height = grid::unit(0.92, "npc"),
    just = c("center", "bottom")
  ))
  on.exit({
    try(grid::popViewport(), silent = TRUE)
  }, add = TRUE)
  draw_fun()
  grid::popViewport()
  draw_outer_title(title)
  invisible(NULL)
}

render_tcgaplot_output <- function(result, title = NULL, add_outer_title = FALSE) {
  # TCGAplot不同函数返回值类型不完全一致：
  # 有的直接返回ggplot，有的返回survminer对象，有的返回ComplexHeatmap或grob。
  # 这里集中做一次类型分派，后续保存图片时只需要调用这一层即可。
  if (is.null(result)) {
    if (add_outer_title) {
      draw_outer_title(title)
    }
    return(invisible(NULL))
  }

  if (inherits(result, "ggplot") || inherits(result, "patchwork")) {
    print(add_title_to_ggplot_if_needed(result, title))
    return(invisible(NULL))
  }

  if (inherits(result, "ggsurvplot")) {
    if (!is.null(title) && nzchar(title) && inherits(result$plot, "ggplot")) {
      result$plot <- add_title_to_ggplot_if_needed(result$plot, title)
    }
    print(result)
    return(invisible(NULL))
  }

  if (inherits(result, "Heatmap") && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    if (add_outer_title) {
      draw_in_title_viewport(title, function() {
        ComplexHeatmap::draw(result, newpage = FALSE)
      })
    } else {
      ComplexHeatmap::draw(result)
    }
    return(invisible(NULL))
  }

  if (inherits(result, "gforge_forestplot")) {
    print(result)
    return(invisible(NULL))
  }

  if (inherits(result, "grob") || inherits(result, "gTree") || inherits(result, "gtable")) {
    if (add_outer_title) {
      draw_in_title_viewport(title, function() {
        grid::grid.draw(result)
      })
    } else {
      grid::grid.draw(result)
    }
    return(invisible(NULL))
  }

  if (is.list(result) && "plot" %in% names(result)) {
    render_tcgaplot_output(result$plot, title = title, add_outer_title = add_outer_title)
    return(invisible(NULL))
  }

  try(print(result), silent = TRUE)
  if (add_outer_title) {
    draw_outer_title(title)
  }
  invisible(NULL)
}

get_pdf_page_count <- function(pdf_file) {
  if (!file.exists(pdf_file) || is.na(file.info(pdf_file)$size) || file.info(pdf_file)$size <= 0) {
    return(NA_integer_)
  }

  if (requireNamespace("pdftools", quietly = TRUE)) {
    page_count <- tryCatch(pdftools::pdf_info(pdf_file)$pages, error = function(error) NA_integer_)
    if (!is.na(page_count)) {
      return(as.integer(page_count))
    }
  }

  pdfinfo <- Sys.which("pdfinfo")
  if (!nzchar(pdfinfo)) {
    return(NA_integer_)
  }

  info <- tryCatch(
    suppressWarnings(system2(pdfinfo, args = pdf_file, stdout = TRUE, stderr = TRUE)),
    error = function(error) character(0)
  )
  page_line <- grep("^Pages:", info, value = TRUE)
  if (length(page_line) == 0L) {
    return(NA_integer_)
  }

  suppressWarnings(as.integer(trimws(sub("^Pages:", "", page_line[1]))))
}

get_png_nonwhite_fraction <- function(png_file, white_tolerance = 0.02) {
  if (!file.exists(png_file) || !requireNamespace("png", quietly = TRUE)) {
    return(NA_real_)
  }

  img <- tryCatch(png::readPNG(png_file), error = function(error) NULL)
  if (is.null(img)) {
    return(NA_real_)
  }

  dims <- dim(img)
  if (length(dims) == 2L) {
    return(mean(abs(img - 1) > white_tolerance))
  }

  channels <- min(3L, dims[3])
  nonwhite <- matrix(FALSE, nrow = dims[1], ncol = dims[2])
  for (channel in seq_len(channels)) {
    nonwhite <- nonwhite | abs(img[, , channel] - 1) > white_tolerance
  }
  if (dims[3] >= 4L) {
    nonwhite <- nonwhite & img[, , 4] > 0.01
  }

  mean(nonwhite)
}

is_png_blank <- function(png_file, nonwhite_fraction_cutoff = 0.001) {
  nonwhite_fraction <- get_png_nonwhite_fraction(png_file)
  is.na(nonwhite_fraction) || nonwhite_fraction < nonwhite_fraction_cutoff
}

render_one_pdf_page_png <- function(pdf_file, page, png_file, dpi = PNG_DPI) {
  dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)
  unlink(png_file, force = TRUE)

  if (requireNamespace("pdftools", quietly = TRUE) && requireNamespace("png", quietly = TRUE)) {
    converted <- tryCatch(
      {
        bitmap <- pdftools::pdf_render_page(pdf_file, page = page, dpi = dpi)
        png::writePNG(bitmap, target = png_file)
        file.exists(png_file)
      },
      error = function(error) FALSE
    )
    if (converted) {
      return(png_file)
    }
  }

  pdftoppm <- Sys.which("pdftoppm")
  if (nzchar(pdftoppm)) {
    tmp_dir <- tempfile("tcgaplot_pdftoppm_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
    tmp_prefix <- file.path(tmp_dir, "page")
    status <- tryCatch(
      suppressWarnings(system2(
        pdftoppm,
        args = c("-png", "-r", as.character(dpi), "-f", as.character(page), "-l", as.character(page), pdf_file, tmp_prefix),
        stdout = TRUE,
        stderr = TRUE
      )),
      error = function(error) structure(character(0), status = 1L)
    )
    rendered <- sort(list.files(tmp_dir, pattern = "^page-.*[.]png$", full.names = TRUE))
    if (length(rendered) > 0L && is.null(attr(status, "status"))) {
      file.copy(rendered[1], png_file, overwrite = TRUE)
      return(png_file)
    }
  }

  sips <- Sys.which("sips")
  if (nzchar(sips) && identical(as.integer(page), 1L)) {
    status <- tryCatch(
      suppressWarnings(system2(
        sips,
        args = c("-s", "format", "png", pdf_file, "--out", png_file),
        stdout = TRUE,
        stderr = TRUE
      )),
      error = function(error) structure(character(0), status = 1L)
    )
    if (file.exists(png_file) && is.null(attr(status, "status"))) {
      return(png_file)
    }
  }

  character(0)
}

detect_nonblank_pdf_pages <- function(pdf_file, page_count = get_pdf_page_count(pdf_file)) {
  if (is.na(page_count) || page_count < 1L) {
    return(integer(0))
  }

  nonblank_pages <- integer(0)
  tmp_dir <- tempfile("tcgaplot_pdf_page_check_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  for (page in seq_len(page_count)) {
    page_png <- file.path(tmp_dir, sprintf("page_%03d.png", page))
    rendered <- render_one_pdf_page_png(pdf_file, page = page, png_file = page_png, dpi = 72)
    if (length(rendered) > 0L && !is_png_blank(page_png)) {
      nonblank_pages <- c(nonblank_pages, page)
    }
  }

  nonblank_pages
}

remove_blank_pdf_pages <- function(pdf_file) {
  page_count <- get_pdf_page_count(pdf_file)
  if (is.na(page_count)) {
    stop("PDF is missing, empty, or unreadable: ", pdf_file)
  }
  if (page_count <= 1L) {
    return(invisible(page_count))
  }

  nonblank_pages <- detect_nonblank_pdf_pages(pdf_file, page_count)
  if (length(nonblank_pages) == 0L) {
    stop("All PDF pages appear blank: ", pdf_file)
  }
  if (identical(nonblank_pages, seq_len(page_count))) {
    return(invisible(page_count))
  }

  pdfseparate <- Sys.which("pdfseparate")
  pdfunite <- Sys.which("pdfunite")
  if (!nzchar(pdfseparate) || !nzchar(pdfunite)) {
    warning("pdfseparate/pdfunite are unavailable; blank PDF pages cannot be removed: ", pdf_file)
    return(invisible(page_count))
  }

  tmp_dir <- tempfile("tcgaplot_pdf_clean_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  page_pattern <- file.path(tmp_dir, "page_%03d.pdf")
  split_status <- tryCatch(
    suppressWarnings(system2(pdfseparate, args = c(pdf_file, page_pattern), stdout = TRUE, stderr = TRUE)),
    error = function(error) structure(character(0), status = 1L)
  )
  if (!is.null(attr(split_status, "status"))) {
    warning("Failed to split PDF for blank-page cleanup: ", pdf_file)
    return(invisible(page_count))
  }

  split_files <- file.path(tmp_dir, sprintf("page_%03d.pdf", nonblank_pages))
  clean_pdf <- file.path(tmp_dir, "clean.pdf")
  unite_status <- tryCatch(
    suppressWarnings(system2(pdfunite, args = c(split_files, clean_pdf), stdout = TRUE, stderr = TRUE)),
    error = function(error) structure(character(0), status = 1L)
  )
  if (!is.null(attr(unite_status, "status")) || !file.exists(clean_pdf) || file.info(clean_pdf)$size <= 0) {
    warning("Failed to rebuild PDF without blank pages: ", pdf_file)
    return(invisible(page_count))
  }

  file.copy(clean_pdf, pdf_file, overwrite = TRUE)
  invisible(length(nonblank_pages))
}

convert_pdf_to_png_outputs <- function(pdf_file, png_file, dpi = PNG_DPI) {
  remove_blank_pdf_pages(pdf_file)
  page_count <- get_pdf_page_count(pdf_file)
  if (is.na(page_count) || page_count < 1L) {
    stop("PDF is missing, empty, or unreadable after cleanup: ", pdf_file)
  }

  output_stem <- tools::file_path_sans_ext(basename(png_file))
  if (page_count == 1L) {
    rendered <- render_one_pdf_page_png(pdf_file, page = 1L, png_file = png_file, dpi = dpi)
    if (length(rendered) == 0L || is_png_blank(png_file)) {
      stop("PNG preview is blank or could not be rendered: ", png_file)
    }
    return(png_file)
  }

  output_dir <- file.path(dirname(png_file), output_stem)
  unlink(output_dir, recursive = TRUE, force = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  png_files <- character(0)
  for (page in seq_len(page_count)) {
    page_png <- file.path(output_dir, sprintf("%s_page_%03d.png", output_stem, page))
    rendered <- render_one_pdf_page_png(pdf_file, page = page, png_file = page_png, dpi = dpi)
    if (length(rendered) == 0L || is_png_blank(page_png)) {
      stop("Rendered PDF page is blank or unavailable: ", pdf_file, " page ", page)
    }
    png_files <- c(png_files, page_png)
  }

  png_files
}

render_pdf_preview_png <- function(pdf_file, png_file, dpi = PNG_DPI) {
  # TCGAplot少数函数只在当前目录直接写PDF，不返回可重绘对象。
  # 这里统一清理空白页并转PNG；多页PDF会拆成同名子目录中的page_001/page_002...
  tryCatch(
    convert_pdf_to_png_outputs(pdf_file = pdf_file, png_file = png_file, dpi = dpi),
    error = function(error) {
      warning(conditionMessage(error))
      character(0)
    }
  )
}

with_task_workspace <- function(task, expr) {
  # 很多TCGAplot函数会在当前工作目录写PDF或中间文件。
  # 因此每个任务都进入独立临时目录运行，再把生成文件复制到规范结果目录。
  task_temp_dir <- file.path(
    TEMP_ROOT,
    "tasks",
    safe_name(task$analysis),
    safe_name(task$target),
    safe_name(task$context)
  )
  dir.create(task_temp_dir, recursive = TRUE, showWarnings = FALSE)
  if (CLEAN_TASK_OUTPUT_DIR) {
    unlink(list.files(task_temp_dir, all.files = FALSE, full.names = TRUE, recursive = FALSE))
  }

  old_wd <- getwd()
  setwd(task_temp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  force(expr)
}

copy_generated_files <- function(task) {
  # 对于TCGAplot内部直接写出的文件，按扩展名分流：
  # 图片进入plots/TCGAplot/pdf或plots/TCGAplot/png，表格进入tables/TCGAplot根目录。
  generated_files <- list.files(
    ".",
    all.files = FALSE,
    recursive = FALSE,
    full.names = TRUE
  )
  generated_files <- generated_files[file.info(generated_files)$isdir == FALSE]

  if (length(generated_files) == 0) {
    return(character(0))
  }

  copied_files <- character(0)

  for (source_file in generated_files) {
    extension <- tolower(tools::file_ext(source_file))
    file_stem <- tools::file_path_sans_ext(basename(source_file))
    output_stem <- task_file_stem(task, file_stem)

    if (extension == "pdf") {
      destination_dir <- PLOT_PDF_DIR
      destination_file <- file.path(destination_dir, paste0(output_stem, ".", extension))
    } else if (extension == "png") {
      destination_dir <- PLOT_PNG_DIR
      destination_file <- file.path(destination_dir, paste0(output_stem, ".", extension))
    } else if (extension %in% c("jpg", "jpeg")) {
      destination_dir <- PLOT_PNG_DIR
      destination_file <- file.path(destination_dir, paste0(output_stem, ".", extension))
    } else {
      destination_dir <- TABLE_ROOT
      destination_file <- file.path(destination_dir, paste0(output_stem, ".", extension))
    }

    dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(source_file, destination_file, overwrite = TRUE)
    copied_files <- c(copied_files, destination_file)

    if (extension == "pdf") {
      png_file <- file.path(PLOT_PNG_DIR, paste0(output_stem, ".png"))
      rendered_png_files <- render_pdf_preview_png(destination_file, png_file)
      if (length(rendered_png_files) > 0L) {
        copied_files <- c(copied_files, rendered_png_files)
      }
    }
  }

  copied_files
}

save_tcgaplot_pdf_png <- function(pdf_file, width, height, draw_fun, png_dpi = PNG_DPI) {
  # TCGAplot的部分函数绘图时有副作用，甚至会联网下载数据。
  # 因此这里只绘制一次PDF，再由PDF生成PNG，避免同一个分析重复运行两遍。
  plot_paths <- get_plot_output_paths(pdf_file)
  dir.create(dirname(plot_paths$pdf), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_paths$png), recursive = TRUE, showWarnings = FALSE)
  unlink(c(plot_paths$pdf, plot_paths$png), force = TRUE)
  unlink(file.path(dirname(plot_paths$png), tools::file_path_sans_ext(basename(plot_paths$png))), recursive = TRUE, force = TRUE)

  device_opened <- FALSE
  Cairo::CairoPDF(
    file = plot_paths$pdf,
    width = width,
    height = height,
    bg = "white"
  )
  device_opened <- TRUE
  on.exit({
    if (device_opened && grDevices::dev.cur() > 1L) {
      try(grDevices::dev.off(), silent = TRUE)
    }
  }, add = TRUE)

  draw_fun()
  invisible(grDevices::dev.off())
  device_opened <- FALSE

  if (!file.exists(plot_paths$pdf) || is.na(file.info(plot_paths$pdf)$size) || file.info(plot_paths$pdf)$size <= 0) {
    stop("PDF output is missing or empty: ", plot_paths$pdf)
  }

  png_files <- convert_pdf_to_png_outputs(
    pdf_file = plot_paths$pdf,
    png_file = plot_paths$png,
    dpi = png_dpi
  )

  invisible(list(pdf_file = plot_paths$pdf, png_file = png_files))
}

validate_plot_output_files <- function(files) {
  files <- unique(files[nzchar(files)])
  if (length(files) == 0L) {
    return(invisible(TRUE))
  }

  invalid <- character(0)
  for (file in files) {
    if (!file.exists(file) || is.na(file.info(file)$size) || file.info(file)$size <= 0) {
      invalid <- c(invalid, paste0(file, " [missing_or_empty]"))
      next
    }

    extension <- tolower(tools::file_ext(file))
    if (extension == "pdf" && is.na(get_pdf_page_count(file))) {
      invalid <- c(invalid, paste0(file, " [unreadable_pdf]"))
    }
    if (extension == "png" && is_png_blank(file)) {
      invalid <- c(invalid, paste0(file, " [blank_png]"))
    }
  }

  if (length(invalid) > 0L) {
    stop("Invalid plot output detected: ", paste(invalid, collapse = "; "))
  }

  invisible(TRUE)
}

cleanup_task_output_files <- function(task) {
  output_stem <- task_file_stem(task)
  cleanup_dir <- function(directory) {
    if (!dir.exists(directory)) {
      return(invisible(FALSE))
    }
    entries <- list.files(directory, all.files = FALSE, full.names = TRUE, recursive = FALSE)
    entries <- entries[startsWith(basename(entries), output_stem)]
    unlink(entries, recursive = TRUE, force = TRUE)
    invisible(TRUE)
  }

  cleanup_dir(PLOT_PDF_DIR)
  cleanup_dir(PLOT_PNG_DIR)
  cleanup_dir(TABLE_ROOT)
  invisible(TRUE)
}

capture_task <- function(expr) {
  # 不让单个任务的warning/message刷屏；统一收集到tcgaplot_task_summary.csv里复查。
  warnings <- character(0)
  messages <- character(0)

  result <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      },
      message = function(message) {
        messages <<- c(messages, conditionMessage(message))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(error) error
  )

  list(
    result = result,
    warnings = unique(warnings),
    messages = unique(messages)
  )
}

make_task_result <- function(
    task,
    status,
    output_files = character(0),
    generated_files = character(0),
    warnings = character(0),
    messages = character(0),
    error = "",
    start_time = NULL,
    end_time = Sys.time(),
    cache_file = "") {
  runtime_seconds <- if (is.null(start_time)) {
    NA_real_
  } else {
    as.numeric(difftime(end_time, start_time, units = "secs"))
  }

  data.frame(
    Analysis = task$analysis,
    Function = task$function_name,
    Target = task$target,
    Context = task$context,
    Status = status,
    Start_Time = format_task_timestamp(start_time),
    End_Time = format_task_timestamp(end_time),
    Runtime_Seconds = runtime_seconds,
    Runtime = ifelse(is.na(runtime_seconds), "", format_runtime_seconds(runtime_seconds)),
    Output_Files = paste(output_files, collapse = ";"),
    Generated_Files = paste(generated_files, collapse = ";"),
    Cache_File = cache_file,
    Warning = paste(warnings, collapse = " | "),
    Message = paste(messages, collapse = " | "),
    Error = error,
    stringsAsFactors = FALSE
  )
}

make_special_task_result <- function(
    task_name,
    status,
    output_files = character(0),
    rows = NA_integer_,
    warnings = character(0),
    messages = character(0),
    error = "",
    start_time = NULL,
    end_time = Sys.time()) {
  runtime_seconds <- if (is.null(start_time)) {
    NA_real_
  } else {
    as.numeric(difftime(end_time, start_time, units = "secs"))
  }

  data.frame(
    Analysis = task_name,
    Status = status,
    Start_Time = format_task_timestamp(start_time),
    End_Time = format_task_timestamp(end_time),
    Runtime_Seconds = runtime_seconds,
    Runtime = ifelse(is.na(runtime_seconds), "", format_runtime_seconds(runtime_seconds)),
    Rows = rows,
    Output_Files = paste(output_files, collapse = ";"),
    Warning = paste(warnings, collapse = " | "),
    Message = paste(messages, collapse = " | "),
    Error = error,
    stringsAsFactors = FALSE
  )
}

CACHEABLE_TCGAPLOT_ANALYSES <- c(
  "gene_network_go",
  "gene_network_kegg",
  "gene_gsea_go",
  "gene_gsea_kegg",
  "gene_coexp_heatmap"
)

split_manifest_files <- function(x) {
  x <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}

manifest_files_exist <- function(manifest) {
  if (!is.data.frame(manifest) || nrow(manifest) == 0L) {
    return(FALSE)
  }

  output_files <- split_manifest_files(manifest$Output_Files[1])
  generated_files <- split_manifest_files(manifest$Generated_Files[1])
  files <- unique(c(output_files, generated_files))
  length(files) > 0L && all(file.exists(files))
}

is_cacheable_tcgaplot_task <- function(task) {
  USE_TCGAPLOT_TASK_CACHE && task$analysis %in% CACHEABLE_TCGAPLOT_ANALYSES
}

get_tcgaplot_task_skip_reason <- function(task) {
  if (
    task$analysis %in% c("gene_network_go", "gene_network_kegg") &&
      length(unique(task$args$gene)) < NETWORK_MIN_GENES
  ) {
    return(paste0(
      "Skipped because ",
      task$analysis,
      " requires at least ",
      NETWORK_MIN_GENES,
      " genes for a stable enrichment network; current NETWORK_GENES=",
      paste(task$args$gene, collapse = ",")
    ))
  }

  ""
}

get_tcgaplot_task_cache_file <- function(task) {
  file.path(
    TCGAPLOT_TASK_CACHE_ROOT,
    safe_name(task$analysis),
    safe_name(task$target),
    paste0(safe_name(task$context), ".rds")
  )
}

make_tcgaplot_task_cache_metadata <- function(task) {
  make_resource_cache_metadata(
    resource_name = paste0("TCGAplot::", task$function_name),
    species = "human",
    extra = list(
      package_version = as.character(utils::packageVersion("TCGAplot")),
      output_layout = "flat_pdf_png_tables_v1",
      task_id = task$task_id %||% NA_integer_,
      analysis = task$analysis,
      target = task$target,
      context = task$context,
      renderer = task$renderer,
      width = task$width,
      height = task$height,
      args = task$args
    )
  )
}

read_tcgaplot_task_cache <- function(task) {
  if (!is_cacheable_tcgaplot_task(task)) {
    return(NULL)
  }

  cache_file <- get_tcgaplot_task_cache_file(task)
  metadata <- make_tcgaplot_task_cache_metadata(task)
  cached <- read_reference_cache(
    cache_file = cache_file,
    expected_metadata = metadata,
    max_age_days = TCGAPLOT_TASK_CACHE_MAX_AGE_DAYS,
    use_cache = USE_TCGAPLOT_TASK_CACHE
  )

  if (!cached$found || !manifest_files_exist(cached$result)) {
    return(NULL)
  }

  cached$result$Status <- "cached"
  cached$result$Cache_File <- cache_file
  cached$result
}

write_tcgaplot_task_cache <- function(task, manifest) {
  if (!is_cacheable_tcgaplot_task(task) || !identical(manifest$Status[1], "success")) {
    return(invisible(FALSE))
  }

  cache_file <- get_tcgaplot_task_cache_file(task)
  metadata <- make_tcgaplot_task_cache_metadata(task)
  write_reference_cache(
    cache_file = cache_file,
    metadata = metadata,
    result = manifest
  )
  invisible(TRUE)
}

get_tcgaplot_data_object <- function(name) {
  get(name, envir = asNamespace("TCGAplot"), inherits = FALSE)
}

make_pan_forest_exprset <- function(cancer, gene) {
  tpm_dat <- get_tcgaplot_data_object("tpm")
  meta_dat <- get_tcgaplot_data_object("meta")

  expr_set <- subset(tpm_dat, Group == "Tumor" & Cancer == cancer) %>%
    tibble::add_column(ID = stringr::str_sub(rownames(.), 1, 12), .before = "Cancer") %>%
    dplyr::filter(!duplicated(ID)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames("ID") %>%
    dplyr::filter(rownames(.) %in% rownames(subset(meta_dat, Cancer == cancer)))

  expr_set <- expr_set[, -(1:2), drop = FALSE]
  as.matrix(t(expr_set))
}

fit_symbol_cox_for_forest <- function(cancer, gene, adjust_variable = NULL) {
  meta_dat <- get_tcgaplot_data_object("meta")
  expr_set <- make_pan_forest_exprset(cancer, gene)
  cl <- meta_dat[colnames(expr_set), , drop = FALSE]
  cl$symbol <- expr_set[gene, ]

  formula_text <- if (is.null(adjust_variable) || !nzchar(adjust_variable)) {
    "survival::Surv(time, event) ~ symbol"
  } else {
    paste0("survival::Surv(time, event) ~ symbol + ", adjust_variable)
  }

  cl <- cl[, unique(c("time", "event", "symbol", adjust_variable)), drop = FALSE]
  cl <- stats::na.omit(cl)
  if (nrow(cl) < 10L || length(unique(cl$event)) < 2L) {
    stop("Insufficient survival records after omitting missing clinical values")
  }

  if (!is.null(adjust_variable) && !is.numeric(cl[[adjust_variable]])) {
    cl[[adjust_variable]] <- as.factor(cl[[adjust_variable]])
    if (nlevels(cl[[adjust_variable]]) < 2L) {
      stop("Clinical covariate has fewer than 2 levels")
    }
  }

  model <- survival::coxph(stats::as.formula(formula_text), data = cl)
  beta <- coef(model)
  se <- sqrt(diag(vcov(model)))
  hr <- exp(beta)
  hrse <- hr * se
  tmp <- round(cbind(
    coef = beta,
    se = se,
    z = beta / se,
    p = 1 - stats::pchisq((beta / se)^2, 1),
    HR = hr,
    HRse = hrse,
    HRz = (hr - 1) / hrse,
    HRp = 1 - stats::pchisq(((hr - 1) / hrse)^2, 1),
    HRCILL = exp(beta - stats::qnorm(0.975, 0, 1) * se),
    HRCIUL = exp(beta + stats::qnorm(0.975, 0, 1) * se)
  ), 3)

  tmp["symbol", ]
}

pan_forest_clinical_adjusted <- function(gene, adjust_variable) {
  cancers <- get_tcgaplot_data_object("cancers")
  cox_results <- list()
  failed <- character(0)

  for (cancer in cancers) {
    fit <- tryCatch(
      fit_symbol_cox_for_forest(
        cancer = cancer,
        gene = gene,
        adjust_variable = adjust_variable
      ),
      error = function(error) {
        failed[[cancer]] <<- conditionMessage(error)
        NULL
      }
    )
    if (!is.null(fit)) {
      cox_results[[cancer]] <- fit
    }
  }

  if (length(cox_results) == 0L) {
    stop("No cancer type could be fitted for ", gene, " adjusted by ", adjust_variable)
  }

  cox_results <- do.call(rbind, cox_results)
  cox_results <- as.data.frame(cox_results[, c(5, 9:10, 4), drop = FALSE])
  np <- paste0(cox_results$HR, " (", cox_results$HRCILL, "-", cox_results$HRCIUL, ")")
  tabletext <- cbind(
    c("Cancer", rownames(cox_results)),
    c("HR (95%CI)", np),
    c("P Value", cox_results$p)
  )

  hrzl_line_last <- nrow(cox_results) + 2L
  hrzl_lines <- list(
    grid::gpar(lwd = 2, col = "black"),
    grid::gpar(lwd = 2, col = "black"),
    grid::gpar(lwd = 2, col = "black")
  )
  names(hrzl_lines) <- as.character(c(1L, 2L, hrzl_line_last))

  forestplot::forestplot(
    labeltext = tabletext,
    graph.pos = 3,
    mean = c(NA, cox_results$HR),
    lower = c(NA, cox_results$HRCILL),
    upper = c(NA, cox_results$HRCIUL),
    title = paste0("Hazard Ratio Plot of ", gene, " adjusted by ", adjust_variable),
    hrzl_lines = hrzl_lines,
    is.summary = c(TRUE, rep(FALSE, nrow(cox_results))),
    col = forestplot::fpColors(box = "#1c61b6", lines = "#1c61b6", zero = "gray50"),
    zero = 1,
    cex = 0.9,
    lineheight = "auto",
    colgap = grid::unit(8, "mm"),
    txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 1)),
    boxsize = 0.5,
    ci.vertices = TRUE,
    ci.vertices.height = 0.3
  )
}


# 4. 样本一致性审计 ------------------------------------------------------------

parse_tcga_barcode <- function(sample_ids) {
  # TCGA条形码关键位点：
  # 前12位是患者ID；第14-15位是样本类型；第16位是vial。
  # COAD差异常见来源是TCGAplot保留01B/01C，而本地SE通常只保留01A。
  sample_ids <- as.character(sample_ids)
  data.frame(
    Sample_ID = sample_ids,
    Patient_ID = substr(sample_ids, 1, 12),
    Sample_Type_Code = substr(sample_ids, 14, 15),
    Vial = substr(sample_ids, 16, 16),
    Short_Sample = substr(sample_ids, 1, 16),
    stringsAsFactors = FALSE
  )
}

explain_sample_difference <- function(status, sample_type_code, vial) {
  if (status == "shared") {
    return("shared_by_exact_sample_id")
  }

  if (sample_type_code == "01" && vial != "A") {
    return("additional_primary_tumor_vial_non_01A; local_SE_keeps_01A_preferred_sample")
  }

  if (sample_type_code == "11" && vial != "A") {
    return("additional_solid_normal_vial_non_11A; local_SE_keeps_11A_preferred_sample")
  }

  if (sample_type_code == "06") {
    return("metastatic_sample")
  }

  if (sample_type_code == "02") {
    return("recurrent_tumor_sample")
  }

  "unclassified_difference"
}

audit_one_cancer_samples <- function(cancer) {
  # 该审计不会修改任何SE对象，只比较样本ID集合并输出差异原因。
  if (!cancer %in% names(LOCAL_TCGA_SE_FILES)) {
    return(data.frame(
      Cancer = cancer,
      Status = "skipped",
      Reason = "No local SE file configured for this cancer",
      stringsAsFactors = FALSE
    ))
  }

  se_file <- LOCAL_TCGA_SE_FILES[[cancer]]
  if (!file.exists(se_file)) {
    return(data.frame(
      Cancer = cancer,
      Status = "skipped",
      Reason = paste("Local SE file does not exist:", se_file),
      stringsAsFactors = FALSE
    ))
  }

  tcgaplot_tpm <- get_tpm(cancer)
  local_se <- readRDS(se_file)

  tcgaplot_samples <- parse_tcga_barcode(rownames(tcgaplot_tpm))
  tcgaplot_samples$Source <- "TCGAplot"
  tcgaplot_samples$Group <- as.character(tcgaplot_tpm$Group)

  local_coldata <- as.data.frame(colData(local_se), stringsAsFactors = FALSE)
  local_samples <- parse_tcga_barcode(colnames(local_se))
  local_samples$Source <- "local_SE"
  local_samples$Group <- as.character(local_coldata$group)
  local_samples$Group_Detail <- as.character(local_coldata$group_detail)
  local_samples$Group_Paired <- as.character(local_coldata$group_paired)

  all_sample_ids <- union(tcgaplot_samples$Sample_ID, local_samples$Sample_ID)
  comparison <- data.frame(Sample_ID = all_sample_ids, stringsAsFactors = FALSE)
  comparison <- merge(
    comparison,
    tcgaplot_samples,
    by = "Sample_ID",
    all.x = TRUE,
    suffixes = c("", "_TCGAplot")
  )
  comparison <- merge(
    comparison,
    local_samples,
    by = "Sample_ID",
    all.x = TRUE,
    suffixes = c("_TCGAplot", "_Local")
  )

  comparison$In_TCGAplot <- comparison$Sample_ID %in% tcgaplot_samples$Sample_ID
  comparison$In_Local_SE <- comparison$Sample_ID %in% local_samples$Sample_ID
  comparison$Comparison_Status <- ifelse(
    comparison$In_TCGAplot & comparison$In_Local_SE,
    "shared",
    ifelse(comparison$In_TCGAplot, "tcgaplot_only", "local_only")
  )

  barcode <- parse_tcga_barcode(comparison$Sample_ID)
  comparison$Patient_ID <- barcode$Patient_ID
  comparison$Sample_Type_Code <- barcode$Sample_Type_Code
  comparison$Vial <- barcode$Vial
  comparison$Likely_Reason <- vapply(
    seq_len(nrow(comparison)),
    function(i) explain_sample_difference(
      status = comparison$Comparison_Status[i],
      sample_type_code = comparison$Sample_Type_Code[i],
      vial = comparison$Vial[i]
    ),
    character(1)
  )

  comparison <- comparison[order(comparison$Comparison_Status, comparison$Sample_ID), ]
  rownames(comparison) <- NULL

  cancer_stem <- paste("000", "sample_audit", safe_name(cancer), sep = "_")
  write_table(comparison, TABLE_ROOT, paste(cancer_stem, "sample_id_comparison", sep = "_"))
  write_table(tcgaplot_samples, TABLE_ROOT, paste(cancer_stem, "tcgaplot_samples", sep = "_"))
  write_table(local_samples, TABLE_ROOT, paste(cancer_stem, "local_se_samples", sep = "_"))

  summary_table <- data.frame(
    Cancer = cancer,
    TCGAplot_Samples = nrow(tcgaplot_samples),
    Local_SE_Samples = nrow(local_samples),
    Shared_Exact_Sample_IDs = sum(comparison$Comparison_Status == "shared"),
    TCGAplot_Only_Samples = sum(comparison$Comparison_Status == "tcgaplot_only"),
    Local_SE_Only_Samples = sum(comparison$Comparison_Status == "local_only"),
    TCGAplot_Tumor = sum(tcgaplot_samples$Group == "Tumor", na.rm = TRUE),
    TCGAplot_Normal = sum(tcgaplot_samples$Group == "Normal", na.rm = TRUE),
    Local_Tumor_Strict_01A = sum(local_samples$Group == "tumor", na.rm = TRUE),
    Local_Normal_Strict_11A = sum(local_samples$Group == "normal", na.rm = TRUE),
    Local_Non_01A_Tumor_Or_Other = sum(
      local_samples$Group %in% c(
        "tumor_primary_non_a",
        "tumor_metastatic",
        "tumor_recurrent_solid"
      ),
      na.rm = TRUE
    ),
    Main_Diagnosis = paste(
      unique(comparison$Likely_Reason[comparison$Comparison_Status != "shared"]),
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )

  write_table(summary_table, TABLE_ROOT, paste(cancer_stem, "summary", sep = "_"))
  summary_table
}

run_sample_audit <- function() {
  audit_summaries <- do.call(
    rbind,
    lapply(TARGET_CANCERS, audit_one_cancer_samples)
  )
  write_table(audit_summaries, SUMMARY_ROOT, "000_sample_audit_summary")
  audit_summaries
}


# 5. TCGAplot内置数据概览与可选提取 -------------------------------------------

summarize_saved_object <- function(label, object, file) {
  data.frame(
    Object = label,
    Class = paste(class(object), collapse = ";"),
    Rows = ifelse(is.null(dim(object)), NA_integer_, dim(object)[1]),
    Columns = ifelse(is.null(dim(object)), NA_integer_, dim(object)[2]),
    Length = length(object),
    File = file,
    Status = "saved",
    Error = "",
    stringsAsFactors = FALSE
  )
}

save_builtin_extract <- function(label, file, expr) {
  # 只有TCGAPLOT_SAVE_BUILTIN_DATA=1时才会调用。
  # 所有提取对象保存到data/tcgaplot，不作为正式结果图表。
  tryCatch(
    {
      object <- force(expr)
      saveRDS(object, file)
      summarize_saved_object(label, object, file)
    },
    error = function(error) {
      data.frame(
        Object = label,
        Class = NA_character_,
        Rows = NA_integer_,
        Columns = NA_integer_,
        Length = NA_integer_,
        File = file,
        Status = "failed",
        Error = conditionMessage(error),
        stringsAsFactors = FALSE
      )
    }
  )
}

run_data_summary <- function() {
  dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
  dir.create(DATA_ROOT, recursive = TRUE, showWarnings = FALSE)

  cancers_summary <- as.data.frame.matrix(get_cancers())
  cancers_summary$Cancer <- rownames(cancers_summary)
  cancers_summary <- cancers_summary[, c("Cancer", setdiff(colnames(cancers_summary), "Cancer"))]
  write_table(cancers_summary, TABLE_ROOT, "000_data_summary_tcgaplot_cancers")

  paired_summary <- as.data.frame.matrix(get_paired_cancers())
  paired_summary$Cancer <- rownames(paired_summary)
  paired_summary <- paired_summary[, c("Cancer", setdiff(colnames(paired_summary), "Cancer"))]
  write_table(paired_summary, TABLE_ROOT, "000_data_summary_tcgaplot_paired_cancers")

  package_info <- data.frame(
    Package = "TCGAplot",
    Version = as.character(utils::packageVersion("TCGAplot")),
    Library_Path = find.package("TCGAplot"),
    Data_Root = DATA_ROOT,
    Temporary_Root = TEMP_ROOT,
    Result_Root = RESULT_ROOT,
    stringsAsFactors = FALSE
  )
  write_table(package_info, TABLE_ROOT, "000_data_summary_tcgaplot_package_info")

  if (SAVE_BUILTIN_DATA_EXTRACTS) {
    saved_extracts <- list()
    add_extract <- function(label, file, expr) {
      saved_extracts[[length(saved_extracts) + 1L]] <<- save_builtin_extract(label, file, expr)
    }

    add_extract("get_tmb", file.path(DATA_ROOT, "tcgaplot_tmb.rds"), get_tmb())
    add_extract("get_msi", file.path(DATA_ROOT, "tcgaplot_msi.rds"), get_msi())
    add_extract("get_immu_ratio", file.path(DATA_ROOT, "tcgaplot_immune_ratio.rds"), get_immu_ratio())
    add_extract("get_immuscore", file.path(DATA_ROOT, "tcgaplot_immune_score.rds"), get_immuscore())
    add_extract("get_geneset", file.path(DATA_ROOT, "tcgaplot_genesets.rds"), get_geneset())
    add_extract("get_methy", file.path(DATA_ROOT, "tcgaplot_methy.rds"), get_methy())
    add_extract("get_all_meta", file.path(DATA_ROOT, "tcgaplot_all_meta.rds"), get_all_meta())
    add_extract("get_all_tpm", file.path(DATA_ROOT, "tcgaplot_all_tpm.rds"), get_all_tpm())
    add_extract(
      "get_all_paired_tpm",
      file.path(DATA_ROOT, "tcgaplot_all_paired_tpm.rds"),
      get_all_paired_tpm()
    )
    add_extract(
      "get_all_promoter_methy",
      file.path(DATA_ROOT, "tcgaplot_all_promoter_methy.rds"),
      get_all_promoter_methy()
    )

    for (cancer in TARGET_CANCERS) {
      add_extract(
        paste0("get_tpm_", cancer),
        file.path(DATA_ROOT, paste0("tcgaplot_tpm_", cancer, ".rds")),
        get_tpm(cancer)
      )
      add_extract(
        paste0("get_meta_", cancer),
        file.path(DATA_ROOT, paste0("tcgaplot_meta_", cancer, ".rds")),
        get_meta(cancer)
      )
      add_extract(
        paste0("get_promoter_methy_", cancer),
        file.path(DATA_ROOT, paste0("tcgaplot_promoter_methy_", cancer, ".rds")),
        get_promoter_methy(cancer)
      )

      if (cancer %in% rownames(get_paired_cancers())) {
        add_extract(
          paste0("get_paired_tpm_", cancer),
          file.path(DATA_ROOT, paste0("tcgaplot_paired_tpm_", cancer, ".rds")),
          get_paired_tpm(cancer)
        )
      }
    }

    saved_extracts <- do.call(rbind, saved_extracts)
    write_table(saved_extracts, TABLE_ROOT, "000_data_summary_tcgaplot_saved_builtin_extracts")
  }

  package_info
}


# 6. TCGAplot任务目录 ----------------------------------------------------------

LITERATURE_QUICK_SCAN_ANALYSES <- c(
  "sample_audit",
  "data_summary",
  "pan_boxplot",
  "pan_paired_boxplot",
  "pan_tumor_boxplot",
  "pan_forest",
  "gene_TMB_radar",
  "gene_MSI_radar",
  "gene_checkpoint_heatmap",
  "gene_chemokine_heatmap",
  "gene_receptor_heatmap",
  "gene_immustimulator_heatmap",
  "gene_immuinhibitor_heatmap",
  "gene_immucell_heatmap",
  "gene_immunescore_heatmap",
  "gene_immunescore_triangle",
  "tcga_boxplot",
  "paired_boxplot",
  "gene_age",
  "gene_3age",
  "gene_gender",
  "gene_stage",
  "tcga_roc",
  "tcga_kmplot",
  "gene_gene_scatter",
  "gene_methylation_scatter",
  "methy_kmplot",
  "gene_deg_heatmap",
  "gene_gsea_go",
  "gene_gsea_kegg",
  "gene_coexp_heatmap",
  "gene_network_go",
  "gene_network_kegg"
)

GENESET_ANALYSES <- c(
  "gs_pan_boxplot",
  "gs_pan_paired_boxplot",
  "gs_pan_tumor_boxplot",
  "gs_pan_forest",
  "gs_TMB_radar",
  "gs_MSI_radar",
  "gs_checkpoint_heatmap",
  "gs_chemokine_heatmap",
  "gs_receptor_heatmap",
  "gs_immustimulator_heatmap",
  "gs_immuinhibitor_heatmap",
  "gs_immucell_heatmap",
  "gs_immunescore_heatmap",
  "gs_boxplot",
  "gs_paired_boxplot",
  "gs_age",
  "gs_3age",
  "gs_gender",
  "gs_stage",
  "gs_roc",
  "gs_kmplot"
)

ALL_WRAPPED_ANALYSES <- unique(c(LITERATURE_QUICK_SCAN_ANALYSES, GENESET_ANALYSES))

DATA_EXTRACTION_FUNCTIONS <- c(
  "get_cancers",
  "get_paired_cancers",
  "get_tpm",
  "get_paired_tpm",
  "get_meta",
  "get_tmb",
  "get_msi",
  "get_immu_ratio",
  "get_immuscore",
  "get_geneset",
  "get_methy",
  "get_promoter_methy",
  "get_all_tpm",
  "get_all_paired_tpm",
  "get_all_meta",
  "get_all_promoter_methy"
)

make_analysis_catalog <- function(selected_analyses = character(0)) {
  # 这张表是脚本内置说明书：
  # - Function为TCGAplot原始函数名；
  # - Script_Analysis为本脚本中用于选择该功能的任务名；
  # - Default_Run说明默认literature_quick_scan是否会运行；
  # - Requirement说明是否需要额外配置基因对、基因集或保存大矩阵。
  rows <- list()
  add_row <- function(
      function_name,
      script_analysis,
      module,
      scope,
      target_type,
      description,
      requirement = "无需额外配置",
      implemented = TRUE) {
    rows[[length(rows) + 1L]] <<- data.frame(
      Function = function_name,
      Script_Analysis = script_analysis,
      Module = module,
      Scope = scope,
      Target_Type = target_type,
      Description_CN = description,
      Requirement = requirement,
      Implemented_In_Script = implemented,
      Default_Run = script_analysis %in% LITERATURE_QUICK_SCAN_ANALYSES,
      Selected_This_Run = script_analysis %in% selected_analyses ||
        (
          function_name %in% DATA_EXTRACTION_FUNCTIONS &&
            "data_summary" %in% selected_analyses
        ),
      stringsAsFactors = FALSE
    )
  }

  add_row("sample_audit", "sample_audit", "样本审计", "COAD本地SE对照", "样本", "比较TCGAplot内置COAD样本与本项目本地SE样本是否一致")
  add_row("get_cancers", "data_summary", "内置数据概览", "泛癌", "数据", "输出TCGAplot内置33癌种tumor/normal样本数")
  add_row("get_paired_cancers", "data_summary", "内置数据概览", "泛癌", "数据", "输出TCGAplot内置可配对癌种样本数")

  add_row("pan_boxplot", "pan_boxplot", "单基因表达", "泛癌", "基因", "泛癌tumor-normal表达差异箱线图")
  add_row("pan_paired_boxplot", "pan_paired_boxplot", "单基因表达", "泛癌", "基因", "泛癌配对tumor-normal表达差异箱线图")
  add_row("pan_tumor_boxplot", "pan_tumor_boxplot", "单基因表达", "泛癌", "基因", "泛癌肿瘤样本表达分布箱线图")
  add_row("pan_forest", "pan_forest", "单基因预后", "泛癌", "基因", "泛癌Cox回归森林图，同时输出未校正和年龄校正版本")
  add_row("gene_TMB_radar", "gene_TMB_radar", "单基因相关", "泛癌", "基因", "目标基因表达与TMB相关性的泛癌雷达图")
  add_row("gene_MSI_radar", "gene_MSI_radar", "单基因相关", "泛癌", "基因", "目标基因表达与MSI相关性的泛癌雷达图")
  add_row("gene_checkpoint_heatmap", "gene_checkpoint_heatmap", "单基因免疫相关", "泛癌", "基因", "目标基因与免疫检查点基因相关性热图")
  add_row("gene_chemokine_heatmap", "gene_chemokine_heatmap", "单基因免疫相关", "泛癌", "基因", "目标基因与趋化因子相关性热图")
  add_row("gene_receptor_heatmap", "gene_receptor_heatmap", "单基因免疫相关", "泛癌", "基因", "目标基因与趋化因子受体相关性热图")
  add_row("gene_immustimulator_heatmap", "gene_immustimulator_heatmap", "单基因免疫相关", "泛癌", "基因", "目标基因与免疫刺激因子相关性热图")
  add_row("gene_immuinhibitor_heatmap", "gene_immuinhibitor_heatmap", "单基因免疫相关", "泛癌", "基因", "目标基因与免疫抑制因子相关性热图")
  add_row("gene_immucell_heatmap", "gene_immucell_heatmap", "单基因免疫相关", "泛癌", "基因", "目标基因与免疫细胞比例相关性热图")
  add_row("gene_immunescore_heatmap", "gene_immunescore_heatmap", "单基因免疫相关", "泛癌", "基因", "目标基因与Stromal/Immune/ESTIMATE score相关性热图")
  add_row("gene_immunescore_triangle", "gene_immunescore_triangle", "单基因免疫相关", "泛癌", "基因", "目标基因与免疫评分相关性的三角热图")

  add_row("tcga_boxplot", "tcga_boxplot", "单基因表达", "指定癌种", "基因", "COAD tumor-normal表达差异箱线图")
  add_row("paired_boxplot", "paired_boxplot", "单基因表达", "指定癌种", "基因", "COAD配对tumor-normal表达差异箱线图")
  add_row("gene_age", "gene_age", "临床分组", "指定癌种", "基因", "按年龄二分组比较目标基因表达")
  add_row("gene_3age", "gene_3age", "临床分组", "指定癌种", "基因", "按年轻/中间/年老三组比较目标基因表达")
  add_row("gene_gender", "gene_gender", "临床分组", "指定癌种", "基因", "按性别比较目标基因表达")
  add_row("gene_stage", "gene_stage", "临床分组", "指定癌种", "基因", "按病理分期比较目标基因表达")
  add_row("tcga_roc", "tcga_roc", "诊断分析", "指定癌种", "基因", "COAD tumor-normal诊断ROC曲线")
  add_row("tcga_kmplot", "tcga_kmplot", "生存分析", "指定癌种", "基因", "按目标基因高低表达分组绘制KM生存曲线")
  add_row("gene_methylation_scatter", "gene_methylation_scatter", "甲基化分析", "指定癌种", "基因", "目标基因表达与启动子甲基化探针相关散点图")
  add_row("methy_kmplot", "methy_kmplot", "甲基化生存分析", "指定癌种", "基因", "按目标基因启动子甲基化高低分组绘制KM生存曲线")
  add_row("gene_deg_heatmap", "gene_deg_heatmap", "高低表达分组机制分析", "指定癌种", "基因", "按目标基因高低表达分组计算DEG并绘制Top DEG热图")
  add_row("gene_gsea_go", "gene_gsea_go", "高低表达分组机制分析", "指定癌种", "基因", "按目标基因高低表达分组做GSEA-GO")
  add_row("gene_gsea_kegg", "gene_gsea_kegg", "高低表达分组机制分析", "指定癌种", "基因", "按目标基因高低表达分组做GSEA-KEGG")
  add_row("gene_coexp_heatmap", "gene_coexp_heatmap", "共表达分析", "指定癌种", "基因", "目标基因正/负共表达基因热图及GO富集")
  add_row("gene_gene_scatter", "gene_gene_scatter", "基因相关", "指定癌种", "基因对", "两个指定基因的表达相关散点图", "需要在GENE_PAIR_DESIGNS配置基因对")
  add_row("gene_network_go", "gene_network_go", "基因网络", "不限癌种", "基因/基因列表", "输入基因列表的GO cnetplot网络图", "默认至少2个NETWORK_GENES；ATF3单基因默认跳过")
  add_row("gene_network_kegg", "gene_network_kegg", "基因网络", "不限癌种", "基因/基因列表", "输入基因列表的KEGG cnetplot网络图", "默认至少2个NETWORK_GENES；ATF3单基因默认跳过")

  add_row("gs_pan_boxplot", "gs_pan_boxplot", "基因集表达", "泛癌", "基因集", "基因集泛癌tumor-normal表达差异箱线图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_pan_paired_boxplot", "gs_pan_paired_boxplot", "基因集表达", "泛癌", "基因集", "基因集泛癌配对tumor-normal表达差异箱线图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_pan_tumor_boxplot", "gs_pan_tumor_boxplot", "基因集表达", "泛癌", "基因集", "基因集泛癌肿瘤样本表达分布箱线图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_pan_forest", "gs_pan_forest", "基因集预后", "泛癌", "基因集", "基因集泛癌Cox回归森林图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_TMB_radar", "gs_TMB_radar", "基因集相关", "泛癌", "基因集", "基因集活性与TMB相关性雷达图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_MSI_radar", "gs_MSI_radar", "基因集相关", "泛癌", "基因集", "基因集活性与MSI相关性雷达图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_checkpoint_heatmap", "gs_checkpoint_heatmap", "基因集免疫相关", "泛癌", "基因集", "基因集与免疫检查点相关性热图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_chemokine_heatmap", "gs_chemokine_heatmap", "基因集免疫相关", "泛癌", "基因集", "基因集与趋化因子相关性热图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_receptor_heatmap", "gs_receptor_heatmap", "基因集免疫相关", "泛癌", "基因集", "基因集与趋化因子受体相关性热图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_immustimulator_heatmap", "gs_immustimulator_heatmap", "基因集免疫相关", "泛癌", "基因集", "基因集与免疫刺激因子相关性热图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_immuinhibitor_heatmap", "gs_immuinhibitor_heatmap", "基因集免疫相关", "泛癌", "基因集", "基因集与免疫抑制因子相关性热图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_immucell_heatmap", "gs_immucell_heatmap", "基因集免疫相关", "泛癌", "基因集", "基因集与免疫细胞比例相关性热图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_immunescore_heatmap", "gs_immunescore_heatmap", "基因集免疫相关", "泛癌", "基因集", "基因集与免疫评分相关性热图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_boxplot", "gs_boxplot", "基因集表达", "指定癌种", "基因集", "COAD基因集tumor-normal表达差异箱线图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_paired_boxplot", "gs_paired_boxplot", "基因集表达", "指定癌种", "基因集", "COAD基因集配对tumor-normal表达差异箱线图", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_age", "gs_age", "基因集临床分组", "指定癌种", "基因集", "按年龄二分组比较基因集活性", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_3age", "gs_3age", "基因集临床分组", "指定癌种", "基因集", "按三年龄组比较基因集活性", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_gender", "gs_gender", "基因集临床分组", "指定癌种", "基因集", "按性别比较基因集活性", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_stage", "gs_stage", "基因集临床分组", "指定癌种", "基因集", "按病理分期比较基因集活性", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_roc", "gs_roc", "基因集诊断分析", "指定癌种", "基因集", "基因集活性诊断ROC曲线", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")
  add_row("gs_kmplot", "gs_kmplot", "基因集生存分析", "指定癌种", "基因集", "按基因集活性高低分组绘制KM生存曲线", "需要RUN_GENESET_ANALYSES=TRUE且配置TARGET_GENESETS")

  add_row("get_tpm", "data_summary", "数据提取", "指定癌种", "数据", "提取指定癌种TPM矩阵", "默认只写概览；如需保存矩阵，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_paired_tpm", "data_summary", "数据提取", "指定癌种", "数据", "提取指定癌种配对TPM矩阵", "默认只写概览；如需保存矩阵，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_meta", "data_summary", "数据提取", "指定癌种", "数据", "提取指定癌种临床信息", "默认只写概览；如需保存矩阵，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_tmb", "data_summary", "数据提取", "泛癌", "数据", "提取泛癌TMB矩阵", "默认只写概览；如需保存矩阵，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_msi", "data_summary", "数据提取", "泛癌", "数据", "提取泛癌MSI矩阵", "默认只写概览；如需保存矩阵，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_immu_ratio", "data_summary", "数据提取", "泛癌", "数据", "提取泛癌免疫细胞比例矩阵", "默认只写概览；如需保存矩阵，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_immuscore", "data_summary", "数据提取", "泛癌", "数据", "提取泛癌免疫评分矩阵", "默认只写概览；如需保存矩阵，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_geneset", "data_summary", "数据提取", "泛癌", "数据", "提取TCGAplot内置基因集列表", "默认只写概览；如需保存对象，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_methy", "data_summary", "数据提取", "泛癌", "数据", "提取或提示TCGAplot内置甲基化数据", "默认只写概览；如需保存对象，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_promoter_methy", "data_summary", "数据提取", "指定癌种", "数据", "提取指定癌种启动子甲基化矩阵", "默认只写概览；如需保存对象，设置TCGAPLOT_SAVE_BUILTIN_DATA=1")
  add_row("get_all_tpm", "data_summary", "数据提取", "泛癌", "数据", "提取全部TPM矩阵", "对象很大；仅在TCGAPLOT_SAVE_BUILTIN_DATA=1时建议保存")
  add_row("get_all_paired_tpm", "data_summary", "数据提取", "泛癌", "数据", "提取全部配对TPM矩阵", "对象较大；仅在TCGAPLOT_SAVE_BUILTIN_DATA=1时建议保存")
  add_row("get_all_meta", "data_summary", "数据提取", "泛癌", "数据", "提取全部临床信息", "仅在TCGAPLOT_SAVE_BUILTIN_DATA=1时保存")
  add_row("get_all_promoter_methy", "data_summary", "数据提取", "泛癌", "数据", "提取全部启动子甲基化矩阵", "对象很大；仅在TCGAPLOT_SAVE_BUILTIN_DATA=1时建议保存")

  catalog <- do.call(rbind, rows)
  catalog <- catalog[order(catalog$Module, catalog$Function), , drop = FALSE]
  rownames(catalog) <- NULL
  catalog
}

write_analysis_catalog <- function(selected_analyses) {
  catalog <- make_analysis_catalog(selected_analyses)
  write_table(catalog, SUMMARY_ROOT, "000_tcgaplot_analysis_catalog")

  exported_functions <- setdiff(sort(ls("package:TCGAplot")), "%>%")
  cataloged_functions <- unique(catalog$Function)
  missing_from_script <- setdiff(exported_functions, cataloged_functions)

  missing_table <- data.frame(
    Function = missing_from_script,
    Status = rep("TCGAplot导出但尚未登记到脚本目录", length(missing_from_script)),
    stringsAsFactors = FALSE
  )
  write_table(missing_table, SUMMARY_ROOT, "000_tcgaplot_functions_not_in_script_catalog")

  invisible(list(catalog = catalog, missing = missing_table))
}

resolve_requested_analyses <- function() {
  requested <- ANALYSES_TO_RUN
  if (length(requested) == 1 && requested == "literature_quick_scan") {
    return(LITERATURE_QUICK_SCAN_ANALYSES)
  }
  if (length(requested) == 1 && requested == "all") {
    return(ALL_WRAPPED_ANALYSES)
  }

  unknown <- setdiff(requested, ALL_WRAPPED_ANALYSES)
  if (length(unknown) > 0) {
    stop("Unknown TCGAplot analyses requested: ", paste(unknown, collapse = ", "))
  }

  requested
}

make_task <- function(
    analysis,
    function_name,
    args,
    target,
    context,
    renderer = "grid",
    width = 7,
    height = 6,
    title = NULL) {
  list(
    analysis = analysis,
    function_name = function_name,
    args = args,
    target = target,
    context = context,
    renderer = renderer,
    width = width,
    height = height,
    title = title
  )
}

make_default_task_title <- function(task) {
  if (!is.null(task$title) && nzchar(task$title)) {
    return(task$title)
  }

  paste(
    gsub("_", " ", task$analysis, fixed = TRUE),
    task$target,
    task$context,
    sep = " | "
  )
}

estimate_title_extra_height <- function(title) {
  if (is.null(title) || !nzchar(title)) {
    return(0)
  }
  line_count <- max(1L, ceiling(nchar(title, type = "width") / 70))
  0.35 + 0.18 * (line_count - 1L)
}

adjust_task_plot_size <- function(task) {
  title <- make_default_task_title(task)
  task$height <- task$height + estimate_title_extra_height(title)

  if (task$analysis %in% c("pan_boxplot", "pan_paired_boxplot", "pan_tumor_boxplot")) {
    task$width <- max(task$width, 12.5)
    task$height <- max(task$height, 6.4)
  }

  if (task$analysis == "pan_forest") {
    task$width <- max(task$width, 9.2)
    task$height <- max(task$height, 10.8)
  }

  if (grepl("heatmap", task$analysis, fixed = TRUE)) {
    task$width <- max(task$width, 7.8)
    task$height <- task$height + 0.25
  }

  if (task$analysis %in% c("gene_gsea_go", "gene_gsea_kegg")) {
    task$width <- max(task$width, 8.8)
    task$height <- max(task$height, 7.2)
  }

  if (task$analysis %in% c("tcga_kmplot", "methy_kmplot")) {
    task$width <- max(task$width, 7.8)
    task$height <- max(task$height, 6.6)
  }

  task
}

assign_task_indices <- function(tasks) {
  if (length(tasks) == 0L) {
    return(tasks)
  }

  for (i in seq_along(tasks)) {
    tasks[[i]]$task_id <- i
    tasks[[i]] <- adjust_task_plot_size(tasks[[i]])
  }

  tasks
}

task_file_stem <- function(task, suffix = NULL) {
  stem_parts <- c(
    sprintf("%03d", task$task_id %||% 0L),
    safe_name(task$analysis),
    safe_name(task$target),
    safe_name(task$context)
  )
  if (!is.null(suffix) && nzchar(suffix)) {
    stem_parts <- c(stem_parts, safe_name(suffix))
  }

  paste(stem_parts, collapse = "_")
}

`%||%` <- function(x, y) {
  if (is.null(x)) {
    return(y)
  }
  x
}

build_tcgaplot_tasks <- function(selected_analyses) {
  tasks <- list()
  add_task <- function(task) {
    tasks[[length(tasks) + 1L]] <<- task
  }

  for (gene in TARGET_GENES) {
    if ("pan_boxplot" %in% selected_analyses) {
      add_task(make_task(
        "pan_boxplot", "pan_boxplot",
        list(gene = gene, palette = BOXPLOT_PALETTE, legend = PAN_LEGEND_POSITION, method = GROUP_TEST_METHOD),
        gene, "pan_cancer", width = 12, height = 5.8
      ))
    }
    if ("pan_paired_boxplot" %in% selected_analyses) {
      add_task(make_task(
        "pan_paired_boxplot", "pan_paired_boxplot",
        list(gene = gene, palette = BOXPLOT_PALETTE, legend = PAN_LEGEND_POSITION, method = GROUP_TEST_METHOD),
        gene, "pan_cancer", width = 12, height = 5.8
      ))
    }
    if ("pan_tumor_boxplot" %in% selected_analyses) {
      add_task(make_task(
        "pan_tumor_boxplot", "pan_tumor_boxplot",
        list(gene = gene),
        gene, "pan_cancer", width = 12, height = 5.8
      ))
    }
    if ("pan_forest" %in% selected_analyses) {
      add_task(make_task(
        "pan_forest", "pan_forest",
        list(gene = gene, adjust = FALSE),
        gene, "pan_cancer_unadjusted", width = 8.5, height = 10,
        title = paste0("Hazard Ratio Plot of ", gene)
      ))
      for (adjust_variable in PAN_FOREST_ADJUST_VARIABLES) {
        if (adjust_variable == "age") {
          add_task(make_task(
            "pan_forest", "pan_forest",
            list(gene = gene, adjust = TRUE),
            gene, "pan_cancer_age_adjusted", width = 8.5, height = 10,
            title = paste0("Hazard Ratio Plot of ", gene, " adjusted by age")
          ))
        } else {
          add_task(make_task(
            "pan_forest", "pan_forest_clinical_adjusted",
            list(gene = gene, adjust_variable = adjust_variable),
            gene, paste0("pan_cancer_", adjust_variable, "_adjusted"),
            width = 8.5,
            height = 10,
            title = paste0("Hazard Ratio Plot of ", gene, " adjusted by ", adjust_variable)
          ))
        }
      }
    }
    if ("gene_TMB_radar" %in% selected_analyses) {
      add_task(make_task(
        "gene_TMB_radar", "gene_TMB_radar",
        list(gene = gene, method = CORRELATION_METHOD),
        gene, "pan_cancer", width = 8, height = 8
      ))
    }
    if ("gene_MSI_radar" %in% selected_analyses) {
      add_task(make_task(
        "gene_MSI_radar", "gene_MSI_radar",
        list(gene = gene, method = CORRELATION_METHOD),
        gene, "pan_cancer", width = 8, height = 8
      ))
    }

    immune_heatmap_names <- c(
      "gene_checkpoint_heatmap",
      "gene_chemokine_heatmap",
      "gene_receptor_heatmap",
      "gene_immustimulator_heatmap",
      "gene_immuinhibitor_heatmap",
      "gene_immucell_heatmap",
      "gene_immunescore_heatmap"
    )
    for (analysis_name in intersect(immune_heatmap_names, selected_analyses)) {
      add_task(make_task(
        analysis_name,
        analysis_name,
        list(
          gene = gene,
          method = CORRELATION_METHOD,
          lowcol = HEATMAP_LOW_COLOR,
          highcol = HEATMAP_HIGH_COLOR,
          cluster_row = HEATMAP_CLUSTER_ROW,
          cluster_col = HEATMAP_CLUSTER_COL,
          legend = HEATMAP_LEGEND
        ),
        gene,
        "pan_cancer",
        width = 12,
        height = 7.5
      ))
    }
    if ("gene_immunescore_triangle" %in% selected_analyses) {
      add_task(make_task(
        "gene_immunescore_triangle", "gene_immunescore_triangle",
        list(gene = gene, method = CORRELATION_METHOD),
        gene, "pan_cancer", width = 8.5, height = 8
      ))
    }

    for (cancer in TARGET_CANCERS) {
      if ("tcga_boxplot" %in% selected_analyses) {
        add_task(make_task(
          "tcga_boxplot", "tcga_boxplot",
          list(
            cancer = cancer, gene = gene, add = BOXPLOT_ADD_LAYER,
            palette = BOXPLOT_PALETTE, legend = LEGEND_POSITION,
            label = P_VALUE_LABEL, method = GROUP_TEST_METHOD
          ),
          gene, cancer, width = 4.8, height = 4.8
        ))
      }
      if ("paired_boxplot" %in% selected_analyses) {
        add_task(make_task(
          "paired_boxplot", "paired_boxplot",
          list(
            cancer = cancer, gene = gene, palette = BOXPLOT_PALETTE,
            legend = LEGEND_POSITION, label = P_VALUE_LABEL,
            method = GROUP_TEST_METHOD
          ),
          gene, cancer, width = 4.8, height = 4.8
        ))
      }
      if ("gene_age" %in% selected_analyses) {
        add_task(make_task(
          "gene_age", "gene_age",
          list(
            cancer = cancer, gene = gene, age = AGE_CUTOFF,
            add = BOXPLOT_ADD_LAYER, palette = BOXPLOT_PALETTE,
            legend = LEGEND_POSITION, label = P_VALUE_LABEL,
            method = GROUP_TEST_METHOD
          ),
          gene, cancer, width = 4.8, height = 4.8
        ))
      }
      if ("gene_3age" %in% selected_analyses) {
        add_task(make_task(
          "gene_3age", "gene_3age",
          list(
            cancer = cancer, gene = gene, age1 = AGE_YOUNG_CUTOFF,
            age2 = AGE_OLD_CUTOFF, add = BOXPLOT_ADD_LAYER,
            palette = BOXPLOT_PALETTE, legend = LEGEND_POSITION,
            label = P_VALUE_LABEL, method = GROUP_TEST_METHOD
          ),
          gene, cancer, width = 5.2, height = 4.8
        ))
      }
      if ("gene_gender" %in% selected_analyses) {
        add_task(make_task(
          "gene_gender", "gene_gender",
          list(
            cancer = cancer, gene = gene, add = BOXPLOT_ADD_LAYER,
            palette = BOXPLOT_PALETTE, legend = LEGEND_POSITION,
            label = P_VALUE_LABEL, method = GROUP_TEST_METHOD
          ),
          gene, cancer, width = 4.8, height = 4.8
        ))
      }
      if ("gene_stage" %in% selected_analyses) {
        add_task(make_task(
          "gene_stage", "gene_stage",
          list(
            cancer = cancer, gene = gene, add = BOXPLOT_ADD_LAYER,
            palette = BOXPLOT_PALETTE, legend = LEGEND_POSITION,
            label = P_VALUE_LABEL, method = GROUP_TEST_METHOD
          ),
          gene, cancer, width = 5.8, height = 4.8
        ))
      }
      if ("tcga_roc" %in% selected_analyses) {
        add_task(make_task(
          "tcga_roc", "tcga_roc",
          list(cancer = cancer, gene = gene),
          gene, cancer, renderer = "direct_once", width = 5.2, height = 5.2
        ))
      }
      if ("tcga_kmplot" %in% selected_analyses) {
        add_task(make_task(
          "tcga_kmplot", "tcga_kmplot",
          list(cancer = cancer, gene = gene, palette = SURVIVAL_PALETTE),
          gene, cancer, width = 7.2, height = 6.2
        ))
      }
      if ("gene_methylation_scatter" %in% selected_analyses) {
        add_task(make_task(
          "gene_methylation_scatter", "gene_methylation_scatter",
          list(cancer = cancer, gene = gene),
          gene, cancer, renderer = "internal_pdf", width = 5, height = 4
        ))
      }
      if ("methy_kmplot" %in% selected_analyses) {
        add_task(make_task(
          "methy_kmplot", "methy_kmplot",
          list(cancer = cancer, gene = gene, palette = SURVIVAL_PALETTE),
          gene, cancer, renderer = "internal_pdf", width = 8, height = 6
        ))
      }
      if ("gene_deg_heatmap" %in% selected_analyses) {
        add_task(make_task(
          "gene_deg_heatmap", "gene_deg_heatmap",
          list(cancer = cancer, gene = gene, top_n = TOP_N_GENES),
          gene, cancer, width = 7.5, height = 8
        ))
      }
      if ("gene_gsea_go" %in% selected_analyses) {
        add_task(make_task(
          "gene_gsea_go", "gene_gsea_go",
          list(cancer = cancer, gene = gene),
          gene, cancer, width = 8, height = 6.5
        ))
      }
      if ("gene_gsea_kegg" %in% selected_analyses) {
        add_task(make_task(
          "gene_gsea_kegg", "gene_gsea_kegg",
          list(cancer = cancer, gene = gene),
          gene, cancer, width = 8, height = 6.5
        ))
      }
      if ("gene_coexp_heatmap" %in% selected_analyses) {
        add_task(make_task(
          "gene_coexp_heatmap", "gene_coexp_heatmap",
          list(cancer = cancer, gene = gene, top_n = TOP_N_GENES, method = CORRELATION_METHOD),
          gene, cancer, width = 11, height = 8
        ))
      }
    }
  }

  if ("gene_gene_scatter" %in% selected_analyses && nrow(GENE_PAIR_DESIGNS) > 0) {
    for (i in seq_len(nrow(GENE_PAIR_DESIGNS))) {
      design <- GENE_PAIR_DESIGNS[i, , drop = FALSE]
      add_task(make_task(
        "gene_gene_scatter", "gene_gene_scatter",
        list(
          cancer = design$Cancer,
          gene1 = design$Gene1,
          gene2 = design$Gene2,
          density = design$Density
        ),
        paste(design$Gene1, design$Gene2, sep = "_vs_"),
        design$Cancer,
        width = 6,
        height = 5.5
      ))
    }
  }

  if ("gene_network_go" %in% selected_analyses) {
    add_task(make_task(
      "gene_network_go", "gene_network_go",
      list(gene = NETWORK_GENES),
      paste_compact(NETWORK_GENES), "gene_network", width = 9, height = 8
    ))
  }
  if ("gene_network_kegg" %in% selected_analyses) {
    add_task(make_task(
      "gene_network_kegg", "gene_network_kegg",
      list(gene = NETWORK_GENES),
      paste_compact(NETWORK_GENES), "gene_network", width = 9, height = 8
    ))
  }

  if (RUN_GENESET_ANALYSES && nrow(TARGET_GENESETS) > 0) {
    tasks <- c(tasks, build_geneset_tasks(selected_analyses))
  }

  tasks
}

build_geneset_tasks <- function(selected_analyses) {
  tasks <- list()
  add_task <- function(task) {
    tasks[[length(tasks) + 1L]] <<- task
  }

  for (i in seq_len(nrow(TARGET_GENESETS))) {
    geneset <- TARGET_GENESETS$GeneSet[i]
    alias <- TARGET_GENESETS$Alias[i]
    target_name <- safe_name(alias, default = safe_name(geneset, default = "geneset"))

    pan_geneset_functions <- c(
      "gs_pan_boxplot",
      "gs_pan_paired_boxplot",
      "gs_pan_tumor_boxplot",
      "gs_pan_forest",
      "gs_TMB_radar",
      "gs_MSI_radar",
      "gs_checkpoint_heatmap",
      "gs_chemokine_heatmap",
      "gs_receptor_heatmap",
      "gs_immustimulator_heatmap",
      "gs_immuinhibitor_heatmap",
      "gs_immucell_heatmap",
      "gs_immunescore_heatmap"
    )

    for (analysis_name in intersect(pan_geneset_functions, selected_analyses)) {
      formals_names <- names(formals(get_tcgaplot_function(analysis_name)))
      args <- list(geneset = geneset)
      if ("geneset_alias" %in% formals_names) {
        args$geneset_alias <- alias
      }
      if ("palette" %in% formals_names) {
        args$palette <- BOXPLOT_PALETTE
      }
      if ("legend" %in% formals_names) {
        args$legend <- PAN_LEGEND_POSITION
      }
      if ("method" %in% formals_names) {
        args$method <- if (analysis_name %in% c("gs_pan_boxplot", "gs_pan_paired_boxplot")) {
          GROUP_TEST_METHOD
        } else {
          CORRELATION_METHOD
        }
      }
      if ("adjust" %in% formals_names) {
        args$adjust <- FALSE
      }
      if ("lowcol" %in% formals_names) {
        args$lowcol <- HEATMAP_LOW_COLOR
        args$highcol <- HEATMAP_HIGH_COLOR
        args$cluster_row <- HEATMAP_CLUSTER_ROW
        args$cluster_col <- HEATMAP_CLUSTER_COL
        args$legend <- HEATMAP_LEGEND
      }

      add_task(make_task(
        analysis_name,
        analysis_name,
        args,
        target_name,
        "pan_cancer",
        width = 12,
        height = 7
      ))
    }

    cancer_geneset_functions <- c(
      "gs_boxplot",
      "gs_paired_boxplot",
      "gs_age",
      "gs_3age",
      "gs_gender",
      "gs_stage",
      "gs_roc",
      "gs_kmplot"
    )
    for (cancer in TARGET_CANCERS) {
      for (analysis_name in intersect(cancer_geneset_functions, selected_analyses)) {
        formals_names <- names(formals(get_tcgaplot_function(analysis_name)))
        args <- list(cancer = cancer, geneset = geneset, geneset_alias = alias)
        if ("add" %in% formals_names) {
          args$add <- BOXPLOT_ADD_LAYER
        }
        if ("palette" %in% formals_names) {
          args$palette <- BOXPLOT_PALETTE
        }
        if ("legend" %in% formals_names) {
          args$legend <- LEGEND_POSITION
        }
        if ("label" %in% formals_names) {
          args$label <- P_VALUE_LABEL
        }
        if ("method" %in% formals_names) {
          args$method <- GROUP_TEST_METHOD
        }
        if ("age" %in% formals_names) {
          args$age <- AGE_CUTOFF
        }
        if ("age1" %in% formals_names) {
          args$age1 <- AGE_YOUNG_CUTOFF
          args$age2 <- AGE_OLD_CUTOFF
        }

        add_task(make_task(
          analysis_name,
          analysis_name,
          args,
          target_name,
          cancer,
          width = 6,
          height = 5.5
        ))
      }
    }
  }

  tasks
}

validate_tcgaplot_tasks <- function(tasks) {
  if (length(tasks) == 0L) {
    validation <- data.frame()
    write_table(validation, SUMMARY_ROOT, "000_tcgaplot_task_argument_validation")
    return(invisible(validation))
  }

  validation <- do.call(
    rbind,
    lapply(seq_along(tasks), function(i) {
      task <- tasks[[i]]
      fn <- get_tcgaplot_function(task$function_name)
      formal_names <- names(formals(fn))
      arg_names <- names(task$args)
      unknown_args <- setdiff(arg_names, formal_names)

      data.frame(
        Task_ID = i,
        Analysis = task$analysis,
        Function = task$function_name,
        Target = task$target,
        Context = task$context,
        Argument_Names = paste(arg_names, collapse = ";"),
        Unknown_Arguments = paste(unknown_args, collapse = ";"),
        Status = ifelse(length(unknown_args) == 0L, "ok", "invalid_args"),
        stringsAsFactors = FALSE
      )
    })
  )

  write_table(validation, SUMMARY_ROOT, "000_tcgaplot_task_argument_validation")

  invalid <- validation[validation$Status != "ok", , drop = FALSE]
  if (nrow(invalid) > 0L) {
    stop(
      "TCGAplot task argument validation failed. See: ",
      file.path(SUMMARY_ROOT, "000_tcgaplot_task_argument_validation.csv")
    )
  }

  invisible(validation)
}


# 7. TCGAplot任务执行器 --------------------------------------------------------

should_run_task_in_subprocess <- function(task) {
  !DISABLE_PROCESS_ISOLATION && task$analysis %in% PROCESS_ISOLATED_ANALYSES
}

read_log_excerpt <- function(log_file, max_lines = 80L) {
  if (!file.exists(log_file) || file.info(log_file)$size <= 0) {
    return(character(0))
  }
  lines <- readLines(log_file, warn = FALSE)
  if (length(lines) > max_lines) {
    lines <- c(lines[seq_len(max_lines)], "... log truncated ...")
  }
  lines
}

run_tcgaplot_task_in_subprocess <- function(task, task_start_time, cache_file = "") {
  worker_root <- file.path(TEMP_ROOT, "process-isolated", task_file_stem(task))
  unlink(worker_root, recursive = TRUE, force = TRUE)
  dir.create(worker_root, recursive = TRUE, showWarnings = FALSE)

  task_file <- file.path(worker_root, "task.rds")
  result_file <- file.path(worker_root, "result.rds")
  stdout_log <- file.path(worker_root, "stdout.log")
  stderr_log <- file.path(worker_root, "stderr.log")
  worker_script <- file.path(worker_root, "worker.R")

  saveRDS(task, task_file)
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "project_root <- args[[1]]",
    "script_file <- args[[2]]",
    "task_file <- args[[3]]",
    "result_file <- args[[4]]",
    "setwd(project_root)",
    "Sys.setenv(TCGAPLOT_DISABLE_PROCESS_ISOLATION = '1')",
    "Sys.setenv(TCGAPLOT_CLEAR_PREVIOUS_OUTPUTS = '0')",
    "Sys.setenv(TCGAPLOT_PARALLEL_WORKERS = '1')",
    "expr <- parse(script_file)",
    "labels <- vapply(expr, function(e) paste(deparse(e, nlines = 1L), collapse = ' '), character(1))",
    "stop_at <- which(grepl('^clear_previous_tcgaplot_run_outputs <-', labels))[1]",
    "if (is.na(stop_at)) stop('Could not find TCGAplot main-run boundary')",
    "for (i in seq_len(stop_at)) eval(expr[[i]], envir = .GlobalEnv)",
    "task <- readRDS(task_file)",
    "manifest <- run_one_tcgaplot_task(task)",
    "saveRDS(manifest, result_file)"
  ), worker_script)

  rscript <- file.path(R.home("bin"), "Rscript")
  status <- tryCatch(
    suppressWarnings(system2(
      rscript,
      args = c("--vanilla", worker_script, PROJECT_ROOT, SCRIPT_FILE, task_file, result_file),
      stdout = stdout_log,
      stderr = stderr_log
    )),
    error = function(error) structure(character(0), status = 1L)
  )
  status_code <- attr(status, "status")
  if (is.null(status_code)) {
    status_code <- 0L
  }

  if (status_code != 0L || !file.exists(result_file)) {
    cleanup_task_output_files(task)
    log_excerpt <- c(read_log_excerpt(stderr_log), read_log_excerpt(stdout_log))
    return(make_task_result(
      task = task,
      status = "failed",
      messages = paste0("Process-isolated worker logs: ", stderr_log, "; ", stdout_log),
      error = paste(
        paste0("Process-isolated worker exited with status ", status_code),
        paste(log_excerpt, collapse = "\n"),
        sep = "\n"
      ),
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  manifest <- readRDS(result_file)
  validation_error <- tryCatch(
    {
      validate_plot_output_files(c(
        split_manifest_files(manifest$Output_Files[1]),
        split_manifest_files(manifest$Generated_Files[1])
      ))
      NULL
    },
    error = function(error) error
  )
  if (inherits(validation_error, "error")) {
    cleanup_task_output_files(task)
    return(make_task_result(
      task = task,
      status = "failed",
      messages = paste0("Process-isolated worker logs: ", stderr_log, "; ", stdout_log),
      error = conditionMessage(validation_error),
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  manifest
}

run_one_tcgaplot_task <- function(task) {
  # 每个TCGAplot函数被抽象成一个task：
  # analysis：脚本任务名；function_name：TCGAplot原始函数名；
  # target/context：用于组织输出目录；renderer：用于处理不同函数的输出方式。
  task_start_time <- Sys.time()
  cache_file <- if (is_cacheable_tcgaplot_task(task)) {
    get_tcgaplot_task_cache_file(task)
  } else {
    ""
  }

  skip_reason <- get_tcgaplot_task_skip_reason(task)
  if (nzchar(skip_reason)) {
    return(make_task_result(
      task = task,
      status = "skipped",
      messages = skip_reason,
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  dir.create(PLOT_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(PLOT_PNG_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)

  cached_manifest <- read_tcgaplot_task_cache(task)
  if (!is.null(cached_manifest)) {
    return(make_task_result(
      task = task,
      status = "cached",
      output_files = split_manifest_files(cached_manifest$Output_Files[1]),
      generated_files = split_manifest_files(cached_manifest$Generated_Files[1]),
      messages = "Reused TCGAplot task cache because output files already exist",
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  if (should_run_task_in_subprocess(task)) {
    return(run_tcgaplot_task_in_subprocess(
      task = task,
      task_start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  task_stem <- task_file_stem(task)
  task_title <- make_default_task_title(task)

  task_call <- function() {
    fn <- get_tcgaplot_function(task$function_name)

    if (task$renderer == "internal_pdf") {
      # 少数TCGAplot函数不返回图形对象，而是自己在工作目录写PDF。
      # 这类任务只运行函数，然后复制它产生的文件。
      do.call(fn, task$args)
      generated_files <- copy_generated_files(task)
      return(list(output_files = character(0), generated_files = generated_files))
    }

    if (task$renderer == "once") {
      # 部分对象需要先在内存中生成，再打开设备保存，避免重复运行函数。
      result <- do.call(fn, task$args)
      output_files <- save_tcgaplot_pdf_png(
        pdf_file = file.path(PLOT_ROOT, paste0(task_stem, ".pdf")),
        width = task$width,
        height = task$height,
        draw_fun = function() {
          render_tcgaplot_output(result, title = task_title, add_outer_title = TRUE)
        }
      )
      generated_files <- copy_generated_files(task)
      return(list(
        output_files = unname(unlist(output_files)),
        generated_files = generated_files
      ))
    }

    output_files <- save_tcgaplot_pdf_png(
      pdf_file = file.path(PLOT_ROOT, paste0(task_stem, ".pdf")),
      width = task$width,
      height = task$height,
      draw_fun = function() {
        result <- do.call(fn, task$args)
        render_tcgaplot_output(result, title = task_title, add_outer_title = TRUE)
      }
    )
    generated_files <- copy_generated_files(task)
    list(
      output_files = unname(unlist(output_files)),
      generated_files = generated_files
    )
  }

  if (task$renderer == "direct_once") {
    # tcga_roc在当前pROC版本下经由do.call()可能触发call-stack判断错误，
    # 因此这里为它保留显式调用分支。
    direct_stage <- new.env(parent = emptyenv())
    direct_stage$value <- "initializing"
    direct_captured <- with_task_workspace(task, capture_task(tryCatch({
      direct_stage$value <- "get_tcgaplot_function"
      fn <- get_tcgaplot_function(task$function_name)
      direct_stage$value <- "call_tcgaplot_function"
      result <- if (task$function_name == "tcga_roc") {
        # pROC::roc_ used inside TCGAplot::tcga_roc is sensitive to do.call()
        # call stacks in current pROC versions. Use an explicit call here.
        fn(cancer = task$args$cancer, gene = task$args$gene)
      } else {
        do.call(fn, task$args)
      }
      direct_stage$value <- "save_plot"
      output_files <- save_tcgaplot_pdf_png(
        pdf_file = file.path(PLOT_ROOT, paste0(task_stem, ".pdf")),
        width = task$width,
        height = task$height,
        draw_fun = function() {
          render_tcgaplot_output(result, title = task_title, add_outer_title = TRUE)
        }
      )
      direct_stage$value <- "copy_generated_files"
      generated_files <- copy_generated_files(task)
      list(
        output_files = unname(unlist(output_files)),
        generated_files = generated_files
      )
    }, error = function(error) {
      attr(error, "task_stage") <- direct_stage$value
      error
    })))

    if (inherits(direct_captured$result, "error")) {
      cleanup_task_output_files(task)
      return(make_task_result(
        task = task,
        status = "failed",
        warnings = direct_captured$warnings,
        messages = direct_captured$messages,
        error = paste(
          paste0("stage=", attr(direct_captured$result, "task_stage")),
          conditionMessage(direct_captured$result),
          sep = "\n"
        ),
        start_time = task_start_time,
        cache_file = cache_file
      ))
    }

    direct_validation_error <- tryCatch(
      {
        validate_plot_output_files(c(direct_captured$result$output_files, direct_captured$result$generated_files))
        NULL
      },
      error = function(error) error
    )
    if (inherits(direct_validation_error, "error")) {
      cleanup_task_output_files(task)
      return(make_task_result(
        task = task,
        status = "failed",
        warnings = direct_captured$warnings,
        messages = direct_captured$messages,
        error = conditionMessage(direct_validation_error),
        start_time = task_start_time,
        cache_file = cache_file
      ))
    }

    manifest <- make_task_result(
      task = task,
      status = "success",
      output_files = direct_captured$result$output_files,
      generated_files = direct_captured$result$generated_files,
      warnings = direct_captured$warnings,
      messages = direct_captured$messages,
      start_time = task_start_time,
      cache_file = cache_file
    )
    write_tcgaplot_task_cache(task, manifest)
    return(manifest)
  }

  captured <- with_task_workspace(task, capture_task(task_call()))

  if (inherits(captured$result, "error")) {
    cleanup_task_output_files(task)
    return(make_task_result(
      task = task,
      status = "failed",
      warnings = captured$warnings,
      messages = captured$messages,
      error = conditionMessage(captured$result),
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  validation_error <- tryCatch(
    {
      validate_plot_output_files(c(captured$result$output_files, captured$result$generated_files))
      NULL
    },
    error = function(error) error
  )
  if (inherits(validation_error, "error")) {
    cleanup_task_output_files(task)
    return(make_task_result(
      task = task,
      status = "failed",
      warnings = captured$warnings,
      messages = captured$messages,
      error = conditionMessage(validation_error),
      start_time = task_start_time,
      cache_file = cache_file
    ))
  }

  manifest <- make_task_result(
    task = task,
    status = "success",
    output_files = captured$result$output_files,
    generated_files = captured$result$generated_files,
    warnings = captured$warnings,
    messages = captured$messages,
    start_time = task_start_time,
    cache_file = cache_file
  )
  write_tcgaplot_task_cache(task, manifest)
  manifest
}

run_one_special_tcgaplot_task <- function(task_name) {
  task_start_time <- Sys.time()
  captured <- capture_task({
    result <- switch(
      task_name,
      sample_audit = run_sample_audit(),
      data_summary = run_data_summary(),
      stop("Unknown special TCGAplot task: ", task_name)
    )

    rows <- if (is.data.frame(result)) {
      nrow(result)
    } else {
      length(result)
    }

    output_files <- switch(
      task_name,
      sample_audit = file.path(SUMMARY_ROOT, "000_sample_audit_summary.csv"),
      data_summary = c(
        file.path(TABLE_ROOT, "000_data_summary_tcgaplot_cancers.csv"),
        file.path(TABLE_ROOT, "000_data_summary_tcgaplot_paired_cancers.csv"),
        file.path(TABLE_ROOT, "000_data_summary_tcgaplot_package_info.csv")
      ),
      character(0)
    )

    list(rows = rows, output_files = output_files)
  })

  if (inherits(captured$result, "error")) {
    return(make_special_task_result(
      task_name = task_name,
      status = "failed",
      warnings = captured$warnings,
      messages = captured$messages,
      error = conditionMessage(captured$result),
      start_time = task_start_time
    ))
  }

  make_special_task_result(
    task_name = task_name,
    status = "success",
    output_files = captured$result$output_files,
    rows = captured$result$rows,
    warnings = captured$warnings,
    messages = captured$messages,
    start_time = task_start_time
  )
}

normalize_parallel_tcgaplot_results <- function(results, tasks) {
  normalized <- vector("list", length(results))

  for (i in seq_along(results)) {
    result <- results[[i]]
    if (inherits(result, "try-error") || is.null(result)) {
      normalized[[i]] <- make_task_result(
        task = tasks[[i]],
        status = "failed",
        error = paste(as.character(result), collapse = "\n")
      )
    } else {
      normalized[[i]] <- result
    }
  }

  normalized
}

normalize_parallel_special_results <- function(results, task_names) {
  normalized <- vector("list", length(results))

  for (i in seq_along(results)) {
    result <- results[[i]]
    if (inherits(result, "try-error") || is.null(result)) {
      normalized[[i]] <- make_special_task_result(
        task_name = task_names[i],
        status = "failed",
        error = paste(as.character(result), collapse = "\n")
      )
    } else {
      normalized[[i]] <- result
    }
  }

  normalized
}

is_success_status <- function(status) {
  status %in% c("success", "cached", "skipped")
}

make_runtime_summary <- function(special_task_summary, task_summary, start_time) {
  total_runtime_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  special_status <- if (nrow(special_task_summary) > 0L) special_task_summary$Status else character(0)
  plot_status <- if (nrow(task_summary) > 0L) task_summary$Status else character(0)

  data.frame(
    Section = c("special_tasks", "tcgaplot_tasks", "total"),
    Tasks = c(nrow(special_task_summary), nrow(task_summary), nrow(special_task_summary) + nrow(task_summary)),
    Completed_Without_Error = c(
      sum(is_success_status(special_status)),
      sum(is_success_status(plot_status)),
      sum(is_success_status(c(special_status, plot_status)))
    ),
    Failed = c(
      sum(!is_success_status(special_status)),
      sum(!is_success_status(plot_status)),
      sum(!is_success_status(c(special_status, plot_status)))
    ),
    Runtime_Seconds = c(
      sum(suppressWarnings(as.numeric(special_task_summary$Runtime_Seconds)), na.rm = TRUE),
      sum(suppressWarnings(as.numeric(task_summary$Runtime_Seconds)), na.rm = TRUE),
      total_runtime_seconds
    ),
    Runtime = c(
      format_runtime_seconds(sum(suppressWarnings(as.numeric(special_task_summary$Runtime_Seconds)), na.rm = TRUE)),
      format_runtime_seconds(sum(suppressWarnings(as.numeric(task_summary$Runtime_Seconds)), na.rm = TRUE)),
      format_runtime_seconds(total_runtime_seconds)
    ),
    stringsAsFactors = FALSE
  )
}

clear_previous_tcgaplot_run_outputs <- function() {
  if (!CLEAR_PREVIOUS_RUN_OUTPUTS) {
    return(invisible(FALSE))
  }

  unlink(RESULT_ROOT, recursive = TRUE, force = TRUE)
  unlink(TEMP_ROOT, recursive = TRUE, force = TRUE)
  unlink(TCGAPLOT_TASK_CACHE_ROOT, recursive = TRUE, force = TRUE)

  invisible(TRUE)
}


# 8. 主运行入口 ----------------------------------------------------------------

clear_previous_tcgaplot_run_outputs()

dir.create(RESULT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_PNG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TCGAPLOT_REFERENCE_CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TCGAPLOT_TASK_CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)
configure_omnipathr_runtime(
  cache_dir = OMNIPATHR_CACHE_DIR,
  log_dir = OMNIPATHR_LOG_DIR
)

selected_analyses <- resolve_requested_analyses()
plot_analyses <- setdiff(selected_analyses, c("sample_audit", "data_summary"))
tcgaplot_tasks <- build_tcgaplot_tasks(plot_analyses)
tcgaplot_tasks <- assign_task_indices(tcgaplot_tasks)
validate_tcgaplot_tasks(tcgaplot_tasks)

special_task_names <- intersect(c("sample_audit", "data_summary"), selected_analyses)
total_runtime_tasks <- length(special_task_names) + length(tcgaplot_tasks)
parallel_strategy <- setup_parallel_strategy(
  total_tasks = max(total_runtime_tasks, 1L),
  max_workers = PARALLEL_WORKERS,
  inner_label = "TCGAplot inner workers",
  nested_label = "Nested workers"
)

if (QUICKANALYSIS_VERBOSE) {
  cat("\nTCGAplot quick analysis configuration:\n")
  cat("TCGAplot version: ", as.character(utils::packageVersion("TCGAplot")), "\n", sep = "")
  cat("Target genes: ", paste(TARGET_GENES, collapse = ", "), "\n", sep = "")
  cat("Target cancers: ", paste(TARGET_CANCERS, collapse = ", "), "\n", sep = "")
  cat("Selected analyses: ", paste(selected_analyses, collapse = ", "), "\n", sep = "")
  cat("Prepared special tasks: ", length(special_task_names), "\n", sep = "")
  cat("Prepared plotting/statistical tasks: ", length(tcgaplot_tasks), "\n", sep = "")
  cat("Result root: ", RESULT_ROOT, "\n", sep = "")
  cat("Temporary root: ", TEMP_ROOT, "\n", sep = "")
  cat("TCGAplot data/cache root: ", DATA_ROOT, "\n", sep = "")
  cat("\nWriting TCGAplot analysis catalog...\n")
}
write_analysis_catalog(selected_analyses)

special_task_summary <- if (length(special_task_names) > 0L) {
  if (QUICKANALYSIS_VERBOSE) {
    cat("\nRunning TCGAplot setup/data tasks...\n")
  }
  special_results <- run_indexed_tasks_with_progress(
    total_tasks = length(special_task_names),
    task_function = function(task_id) {
      run_one_special_tcgaplot_task(special_task_names[task_id])
    },
    workers = min(parallel_strategy$task_workers, length(special_task_names)),
    progress_label = "TCGAplot setup"
  )
  special_results <- normalize_parallel_special_results(
    results = special_results,
    task_names = special_task_names
  )
  do.call(rbind, special_results)
} else {
  data.frame()
}

if (nrow(special_task_summary) > 0L) {
  write_table(special_task_summary, SUMMARY_ROOT, "000_tcgaplot_special_task_summary")
}

if (QUICKANALYSIS_VERBOSE) {
  cat("\nTCGAplot plotting/statistical tasks: ", length(tcgaplot_tasks), "\n", sep = "")
}

task_summary <- if (length(tcgaplot_tasks) > 0L) {
  raw_task_results <- run_indexed_tasks_with_progress(
    total_tasks = length(tcgaplot_tasks),
    task_function = function(task_id) {
      run_one_tcgaplot_task(tcgaplot_tasks[[task_id]])
    },
    workers = min(parallel_strategy$task_workers, length(tcgaplot_tasks)),
    progress_label = "TCGAplot tasks"
  )
  normalized_task_results <- normalize_parallel_tcgaplot_results(
    results = raw_task_results,
    tasks = tcgaplot_tasks
  )
  do.call(rbind, normalized_task_results)
} else {
  data.frame()
}

if (nrow(task_summary) > 0) {
  write_table(task_summary, SUMMARY_ROOT, "000_tcgaplot_task_summary")
  failed_tasks <- task_summary[!is_success_status(task_summary$Status), , drop = FALSE]
  write_table(failed_tasks, SUMMARY_ROOT, "000_tcgaplot_failed_tasks")
}

runtime_summary <- make_runtime_summary(
  special_task_summary = special_task_summary,
  task_summary = task_summary,
  start_time = SCRIPT_START_TIME
)
write_table(runtime_summary, SUMMARY_ROOT, "000_tcgaplot_runtime_summary")

if (QUICKANALYSIS_VERBOSE) {
  cat("\nTCGAplot quick analysis finished.\n")
  if (nrow(special_task_summary) > 0) {
    cat("Special tasks successful: ", sum(is_success_status(special_task_summary$Status)), "\n", sep = "")
    cat("Special tasks failed: ", sum(!is_success_status(special_task_summary$Status)), "\n", sep = "")
    cat("Special task summary: ", file.path(SUMMARY_ROOT, "000_tcgaplot_special_task_summary.csv"), "\n", sep = "")
  }
  if (nrow(task_summary) > 0) {
    cat("Completed without error: ", sum(is_success_status(task_summary$Status)), "\n", sep = "")
    cat("Failed tasks: ", sum(!is_success_status(task_summary$Status)), "\n", sep = "")
    cat("Task summary: ", file.path(SUMMARY_ROOT, "000_tcgaplot_task_summary.csv"), "\n", sep = "")
    cat("Failed task summary: ", file.path(SUMMARY_ROOT, "000_tcgaplot_failed_tasks.csv"), "\n", sep = "")
  }
  if ("sample_audit" %in% selected_analyses) {
    cat("Sample audit summary: ", file.path(SUMMARY_ROOT, "000_sample_audit_summary.csv"), "\n", sep = "")
  }
  cat("Runtime summary: ", file.path(SUMMARY_ROOT, "000_tcgaplot_runtime_summary.csv"), "\n", sep = "")
  print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
}

cat("\n01 TCGAplot quick analysis finished: ", RESULT_ROOT, "\n", sep = "")
