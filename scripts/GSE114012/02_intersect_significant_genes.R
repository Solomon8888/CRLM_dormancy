# GSE114012显著差异基因交集分析
#
# 读取01号limma脚本生成的DEG结果，按配置计算一套或多套显著基因交集。
# 交集结果统一保存到：
# results/ngs/GSE114012/intersect/<intersection_scheme>/
#
# 推荐用法：先修改第0节配置，再按需整段运行或逐节交互式运行。


# 0. 可修改配置 ---------------------------------------------------------------

DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"
PARALLEL_FUNCTION_FILE <- "scripts/functions/parallel_runtime_functions.R"

TABLE_ROOT <- file.path("results", DATA_TYPE, DATASET_ID, "tables")
DEG_DIR_NAME <- "DEG"
INTERSECT_DIR_NAME <- "intersect"

RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
INTERSECT_ROOT <- file.path(RESULT_ROOT, INTERSECT_DIR_NAME)

# 交集使用的基因ID列和上下调判断列。
# 01号脚本的最终CSV不再保存Feature_ID，因此这里使用GeneID。
GENE_ID_COLUMN <- "GeneID"
LOGFC_COLUMN <- "logFC"

# 交集基因列表保留的注释列；不包含logFC、P值等差异分析统计量。
GENE_ANNOTATION_COLUMNS <- c(
  "GeneID", "Symbol", "Ensembl", "Entrez"
)

# 最终输出表中不保存的辅助注释列。
OUTPUT_DROP_COLUMNS <- c("Feature_ID", "Biotype", "Length")

# 多套交集方案。
# 新增方案时，只需要按下面格式增加一行；顺序会影响交集基因的输出排序。
INTERSECTION_SCHEMES <- list(
  DLD1_HCT15_SW48 = c("DLD1", "HCT15", "SW48"),
  DLD1_HCT15 = c("DLD1", "HCT15"),
  HT55_SW948 = c("HT55", "SW948"),
  SW948_RKO = c("SW948", "RKO")
)

# 运行哪些交集方案。
# 可设为names(INTERSECTION_SCHEMES)运行全部；也可只写部分方案名。
SCHEMES_TO_RUN <- names(INTERSECTION_SCHEMES)
# SCHEMES_TO_RUN <- c("DLD1_HCT15_SW48")

# 重跑时清空对应交集方案目录，避免旧文件与新结果混在一起。
OVERWRITE_SCHEME_OUTPUT <- TRUE

# 清理旧版脚本生成在tables/<analysis_name>/intersect/下的分散结果。
REMOVE_LEGACY_ANALYSIS_INTERSECT_DIRS <- TRUE

# 清理旧版脚本生成在tables/intersect/下的集中结果。
REMOVE_LEGACY_TABLE_INTERSECT_DIR <- TRUE

# 交集后的基因在单个差异分析结果中展示。
# TRUE时，DLD1_HCT15_SW48会展开为DLD1、HCT15、SW48，不输出带"_"的组合分析目录。
OUTPUT_SINGLE_ANALYSIS_ONLY <- TRUE

options(width = 200)


# 1. 常用函数 -----------------------------------------------------------------

source(REPORT_TABLE_FUNCTION_FILE)
source(PARALLEL_FUNCTION_FILE)

SCRIPT_START_TIME <- start_runtime_timer()

sanitize_file_name <- function(x) {
  # 文件夹名保留字母、数字、下划线、点和短横线，避免不同系统下路径出错。
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

get_direction <- function(logfc) {
  # 交集后的方向统计只看logFC正负；显著性阈值已在01号脚本中完成。
  ifelse(logfc > 0, "Up", ifelse(logfc < 0, "Down", "No_change"))
}

get_analysis_name_from_deg_file <- function(file_name) {
  # 兼容tables/<analysis_name>/DEG/*.csv和tables/<analysis_name>/DEG/csv/*.csv。
  if (basename(dirname(file_name)) == "csv") {
    return(basename(dirname(dirname(dirname(file_name)))))
  }

  basename(dirname(dirname(file_name)))
}

get_deg_file_info <- function(table_root) {
  # 同时定位显著基因表和全基因表；交集基因的DEG展示要从all_genes中截取。
  significant_files <- list.files(
    table_root,
    pattern = "significant_genes[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  all_genes_files <- list.files(
    table_root,
    pattern = "all_genes[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  significant_files <- significant_files[
    basename(dirname(significant_files)) == DEG_DIR_NAME |
      (
        basename(dirname(significant_files)) == "csv" &
          basename(dirname(dirname(significant_files))) == DEG_DIR_NAME
      )
  ]
  all_genes_files <- all_genes_files[
    basename(dirname(all_genes_files)) == DEG_DIR_NAME |
      (
        basename(dirname(all_genes_files)) == "csv" &
          basename(dirname(dirname(all_genes_files))) == DEG_DIR_NAME
      )
  ]

  significant_files <- prefer_report_csv_files(
    significant_files,
    get_analysis_name_from_deg_file
  )
  all_genes_files <- prefer_report_csv_files(
    all_genes_files,
    get_analysis_name_from_deg_file
  )

  significant_files <- significant_files[
    !grepl("_intersect_", basename(significant_files))
  ]
  all_genes_files <- all_genes_files[
    !grepl("_intersect_", basename(all_genes_files))
  ]

  stopifnot(length(significant_files) > 0)
  stopifnot(length(all_genes_files) > 0)

  significant_info <- data.frame(
    Analysis_Name = vapply(
      significant_files,
      get_analysis_name_from_deg_file,
      character(1)
    ),
    Significant_File = significant_files,
    stringsAsFactors = FALSE
  )

  all_genes_info <- data.frame(
    Analysis_Name = vapply(
      all_genes_files,
      get_analysis_name_from_deg_file,
      character(1)
    ),
    All_Genes_File = all_genes_files,
    stringsAsFactors = FALSE
  )

  if (any(duplicated(significant_info$Analysis_Name))) {
    duplicated_names <- unique(significant_info$Analysis_Name[
      duplicated(significant_info$Analysis_Name)
    ])
    stop(
      "More than one significant_genes file was found for: ",
      paste(duplicated_names, collapse = ", ")
    )
  }

  if (any(duplicated(all_genes_info$Analysis_Name))) {
    duplicated_names <- unique(all_genes_info$Analysis_Name[
      duplicated(all_genes_info$Analysis_Name)
    ])
    stop(
      "More than one all_genes file was found for: ",
      paste(duplicated_names, collapse = ", ")
    )
  }

  file_info <- merge(
    significant_info,
    all_genes_info,
    by = "Analysis_Name",
    all = FALSE
  )

  file_info <- file_info[order(file_info$Analysis_Name), , drop = FALSE]
  rownames(file_info) <- NULL

  file_info
}

read_deg_table <- function(gene_file) {
  gene_file <- resolve_report_csv_file(gene_file)
  dat <- read.csv(
    gene_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  stopifnot(GENE_ID_COLUMN %in% colnames(dat))
  stopifnot(LOGFC_COLUMN %in% colnames(dat))

  dat <- dat[
    !is.na(dat[[GENE_ID_COLUMN]]) & dat[[GENE_ID_COLUMN]] != "",
    ,
    drop = FALSE
  ]

  dat$Direction <- get_direction(dat[[LOGFC_COLUMN]])
  dat
}

prepare_output_table <- function(dat) {
  keep_columns <- setdiff(colnames(dat), OUTPUT_DROP_COLUMNS)
  dat <- dat[, keep_columns, drop = FALSE]

  # Entrez是基因ID，不作为连续数值展示，避免写出为9.62e+03一类科学计数法。
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

make_gene_annotation_table <- function(intersect_gene_ids, selected_analyses, significant_tables) {
  first_table <- significant_tables[[selected_analyses[1]]]

  annotation_columns <- intersect(GENE_ANNOTATION_COLUMNS, colnames(first_table))
  annotation_columns <- unique(c(GENE_ID_COLUMN, annotation_columns))

  first_table[
    match(intersect_gene_ids, first_table[[GENE_ID_COLUMN]]),
    annotation_columns,
    drop = FALSE
  ]
}

get_deg_output_analyses <- function(selected_analyses) {
  if (OUTPUT_SINGLE_ANALYSIS_ONLY) {
    selected_analyses <- unlist(strsplit(selected_analyses, "_", fixed = TRUE))
  }

  unique(selected_analyses)
}


# 2. 查找差异分析结果 ----------------------------------------------------------

dir.create(INTERSECT_ROOT, recursive = TRUE, showWarnings = FALSE)

if (REMOVE_LEGACY_TABLE_INTERSECT_DIR) {
  legacy_table_intersect_dir <- file.path(TABLE_ROOT, INTERSECT_DIR_NAME)
  if (dir.exists(legacy_table_intersect_dir)) {
    unlink(legacy_table_intersect_dir, recursive = TRUE)
  }
}

if (REMOVE_LEGACY_ANALYSIS_INTERSECT_DIRS) {
  legacy_intersect_dirs <- file.path(
    list.dirs(TABLE_ROOT, recursive = FALSE, full.names = TRUE),
    INTERSECT_DIR_NAME
  )
  legacy_intersect_dirs <- legacy_intersect_dirs[
    dir.exists(legacy_intersect_dirs)
  ]

  if (length(legacy_intersect_dirs) > 0) {
    unlink(legacy_intersect_dirs, recursive = TRUE)
  }
}

file_info <- get_deg_file_info(TABLE_ROOT)

stopifnot(length(INTERSECTION_SCHEMES) > 0)
stopifnot(!is.null(names(INTERSECTION_SCHEMES)))
stopifnot(!any(names(INTERSECTION_SCHEMES) == ""))
stopifnot(all(SCHEMES_TO_RUN %in% names(INTERSECTION_SCHEMES)))

active_schemes <- INTERSECTION_SCHEMES[SCHEMES_TO_RUN]
intersection_analyses <- unique(unlist(active_schemes))

deg_output_analyses_by_scheme <- lapply(active_schemes, get_deg_output_analyses)
required_analyses <- unique(c(
  intersection_analyses,
  unique(unlist(deg_output_analyses_by_scheme))
))

missing_analyses <- setdiff(required_analyses, file_info$Analysis_Name)
if (length(missing_analyses) > 0) {
  stop(
    "No DEG file was found for: ",
    paste(missing_analyses, collapse = ", ")
  )
}

selected_file_info <- file_info[
  match(required_analyses, file_info$Analysis_Name),
  ,
  drop = FALSE
]


# 3. 读取差异分析结果 ----------------------------------------------------------

significant_tables <- vector("list", length(intersection_analyses))
names(significant_tables) <- intersection_analyses

all_gene_tables <- vector("list", length(required_analyses))
names(all_gene_tables) <- required_analyses

all_gene_columns <- vector("list", length(required_analyses))
names(all_gene_columns) <- required_analyses

for (analysis_name in intersection_analyses) {
  gene_file <- selected_file_info$Significant_File[
    selected_file_info$Analysis_Name == analysis_name
  ]

  significant_tables[[analysis_name]] <- read_deg_table(gene_file)
}

for (analysis_name in required_analyses) {
  gene_file <- selected_file_info$All_Genes_File[
    selected_file_info$Analysis_Name == analysis_name
  ]

  dat <- read_deg_table(gene_file)

  all_gene_tables[[analysis_name]] <- dat
  all_gene_columns[[analysis_name]] <- setdiff(colnames(dat), "Direction")
}


# 4. 输入结果概览 --------------------------------------------------------------

input_summary <- do.call(rbind, lapply(intersection_analyses, function(analysis_name) {
  dat <- significant_tables[[analysis_name]]

  data.frame(
    Analysis_Name = analysis_name,
    Total_Significant_Genes = nrow(dat),
    Up = sum(dat$Direction == "Up"),
    Down = sum(dat$Direction == "Down"),
    No_change = sum(dat$Direction == "No_change"),
    stringsAsFactors = FALSE
  )
}))

cat("\nInput significant gene summary:\n")
print(input_summary, row.names = FALSE)


# 5. 按方案计算交集并保存 ------------------------------------------------------

run_one_intersection_scheme <- function(scheme_index) {
  intersection_name <- names(active_schemes)[scheme_index]
  selected_analyses <- unique(active_schemes[[scheme_index]])
  deg_output_analyses <- get_deg_output_analyses(selected_analyses)

  stopifnot(length(selected_analyses) >= 2)

  safe_intersection_name <- sanitize_file_name(intersection_name)
  intersection_dir <- file.path(INTERSECT_ROOT, safe_intersection_name)

  if (OVERWRITE_SCHEME_OUTPUT && dir.exists(intersection_dir)) {
    unlink(intersection_dir, recursive = TRUE)
  }
  dir.create(intersection_dir, recursive = TRUE, showWarnings = FALSE)

  gene_id_list <- lapply(selected_analyses, function(analysis_name) {
    significant_tables[[analysis_name]][[GENE_ID_COLUMN]]
  })

  intersect_gene_set <- Reduce(intersect, gene_id_list)
  first_gene_ids <- significant_tables[[selected_analyses[1]]][[GENE_ID_COLUMN]]
  intersect_gene_ids <- first_gene_ids[first_gene_ids %in% intersect_gene_set]

  if (length(intersect_gene_ids) > 0) {
    direction_matrix <- do.call(cbind, lapply(selected_analyses, function(analysis_name) {
      dat <- significant_tables[[analysis_name]]
      dat$Direction[match(intersect_gene_ids, dat[[GENE_ID_COLUMN]])]
    }))
    colnames(direction_matrix) <- selected_analyses

    common_up <- apply(direction_matrix, 1, function(x) all(x == "Up"))
    common_down <- apply(direction_matrix, 1, function(x) all(x == "Down"))
  } else {
    common_up <- logical(0)
    common_down <- logical(0)
  }

  gene_annotation_table <- make_gene_annotation_table(
    intersect_gene_ids = intersect_gene_ids,
    selected_analyses = selected_analyses,
    significant_tables = significant_tables
  )
  gene_annotation_table <- prepare_output_table(gene_annotation_table)

  gene_list_file <- file.path(
    intersection_dir,
    "gene_list.csv"
  )

  gene_list_file <- write_csv_with_report_previews(gene_annotation_table, gene_list_file)

  saved_gene_annotation_table <- read.csv(
    gene_list_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  stopifnot(nrow(saved_gene_annotation_table) == length(intersect_gene_ids))
  stopifnot(!any(OUTPUT_DROP_COLUMNS %in% colnames(saved_gene_annotation_table)))

  scheme_summary <- data.frame(
    Selected_Analyses = paste(selected_analyses, collapse = ";"),
    DEG_Result_Analyses = paste(deg_output_analyses, collapse = ";"),
    Total_Intersected_Genes = length(intersect_gene_ids),
    Common_Up = sum(common_up),
    Common_Down = sum(common_down),
    Mixed_Direction = length(intersect_gene_ids) - sum(common_up) - sum(common_down),
    stringsAsFactors = FALSE
  )

  for (analysis_name in deg_output_analyses) {
    dat <- all_gene_tables[[analysis_name]]
    match_index <- match(intersect_gene_ids, dat[[GENE_ID_COLUMN]])
    stopifnot(!any(is.na(match_index)))

    intersect_deg <- dat[
      match_index,
      all_gene_columns[[analysis_name]],
      drop = FALSE
    ]
    intersect_deg <- prepare_output_table(intersect_deg)

    analysis_output_dir <- file.path(intersection_dir, sanitize_file_name(analysis_name))
    dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

    deg_output_file <- file.path(
      analysis_output_dir,
      "deg_results.csv"
    )

    deg_output_file <- write_csv_with_report_previews(intersect_deg, deg_output_file)

    saved_intersect_deg <- read.csv(
      deg_output_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    stopifnot(nrow(saved_intersect_deg) == length(intersect_gene_ids))
    stopifnot(!any(OUTPUT_DROP_COLUMNS %in% colnames(saved_intersect_deg)))

  }

  scheme_summary_file <- file.path(
    intersection_dir,
    "summary.csv"
  )

  scheme_summary_file <- write_csv_with_report_previews(scheme_summary, scheme_summary_file)

  saved_scheme_summary <- read.csv(
    scheme_summary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  stopifnot(nrow(saved_scheme_summary) == 1)
  stopifnot(saved_scheme_summary$Total_Intersected_Genes == length(intersect_gene_ids))

  scheme_summary
}

cat("\nRunning intersection schemes...\n")
parallel_strategy <- setup_parallel_strategy(
  total_tasks = length(active_schemes),
  inner_label = "Intersection inner workers",
  nested_label = "Nested workers"
)

scheme_summary_list <- run_indexed_tasks_with_progress(
  total_tasks = length(active_schemes),
  workers = parallel_strategy$task_workers,
  task_function = run_one_intersection_scheme
)
stop_on_parallel_errors(
  scheme_summary_list,
  task_ids = names(active_schemes),
  label = "intersection schemes"
)
names(scheme_summary_list) <- names(active_schemes)


# 6. 终端快速汇总 --------------------------------------------------------------

scheme_summary_table <- do.call(rbind, scheme_summary_list)
rownames(scheme_summary_table) <- names(scheme_summary_list)

cat("\nIntersection scheme summary:\n")
print(
  scheme_summary_table[
    ,
    c(
      "Total_Intersected_Genes", "Common_Up", "Common_Down", "Mixed_Direction",
      "Selected_Analyses", "DEG_Result_Analyses"
    )
  ],
  row.names = TRUE
)

cat("\nIntersected gene lists and DEG tables were saved in:\n")
cat(file.path(INTERSECT_ROOT, "<INTERSECTION_SCHEME>"), "\n", sep = "")
cat("\nSignificant DEG intersection analysis finished.\n")
print_runtime_summary(SCRIPT_START_TIME, label = "Total runtime")
