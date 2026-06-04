# GSE114012批量GSEA分析
#
# 自动读取01号脚本输出的每套DEG/all_genes.csv，
# 使用clusterProfiler::GSEA和msigdbr基因集批量运行GSEA分析。
# GSEA结果表保存为csv/md/tex三种完整格式，Top通路气泡图使用GseaVis绘制。

SCRIPT_START_TIME <- Sys.time()


# 0. 可修改配置 ---------------------------------------------------------------

# 当前脚本服务于GSE114012这个NGS数据集；目录结构依赖这两个字段。
DATASET_ID <- "GSE114012"
DATA_TYPE <- "ngs"

# 输入为tables/<analysis_name>/DEG/all_genes.csv；
# 单通路GSEA图会额外读取SE对象中的TPM矩阵，用于绘制核心富集基因表达热图。
# GSEA表格输出到tables/<analysis_name>/GSEA/，图片输出到plots/GSEA/<analysis_name>/。
SE_RDS_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_se_raw.rds"
CLINICAL_FILE <- "data/ngs/GSE114012/data_prepare/GSE114012_clinical_edit.csv"
FUNCTION_FILE <- "scripts/functions/limma_de_functions.R"
RESULT_ROOT <- file.path("results", DATA_TYPE, DATASET_ID)
TABLE_ROOT <- file.path(RESULT_ROOT, "tables")

# 默认输出到results目录。测试脚本性能时可临时设置环境变量GSEA_OUTPUT_ROOT，
# 将GSEA表格和图片输出到temporary目录，避免覆盖正式全量结果。
OUTPUT_ROOT <- Sys.getenv("GSEA_OUTPUT_ROOT", unset = RESULT_ROOT)
TABLE_OUTPUT_ROOT <- file.path(OUTPUT_ROOT, "tables")
PLOT_ROOT <- file.path(OUTPUT_ROOT, "plots", "GSEA")

# 公共函数脚本：差异分析文件定位、绘图保存、统一绘图风格、文本表格转义。
PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
REPORT_TABLE_FUNCTION_FILE <- "scripts/functions/report_table_functions.R"

# 需要运行哪些差异分析设计。设为"all"时自动运行全部DEG/all_genes.csv。
ANALYSES_TO_RUN <- "all"
# ANALYSES_TO_RUN <- c("DLD1", "DLD1_HCT15_SW48")

# 物种配置。常用可选项："human"、"Mus musculus"、"Rattus norvegicus"。
SPECIES <- "human"

# GSEA排序基因列表配置。
# 默认使用limma的t统计量排序；若希望更直接展示表达变化幅度，可改为"logFC"。
GENE_ID_TYPE <- "ENTREZ"  # 可选："ENTREZ", "SYMBOL", "ENSEMBL"
RANK_METRIC_COLUMN <- "t"

# clusterProfiler::GSEA官方运算参数。
# 这些参数会直接传入clusterProfiler::GSEA；不同版本的clusterProfiler对method等参数
# 有明确默认值，这里保留官方推荐默认口径，并在每个参数旁说明其含义和常见调整方向。
GSEA_PARAMS <- list(
  # exponent：加权GSEA的权重指数。1是经典加权GSEA；
  # 设为0时接近不加权KS统计，设为>1会更强调排序靠前且统计量绝对值更大的基因。
  exponent = 1,

  # minGSSize：参与检验的基因集最小基因数。调小可保留更小通路，但稳定性下降；
  # 调大可减少噪音，但可能丢掉精细生物过程。
  minGSSize = 5,

  # maxGSSize：参与检验的基因集最大基因数。过大的基因集通常过于宽泛；
  # 500是GSEA常用上限，适合GO/MSigDB批量分析。
  maxGSSize = 500,

  # pvalueCutoff：clusterProfiler官方手册中说明为adjusted pvalue cutoff；
  # 当前0.05只返回校正后显著通路。若想完整保留全部通路用于二次筛选，可改为1.0。
  pvalueCutoff = 0.05,

  # pAdjustMethod：多重检验校正方法。"BH"控制FDR，是转录组富集分析常用选择；
  # 也可改为"bonferroni"等更严格方法。
  pAdjustMethod = "BH",

  # verbose：是否打印clusterProfiler内部运行信息。批量脚本中通常TRUE便于追踪进度。
  verbose = TRUE,

  # nPerm：传统置换法的置换次数。当前method="multilevel"时主要保留兼容官方接口；
  # 若改用置换法，可适当增大以提高P值精度，但运行时间会显著增加。
  nPerm = 1000,

  # method：GSEA P值估计方法。clusterProfiler官方可选"multilevel"、"monte carlo"、"fgsea"。
  # "multilevel"为当前官方默认，适合批量分析且对较小P值估计更稳定；
  # "monte carlo"偏传统随机置换思路，结果精度更依赖nPerm；
  # "fgsea"用于兼容旧fgsea风格流程，适合需要与旧分析复现保持一致时使用。
  method = "multilevel",

  # adaptive：fgsea相关的自适应置换开关。官方默认FALSE；通常保持默认即可。
  adaptive = FALSE,

  # minPerm/maxPerm：自适应/多层算法允许的置换下限和上限。
  # 增大maxPerm可提升极小P值精度，但会明显增加耗时。
  minPerm = 101,
  maxPerm = 1e5,

  # pvalThreshold：fgsea内部进一步精化P值估计的阈值；
  # 较小值会减少精化范围，较大值可能增加计算量。
  pvalThreshold = 0.1
)

# 当前批量运行哪些MSigDB基因集。
# 常用写法示例："H"、"C1"、"C2:CP:REACTOME"、"C5:GO:BP"、"C5:GO:CC"、
# "C5:GO:MF"、"C6"、"C7:IMMUNESIGDB"、"C8"、"C9"。
# 设为"all"时，会自动运行msigdbr当前数据库里全部可用基因集类别。
GSEA_GENESETS_TO_RUN <- "all"
# GSEA_GENESETS_TO_RUN <- c("H", "C5:GO:BP", "C6")

# GseaVis::dotplotGsea气泡图配置。
# SIMPLIFY_PATHWAY_PREFIX_IN_PLOT只影响图片标签，不改变GSEA官方结果表。
SIMPLIFY_PATHWAY_PREFIX_IN_PLOT <- TRUE
REPLACE_UNDERSCORE_WITH_SPACE_IN_PLOT <- TRUE

GSEAVIS_DOTPLOT_PARAMS <- list(
  # topn：每个富集方向展示的Top通路数量；NULL则展示通过阈值的全部通路。
  topn = 20,

  # pval：按原始P值过滤；NULL表示不用原始P值过滤，只按pajust等条件绘图。
  pval = NULL,

  # pajust：按校正P值过滤；只影响图片中展示哪些通路，不改变CSV/MD/TEX结果表。
  pajust = 0.05,

  # order.by：气泡横坐标和Top排序字段。常用"GeneRatio"或"NES"。
  # GeneRatio更直观显示核心富集基因比例，NES更强调富集强度。
  order.by = "GeneRatio",

  # str.width：通路名称单行最大字符数。超过该长度自动换行；45较适合汇报图。
  str.width = 45,

  # base_size：GseaVis基础字号；最终仍会叠加本项目Helvetica黑色粗体风格。
  base_size = 10,

  # scales：分面横坐标缩放。"free_x"让activated/suppressed各自适配；
  # "fixed"则强制两侧统一横坐标，便于严格比较。
  scales = "free_x",

  # add.seg：是否给气泡增加辅助线。FALSE更清爽；TRUE适合通路较少时增强定位。
  add.seg = FALSE,

  # line.col/line.size/line.type：辅助线颜色、粗细、线型；仅add.seg=TRUE时明显生效。
  line.col = "grey80",
  line.size = 1.5,
  line.type = "solid"
)

# GseaVis气泡大小范围。这里在GseaVis默认基础上整体放大一点。
GSEAVIS_POINT_SIZE_RANGE <- c(4.2, 9.8)

# 气泡图导出尺寸。
# 通路名会按str.width换行；图片高度按实际通路名总行数动态增加，避免多行通路名重叠。
DOTPLOT_BODY_BASE_SIZE <- 4.8
DOTPLOT_LABEL_LINE_HEIGHT <- 0.34
DOTPLOT_TERM_GAP_HEIGHT <- 0.14
DOTPLOT_BODY_MIN_SIZE <- 5.2
DOTPLOT_BODY_MAX_SIZE <- 18.0
DOTPLOT_LABEL_BASE_WIDTH <- 1.8
DOTPLOT_LABEL_WIDTH_PER_CHARACTER <- 0.045
DOTPLOT_LABEL_MIN_WIDTH <- 2.4
DOTPLOT_LABEL_MAX_WIDTH <- 6.2
DOTPLOT_LEGEND_WIDTH <- 1.4
DOTPLOT_VERTICAL_PADDING <- 0.45

# 单通路GSEA图配置。默认绘制每个分析设计、每个基因集中的Top20显著通路。
# 这里与上方GseaVis气泡图topn保持一致，便于一张dotplot对应一批单通路细节图。
DRAW_SINGLE_PATHWAY_GSEA <- TRUE
SINGLE_PATHWAY_TOP_N <- GSEAVIS_DOTPLOT_PARAMS$topn

# 若填写关键词，则不再按Top20，而是在当前基因集的显著GSEA结果中按ID/Description搜索。
# 例如：SINGLE_PATHWAY_KEYWORDS <- c("TGF", "IL6", "WNT")
SINGLE_PATHWAY_KEYWORDS <- character(0)
SINGLE_PATHWAY_MAX_KEYWORD_TERMS <- 20
SINGLE_PATHWAY_PVALUE_COLUMN <- "p.adjust"
SINGLE_PATHWAY_PVALUE_CUTOFF <- 0.05

# 单通路图表达热图使用TPM矩阵；默认log2(TPM + 1)后交给GseaVis按行Z-score缩放。
TPM_ASSAY_NAME <- "tpm"
SINGLE_PATHWAY_SAMPLE_LABEL_COLUMN <- "Title"

# GseaVis::gseaNb单通路图配置。
GSEAVIS_SINGLE_PATHWAY_PARAMS <- list(
  # subPlot：传统GSEA子图数量。2表示富集曲线+hit分布；3会额外展示ranked list。
  # 这里配合newGsea=TRUE使用新式GSEA曲线，画面更简洁。
  subPlot = 2,

  # lineSize：富集曲线、基准线等线条粗细。
  lineSize = 1.0,

  # rmSegment：是否去掉被标注基因到曲线的线段；FALSE保留线段。
  rmSegment = FALSE,

  # termWidth：单通路标题换行宽度；过长通路名会按该宽度换行。
  termWidth = 45,

  # segCol：标注基因线段颜色。
  segCol = "black",

  # curveCol：旧式GSEA曲线颜色；newGsea=TRUE时主要使用newCurveCol。
  curveCol = c("#2166AC", "#D73027", "#762A83"),

  # htCol：hit分布条颜色；newGsea=TRUE时主要使用newHtCol。
  htCol = c("#2166AC", "#D73027"),

  # rankCol：ranked list柱图渐变色；subPlot=3时更明显。
  rankCol = c("#2166AC", "white", "#D73027"),

  # rankSeq：ranked list横坐标刻度间隔。
  rankSeq = 5000,

  # htHeight：hit分布条高度。
  htHeight = 0.30,

  # force/max.overlaps/geneSize：ggrepel标注基因时的排斥力、最大重叠、字号。
  force = 20,
  max.overlaps = 50,
  geneSize = 4,

  # newGsea：使用GseaVis新式GSEA曲线样式。
  newGsea = TRUE,

  # addPoint：在新式GSEA曲线上叠加点。
  addPoint = TRUE,

  # newCurveCol/newHtCol：新式GSEA曲线和hit条渐变配色，沿用本项目红蓝风格。
  newCurveCol = c("#2166AC", "white", "#D73027"),
  newHtCol = c("#2166AC", "white", "#D73027"),

  # rmHt：是否移除hit分布条。FALSE保留。
  rmHt = FALSE,

  # addPval：是否在图中写NES/P值。TRUE更适合汇报单图。
  addPval = TRUE,
  pvalX = 0.96,
  pvalY = 0.92,
  pvalSize = 4,
  pCol = "black",
  pHjust = 1,

  # rmPrefix：标题是否去掉HALLMARK/GO等前缀；只影响图片标题。
  rmPrefix = TRUE,
  nesDigit = 2,
  pDigit = 3,

  # markTopgene/topGeneN：是否自动标注曲线附近Top基因。
  markTopgene = FALSE,
  topGeneN = 5,

  # legend.position：多通路曲线图图例位置。当前单通路图通常不显示曲线图例。
  legend.position = "right",

  # add.geneExpHt：增加核心富集基因的表达热图注释。
  add.geneExpHt = TRUE,

  # scale.exp：表达热图是否按基因做Z-score缩放。
  scale.exp = TRUE,

  # exp.col：表达热图颜色，保持本项目红蓝风格。
  exp.col = c("#2166AC", "white", "#D73027"),

  # ht.legend：是否展示表达热图图例。
  ht.legend = TRUE,

  # ght.relHight：表达热图占单通路图总高度的比例。
  ght.relHight = 0.42,

  # ght.geneText.size：表达热图核心基因名字号。
  ght.geneText.size = 7,

  # ght.facet/ght.facet.scale：多个通路合并画热图时的分面设置；当前逐通路输出，默认FALSE。
  ght.facet = FALSE,
  ght.facet.scale = "free"
)

SINGLE_PATHWAY_PLOT_BASE_SIZE <- 7.2
SINGLE_PATHWAY_PLOT_SIZE_PER_GENE <- 0.090
SINGLE_PATHWAY_PLOT_SIZE_PER_SAMPLE <- 0.115
SINGLE_PATHWAY_PLOT_SIZE_PER_TITLE_LINE <- 0.32
SINGLE_PATHWAY_PLOT_MIN_SIZE <- 7.2
SINGLE_PATHWAY_PLOT_MAX_SIZE <- 34.0
SINGLE_PATHWAY_GENE_TEXT_MIN_SIZE <- 2.2
SINGLE_PATHWAY_GENE_TEXT_MAX_SIZE <- 7.0
SINGLE_PATHWAY_GENE_TEXT_WIDTH_FACTOR <- 0.62
SINGLE_PATHWAY_SAMPLE_TEXT_MIN_SIZE <- 2.6
SINGLE_PATHWAY_SAMPLE_TEXT_MAX_SIZE <- 6.6
SINGLE_PATHWAY_SAMPLE_TEXT_HEIGHT_FACTOR <- 0.74

options(width = 200)
options(lifecycle_verbosity = "quiet")


# 1. 加载包、公共函数和内部声明 ------------------------------------------------

required_packages <- c(
  "SummarizedExperiment",
  "clusterProfiler",
  "msigdbr",
  "GseaVis",
  "ggplot2",
  "qs2",
  "parallel"
)

is_package_available <- function(package_name) {
  suppressWarnings(
    suppressPackageStartupMessages(
      requireNamespace(package_name, quietly = TRUE)
    )
  )
}

missing_packages <- required_packages[
  !vapply(required_packages, is_package_available, logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Please install required R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

source(FUNCTION_FILE)
source(PLOTTING_FUNCTION_FILE)
source(REPORT_TABLE_FUNCTION_FILE)

# 以下为脚本内部固定声明，通常不需要在日常分析中修改。
CLEAN_GSEA_OUTPUT_DIR <- TRUE
READABLE_GENE_SYMBOLS <- TRUE
USE_QS2_CACHE <- TRUE

# 默认重新计算核心GSEA缓存，适合每次调整参数后重新运行。
# 若确认参数完全不变且只想复用缓存，可在终端临时设置：
# GSEA_REFRESH_QS2_CACHE=0 Rscript scripts/GSE114012/06_gsea_analysis.R
REFRESH_QS2_CACHE <- Sys.getenv("GSEA_REFRESH_QS2_CACHE", unset = "1") == "1"
QS2_CACHE_DIR <- file.path("temporary", DATA_TYPE, DATASET_ID, "GSEA_qs2_cache")
MSIGDB_REFERENCE_DIR <- file.path("data", "reference", "msigdb")
MSIGDB_REFERENCE_MAX_AGE_DAYS <- 7
QS2_NTHREADS <- max(1L, parallel::detectCores(logical = TRUE))
GSEA_TASK_WORKERS <- QS2_NTHREADS
GSEA_INNER_NPROC <- 1L
SINGLE_PATHWAY_PLOT_WORKERS <- 1L

qs2::qopt("nthreads", QS2_NTHREADS)
options(mc.cores = GSEA_TASK_WORKERS)

if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(
    BiocParallel::MulticoreParam(workers = GSEA_TASK_WORKERS),
    default = TRUE
  )
}

GENE_ID_COLUMN_BY_TYPE <- c(
  ENTREZ = "Entrez",
  SYMBOL = "Symbol",
  ENSEMBL = "Ensembl"
)

READABLE_KEY_TYPE_BY_GENE_ID_TYPE <- c(
  ENTREZ = "ENTREZID",
  SYMBOL = "SYMBOL",
  ENSEMBL = "ENSEMBL"
)

ORGDB_PACKAGE_BY_SPECIES <- c(
  human = "org.Hs.eg.db",
  `homo sapiens` = "org.Hs.eg.db",
  mouse = "org.Mm.eg.db",
  `mus musculus` = "org.Mm.eg.db",
  rat = "org.Rn.eg.db",
  `rattus norvegicus` = "org.Rn.eg.db"
)

MSIGDB_GENE_COLUMN_BY_TYPE <- c(
  ENTREZ = "ncbi_gene",
  SYMBOL = "gene_symbol",
  ENSEMBL = "ensembl_gene"
)

BREAK_RANK_TIES <- TRUE
RANK_TIE_EPSILON <- 1e-12
GSEA_RANDOM_SEED <- 20260604

GSEA_WARNING_PATTERNS_TO_HIDE <- c(
  "P-values are less than",
  "P-values were not calculated properly",
  "Invalid p-values detected",
  "NA values detected in gene set IDs",
  "Duplicate gene set IDs detected",
  "aes_string() was deprecated",
  "aes_() was deprecated",
  "Using `size` aesthetic for lines was deprecated"
)

GSEA_RESULT_COLUMNS <- c(
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment"
)


# 2. 常用函数 -----------------------------------------------------------------

get_msigdb_db_species <- function(species) {
  # msigdbr区分MSigDB数据库物种(db_species)和目标物种(species)。
  # 为减少头部配置，人类/大鼠默认使用人类MSigDB库，小鼠使用小鼠MSigDB库。
  species_lower <- tolower(trimws(species))

  if (species_lower %in% c("mouse", "mus musculus", "mm", "mmusculus")) {
    return("MM")
  }

  "HS"
}

get_orgdb_package <- function(species) {
  species_key <- tolower(trimws(species))
  orgdb_package <- ORGDB_PACKAGE_BY_SPECIES[[species_key]]

  if (is.null(orgdb_package)) {
    stop("No OrgDb package is configured for SPECIES = ", species)
  }

  orgdb_package
}

get_orgdb_object <- function(species) {
  orgdb_package <- get_orgdb_package(species)
  if (!requireNamespace(orgdb_package, quietly = TRUE)) {
    stop(
      "READABLE_GENE_SYMBOLS requires package ",
      orgdb_package,
      ". Please install it before running GSEA."
    )
  }

  get(orgdb_package, envir = asNamespace(orgdb_package))
}

get_readable_key_type <- function(gene_id_type) {
  key_type <- READABLE_KEY_TYPE_BY_GENE_ID_TYPE[[gene_id_type]]
  if (is.null(key_type)) {
    stop("No readable keyType is configured for GENE_ID_TYPE = ", gene_id_type)
  }

  key_type
}

make_msigdb_key <- function(collection, subcollection) {
  if (is.na(subcollection) || subcollection == "") {
    return(collection)
  }

  paste(collection, subcollection, sep = ":")
}

make_msigdb_output_name <- function(collection, subcollection) {
  if (is.na(subcollection) || subcollection == "") {
    if (collection == "H") {
      return("hallmark")
    }
    return(collection)
  }

  # 有子类别时优先用子类别命名输出目录，例如GO:BP保存为GO_BP。
  # 目录结构已经位于GSEA下，不再重复写入C5等大类前缀。
  gsub("[:/ -]+", "_", subcollection)
}

build_msigdb_geneset_catalog <- function() {
  catalog_table <- load_msigdb_catalog_table()

  catalog_list <- lapply(seq_len(nrow(catalog_table)), function(i) {
    collection <- as.character(catalog_table$gs_collection[i])
    subcollection <- as.character(catalog_table$gs_subcollection[i])
    if (is.na(subcollection)) {
      subcollection <- ""
    }

    list(
      key = make_msigdb_key(collection, subcollection),
      collection = collection,
      subcollection = if (subcollection == "") NULL else subcollection,
      output_name = make_msigdb_output_name(collection, subcollection),
      description = as.character(catalog_table$gs_collection_name[i])
    )
  })

  names(catalog_list) <- vapply(catalog_list, `[[`, character(1), "key")
  catalog_list
}

select_msigdb_genesets <- function(catalog, genesets_to_run) {
  genesets_to_run <- unique(trimws(as.character(genesets_to_run)))

  if (length(genesets_to_run) == 1 && tolower(genesets_to_run) == "all") {
    selected_catalog <- catalog
  } else {
    missing_genesets <- setdiff(genesets_to_run, names(catalog))
    if (length(missing_genesets) > 0) {
      stop(
        "Undefined MSigDB gene set keys in GSEA_GENESETS_TO_RUN: ",
        paste(missing_genesets, collapse = ", ")
      )
    }

    selected_catalog <- catalog[genesets_to_run]
  }

  names(selected_catalog) <- vapply(
    selected_catalog,
    `[[`,
    character(1),
    "output_name"
  )

  if (any(duplicated(names(selected_catalog)))) {
    names(selected_catalog) <- make.unique(names(selected_catalog), sep = "_")
  }

  selected_catalog
}

get_runtime_genesets_to_run <- function() {
  # 日常运行使用脚本头部GSEA_GENESETS_TO_RUN。
  # 仅测试脚本时，可在命令行临时设置环境变量：
  # GSEA_TEST_GENESETS="H,C6" Rscript scripts/GSE114012/06_gsea_analysis.R
  # 这样无需修改脚本即可只跑少数基因集。
  test_genesets <- Sys.getenv("GSEA_TEST_GENESETS", unset = "")
  test_genesets <- trimws(test_genesets)

  if (test_genesets == "") {
    return(GSEA_GENESETS_TO_RUN)
  }

  unique(trimws(strsplit(test_genesets, ",", fixed = TRUE)[[1]]))
}

clean_previous_gsea_outputs <- function(selected_analyses) {
  # 每次重跑都清空所选分析设计的GSEA表格和图片结果，避免旧结果残留。
  # ANALYSES_TO_RUN为"all"时，selected_analyses即全部差异分析设计。
  if (!CLEAN_GSEA_OUTPUT_DIR) {
    return(invisible(FALSE))
  }

  table_dirs <- file.path(TABLE_OUTPUT_ROOT, selected_analyses, "GSEA")
  plot_dirs <- file.path(PLOT_ROOT, sanitize_file_name(selected_analyses))

  unlink(table_dirs, recursive = TRUE, force = TRUE)
  unlink(plot_dirs, recursive = TRUE, force = TRUE)

  # 若当前运行覆盖全部分析，则直接清理GSEA图片根目录下的陈旧散落文件/目录。
  if (identical(ANALYSES_TO_RUN, "all") && dir.exists(PLOT_ROOT)) {
    unlink(PLOT_ROOT, recursive = TRUE, force = TRUE)
  }

  invisible(TRUE)
}

with_gsea_warnings_suppressed <- function(expr) {
  # 定向静默批量运行中反复出现、但不影响结果的第三方包提示。
  withCallingHandlers(
    expr,
    warning = function(w) {
      warning_text <- conditionMessage(w)
      if (any(vapply(
        GSEA_WARNING_PATTERNS_TO_HIDE,
        grepl,
        logical(1),
        x = warning_text,
        fixed = TRUE
      ))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

object_signature <- function(x) {
  paste(capture.output(dput(x)), collapse = "")
}

cache_hash <- function(...) {
  # 使用md5短哈希避免缓存文件名过长，同时让配置变化能自动生成新缓存。
  text <- paste(..., collapse = "|")
  tmp_file <- tempfile("gsea_cache_hash_")
  writeLines(text, tmp_file, useBytes = TRUE)
  hash_value <- unname(tools::md5sum(tmp_file))
  unlink(tmp_file)
  hash_value
}

get_file_stamp <- function(file) {
  file_info <- file.info(file)
  paste(file_info$size, as.integer(file_info$mtime), sep = "_")
}

read_qs2_cache <- function(cache_file) {
  if (!USE_QS2_CACHE || REFRESH_QS2_CACHE || !file.exists(cache_file)) {
    return(NULL)
  }

  qs2::qs_read(cache_file)
}

write_qs2_cache <- function(object, cache_file) {
  if (!USE_QS2_CACHE) {
    return(invisible(FALSE))
  }

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(object, cache_file)
  invisible(TRUE)
}

is_reference_cache_fresh <- function(cache_file) {
  if (!file.exists(cache_file)) {
    return(FALSE)
  }

  cache_age_days <- as.numeric(
    difftime(Sys.time(), file.info(cache_file)$mtime, units = "days")
  )
  cache_age_days <= MSIGDB_REFERENCE_MAX_AGE_DAYS
}

read_reference_qs2_cache <- function(cache_file) {
  if (!is_reference_cache_fresh(cache_file)) {
    return(NULL)
  }

  qs2::qs_read(cache_file)
}

write_reference_qs2_cache <- function(object, cache_file) {
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  qs2::qs_save(object, cache_file)
  invisible(TRUE)
}

get_msigdb_catalog_cache_file <- function() {
  file.path(
    MSIGDB_REFERENCE_DIR,
    paste0("catalog_", get_msigdb_db_species(SPECIES), ".qs2")
  )
}

load_msigdb_catalog_table <- function() {
  cache_file <- get_msigdb_catalog_cache_file()
  cached_catalog <- read_reference_qs2_cache(cache_file)
  if (!is.null(cached_catalog)) {
    return(cached_catalog)
  }

  catalog_table <- msigdbr::msigdbr_collections(
    db_species = get_msigdb_db_species(SPECIES)
  )
  write_reference_qs2_cache(catalog_table, cache_file)
  catalog_table
}

get_msigdb_data_cache_file <- function(config) {
  subcollection <- ifelse(
    is.null(config$subcollection),
    "all",
    sanitize_file_name(config$subcollection)
  )

  file.path(
    MSIGDB_REFERENCE_DIR,
    get_msigdb_db_species(SPECIES),
    sanitize_file_name(SPECIES),
    paste0(
      sanitize_file_name(config$collection),
      "_",
      subcollection,
      ".qs2"
    )
  )
}

get_gsea_cache_file <- function(analysis_name, deg_file, geneset_name, config) {
  hash_value <- cache_hash(
    "gsea",
    SPECIES,
    GENE_ID_TYPE,
    RANK_METRIC_COLUMN,
    READABLE_GENE_SYMBOLS,
    if (READABLE_GENE_SYMBOLS) get_orgdb_package(SPECIES) else "not_readable",
    get_file_stamp(deg_file),
    geneset_name,
    object_signature(config),
    object_signature(GSEA_PARAMS)
  )

  file.path(
    QS2_CACHE_DIR,
    "gsea_result",
    sanitize_file_name(analysis_name),
    paste0(sanitize_file_name(geneset_name), "_", hash_value, ".qs2")
  )
}

get_msigdb_collection <- function(config) {
  # 根据脚本头部配置获取一套msigdbr基因集。
  cache_file <- get_msigdb_data_cache_file(config)
  cached_msigdb_data <- read_reference_qs2_cache(cache_file)
  if (!is.null(cached_msigdb_data)) {
    return(list(
      data = cached_msigdb_data,
      cache_source = "reference_cache",
      cache_file = cache_file
    ))
  }

  args <- list(
    db_species = get_msigdb_db_species(SPECIES),
    species = SPECIES,
    collection = config$collection
  )

  if (!is.null(config$subcollection)) {
    args$subcollection <- config$subcollection
  }

  msigdb_data <- do.call(msigdbr::msigdbr, args)
  write_reference_qs2_cache(msigdb_data, cache_file)

  list(
    data = msigdb_data,
    cache_source = "computed",
    cache_file = cache_file
  )
}

prepare_msigdb_terms <- function(msigdb_data, gene_id_type) {
  # 将msigdbr结果整理为clusterProfiler::GSEA需要的TERM2GENE/TERM2NAME格式。
  gene_column <- MSIGDB_GENE_COLUMN_BY_TYPE[[gene_id_type]]
  stopifnot(!is.null(gene_column))
  stopifnot(gene_column %in% colnames(msigdb_data))
  stopifnot("gs_name" %in% colnames(msigdb_data))

  term2gene <- msigdb_data[, c("gs_name", gene_column), drop = FALSE]
  colnames(term2gene) <- c("term", "gene")

  term2gene$term <- as.character(term2gene$term)
  term2gene$gene <- as.character(term2gene$gene)
  term2gene <- term2gene[
    !is.na(term2gene$term) & term2gene$term != "" &
      !is.na(term2gene$gene) & term2gene$gene != "",
    ,
    drop = FALSE
  ]
  term2gene <- unique(term2gene)

  term2name <- unique(data.frame(
    term = as.character(msigdb_data$gs_name),
    name = as.character(msigdb_data$gs_name),
    stringsAsFactors = FALSE
  ))

  list(
    term2gene = term2gene,
    term2name = term2name
  )
}

load_msigdb_terms <- function(geneset_name, config) {
  msigdb_cache <- get_msigdb_collection(config)
  terms <- prepare_msigdb_terms(msigdb_cache$data, GENE_ID_TYPE)
  terms$config <- config
  terms$Cache_Source <- msigdb_cache$cache_source
  terms$Reference_Cache_File <- msigdb_cache$cache_file
  terms
}

count_nes_direction <- function(result_table, direction = c("positive", "negative")) {
  # 终端summary用；不参与官方GSEA结果表保存。
  direction <- match.arg(direction)
  if (!"NES" %in% colnames(result_table)) {
    return(0L)
  }

  nes <- as.numeric(result_table$NES)
  if (direction == "positive") {
    return(sum(nes > 0, na.rm = TRUE))
  }

  sum(nes < 0, na.rm = TRUE)
}

prepare_gsea_table_for_text_output <- function(dat) {
  # md/tex属于同一官方结果表的不同文本格式；NA按R输出习惯保留为"NA"。
  dat <- as.data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
  dat <- data.frame(
    lapply(dat, function(column) {
      column <- as.character(column)
      column[is.na(column)] <- "NA"
      column
    }),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  dat
}

write_gsea_markdown_table <- function(dat, md_file) {
  dat <- prepare_gsea_table_for_text_output(dat)

  header <- paste(escape_markdown_vector(colnames(dat)), collapse = " | ")
  separator <- paste(rep("---", ncol(dat)), collapse = " | ")
  rows <- if (nrow(dat) > 0) {
    apply(dat, 1, function(row) {
      paste(escape_markdown_vector(row), collapse = " | ")
    })
  } else {
    character(0)
  }

  writeLines(
    c(
      paste0("| ", header, " |"),
      paste0("| ", separator, " |"),
      paste0("| ", rows, " |")
    ),
    md_file,
    useBytes = TRUE
  )
}

write_gsea_latex_table <- function(dat, tex_file) {
  dat <- prepare_gsea_table_for_text_output(dat)

  col_spec <- paste(rep("l", ncol(dat)), collapse = "")
  header <- paste(escape_latex_vector(colnames(dat)), collapse = " & ")
  rows <- if (nrow(dat) > 0) {
    apply(dat, 1, function(row) {
      paste(escape_latex_vector(row), collapse = " & ")
    })
  } else {
    character(0)
  }

  writeLines(
    c(
      "\\begingroup",
      "\\tiny",
      "\\setlength{\\tabcolsep}{2pt}",
      "\\renewcommand{\\arraystretch}{1.18}",
      sprintf("\\begin{tabular}{@{}%s@{}}", col_spec),
      "\\toprule",
      paste0(header, " \\\\"),
      "\\midrule",
      if (length(rows) > 0) paste0(rows, " \\\\"),
      "\\bottomrule",
      "\\end{tabular}",
      "\\endgroup"
    ),
    tex_file,
    useBytes = TRUE
  )
}

prepare_gene_list <- function(deg_table, gene_id_type, rank_metric_column) {
  # 从all_genes.csv构建GSEA排序向量：names为基因ID，values为排序统计量。
  gene_column <- GENE_ID_COLUMN_BY_TYPE[[gene_id_type]]
  stopifnot(!is.null(gene_column))
  stopifnot(gene_column %in% colnames(deg_table))
  stopifnot(rank_metric_column %in% colnames(deg_table))

  dat <- data.frame(
    Gene_ID = trimws(as.character(deg_table[[gene_column]])),
    Rank_Metric = as.numeric(deg_table[[rank_metric_column]]),
    stringsAsFactors = FALSE
  )

  dat <- dat[
    !is.na(dat$Gene_ID) & dat$Gene_ID != "" &
      is.finite(dat$Rank_Metric),
    ,
    drop = FALSE
  ]

  stopifnot(nrow(dat) > 0)

  # 同一基因ID出现多次时，保留绝对统计量最大的记录。
  dat <- dat[order(-abs(dat$Rank_Metric), dat$Gene_ID), , drop = FALSE]
  dat <- dat[!duplicated(dat$Gene_ID), , drop = FALSE]

  if (BREAK_RANK_TIES && any(duplicated(dat$Rank_Metric))) {
    tie_breaker <- rank(dat$Gene_ID, ties.method = "first")
    dat$Rank_Metric <- dat$Rank_Metric + tie_breaker * RANK_TIE_EPSILON
  }

  gene_list <- dat$Rank_Metric
  names(gene_list) <- dat$Gene_ID
  gene_list <- sort(gene_list, decreasing = TRUE)

  stopifnot(length(gene_list) > 1)
  stopifnot(!any(duplicated(names(gene_list))))
  gene_list
}

run_clusterprofiler_gsea <- function(gene_list, term2gene, term2name, params) {
  # 只把官方GSEA支持的参数传入clusterProfiler::GSEA，保持结果字段为官方格式。
  # nproc经由...传递给fgsea/multilevel底层函数。
  # 当前脚本优先做任务级并行；若任务数少，再自动把更多核心分配给单个GSEA。
  set.seed(GSEA_RANDOM_SEED)
  params_for_run <- params
  if (!"nproc" %in% names(params_for_run)) {
    params_for_run$nproc <- GSEA_INNER_NPROC
  }

  gsea_call <- function() {
    do.call(
      clusterProfiler::GSEA,
      c(
        list(
          geneList = gene_list,
          TERM2GENE = term2gene,
          TERM2NAME = term2name
        ),
        params_for_run
      )
    )
  }

  with_gsea_warnings_suppressed(gsea_call())
}

make_gsea_result_readable <- function(gsea_result) {
  # 将core_enrichment里的基因ID转换为Symbol，方便直接阅读富集到的基因。
  if (!READABLE_GENE_SYMBOLS || GENE_ID_TYPE == "SYMBOL") {
    return(gsea_result)
  }
  if (!inherits(gsea_result, c("gseaResult", "enrichResult", "compareClusterResult"))) {
    return(gsea_result)
  }

  clusterProfiler::setReadable(
    x = gsea_result,
    OrgDb = get_orgdb_object(SPECIES),
    keyType = get_readable_key_type(GENE_ID_TYPE),
    toType = "SYMBOL"
  )
}

load_or_run_gsea <- function(
    analysis_name,
    deg_file,
    geneset_name,
    config,
    gene_list,
    term2gene,
    term2name) {
  cache_file <- get_gsea_cache_file(
    analysis_name = analysis_name,
    deg_file = deg_file,
    geneset_name = geneset_name,
    config = config
  )

  cached_result <- read_qs2_cache(cache_file)
  if (!is.null(cached_result)) {
    return(list(result = cached_result, source = "cache"))
  }

  gsea_result <- run_clusterprofiler_gsea(
    gene_list = gene_list,
    term2gene = term2gene,
    term2name = term2name,
    params = GSEA_PARAMS
  )
  gsea_result <- make_gsea_result_readable(gsea_result)
  write_qs2_cache(gsea_result, cache_file)

  list(result = gsea_result, source = "computed")
}

is_gsea_result_object <- function(gsea_result) {
  inherits(gsea_result, c("gseaResult", "enrichResult", "compareClusterResult"))
}

make_empty_gsea_result_table <- function() {
  empty_table <- as.data.frame(
    matrix(nrow = 0, ncol = length(GSEA_RESULT_COLUMNS)),
    stringsAsFactors = FALSE
  )
  colnames(empty_table) <- GSEA_RESULT_COLUMNS
  empty_table
}

get_gsea_result_table <- function(gsea_result) {
  # as.data.frame(gseaResult)保留clusterProfiler官方结果字段。
  if (is_gsea_result_object(gsea_result)) {
    return(as.data.frame(gsea_result))
  }

  if (is.data.frame(gsea_result)) {
    return(gsea_result)
  }

  make_empty_gsea_result_table()
}

write_gsea_result_tables <- function(gsea_result, csv_file) {
  result_table <- get_gsea_result_table(gsea_result)
  write.csv(result_table, csv_file, row.names = FALSE, na = "NA")
  write_gsea_markdown_table(result_table, sub("[.]csv$", ".md", csv_file))
  write_gsea_latex_table(result_table, sub("[.]csv$", ".tex", csv_file))
  result_table
}

format_gsea_description_for_plot <- function(description) {
  # 只优化绘图标签，不改变官方GSEA结果表。
  description <- as.character(description)

  if (SIMPLIFY_PATHWAY_PREFIX_IN_PLOT) {
    description <- gsub("^HALLMARKS?[_ -]+", "", description)
    description <- gsub("^GO[_ -]*BP[_ -]+", "", description)
    description <- gsub("^GOBP[_ -]+", "", description)
    description <- gsub("^GO[_ -]*CC[_ -]+", "", description)
    description <- gsub("^GOCC[_ -]+", "", description)
    description <- gsub("^GO[_ -]*MF[_ -]+", "", description)
    description <- gsub("^GOMF[_ -]+", "", description)
    description <- gsub("^KEGG[_ -]+", "", description)
    description <- gsub("^REACTOME[_ -]+", "", description)
    description <- gsub("^BIOCARTA[_ -]+", "", description)
    description <- gsub("^WP[_ -]+", "", description)
    description <- gsub("^C[0-9]+[_ -]+", "", description)
  }

  if (REPLACE_UNDERSCORE_WITH_SPACE_IN_PLOT) {
    description <- gsub("_", " ", description)
  }

  description <- gsub("\\s+", " ", description)
  trimws(description)
}

prepare_gsea_result_for_plot <- function(gsea_result) {
  # dotplot无法稳定处理ID/NES/P值为NA的通路；这里只过滤绘图对象，不改变CSV结果。
  plot_result <- gsea_result
  result_table <- as.data.frame(plot_result)

  required_columns <- c("ID", "Description", "NES", "pvalue", "p.adjust")
  missing_columns <- setdiff(required_columns, colnames(result_table))
  if (length(missing_columns) > 0) {
    return(plot_result)
  }

  valid_index <- !is.na(result_table$ID) &
    result_table$ID != "" &
    !is.na(result_table$Description) &
    result_table$Description != "" &
    is.finite(as.numeric(result_table$NES)) &
    is.finite(as.numeric(result_table$pvalue)) &
    is.finite(as.numeric(result_table$p.adjust))

  plot_result@result <- result_table[valid_index, , drop = FALSE]
  plot_result@result$Description <- format_gsea_description_for_plot(
    plot_result@result$Description
  )
  plot_result
}

get_longest_plot_label_line <- function(plot_labels) {
  if (length(plot_labels) == 0) {
    return(20L)
  }

  label_lines <- unlist(strsplit(as.character(plot_labels), "\n", fixed = TRUE))
  max(nchar(label_lines), na.rm = TRUE)
}

get_plot_label_line_counts <- function(plot_labels, shown_terms) {
  if (length(plot_labels) == 0) {
    return(rep(1L, max(as.integer(shown_terms), 1L)))
  }

  vapply(strsplit(as.character(plot_labels), "\n", fixed = TRUE), length, integer(1))
}

get_gsea_dotplot_size <- function(shown_terms, plot_labels) {
  shown_terms <- max(as.integer(shown_terms), 1L)
  label_line_counts <- get_plot_label_line_counts(plot_labels, shown_terms)
  total_label_lines <- sum(label_line_counts, na.rm = TRUE)
  label_height <- total_label_lines * DOTPLOT_LABEL_LINE_HEIGHT +
    shown_terms * DOTPLOT_TERM_GAP_HEIGHT

  body_size <- max(DOTPLOT_BODY_BASE_SIZE, label_height)
  body_size <- min(max(body_size, DOTPLOT_BODY_MIN_SIZE), DOTPLOT_BODY_MAX_SIZE)

  longest_label <- get_longest_plot_label_line(plot_labels)
  label_width <- DOTPLOT_LABEL_BASE_WIDTH +
    longest_label * DOTPLOT_LABEL_WIDTH_PER_CHARACTER
  label_width <- min(
    max(label_width, DOTPLOT_LABEL_MIN_WIDTH),
    DOTPLOT_LABEL_MAX_WIDTH
  )

  width <- label_width + body_size + DOTPLOT_LEGEND_WIDTH
  height <- body_size + DOTPLOT_VERTICAL_PADDING

  list(
    width = width,
    height = height
  )
}

make_empty_gsea_plot <- function() {
  # 没有可展示通路时仍输出空白占位图，保证批量结果目录完整。
  ggplot2::ggplot() +
    ggplot2::geom_text(
      ggplot2::aes(x = 0, y = 0, label = "No GSEA terms to display"),
      family = TEXT_FONT_FAMILY,
      fontface = TEXT_FONT_FACE,
      size = 4.2
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void(base_family = TEXT_FONT_FAMILY) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

apply_common_gsea_theme <- function(plot) {
  suppressMessages(plot + ggplot2::scale_size(range = GSEAVIS_POINT_SIZE_RANGE)) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = TEXT_FONT_FAMILY,
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      plot.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      axis.text.y = ggplot2::element_text(
        family = TEXT_FONT_FAMILY,
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      axis.title = ggplot2::element_text(
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      legend.title = ggplot2::element_text(
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      legend.text = ggplot2::element_text(
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      panel.border = ggplot2::element_rect(
        color = TEXT_COLOR,
        linewidth = AXIS_LINE_WIDTH,
        fill = NA
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey90",
        color = TEXT_COLOR,
        linewidth = AXIS_LINE_WIDTH
      ),
      strip.text = ggplot2::element_text(
        family = TEXT_FONT_FAMILY,
        face = TEXT_FONT_FACE,
        color = TEXT_COLOR
      ),
      axis.ticks = ggplot2::element_line(
        color = TEXT_COLOR,
        linewidth = AXIS_LINE_WIDTH
      ),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

has_dotplot_display_terms <- function(result_table) {
  # GseaVis::dotplotGsea会先按p值阈值筛选通路。
  # 若筛选后没有任何通路，ggplot分面在导出PDF/PNG时会报错；
  # 这里提前判断并改为输出空白占位图，不改变GSEA官方结果表。
  if (nrow(result_table) == 0) {
    return(FALSE)
  }

  display_table <- result_table

  if (!is.null(GSEAVIS_DOTPLOT_PARAMS$pval) && "pvalue" %in% colnames(display_table)) {
    display_table <- display_table[
      !is.na(display_table$pvalue) &
        display_table$pvalue <= GSEAVIS_DOTPLOT_PARAMS$pval,
      ,
      drop = FALSE
    ]
  }

  if (!is.null(GSEAVIS_DOTPLOT_PARAMS$pajust) && "p.adjust" %in% colnames(display_table)) {
    display_table <- display_table[
      !is.na(display_table$p.adjust) &
        display_table$p.adjust <= GSEAVIS_DOTPLOT_PARAMS$pajust,
      ,
      drop = FALSE
    ]
  }

  nrow(display_table) > 0
}

make_gsea_dotplot <- function(gsea_result, result_table, analysis_name, geneset_name) {
  if (!has_dotplot_display_terms(result_table)) {
    return(list(
      plot = make_empty_gsea_plot(),
      shown_terms = 0L,
      plot_labels = character(0)
    ))
  }

  dotplot_object <- with_gsea_warnings_suppressed(
    do.call(
      GseaVis::dotplotGsea,
      c(list(data = gsea_result), GSEAVIS_DOTPLOT_PARAMS)
    )
  )

  if (is.list(dotplot_object) && "plot" %in% names(dotplot_object)) {
    plot_object <- dotplot_object$plot
    shown_terms <- if ("df" %in% names(dotplot_object)) {
      nrow(dotplot_object$df)
    } else {
      min(nrow(result_table), GSEAVIS_DOTPLOT_PARAMS$topn * 2L)
    }
  } else {
    plot_object <- dotplot_object
    shown_terms <- min(nrow(result_table), GSEAVIS_DOTPLOT_PARAMS$topn * 2L)
  }

  if (is.list(dotplot_object) && "df" %in% names(dotplot_object) &&
      nrow(dotplot_object$df) == 0) {
    return(list(
      plot = make_empty_gsea_plot(),
      shown_terms = 0L,
      plot_labels = character(0)
    ))
  }

  plot_labels <- if (is.list(dotplot_object) && "df" %in% names(dotplot_object)) {
    levels(dotplot_object$df$Description)
  } else {
    character(0)
  }

  list(
    plot = apply_common_gsea_theme(plot_object + ggplot2::labs(title = NULL)),
    shown_terms = shown_terms,
    plot_labels = plot_labels
  )
}

gsea_scores_compat <- function(gene_list, gene_set, exponent = 1, fortify = TRUE) {
  # 兼容GseaVis旧接口需要的DOSE::gseaScores。
  # 该函数按经典加权GSEA公式计算running enrichment score，只用于绘图，不改变GSEA结果。
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_set <- intersect(as.character(gene_set), names(gene_list))

  hit_index <- names(gene_list) %in% gene_set
  hit_count <- sum(hit_index)
  gene_count <- length(gene_list)

  if (hit_count == 0 || hit_count == gene_count) {
    stop("Invalid gene set for GseaVis plotting.")
  }

  rank_weight <- abs(gene_list)^exponent
  hit_score <- cumsum(ifelse(
    hit_index,
    rank_weight / sum(rank_weight[hit_index]),
    0
  ))
  miss_score <- cumsum(ifelse(!hit_index, 1 / (gene_count - hit_count), 0))
  running_score <- hit_score - miss_score

  data.frame(
    x = seq_along(gene_list),
    runningScore = running_score,
    position = as.integer(hit_index)
  )
}

patch_gseavis_gsinfo_if_needed <- function() {
  # 当前GseaVis版本的gsInfo会调用旧版DOSE内部函数gseaScores；
  # 若本机DOSE已移除该函数，则在当前R会话内替换GseaVis::gsInfo以保证gseaNb可用。
  if (exists("gseaScores", asNamespace("DOSE"), inherits = FALSE)) {
    return(invisible(FALSE))
  }

  gsinfo_compat <- function(object, geneSetID) {
    geneList <- object@geneList
    if (is.numeric(geneSetID)) {
      geneSetID <- object@result[geneSetID, "ID"]
    }

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gsea_scores_compat(geneList, geneSet, exponent, fortify = TRUE)
    df$ymin <- 0
    df$ymax <- 0

    pos <- df$position == 1
    h <- diff(range(df$runningScore)) / 20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    df$Description <- object@result[geneSetID, "Description"]
    df
  }

  gseavis_namespace <- asNamespace("GseaVis")
  was_locked <- bindingIsLocked("gsInfo", gseavis_namespace)
  if (was_locked) {
    unlockBinding("gsInfo", gseavis_namespace)
  }
  assign("gsInfo", gsinfo_compat, envir = gseavis_namespace)
  if (was_locked) {
    lockBinding("gsInfo", gseavis_namespace)
  }
  invisible(TRUE)
}

get_analysis_design_row <- function(analysis_designs, analysis_name) {
  design_index <- match(analysis_name, analysis_designs$Analysis_Name)
  if (is.na(design_index)) {
    stop("Cannot find analysis design for ", analysis_name)
  }

  analysis_designs[design_index, , drop = FALSE]
}

get_single_pathway_sample_info <- function(
    analysis_designs,
    sample_info_all,
    analysis_name) {
  design_row <- get_analysis_design_row(analysis_designs, analysis_name)
  design_samples <- prepare_design_samples(
    sample_info = sample_info_all,
    group_column_index = design_row$Column_Index,
    experiment_group = design_row$Experiment_Group
  )

  sample_info <- design_samples$sample_info
  group_value <- as.character(design_samples$group_list)

  # 表达热图中优先展示实验组(LRC)，再展示对照组(BULK)，便于解读方向。
  sample_order <- order(
    group_value != design_row$Experiment_Group,
    group_value,
    sample_info$Sample_ID
  )
  sample_info <- sample_info[sample_order, , drop = FALSE]

  label_column <- SINGLE_PATHWAY_SAMPLE_LABEL_COLUMN
  if (!label_column %in% colnames(sample_info)) {
    label_column <- "Sample_ID"
  }

  sample_labels <- trimws(as.character(sample_info[[label_column]]))
  empty_label <- is.na(sample_labels) | sample_labels == ""
  sample_labels[empty_label] <- sample_info$Sample_ID[empty_label]
  sample_labels <- make.unique(sample_labels, sep = "_")

  list(
    sample_info = sample_info,
    sample_labels = sample_labels
  )
}

get_expression_cache_file <- function(analysis_name, sample_ids, sample_labels) {
  hash_value <- cache_hash(
    "single_pathway_expression",
    analysis_name,
    get_file_stamp(SE_RDS_FILE),
    TPM_ASSAY_NAME,
    SINGLE_PATHWAY_SAMPLE_LABEL_COLUMN,
    paste(sample_ids, collapse = ","),
    paste(sample_labels, collapse = ",")
  )

  file.path(
    QS2_CACHE_DIR,
    "single_pathway_expression",
    paste0(sanitize_file_name(analysis_name), "_", hash_value, ".qs2")
  )
}

prepare_single_pathway_expression_table <- function(
    analysis_name,
    expression_matrix_all,
    gene_annotation_all,
    analysis_designs,
    sample_info_all) {
  sample_data <- get_single_pathway_sample_info(
    analysis_designs = analysis_designs,
    sample_info_all = sample_info_all,
    analysis_name = analysis_name
  )

  sample_ids <- sample_data$sample_info$Sample_ID
  sample_labels <- sample_data$sample_labels
  cache_file <- get_expression_cache_file(analysis_name, sample_ids, sample_labels)
  cached_expression <- read_qs2_cache(cache_file)
  if (!is.null(cached_expression)) {
    return(cached_expression)
  }

  stopifnot(all(sample_ids %in% colnames(expression_matrix_all)))
  stopifnot("Symbol" %in% colnames(gene_annotation_all))

  expr_matrix <- expression_matrix_all[, sample_ids, drop = FALSE]
  colnames(expr_matrix) <- sample_labels

  gene_symbol <- trimws(as.character(gene_annotation_all$Symbol))
  keep_gene <- !is.na(gene_symbol) & gene_symbol != ""
  expr_matrix <- expr_matrix[keep_gene, , drop = FALSE]
  gene_symbol <- gene_symbol[keep_gene]

  # 同名Symbol只保留平均表达最高的一条记录，避免GseaVis热图中重复基因名。
  mean_expr <- rowMeans(expr_matrix, na.rm = TRUE)
  keep_order <- order(-mean_expr, gene_symbol)
  expr_matrix <- expr_matrix[keep_order, , drop = FALSE]
  gene_symbol <- gene_symbol[keep_order]
  unique_index <- !duplicated(gene_symbol)
  expr_matrix <- expr_matrix[unique_index, , drop = FALSE]
  gene_symbol <- gene_symbol[unique_index]

  expression_table <- data.frame(
    gene_name = gene_symbol,
    expr_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expression_data <- list(
    expression_table = expression_table,
    sample_order = sample_labels,
    sample_info = sample_data$sample_info
  )

  write_qs2_cache(expression_data, cache_file)
  expression_data
}

select_single_pathway_ids <- function(result_table) {
  if (!DRAW_SINGLE_PATHWAY_GSEA || nrow(result_table) == 0) {
    return(character(0))
  }

  stopifnot(SINGLE_PATHWAY_PVALUE_COLUMN %in% colnames(result_table))
  plot_table <- result_table[
    is.finite(as.numeric(result_table[[SINGLE_PATHWAY_PVALUE_COLUMN]])) &
      as.numeric(result_table[[SINGLE_PATHWAY_PVALUE_COLUMN]]) <=
        SINGLE_PATHWAY_PVALUE_CUTOFF,
    ,
    drop = FALSE
  ]

  if (nrow(plot_table) == 0) {
    return(character(0))
  }

  if (length(SINGLE_PATHWAY_KEYWORDS) > 0) {
    search_text <- paste(plot_table$ID, plot_table$Description, sep = " ")
    keyword_pattern <- paste(SINGLE_PATHWAY_KEYWORDS, collapse = "|")
    matched_index <- grepl(keyword_pattern, search_text, ignore.case = TRUE)
    plot_table <- plot_table[matched_index, , drop = FALSE]

    if (nrow(plot_table) == 0) {
      return(character(0))
    }

    plot_table <- plot_table[order(
      as.numeric(plot_table[[SINGLE_PATHWAY_PVALUE_COLUMN]]),
      -abs(as.numeric(plot_table$NES))
    ), , drop = FALSE]
    return(head(unique(plot_table$ID), SINGLE_PATHWAY_MAX_KEYWORD_TERMS))
  }

  plot_table <- plot_table[order(
    as.numeric(plot_table[[SINGLE_PATHWAY_PVALUE_COLUMN]]),
    -abs(as.numeric(plot_table$NES))
  ), , drop = FALSE]
  head(unique(plot_table$ID), SINGLE_PATHWAY_TOP_N)
}

get_core_gene_count <- function(result_table, term_id) {
  term_index <- match(term_id, result_table$ID)
  if (is.na(term_index) || !"core_enrichment" %in% colnames(result_table)) {
    return(20L)
  }

  core_genes <- unique(unlist(strsplit(result_table$core_enrichment[term_index], "/")))
  length(core_genes[core_genes != "" & !is.na(core_genes)])
}

get_single_pathway_plot_size <- function(result_table, term_id, sample_count) {
  core_gene_count <- get_core_gene_count(result_table, term_id)
  term_index <- match(term_id, result_table$ID)
  title_text <- if (is.na(term_index)) term_id else result_table$Description[term_index]
  title_text <- format_gsea_description_for_plot(title_text)
  title_lines <- length(strsplit(
    paste(strwrap(title_text, width = GSEAVIS_SINGLE_PATHWAY_PARAMS$termWidth), collapse = "\n"),
    "\n",
    fixed = TRUE
  )[[1]])

  plot_size <- SINGLE_PATHWAY_PLOT_BASE_SIZE +
    core_gene_count * SINGLE_PATHWAY_PLOT_SIZE_PER_GENE +
    sample_count * SINGLE_PATHWAY_PLOT_SIZE_PER_SAMPLE +
    title_lines * SINGLE_PATHWAY_PLOT_SIZE_PER_TITLE_LINE

  # 单通路图要求整体保持正方形。这里额外根据热图样本名和核心基因名数量
  # 估算最低画布大小，避免样本名或基因名在GseaVis热图区域发生重叠。
  heatmap_ratio <- GSEAVIS_SINGLE_PATHWAY_PARAMS$ght.relHight
  sample_required_size <- sample_count *
    SINGLE_PATHWAY_SAMPLE_TEXT_MIN_SIZE / 72 /
    max(heatmap_ratio * SINGLE_PATHWAY_SAMPLE_TEXT_HEIGHT_FACTOR, 0.05)
  gene_required_size <- core_gene_count *
    SINGLE_PATHWAY_GENE_TEXT_MIN_SIZE / 72 /
    max(SINGLE_PATHWAY_GENE_TEXT_WIDTH_FACTOR, 0.05)

  plot_size <- max(plot_size, sample_required_size, gene_required_size)
  plot_size <- min(
    max(plot_size, SINGLE_PATHWAY_PLOT_MIN_SIZE),
    SINGLE_PATHWAY_PLOT_MAX_SIZE
  )

  # 主图区域按正方形导出，避免上方GSEA曲线和下方表达热图比例失衡。
  list(width = plot_size, height = plot_size)
}

get_single_pathway_gene_text_size <- function(plot_width, core_gene_count) {
  # GseaVis表达热图中基因名为90度竖排；
  # 这里按图宽和核心基因数动态缩小字号，尽量保证相邻基因名不重叠。
  core_gene_count <- max(as.integer(core_gene_count), 1L)
  available_points <- plot_width * 72 * SINGLE_PATHWAY_GENE_TEXT_WIDTH_FACTOR /
    core_gene_count

  min(
    max(available_points, SINGLE_PATHWAY_GENE_TEXT_MIN_SIZE),
    SINGLE_PATHWAY_GENE_TEXT_MAX_SIZE
  )
}

get_single_pathway_sample_text_size <- function(plot_height, sample_count) {
  # GseaVis没有单独暴露表达热图样本名字号参数；
  # 这里通过patchwork全局主题设置axis.text.y，并按热图高度动态计算字号。
  sample_count <- max(as.integer(sample_count), 1L)
  heatmap_height <- plot_height * GSEAVIS_SINGLE_PATHWAY_PARAMS$ght.relHight
  available_points <- heatmap_height * 72 *
    SINGLE_PATHWAY_SAMPLE_TEXT_HEIGHT_FACTOR / sample_count

  min(
    max(available_points, SINGLE_PATHWAY_SAMPLE_TEXT_MIN_SIZE),
    SINGLE_PATHWAY_SAMPLE_TEXT_MAX_SIZE
  )
}

apply_common_single_gsea_theme <- function(plot, sample_text_size) {
  common_theme <- ggplot2::theme(
    text = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    plot.title = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR,
      hjust = 0.5
    ),
    axis.text = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    axis.text.y = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR,
      size = sample_text_size
    ),
    axis.title = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    legend.title = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    legend.text = ggplot2::element_text(
      family = TEXT_FONT_FAMILY,
      face = TEXT_FONT_FACE,
      color = TEXT_COLOR
    ),
    panel.border = ggplot2::element_rect(
      color = TEXT_COLOR,
      linewidth = AXIS_LINE_WIDTH,
      fill = NA
    )
  )

  # gseaNb返回patchwork对象时，&可将主题应用到所有子图；若失败则回退到普通ggplot叠加。
  tryCatch(
    plot & common_theme,
    error = function(e) plot + common_theme
  )
}

make_single_pathway_gsea_plot <- function(
    gsea_result,
    term_id,
    expression_data,
    result_table,
    plot_size) {
  params <- GSEAVIS_SINGLE_PATHWAY_PARAMS
  params$object <- gsea_result
  params$geneSetID <- term_id
  params$ght.geneText.size <- get_single_pathway_gene_text_size(
    plot_width = plot_size$width,
    core_gene_count = get_core_gene_count(result_table, term_id)
  )
  sample_text_size <- get_single_pathway_sample_text_size(
    plot_height = plot_size$height,
    sample_count = length(expression_data$sample_order)
  )

  if (isTRUE(params$add.geneExpHt)) {
    params$exp <- expression_data$expression_table
    params$sample.order <- expression_data$sample_order
  }

  # setReadable后core_enrichment是Symbol，但geneList仍保留原始Entrez。
  # GseaVis在kegg=TRUE时会使用gseaResult@gene2Symbol匹配表达热图，适合当前readable结果。
  params$kegg <- READABLE_GENE_SYMBOLS && GENE_ID_TYPE != "SYMBOL"

  plot <- with_gsea_warnings_suppressed(
    do.call(GseaVis::gseaNb, params)
  )

  apply_common_single_gsea_theme(plot, sample_text_size = sample_text_size)
}

save_single_pathway_gsea_plots <- function(
    gsea_result,
    result_table,
    analysis_name,
    geneset_name,
    plot_output_dir,
    expression_data) {
  selected_term_ids <- select_single_pathway_ids(result_table)
  if (length(selected_term_ids) == 0) {
    return(0L)
  }

  patch_gseavis_gsinfo_if_needed()
  single_pathway_dir <- file.path(plot_output_dir, "single_pathway")
  dir.create(single_pathway_dir, recursive = TRUE, showWarnings = FALSE)

  save_one_single_pathway_plot <- function(i) {
    patch_gseavis_gsinfo_if_needed()
    term_id <- selected_term_ids[i]
    term_dir <- file.path(
      single_pathway_dir,
      paste0(sprintf("%02d_", i), sanitize_file_name(term_id))
    )
    dir.create(term_dir, recursive = TRUE, showWarnings = FALSE)

    plot_size <- get_single_pathway_plot_size(
      result_table = result_table,
      term_id = term_id,
      sample_count = length(expression_data$sample_order)
    )
    single_plot <- make_single_pathway_gsea_plot(
      gsea_result = gsea_result,
      term_id = term_id,
      expression_data = expression_data,
      result_table = result_table,
      plot_size = plot_size
    )

    save_ggplot_pdf_png(
      plot = single_plot,
      pdf_file = file.path(term_dir, "gsea_plot.pdf"),
      width = plot_size$width,
      height = plot_size$height
    )

    TRUE
  }

  plot_index <- seq_along(selected_term_ids)
  if (.Platform$OS.type != "windows" &&
      SINGLE_PATHWAY_PLOT_WORKERS > 1 &&
      length(plot_index) > 1) {
    saved_flags <- parallel::mclapply(
      plot_index,
      save_one_single_pathway_plot,
      mc.cores = min(SINGLE_PATHWAY_PLOT_WORKERS, length(plot_index)),
      mc.preschedule = FALSE
    )
  } else {
    saved_flags <- lapply(plot_index, save_one_single_pathway_plot)
  }

  sum(vapply(saved_flags, isTRUE, logical(1)))
}


# 3. 准备输入文件和基因集 ------------------------------------------------------

stopifnot(GENE_ID_TYPE %in% names(GENE_ID_COLUMN_BY_TYPE))
stopifnot(GENE_ID_TYPE %in% names(MSIGDB_GENE_COLUMN_BY_TYPE))

MSIGDB_GENESET_CATALOG <- build_msigdb_geneset_catalog()
RUNTIME_GENESETS_TO_RUN <- get_runtime_genesets_to_run()
GSEA_GENESET_CONFIG <- select_msigdb_genesets(
  catalog = MSIGDB_GENESET_CATALOG,
  genesets_to_run = RUNTIME_GENESETS_TO_RUN
)

file_info <- get_deg_file_info(TABLE_ROOT)
selected_analyses <- get_selected_analysis_names(
  file_info = file_info,
  analyses_to_plot = ANALYSES_TO_RUN
)

clean_previous_gsea_outputs(selected_analyses)

if (DRAW_SINGLE_PATHWAY_GSEA) {
  se <- readRDS(SE_RDS_FILE)
  clinical_data <- read.csv(
    CLINICAL_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot("Sample_ID" %in% colnames(clinical_data))
  stopifnot(!any(duplicated(clinical_data$Sample_ID)))
  stopifnot(TPM_ASSAY_NAME %in% names(SummarizedExperiment::assays(se)))

  analysis_designs <- get_analysis_designs(clinical_data)
  expression_matrix_all <- as.matrix(SummarizedExperiment::assay(se, TPM_ASSAY_NAME))
  expression_matrix_all <- log2(expression_matrix_all + 1)
  gene_annotation_all <- as.data.frame(
    SummarizedExperiment::rowData(se),
    stringsAsFactors = FALSE
  )

  sample_info_all <- clinical_data[
    match(colnames(expression_matrix_all), clinical_data$Sample_ID),
    ,
    drop = FALSE
  ]
  rownames(sample_info_all) <- sample_info_all$Sample_ID
  stopifnot(all(sample_info_all$Sample_ID == colnames(expression_matrix_all)))
}

cat("\nGSEA runtime configuration:\n")
cat("Selected analyses: ", length(selected_analyses), "\n", sep = "")
cat(
  "Selected MSigDB gene set categories: ",
  length(GSEA_GENESET_CONFIG),
  if (length(RUNTIME_GENESETS_TO_RUN) == 1 &&
      tolower(RUNTIME_GENESETS_TO_RUN) == "all") {
    " (all available categories from msigdbr)"
  } else {
    ""
  },
  "\n",
  sep = ""
)
cat("qs2 threads: ", QS2_NTHREADS, "\n", sep = "")
cat("Output root:  ", OUTPUT_ROOT, "\n", sep = "")
cat("Refresh qs2 cache: ", REFRESH_QS2_CACHE, "\n", sep = "")
cat("Previous GSEA outputs were cleaned before this run.\n")

cat("\nLoading MSigDB gene sets...\n")
geneset_cache <- lapply(names(GSEA_GENESET_CONFIG), function(geneset_name) {
  config <- GSEA_GENESET_CONFIG[[geneset_name]]
  terms <- load_msigdb_terms(geneset_name, config)

  list(
    config = config,
    term2gene = terms$term2gene,
    term2name = terms$term2name,
    cache_source = terms$Cache_Source
  )
})
names(geneset_cache) <- names(GSEA_GENESET_CONFIG)

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


# 4. 批量运行GSEA并绘图 --------------------------------------------------------

cat("\nPreparing analysis-level inputs...\n")
analysis_cache <- lapply(selected_analyses, function(analysis_name) {
  deg_index <- match(analysis_name, file_info$Analysis_Name)
  deg_file <- file_info$All_Genes_File[deg_index]
  deg_result <- read_deg_result(file_info, analysis_name)
  gene_list <- prepare_gene_list(
    deg_table = deg_result,
    gene_id_type = GENE_ID_TYPE,
    rank_metric_column = RANK_METRIC_COLUMN
  )
  expression_data <- if (DRAW_SINGLE_PATHWAY_GSEA) {
    prepare_single_pathway_expression_table(
      analysis_name = analysis_name,
      expression_matrix_all = expression_matrix_all,
      gene_annotation_all = gene_annotation_all,
      analysis_designs = analysis_designs,
      sample_info_all = sample_info_all
    )
  } else {
    NULL
  }

  list(
    deg_file = deg_file,
    gene_list = gene_list,
    expression_data = expression_data
  )
})
names(analysis_cache) <- selected_analyses

task_table <- expand.grid(
  Analysis_Name = selected_analyses,
  GeneSet_Name = names(geneset_cache),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_tasks <- nrow(task_table)

# 任务级并行策略：
# 1. 任务很多时，外层同时运行多个analysis x geneset任务，使CPU长期保持高负荷；
# 2. 任务很少时，单个GSEA任务自动获得更多nproc核心；
# 3. 外层并行开启时，单通路图在任务内部顺序保存，避免嵌套并行导致CPU过度抢占。
GSEA_TASK_WORKERS <- if (.Platform$OS.type != "windows") {
  min(QS2_NTHREADS, total_tasks)
} else {
  1L
}
GSEA_INNER_NPROC <- if (GSEA_TASK_WORKERS > 0) {
  max(1L, floor(QS2_NTHREADS / GSEA_TASK_WORKERS))
} else {
  QS2_NTHREADS
}
QS2_NTHREADS_PER_TASK <- GSEA_INNER_NPROC
SINGLE_PATHWAY_PLOT_WORKERS <- if (GSEA_TASK_WORKERS > 1) {
  1L
} else {
  QS2_NTHREADS
}

qs2::qopt("nthreads", QS2_NTHREADS_PER_TASK)
options(mc.cores = GSEA_TASK_WORKERS)
if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(
    BiocParallel::MulticoreParam(workers = GSEA_INNER_NPROC),
    default = TRUE
  )
}

cat("\nParallel execution strategy:\n")
cat("Total tasks:              ", total_tasks, "\n", sep = "")
cat("Task-level workers:       ", GSEA_TASK_WORKERS, "\n", sep = "")
cat("GSEA nproc per task:      ", GSEA_INNER_NPROC, "\n", sep = "")
cat("qs2 nthreads per task:    ", QS2_NTHREADS_PER_TASK, "\n", sep = "")
cat("Single-pathway workers:   ", SINGLE_PATHWAY_PLOT_WORKERS, "\n", sep = "")

run_gsea_task <- function(task_id) {
  analysis_name <- task_table$Analysis_Name[task_id]
  geneset_name <- task_table$GeneSet_Name[task_id]
  analysis_input <- analysis_cache[[analysis_name]]
  cache <- geneset_cache[[geneset_name]]
  output_name <- cache$config$output_name
  output_dir_name <- sanitize_file_name(output_name)

  table_output_dir <- file.path(
    TABLE_OUTPUT_ROOT,
    analysis_name,
    "GSEA",
    output_dir_name
  )
  plot_output_dir <- file.path(
    PLOT_ROOT,
    sanitize_file_name(analysis_name),
    output_dir_name
  )
  dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)

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
    plot_result_table <- make_empty_gsea_result_table()
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
  pdf_file <- file.path(plot_output_dir, "dotplot.pdf")
  plot_files <- with_gsea_warnings_suppressed(
    save_ggplot_pdf_png(
      plot = dotplot_result$plot,
      pdf_file = pdf_file,
      width = plot_size$width,
      height = plot_size$height
    )
  )

  single_pathway_count <- if (DRAW_SINGLE_PATHWAY_GSEA &&
                              is_gsea_result_object(gsea_result)) {
    with_gsea_warnings_suppressed(
      save_single_pathway_gsea_plots(
        gsea_result = gsea_result,
        result_table = result_table,
        analysis_name = analysis_name,
        geneset_name = geneset_name,
        plot_output_dir = plot_output_dir,
        expression_data = analysis_input$expression_data
      )
    )
  } else {
    0L
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

run_gsea_tasks_with_progress <- function(task_ids, task_function, workers) {
  # 主进程维护进度条，子进程负责实际GSEA与绘图任务。
  # 与parallel::mclapply相比，这里可以实时回收已完成任务并刷新终端进度。
  total_task_count <- length(task_ids)
  progress_bar <- utils::txtProgressBar(
    min = 0,
    max = total_task_count,
    style = 3
  )
  on.exit(close(progress_bar), add = TRUE)

  results <- vector("list", total_task_count)
  names(results) <- as.character(task_ids)

  if (workers <= 1L || .Platform$OS.type == "windows") {
    for (task_position in seq_along(task_ids)) {
      task_id <- task_ids[task_position]
      results[[as.character(task_id)]] <- try(task_function(task_id), silent = TRUE)
      utils::setTxtProgressBar(progress_bar, task_position)
    }

    return(results)
  }

  workers <- max(1L, min(as.integer(workers), total_task_count))
  next_task_position <- 1L
  completed_tasks <- 0L
  active_jobs <- list()
  active_task_by_pid <- list()

  launch_next_task <- function() {
    task_id <- task_ids[next_task_position]
    job <- parallel::mcparallel(
      try(task_function(task_id), silent = TRUE),
      silent = TRUE
    )
    pid <- as.character(job$pid)

    active_jobs[[pid]] <<- job
    active_task_by_pid[[pid]] <<- task_id
    next_task_position <<- next_task_position + 1L
  }

  while (
    next_task_position <= total_task_count &&
      length(active_jobs) < workers
  ) {
    launch_next_task()
  }

  while (completed_tasks < total_task_count) {
    ready_results <- parallel::mccollect(
      active_jobs,
      wait = FALSE,
      timeout = 0.5
    )

    if (is.null(ready_results) || length(ready_results) == 0) {
      Sys.sleep(0.2)
      next
    }

    for (pid in names(ready_results)) {
      task_id <- active_task_by_pid[[pid]]
      results[[as.character(task_id)]] <- ready_results[[pid]]
      active_jobs[[pid]] <- NULL
      active_task_by_pid[[pid]] <- NULL
      completed_tasks <- completed_tasks + 1L
      utils::setTxtProgressBar(progress_bar, completed_tasks)

      while (
        next_task_position <= total_task_count &&
          length(active_jobs) < workers
      ) {
        launch_next_task()
      }
    }
  }

  results
}

cat("\nRunning batch GSEA analyses...\n")
task_ids <- seq_len(total_tasks)
summary_records <- run_gsea_tasks_with_progress(
  task_ids = task_ids,
  task_function = run_gsea_task,
  workers = GSEA_TASK_WORKERS
)

task_errors <- vapply(summary_records, inherits, logical(1), "try-error")
if (any(task_errors)) {
  stop(
    "Some GSEA tasks failed: ",
    paste(which(task_errors), collapse = ", ")
  )
}


# 5. 终端快速汇总 --------------------------------------------------------------

summary_table <- do.call(rbind, summary_records)
rownames(summary_table) <- NULL

summary_output_dir <- file.path(OUTPUT_ROOT, "tables", "GSEA_summary")
dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv_with_report_previews(
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

cat("\nOutput summary:\n")
cat("GSEA table directory: ", file.path(TABLE_OUTPUT_ROOT, "<analysis_name>", "GSEA"), "\n", sep = "")
cat("GSEA plot directory:  ", file.path(PLOT_ROOT, "<analysis_name>"), "\n", sep = "")
cat("GSEA summary table:   ", file.path(summary_output_dir, "summary.csv"), "\n", sep = "")
cat("CSV/MD/TEX result sets: ", nrow(summary_table), " each\n", sep = "")
cat("PDF/PNG dotplots: ", nrow(summary_table), " each\n", sep = "")
cat("Single-pathway GSEA PDF/PNG plots: ", sum(summary_table$Single_Pathway_Plots), " each\n", sep = "")

SCRIPT_END_TIME <- Sys.time()
SCRIPT_RUNTIME_SECONDS <- as.numeric(
  difftime(SCRIPT_END_TIME, SCRIPT_START_TIME, units = "secs")
)
cat(
  "Total runtime: ",
  sprintf("%.2f min (%.1f sec)", SCRIPT_RUNTIME_SECONDS / 60, SCRIPT_RUNTIME_SECONDS),
  "\n",
  sep = ""
)

cat("\nBatch GSEA analysis finished.\n")
