# GSE282081分析流水线一键运行脚本
#
# 默认顺序运行：
# 00 样本设计 -> 01 ATF3表达概览 -> 02 差异分析 ->
# 03 ATF3相关性功能GSEA -> 04 差异分析GSEA


# 0. 可修改配置 ---------------------------------------------------------------

SCRIPTS_TO_RUN <- c(
  "scripts/GSE282081/00_prepare_sample_design.R",
  "scripts/GSE282081/01_atf3_expression_overview.R",
  "scripts/GSE282081/02_limma_differential_expression.R",
  "scripts/GSE282081/03_atf3_correlation_gsea.R",
  "scripts/GSE282081/04_gsea_analysis.R"
)

options(width = 200)


# 1. 运行 ----------------------------------------------------------------------

cat("\nRunning GSE282081 analysis pipeline...\n")
cat("Scripts: ", length(SCRIPTS_TO_RUN), "\n\n", sep = "")

rscript_bin <- file.path(R.home("bin"), "Rscript")

for (script_file in SCRIPTS_TO_RUN) {
  if (!file.exists(script_file)) {
    stop("Missing script: ", script_file)
  }

  cat("\n============================================================\n")
  cat("Running: ", script_file, "\n", sep = "")
  cat("============================================================\n")

  exit_status <- system2(rscript_bin, args = script_file)

  if (!identical(exit_status, 0L)) {
    stop("Script failed: ", script_file, " (exit status ", exit_status, ")")
  }
}

cat("\nGSE282081 analysis pipeline finished.\n")
