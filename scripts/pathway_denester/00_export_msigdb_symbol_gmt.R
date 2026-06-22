# Export MSigDB SYMBOL GMT files for PathwayDenester.
#
# PathwayDenester needs the same term IDs used in the enrichment result and a
# GMT file with full gene membership. Our GSEA tables are readable SYMBOL
# outputs, so these GMT files are exported with SYMBOL IDs.


# 0. Config -------------------------------------------------------------------

DATA_TYPE <- "ngs"
SPECIES <- "human"
GENE_ID_TYPE <- "SYMBOL"

PLOTTING_FUNCTION_FILE <- "scripts/functions/plotting_common_functions.R"
GSEA_FUNCTION_FILE <- "scripts/functions/gsea_common_functions.R"

MSIGDB_REFERENCE_DIR <- file.path("data", "reference", "msigdb")
MSIGDB_REFERENCE_MAX_AGE_DAYS <- 30L
OUTPUT_ROOT <- file.path("data", "reference", "pathway_denester", "msigdb_symbol_gmt")

GSEA_GENESETS_TO_EXPORT <- c(
  "H",
  "C2:CP:BIOCARTA",
  "C2:CP:KEGG_MEDICUS",
  "C2:CP:KEGG_LEGACY",
  "C2:CP:REACTOME",
  "C2:CP:WIKIPATHWAYS",
  "C3:TFT:TFT_LEGACY",
  "C3:TFT:GTRD",
  "C5:GO:BP",
  "C5:GO:CC",
  "C5:GO:MF",
  "C5:HPO",
  "C6",
  "C7:IMMUNESIGDB"
)

options(width = 200)


# 1. Packages and helpers ------------------------------------------------------

suppressPackageStartupMessages({
  library(qs2)
})

source(PLOTTING_FUNCTION_FILE)
source(GSEA_FUNCTION_FILE)

write_gmt_file <- function(term2gene, term2name, gmt_file) {
  term_ids <- unique(as.character(term2gene$term))
  term_ids <- term_ids[!is.na(term_ids) & term_ids != ""]
  term_ids <- sort(term_ids)

  name_map <- setNames(as.character(term2name$name), as.character(term2name$term))

  dir.create(dirname(gmt_file), recursive = TRUE, showWarnings = FALSE)
  con <- file(gmt_file, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  for (term_id in term_ids) {
    genes <- unique(as.character(term2gene$gene[term2gene$term == term_id]))
    genes <- genes[!is.na(genes) & genes != ""]
    if (length(genes) == 0) {
      next
    }

    term_name <- name_map[[term_id]]
    if (is.null(term_name) || is.na(term_name) || term_name == "") {
      term_name <- term_id
    }

    writeLines(
      paste(c(term_id, term_name, genes), collapse = "\t"),
      con = con,
      useBytes = TRUE
    )
  }

  invisible(TRUE)
}


# 2. Export --------------------------------------------------------------------

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

MSIGDB_GENESET_CATALOG <- build_msigdb_geneset_catalog()
GSEA_GENESET_CONFIG <- select_msigdb_genesets(
  catalog = MSIGDB_GENESET_CATALOG,
  genesets_to_run = GSEA_GENESETS_TO_EXPORT
)

summary_table <- do.call(
  rbind,
  lapply(names(GSEA_GENESET_CONFIG), function(geneset_name) {
    config <- GSEA_GENESET_CONFIG[[geneset_name]]
    terms <- load_msigdb_terms(geneset_name, config)
    gmt_file <- file.path(OUTPUT_ROOT, paste0(sanitize_file_name(geneset_name), ".gmt"))

    write_gmt_file(
      term2gene = terms$term2gene,
      term2name = terms$term2name,
      gmt_file = gmt_file
    )

    data.frame(
      GeneSet_Name = geneset_name,
      GMT_File = gmt_file,
      Terms = length(unique(terms$term2gene$term)),
      Term_Gene_Links = nrow(terms$term2gene),
      Source = terms$Cache_Source,
      stringsAsFactors = FALSE
    )
  })
)

summary_file <- file.path(OUTPUT_ROOT, "gmt_export_summary.csv")
write.csv(summary_table, summary_file, row.names = FALSE, na = "NA")

cat("\nPathwayDenester GMT export finished.\n")
cat("GMT root: ", OUTPUT_ROOT, "\n", sep = "")
cat("Summary:  ", summary_file, "\n\n", sep = "")
print(summary_table, row.names = FALSE)
