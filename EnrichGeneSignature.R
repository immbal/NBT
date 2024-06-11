library(clusterProfiler)
library(cogena)
library(org.Hs.eg.db)
library(ggpubr)

#' Enrich Gene Signature
#'
#' This function performs gene set enrichment analysis using a specified database. It supports Reactome and GO databases.
#'
#' @param geneset A character vector of gene symbols to be analyzed.
#' @param database A character string specifying the database to use for enrichment analysis. Default is "Reactome". Supported values are "Reactome" and "GO".
#' @return A data frame containing the enrichment analysis results.
#' @examples
#' geneset <- c("TP53", "BRCA1", "EGFR")
#' result <- enrichGeneSignature(geneset, database = "Reactome")
#' print(result)
#' @export
enrichGeneSignature <- function(geneset,database="Reactome"){
  result <- data.frame()
  if (database=="Reactome") {
    gmt <- read.gmt("d:/project/sparseNN/data/ReactomePathways.gmt")
    result <- enricher(geneset,TERM2GENE = gmt)
  }
  if (database=="GO") {
    result <-  enrichGO(gene = geneset, keyType='SYMBOL',OrgDb = org.Hs.eg.db,ont = "BP",pvalueCutoff = 0.05)
  }
  return(result)
}

#' Visualize Enrichment Results
#'
#' This function creates a dot plot visualization of the enrichment analysis results.
#'
#' @param enrich An enrichment analysis result object, typically returned by functions like `enrichGO` or `enricher`.
#' @return A ggplot object representing the dot plot of the enrichment analysis results.
#' @examples
#' library(clusterProfiler)
#' geneset <- c("TP53", "BRCA1", "EGFR")
#' result <- enrichGO(gene = geneset, keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
#' p <- visualizeEnrichment(result)
#' print(p)
#' @export
visualizeEnrichment <- function(enrich){
  p <- dotplot(enrich ,font.size=7 )+theme_cleveland() 
  return(p)
}

 
