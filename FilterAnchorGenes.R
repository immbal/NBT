library(survival)

source("utils.R")
#' Filter Anchor Genes Based on Concordance Index
#'
#' This function filters anchor genes based on the concordance index (C-index) calculated from the expression data
#' and phenotype data. It requires the 'Hmisc' and 'progress' packages.
#'
#' @param expression A data frame containing gene expression data with columns for genes and rows for samples. It must include an 'ID' column.
#' @param phenotype A data frame containing phenotype data with an 'ID' column and a 'phenotype' column.
#' @param cutoff A numeric value specifying the cutoff for the C-index (default: 0.6).
#' @return A character vector of anchor gene names that have a C-index greater than the specified cutoff.
#' @examples
#' expression <- data.frame(ID = 1:10, Gene1 = rnorm(10), Gene2 = rnorm(10))
#' phenotype <- data.frame(ID = 1:10, phenotype = sample(0:1, 10, replace = TRUE))
#' filterAnchorGenes(expression, phenotype)
#' @export
filterAnchorGenes <- function(expression,phenotype,cutoff=0.6){

  check_columns(expression, c("ID"))
  check_columns(phenotype, c("ID", "phenotype"))

  check_required_packages(c("Hmisc", "progress","survival"))

  genes <- setdiff(colnames(expression), "ID")
  pdata <- merge(expression,phenotype,by="ID")
  gene_ci <- list()

  pb <- progress::progress_bar$new(
    format = "  Processing genes [:bar] :percent in :elapsed",
    total = length(genes), clear = FALSE, width = 60
  )

  for (gene in genes) {
    s <- Surv(pdata[["phenotype"]],rep(1,nrow(pdata )))
    ci <- Hmisc::rcorr.cens(pdata[[gene]],s)[["C Index"]]
    gene_ci[[gene]] <- ci
    pb$tick()
  }
  anchor_genes <- gene_ci[gene_ci>cutoff]
  return( names(anchor_genes))
}

