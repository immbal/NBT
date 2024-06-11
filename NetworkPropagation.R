library(igraph)
library(Matrix)
library(expm)
source("utils.R")
#' Perform Network Propagation to Identify Top Genes
#'
#' This function performs network propagation to identify top genes based on the given anchor genes and protein-protein interactions (PPIs).
#' It requires the 'igraph', 'Matrix', and 'expm' packages.
#'
#' @param anchorGenes A character vector of anchor gene names.
#' @param PPIs A data frame containing protein-protein interaction data with columns for interacting proteins.
#' @param TOPN An integer specifying the number of top genes to return (default: 150).
#' @return A character vector of top gene names based on the network propagation.
#' @examples
#' anchorGenes <- c("Gene1", "Gene2", "Gene3")
#' PPIs <- data.frame(symbol_x = c("Gene1", "Gene2", "Gene4"), symbol_y = c("Gene2", "Gene3", "Gene5"))
#' top_genes <- networkPropagation(anchorGenes, PPIs)
#' @export
networkPropagation <- function(anchorGenes,PPIs,TOPN=150){

    check_required_packages(c("igraph", "Matrix","expm"))
    G <- graph_from_data_frame(PPIs, directed = FALSE)
    genes <- V(G)$name
    heat <- numeric(length(genes))
    for (gene in anchorGenes) {
      if (gene %in% genes) {
        heat[which(genes == gene)] <- 1
      }
    }
    L <- as.matrix(laplacian_matrix(G))
    t <- 0.1
    heat_diffused <- expm(-t * L) %*% heat
    df_outTab <- data.frame(gene = genes, heat = heat_diffused)
    df_outTab <- df_outTab[order(df_outTab$heat, decreasing = TRUE), ]
    top_genes <- head(df_outTab$gene, TOPN)
    return(top_genes)
}
