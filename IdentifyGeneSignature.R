library(survival)
library(caret)
library(ggplot2)
library(ggpubr)
source("utils.R")
#' Identify Gene Signature Based on Concordance Index
#'
#' This function identifies a gene signature from a ranked list of genes based on the concordance index (C-index).
#' It requires the 'caret' and 'survival' packages.
#'
#' @param rankedGeneList A character vector of gene names, ranked by node heat value
#' @param expression A data frame containing gene expression data with columns for genes and rows for samples. It must include an 'ID' column.
#' @param survival A data frame containing survival data with 'ID', 'TIME', and 'STATUS' columns.
#' @return A list containing the identified gene signature and iteration details.
#' @examples
#' rankedGeneList <- c("Gene1", "Gene2", "Gene3")
#' expression <- data.frame(ID = 1:10, Gene1 = rnorm(10), Gene2 = rnorm(10), Gene3 = rnorm(10))
#' survival <- data.frame(ID = 1:10, TIME = rnorm(10, 100, 50), STATUS = sample(0:1, 10, replace = TRUE))
#' result <- identifyGeneSignature(rankedGeneList, expression, survival)
#' @export
identifyGeneSignature <- function(rankedGeneList,expression,survival){
    check_columns(expression, c("ID"))
    check_columns(survival, c("ID", "TIME","STATUS"))
    check_required_packages(c("caret", "survival" ))
    pdata <- merge(expression,survival,by="ID")
    exprs <- pdata[,colnames(expression)]
    survival <- pdata[,colnames(survival)]

    exprs <- column_to_rownames(exprs,var="ID")
    survival <- column_to_rownames(survival,var="ID")

    candidategenes <- rankedGeneList
    outputs <- list()
    convergent <- F
    outTab <- data.frame()
    i<- 1
    while (!convergent) {
      before_count <- length(candidategenes)
      selectedgenes <- c()
      x <- c();y <- c();label <- c()
      ci_prev <- 0
      for (n in 1:length(candidategenes)) {
        genes<- candidategenes[1:n]
        genes <- gsub(" ","",genes)
        genes <- intersect(genes,colnames(exprs))
        clinical <- survival[,1:2]
        colnames(clinical) <- c("status","time")
        data <- cbind(clinical,  exprs[,genes]  )

        result <- coxphcv(data)
        ci <-  result[["ci"]]
        if (ci>ci_prev) {
          selectedgenes <- c(selectedgenes,candidategenes[n])
          label <- c(label,"increase")
        }else{
          label <- c(label,"decrease")

        }
        ci_prev <- ci
        x <- c(x,n);y <- c(y,ci)

        # print(paste0(n,": ",ci))
      }
      outTab <- rbind(outTab,data.frame(iter=rep( paste0("iteration ",i),length(x)),x=x,ci=y,label=label))
      i <- i+1
      candidategenes <- selectedgenes
      after_count <- length(selectedgenes)
      if (after_count==before_count) {
        convergent=TRUE
        outputs[["genesignature"]] <- selectedgenes
        outputs[["interation"]] <- outTab
      }
    }
    return(outputs)
}

#' Perform Cox Proportional Hazards Model Cross-Validation
#'
#' This function performs k-fold cross-validation for the Cox proportional hazards model and calculates the concordance index (C-index)
#' and hazard ratio (HR) for each fold. It requires the 'caret' and 'Hmisc' packages.
#'
#' @param data A data frame containing survival data with 'time' and 'status' columns, along with other predictor variables.
#' @param seed An integer specifying the random seed for reproducibility (default: 123).
#' @param kfold An integer specifying the number of folds for cross-validation (default: 10).
#' @return A list containing the mean C-index and hazard ratio (HR) across all folds.
#' @examples
#' data <- data.frame(time = rnorm(100, 100, 50), status = sample(0:1, 100, replace = TRUE), Gene1 = rnorm(100), Gene2 = rnorm(100))
#' result <- coxphcv(data)
#' @export
coxphcv <- function(data,seed=123,kfold= 10 ) {
  library(caret)
  library(Hmisc)
  set.seed(seed)
  folds <- createFolds(data$time, k = kfold, list = TRUE)

  results <- sapply(folds, function(indices) {
    train_data <- data[-indices, ]
    test_data <- data[indices, ]
    model <- coxph(Surv(time, status) ~., data = train_data)
    predictions <- predict(model, newdata = test_data, type = "risk")
    s <-  Surv(test_data$time, test_data$status)
    ci <- Hmisc::rcorr.cens(-1 * predictions, s)[["C Index"]]
    mod.ph <- coxph(s  ~ predictions)

    return(list(ci=ci,hr=exp(mod.ph$coef)))
  })
  results <- as.data.frame(t(results))
  ci <- mean(unlist(results[["ci"]]),na.rm = T)
  hr <- mean(unlist(results[["hr"]]),na.rm = T)
  return(list(ci=ci,hr=hr))
}


#' Visualize Iteration Data
#'
#' This function visualizes the iteration data by creating a scatter plot of Concordance index (ci) vs gene list size (x),
#' with different shapes and colors for the labels. It also ensures that the required columns are present in the data frame.
#'
#' @param iteration A data frame containing iteration data. It must include the columns: "iter", "x", "ci", and "label".
#' @return A ggplot object representing the visualization.
#' @examples
#' iteration <- data.frame(iter = rep(1:3, each = 10),
#'                         x = rep(1:10, times = 3),
#'                         ci = runif(30),
#'                         label = rep(c("increase", "decrease"), each = 15))
#' p <- visualizeIteration(iteration)
#' print(p)
#' @export
visualizeIteration <- function(iteration){
  required_cols <- c("iter", "x", "ci", "label")
  missing_cols <- setdiff(required_cols, colnames(iteration))
  if (length(missing_cols) > 0) {
    stop("The following columns are missing from the data frame: ", paste(missing_cols, collapse = ", "))
  }

  p <- ggplot(iteration,aes(x,ci))+
    geom_point(shape=ifelse(iteration$label=="increase",24,16),color=ifelse(iteration$label=="increase","chocolate3","cyan4"))+
    geom_line()+
    facet_wrap(iter~.,scales = "free")+
    theme_pubr()+
    labs(x="gene list size",y="Concordance index")+
    theme(panel.border = element_rect(fill = NA,color="gray30") )


  return(p)

}







