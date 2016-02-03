#' Apply edgeR on dataset.
#' This is a general function which applies edgeR on a given dataset.
#' @param dataset Numeric count matrix of read counts.
#' @param filename Filename of the dataset; edgeR results on that dataset will be named accordingly
#' @param group Vector or factor giving the experimental group/condition for each sample/library in the pilot dataset. If not provided, groups are assumed equally divided.
#' @export
#' @examples
#'
#' @return Applies edgeR to given dataset, and returns results table

apply.edgeR <- function(dataset, group=NULL){

  # Note I'm in the ExpDesignApp/Code directory

  G <- dim(dataset)[1]
  N <- dim(dataset)[2]
  libsizes <- colSums(dataset)

  if (is.null(group)){
    group = as.factor(rep(c("A", "B"), each=N/2))
  }
  group = as.factor(group)

  d <- DGEList(counts=dataset,
               group=group)
  design <- model.matrix(~group)

  d <- calcNormFactors(d)
  d <- estimateGLMCommonDisp(d, design)
  d <- estimateGLMTrendedDisp(d, design)
  d <- estimateGLMTagwiseDisp(d, design)

  # Fit genewise GLMs and test treatment effect
  fit <- glmFit(d, design)
  lrt <- glmLRT(fit)
  res.table <- lrt$table
  res.table
}
