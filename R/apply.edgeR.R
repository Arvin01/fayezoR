#' Apply edgeR on dataset. 
#' This is a general function which applies edgeR on a given dataset. 
#' @param dataset Numeric count matrix of read counts. 
#' @param filename Filename of the dataset; edgeR results on that dataset will be named accordingly
#' @param group Vector or factor giving the experimental group/condition for each sample/library in the pilot dataset. If not provided, groups are assumed equally divided. 
#' @export
#' @examples
#' 
#' @return Doesn't return anything, just saves the reults table as RData file.

apply.edgeR <- function(dataset, filename, group=NULL){
  
  # Note: I'm in the extrapolated data directory
  
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
  
  # Things to save: adjusted p-values, indices of DE genes
  res.table$FDR <- p.adjust(res.table$PValue, method="BH")
  
  # Save edgeR results with same filename as pilot data; just make sure it goes in appropriate directory
  filename <- paste("../edgeR_Results", filename, sep="/")
  save(res.table, file=filename)
}
