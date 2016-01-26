#' Impose differential expression. 
#' Alters the re-sampled parameters to reflect differential expression
#' @param object DGEList object.
#' @param logFC Desired log-fold changes for differentially expressed genes. If single number, all DE genes will exhibit constant logFC. If vector of two numbers provided (default \code{c(0.25, 1.5)}), logFC of DE genes will be drawn uniformly within that range. 
#' @param seed Optional seed, if it desired to replicate the simulation. 
#' @keywords internal
#' @export
#' @examples
#' 
#' @return Saves differential expression attributes for simulation to DGEList object: DE-adjusted lambdas, indices of DE and non-DE genes, indices of up- vs. down-regulation
#' 

diff.fun <- function(object, logFC, seed)
{   	 
  
  group <- object$group
  pDiff <- object$pDiff
  pUp <-  object$pUp 
  Lambda <- object$Lambda
  nTags <- object$nTags
  
  # TRUE/FALSE to identify group 1
  g <- group == levels(group)[1]
  
  # Indices for differentially expressed genes
  set.seed(seed)
  ind <- sample(nTags, floor(pDiff*nTags))
  
  # If there are at least some DE genes and their logFC's aren't trivial
  if(length(ind)>0 & logFC!=0) {
    
    # Get logFC for each DE gene
    set.seed(seed)
    if(length(logFC)==1){
      tagwise.logFC <- rep(logFC, length(ind))
    }
    if(length(logFC)==2){
      tagwise.logFC <- runif(length(ind), logFC[1], logFC[2])
    }
    
    # Direction of fold change for those differentially expressed genes
    # -1 if downregulated (group1<group2), 1 if upregulated (group1>group2)
    set.seed(seed)
    fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
    
    # Changing lambdas for the 2 groups to reflect differential expression
    Lambda[ind,g] <- Lambda[ind,g]*exp(tagwise.logFC/2*fcDir)
    Lambda[ind,!g] <- Lambda[ind,!g]*exp(tagwise.logFC/2*(-fcDir)) 
    
    # Saving new attributes to the object
    object$Lambda <- Lambda
    object$indDE <- ind #indices of DE genes
    object$indNonDE <- (1:nTags)[-ind] #indices of non-DE genes
    object$indUp <- ind[fcDir==1] #indices of DE genes for which group1>group2
    object$indDown <- ind[fcDir==-1] #indices of DE genes for which group1<group2
    object$logFC <- tagwise.logFC #log-fold changes of DE genes
    object
  }
  
  else{
    object$indDE <- NA
    object
  }
}
