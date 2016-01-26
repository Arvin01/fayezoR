#' Impose differential expression. 
#' Alters the re-sampled parameters to reflect differential expression
#' @param object DGEList object.
#' @param seed Optional seed, if it desired to replicate the simulation. 
#' @keywords internal
#' @export
#' @examples
#' 
#' @return Saves differential expression attributes for simulation to DGEList object: DE-adjusted lambdas, indices of DE and non-DE genes, indices of up- vs. down-regulation
#' 

diff.fun <- function(object, seed)
{   	 
  
  group <- object$group
  pDiff <- object$pDiff
  pUp <-  object$pUp 
  foldDiff <- object$foldDiff
  Lambda <- object$Lambda
  nTags <- object$nTags
  
  # TRUE/FALSE to identify group 1
  g <- group == levels(group)[1]
  
  # Indices for differentially expressed genes
  set.seed(seed)
  ind <- sample(nTags, floor(pDiff*nTags))
  
  if(length(ind)>0 & !foldDiff == 1 ) {
    
    # Direction of fold change for those differentially expressed genes
    # -1 if downregulated (group1<group2), 1 if upregulated (group1>group)
    set.seed(seed)
    fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
    
    # Changing lambdas for the 2 groups to reflect differential expression
    Lambda[ind,g] <- Lambda[ind,g]*exp(log(foldDiff)/2*fcDir)
    Lambda[ind,!g] <- Lambda[ind,!g]*exp(log(foldDiff)/2*(-fcDir)) 
    
    # Saving new attributes to the object
    object$Lambda <- Lambda
    object$indDE <- ind #indices of DE genes
    object$indNonDE <- (1:nTags)[-ind] #indices of non-DE genes
    object$mask_DEup <- object$mask_DEdown <- object$mask_nonDE <- expandAsMatrix(FALSE, dim = dim(Lambda)) 
    object$mask_DEup[ind[fcDir == 1], g] <- TRUE #group1=TRUE when group1 bigger
    object$mask_DEup[ind[fcDir == -1], !g] <- TRUE #group2=TRUE when group2 bigger
    object$mask_DEdown[ind[fcDir == -1], g] <- TRUE #group1=TRUE when group2 bigger
    object$mask_DEdown[ind[fcDir == 1], !g] <- TRUE #group2=TRUE when group1 bigger
    object$mask_nonDE[-ind,] <- TRUE #TRUE for nonDE
    object$mask_DE <- object$mask_DEup | object$mask_DEdown
  }
  
  if(foldDiff == 1 | pDiff == 0)
    object$indDE <- NA
  object
}
