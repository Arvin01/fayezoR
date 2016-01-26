#' Impose dropouts. 
#' Impose dropouts 
#' @param object DGEList object.
#' @param pDropout Numeric vector of dropout probabilities to consider.
#' @param drop.method Dropout method, with two options: zero-inflation, or low-rate Poisson.
#' @param drop.lambda Mean of low-rate Poisson for the Poisson-dropout scenario.
#' @param seed Optional seed, if it desired to replicate the simulation. 
#' @keywords internal
#' @export
#' @examples
#' 
#' @return Saves attributes of dropout imposition on simulated data to DGEList object: dropout-imposed count matrix, gene-wise probability of dropout, logical matrix indicating dropout observations.
#' 

dropout.fun <- function(object, pDropout, drop.method, drop.lambda, seed){
  
  dim <- dim(object$counts)
  countAddDrop <- object$counts
  nTags <- object$nTags
  nlibs <- object$nlibs
  
  # Get dropout rate for each tag  
  set.seed(seed)
  tagwise.pDropout <- pDropout[sample(length(pDropout), nTags, replace=T)]
  
  # Given p, function that gets indices of dropouts for one row and changes to TRUE
  idxDropout <- function(p){
    vec <- rep(FALSE,nlibs)
    idx <- which(runif(nlibs) < p) #indices of dropouts
    vec[idx] <- TRUE
    vec
  }
  
  # Matrix of TRUE where dropout and FALSE otherwise
  mask_dropout <- do.call(rbind, lapply(tagwise.pDropout, idxDropout))
  
  # Alter dropout counts to be low-level Poisson
  set.seed(seed)
  
  countAddDrop[mask_dropout] <- switch(drop.method,
                                       "zero" = 0,
                                       "poisson" = rpois(sum(mask_dropout), drop.lambda)
  )
  
  # Saving some things
  object$counts.dropouts <- countAddDrop
  object$tagwise.pDropout <- tagwise.pDropout
  object$mask_dropout <- mask_dropout
  object
}
