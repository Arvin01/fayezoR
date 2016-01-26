#' Re-sample real dataset parameters. 
#' Re-sample parameters extracted from real dataset, to be used for the simulation dataset.
#' @param object DGEList object. 
#' @param seed Optional seed, if it desired to replicate the simulation. 
#' @keywords internal
#' @export
#' @examples
#' 
#' @return Saves re-sampled gene-wise lambdas and dispersions to DGEList object. 
#' 

sample.fun <- function(object, seed)
{
  nlibs <- object$nlibs
  nTags <- object$nTags
  
  # Parameters extracted from real dataset
  AveLogCPM <- object$dataset$dataset.AveLogCPM
  dispersion <- object$dataset$dataset.dispersion
  
  # Indices from sampling a subset of those real dataset parameters
  set.seed(seed)
  id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
  
  # Lambda, Dispersion = genewise parameters to be used for simulation dataset
  Lambda <- 2^(AveLogCPM[id_r]) #lambda for each gene
  Lambda <- Lambda/sum(Lambda)
  Dispersion <- dispersion[id_r]
  
  # Get ids for which lambdas are 0
  id_0 <- Lambda == 0
  
  # Keep only the parameters for which lambdas are not 0
  Lambda <- Lambda[!id_0]
  Dispersion <- Dispersion[!id_0]
  Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
  Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
  
  object$Lambda <- Lambda
  object$Dispersion <- Dispersion
  object  	    
}