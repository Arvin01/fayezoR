#' Simulate scRNA-seq dataset.
#' Uses previously re-sampled parameters to simulate data according to Negative Binomial distribution.
#' @param object DGEList object.
#' @param seed Optional seed, if it desired to replicate the simulation.
#' @keywords internal
#' @export
#' @examples
#'
#' @return Saves simulated DE data to DGEList object.
#'

sim.fun <- function(object, seed)
{
  Lambda <- object$Lambda
  Dispersion <- object$Dispersion
  pZero <- object$pZero
  nTags <- object$nTags
  nlibs <- object$nlibs
  lib.size <- object$lib.size

  set.seed(seed)
  counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = nTags, ncol = nlibs)
  I.count <- t(sapply(pZero, function(p) rbinom(nlibs, prob=1-p, size=1)))
  counts <- I.count*counts

  rownames(counts) <- paste("ids", 1:nTags, sep = "")
  object$counts <- counts
  return(object)
}
