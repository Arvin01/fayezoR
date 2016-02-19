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
  mu <- object$mu
  dispersion <- object$dispersion
  pZero <- object$pZero
  nTags <- object$nTags
  nlibs <- object$nlibs
  lib.size <- object$lib.size
  design <- object$design
  beta <- object$beta

  # If either design or coefficients aren't provided
  if(is.null(design) | is.null(beta)){
    cat("Need both design and coefficients to impose differential expression...\n Simulating non-differentially expressed data.")
    lambda <- mu
    lambda <- expandAsMatrix(lambda, dim=c(nTags, nlibs))
  }

  else{
    cat("Simulating differentially expressed data.\n")
    # recall that mu = exp(beta0)
    lambda <- exp(log(mu) + beta %*% t(design))
  }

  dispersion <- expandAsMatrix(dispersion, dim = c(nTags, nlibs))

  set.seed(seed)
  counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(lambda)*lib.size), size = 1/dispersion), nrow = nTags, ncol = nlibs)
  I.count <- t(sapply(pZero, function(p) rbinom(nlibs, prob=1-p, size=1)))
  counts <- I.count*counts

  rownames(counts) <- paste("ids", 1:nTags, sep = "")
  colnames(counts) = paste("sample", 1:nlibs, sep = "")

  object$counts <- counts
  return(object)
}
