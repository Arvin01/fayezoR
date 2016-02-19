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
  aveLogCPM <- object$dataset.params$aveLogCPM
  dispersion <- object$dataset.params$dispersion
  pZero <- object$dataset.params$pZero

  # Indices from sampling a subset of those real dataset parameters
  set.seed(seed)
  id_r <- sample(length(aveLogCPM), nTags, replace = TRUE)

  # mu, dispersion = genewise parameters to be used for simulation dataset
  # mu is the baseline rate per gene
  mu <- 2^(aveLogCPM[id_r]) #lambda for each gene
  mu <- mu/sum(mu)
  dispersion <- dispersion[id_r]
  pZero <- pZero[id_r]

  # Get ids for which gene-wise rates are 0
  id_0 <- mu == 0

  # Keep only the parameters for which lambdas are not 0
  mu <- mu[!id_0]
  dispersion <- dispersion[!id_0]
  pZero <- pZero[!id_0]

  object$mu <- mu
  object$dispersion <- dispersion
  object$pZero <- pZero
  object
}
