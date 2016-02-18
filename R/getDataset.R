#' Extract real dataset parameters.
#' Extract Negative Binomial parameters from real dataset.
#' @param counts Numeric matrix of read counts from which to extract parameters.
#' @param drop.extreme.dispersion Proportion of extreme dispersions to drop.
#' @param drop.low.lambda Logical, whether to drop low lambdas.
#' @keywords internal
#' @export
#' @examples
#'
#' @return List of parameters sampled from the real dataset: vector of gene-wise average log-CPM, vector of gene-wise dispersions, library sizes, number of genes.
#'

getDataset <- function(counts, drop.extreme.dispersion=0.1, drop.low.lambda=TRUE) {

  # Get nn0 = number of nonzero counts per row
  # Filter out rows with only 1 nonzero count
  nsamples = dim(counts)[2]
  counts0 = counts==0
  nn0 = rowSums(!counts0)
  if(any(nn0 == 1)){
    counts = counts[nn0 > 1, ]
    nn0 = nn0[nn0 > 1]
    counts0 = counts==0
  }

  # Turn count matrix input into DGElist object
  # Calculate normalization factors
  # Get the counts per million of normalized library sizes
  d <- DGEList(counts)
  d <- calcNormFactors(d)
  cp <- round(cpm(d, normalized.lib.sizes=TRUE),1)

  # If we decided to drop low lambdas: take rowsums using only cpm > 1, and keep only the rows whose sums are greater than 2
  if(drop.low.lambda) d <- d[rowSums(cp>1) >= 2, ]

  # Save the tagwise average logCPM
  d$AveLogCPM <- log2(rowMeans(cpm(d, prior.count = 1e-5)))

  mu = rowSums(counts)/nsamples
  var = rowSums((counts - mu)^2)/(nsamples-1)
  disp = (var-mu + 0.0001)/mu^2
  disp = ifelse(disp > 0, disp, min(disp[disp > 0]))

  d$tagwise.dispersion <- disp

  # If we specified a percentage of extreme dispersions to drop, then drop those genes with extreme dispersions past the cutoff
  if(is.numeric(drop.extreme.dispersion))
  {
    # bad = cutoff dispersion value
    bad <- quantile(disp, 1-drop.extreme.dispersion, names = FALSE)
    # Index the tags whose dispersions are below the cutoff dispersion, i.e. the ones to keep
    ids <- disp <= bad
    # Subset avelogcpm and dispersions for the tags with ok dispersions
    d$AveLogCPM <- d$AveLogCPM[ids]
    disp <- disp[ids]
  }

  # Return list of parameters sampled from the real dataset
  list(dataset.AveLogCPM = as.numeric(d$AveLogCPM), dataset.dispersion = disp, dataset.lib.size = d$samples$lib.size, dataset.nTags = nrow(d))
}
