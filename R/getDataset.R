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
  
  # Estimate dispersions via edgeR
  d <- estimateGLMCommonDisp(d)
  d <- estimateGLMTrendedDisp(d)
  d <- estimateGLMTagwiseDisp(d)
  
  dispersion <- d$tagwise.dispersion
  AveLogCPM <- d$AveLogCPM
  
  # If we specified a percentage of extreme dispersions to drop, then drop those genes with extreme dispersions past the cutoff
  if(is.numeric(drop.extreme.dispersion))
  {  
    # bad = cutoff dispersion value
    bad <- quantile(dispersion, 1-drop.extreme.dispersion, names = FALSE)
    # Index the tags whose dispersions are below the cutoff dispersion, i.e. the ones to keep
    ids <- dispersion <= bad
    # Subset avelogcpm and dispersions for the tags with ok dispersions
    AveLogCPM <- AveLogCPM[ids]
    dispersion <- dispersion[ids]
  }
  
  dataset.AveLogCPM <- AveLogCPM
  dataset.dispersion <- dispersion
  dataset.lib.size <- d$samples$lib.size
  dataset.nTags <- nrow(d)
  
  # Return list of parameters sampled from the real dataset
  list(dataset.AveLogCPM = AveLogCPM, dataset.dispersion = dispersion, dataset.lib.size = d$samples$lib.size, dataset.nTags = nrow(d))
}
