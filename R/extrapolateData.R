#' Extrapolate pilot data.
#' Takes an input pilot dataset and extrapolates it to a desired size.
#' @param pilotdata Numeric count matrix; pilot dataset to be extrapolated to a desired larger size.
#' @param design Model matrix (without an intercept) that you would like to simulate from
#' @param groups.pilot Vector or factor giving the experimental group/condition for each sample/library in the pilot dataset. If not provided, groups are assumed equally divided.
#' @param nreps Desired number of replicates per group in extrapolated data.
#' @param depth Desired depth (in millions) for each sample/replicate in extrapolated data.
#' @export
#' @examples
#'
#' @return Extrapolated dataset, given pilot data

extrapolateData <- function(pilotdata, design=NULL, nreps, depth){

  counts <- pilotdata

  # Get nn0 = number of nonzero counts per row
  # Filter out rows with only 1 nonzero count (otherwise can't get disp)
  counts0 = counts==0
  nn0 = rowSums(!counts0)
  if(any(nn0 == 1)){
    counts = counts[nn0 > 1, ]
    nn0 = nn0[nn0 > 1]
    counts0 = counts==0
  }

  nsamples <- dim(counts)[2]
  G <- dim(counts)[1]

  ### Estimate dispersion from nonzero data

  mu = rowSums((!counts0)*counts)/nsamples
  var = rowSums((!counts0)*(counts - mu)^2)/(nn0-1)
  disp = (var-mu + 0.0001)/mu^2
  disp = ifelse(disp > 0, disp, min(disp[disp > 0]))

  ### Get proportion of nonzeros for each row

  pZero = 1 - nn0/nsamples

  ### Estimate row means from nonzero data, separately for each group

  # If design matrix indicating groups is not provided, set 2 groups with equal number of replicates
  if (is.null(design)){
    group <- rep(c(-1, 1), each=nsamples/2)
    design <- model.matrix(~-1 + group)
  }

  # Subset the data by groups
  group <- as.factor(group)
  d.grp1 <- counts[,group==levels(group)[1]]
  d.grp2 <- counts[,group==levels(group)[2]]
  counts0.grp1 = d.grp1==0
  counts0.grp2 = d.grp2==0

  # Save the average log CPM, for each group
  # Note that aveLogCPM(x), rowMeans(cpm(x,log=TRUE)) and log2(rowMeans(cpm(x)) all give slightly different results.
  # Use aveLogCPM(x) here as it doesn't return any Inf values and its values are close to log2(rowMeans(cpm(x))) which I had used previously
  grp1.aveLogCPM <- aveLogCPM((!counts0.grp1)*d.grp1)
  grp2.aveLogCPM <- aveLogCPM((!counts0.grp2)*d.grp2)

  # Get genewise mean (for each group) to be used for extrapolated dataset
  mu1 <- 2^(grp1.aveLogCPM)
  mu2 <- 2^(grp2.aveLogCPM)
  mu1 <- mu1/sum(mu1)
  mu2 <- mu2/sum(mu2)

  # Expand genewise lambdas and dispersions
  mu1 <- expandAsMatrix(mu1, dim=c(G, nreps))
  mu2 <- expandAsMatrix(mu2, dim=c(G, nreps))
  mu <- cbind(mu1, mu2)
  dispersion <- expandAsMatrix(disp, dim = c(G, nreps*2))

  lib.size <- rep(depth*1e6, nreps*2)

  ### Extrapolate the data using NB sampling

  extrapolatedData <- matrix(rnbinom(G*nreps*2, mu = t(t(mu)*lib.size), size = 1/dispersion),
                             nrow = G, ncol = nreps*2)
  rownames(extrapolatedData) <- paste("ids", 1:G, sep = "")
  colnames(extrapolatedData) <- c(paste("A", 1:nreps, sep=""),
                                  paste("B", 1:nreps, sep=""))

  extrapolatedData
}

