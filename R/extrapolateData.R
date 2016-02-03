#' Extrapolate pilot data.
#' Takes an input pilot dataset and extrapolates it to a desired size.
#' @param pilotdata Numeric count matrix; pilot dataset to be extrapolated to a desired larger size.
#' @param filename.pilot Filename of pilot dataset; extrapolated dataset will be named accordingly
#' @param groups.pilot Vector or factor giving the experimental group/condition for each sample/library in the pilot dataset. If not provided, groups are assumed equally divided.
#' @param nreps Desired number of replicates per group in extrapolated data.
#' @param depth Desired depth (in millions) for each sample/replicate in extrapolated data.
#' @export
#' @examples
#'
#' @return Doesn't return anything, just saves the down-sampled datasets as RData.

extrapolateData <- function(pilotdata, filename.pilot, group.pilot=NULL, nreps, depth){

  N.pilot <- dim(pilotdata)[2]
  G <- dim(pilotdata)[1]

  ### First, multinomial sampling to increase the depths
  ### Why? Higher depths better for estimating means and dispersions

  # Initialize the matrix
  d.augdepths <- data.frame(matrix(ncol=N.pilot, nrow=G))
  # colnames(d.augdepths) <- colnames(pilotdata)

  # Increase to desired depths
  for(n in 1:N.pilot){
    d.augdepths[,n] <- rmultinom(1, size=depth*1e6, prob=pilotdata[,n])
  }

  ### Extract genewise dispersions from pilot data

  # Turn count matrix input into DGElist object
  d <- DGEList(d.augdepths)
  # Calculate normalization factors
  d <- calcNormFactors(d)
  # Estimate dispersions via edgeR
  d <- estimateGLMCommonDisp(d)
  d <- estimateGLMTrendedDisp(d)
  d <- estimateGLMTagwiseDisp(d)

  dispersion <- d$tagwise.dispersion

  ### Get separate means per group from pilot data

  # If group is not provided, set 2 groups with equal number of replicates
  if (is.null(group.pilot)){
    group.pilot <- as.factor(rep(c("A", "B"), each=N.pilot/2))
  }
  group.pilot <- as.factor(group.pilot)
  d$group <- group.pilot

  # Subset the groups
  d.grp1 <- d[,d$group==levels(d$group)[1]]
  d.grp2 <- d[,d$group==levels(d$group)[2]]

  # Save the average log CPM, for each group
  # Note that aveLogCPM(x), rowMeans(cpm(x,log=TRUE)) and log2(rowMeans(cpm(x)) all give slightly different results.
  # Use aveLogCPM(x) here as it doesn't return any Inf values and its values are close to log2(rowMeans(cpm(x))) which I had used previously
  grp1.AveLogCPM <- aveLogCPM(d.grp1)
  grp2.AveLogCPM <- aveLogCPM(d.grp2)

  # Get genewise lambdas (for each roup) to be used for extrapolated dataset
  lambda1 <- 2^(grp1.AveLogCPM)
  lambda2 <- 2^(grp2.AveLogCPM)
  lambda1 <- lambda1/sum(lambda1)
  lambda2 <- lambda2/sum(lambda2)

  # Expand genewise lambdas and dispersions
  Lambda1 <- expandAsMatrix(lambda1, dim=c(G, nreps))
  Lambda2 <- expandAsMatrix(lambda2, dim=c(G, nreps))
  Lambda <- cbind(Lambda1, Lambda2)
  Dispersion <- expandAsMatrix(dispersion, dim = c(G, nreps*2))

  ### NB sampling, given means, dispersions, library sizes

  lib.size <- rep(depth*1e6, nreps*2)
  extrapolatedData <- matrix(rnbinom(G*nreps*2, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = G, ncol = nreps*2)
  rownames(extrapolatedData) <- paste("ids", 1:G, sep = "")
  colnames(extrapolatedData) <- c(paste(levels(group.pilot)[1], 1:nreps, sep=""),
                                  paste(levels(group.pilot)[2], 1:nreps, sep=""))

  # Save extrapolated data with same filename as pilot data; just make sure it goes in appropriate directory
  filename <- paste("../Datasets/ExtrapolatedData", filename.pilot, sep="/")
  save(extrapolatedData, file=filename)
}

