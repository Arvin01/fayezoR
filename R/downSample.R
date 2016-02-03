#' Down-sample RNA-seq data.
#' Down-sample full dataset to a smaller dataset of specified size.
#' @param dataset Numeric count matrix; full dataset from which to down-sample. Note that for my purposes, this full dataset is simulated.
#' @param group Vector or factor giving the experimental group/condition for each sample/library in the full dataset.
#' @param nreps Desired number of replicates per group.
#' @param depth Desired depth (library size) in millions for each sample.
#' @param nsims Number of simulations to run, i.e. number of datasets to generate of this size.
#' @export
#' @examples
#'
#' @return Doesn't return anything, just saves the down-sampled datasets as RData.

downSample <- function(dataset, group=NULL, nreps, depth, nsims=20){

  N <- nreps*2
  G <- dim(dataset)[1]

  # If group isn't provided, just assume the groups are equally divided
  if (is.null(group)){
    group <- as.factor(rep(c("A","B"), each=dim(dataset)[2]/2))
  }
  group <- as.factor(group)
  samplenames <- paste(group, c(1:length(which(group==levels(group)[1])), 1:length(which(group==levels(group)[2]))), sep="")

  # Downsample nims times
  for (k in 1:nsims){

    # File-naming convention: N20_D1.2_sim5.RData
    filename <- paste("N", nreps, "_D", depth, "_sim", k, ".RData", sep="")

    ### First, select nreps from each group in full dataset
    n1 <- length(which(group=="A"))
    n2 <- length(which(group=="B"))
    # indices for randomly selected control and treatment samples
    control.idx <- sample(1:n1, nreps)
    trt.idx <- sample(1:n2, nreps) + n1
    # dataset subsetted by randomly selected samples
    subReps <- as.matrix(dataset[,c(control.idx, trt.idx)], ncol=N)
    subReps.depths <- colSums(subReps)

    ### Next, down-sample each sample in the subReps and put it in a new matrix
    ### Two ways to down-sample, code both ways for now

    # Initialize the matrix
    subData <- data.frame(matrix(ncol=N, nrow=G))
    colnames(subData) <- samplenames[c(control.idx, trt.idx)]

    # Multinomial sampling at desired depth
    for(n in 1:N){
      subData[,n] <- rmultinom(1, size=depthx, prob=subReps[,n])
    }

    # Binomial sampling (Li and Tibshirani); this is a special case of multinomial sampling and won't be used
    # for(i in 1:N){
    #   for(j in 1:G){
    #     subData[j,i] <- rbinom(1, size=subReps[j,i], prob=depth/subReps.depths[i])
    #   }
    # }

    # Save this dataset as an RData object to be used again later

    save(subData, file=filename)
  }
}

