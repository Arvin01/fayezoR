#' Simulating scRNA-seq expression data This function blah blah blah description
#' @param dataset Numeric matrix of read counts from which to extract
#'   parameters.
#' @param nreps Desired number of replicates per group.
#' @param group Vector or factor giving the experimental group/condition for
#'   each sample/library.
#' @param nTags Desired number of genes. If not provided, simulated dataset will
#'   contain as many genes as the real dataset.
#' @param lib.size Numeric vector giving the total count (sequence depth) for
#'   each library. If not provided, library sizes are extracted from the real
#'   dataset.
#' @param flibs Some wiggle room for the library sizes extracted from real
#'   dataset.
#' @param pDiff Proportion of differentially expressed genes.
#' @param pUp Proportion of differentially expressed genes that are upregulated.
#' @param foldDiff Fold-change difference of differentially expressed genes.
#' @param drop.low.lambda Logical, whether to drop low lambdas.
#' @param drop.extreme.dispersion Proportion of extreme dispersions to drop. 
#' @param add.dropout Logical, whether to impose dropouts.
#' @param pDropout Numeric vector of dropout probabilities to consider. Will only be used if \code{add.dropout=TRUE}.
#' @param drop.method Dropout method, with two options: zero-inflation, or low-rate Poisson. Will only be used if \code{add.dropout=TRUE}.
#' @param drop.lambda Mean of low-rate Poisson for the Poisson-dropout scenario. Will only be used if \code{add.dropout=TRUE} and \code{drop.method="poisson"}.
#' @param verbose Logical, whether to print simulation progress to screen. 
#' @param seed Optional seed, if it desired to replicate the simulation.
#' @export
#' @examples
#' # Simulate toy dataset with 3 reps per group and 20 genes
#' # data(prostatedata)
#' #simData <- sim.scRNAseq(dataset=prostatedata, nreps=3, nTags=20, drop.method=c("zero"))
#' NULL
#' @return A numeric count matrix of desired size, simulated with parameters extracted from a real dataset, containing differentially expressed genes and (optionally imposed) dropouts.

sim.scRNAseq <- function(dataset, nreps, group=NULL, nTags=NULL, lib.size=NULL, flibs=c(0.7,1.3), pDiff=0.1, pUp=0.5, foldDiff=3, drop.low.lambda=TRUE, drop.extreme.dispersion=0.1, add.dropout=TRUE, pDropout=c(0.1,0.2,0.5,0.8), drop.method=c("zero", "poisson"), drop.lambda=1, verbose=TRUE, seed=NULL){
  
  require(edgeR)
  
  # If group is not provided, set 2 groups with equal number of replicates
  if (is.null(group)){
    group = as.factor(rep(0:1, each=nreps)) 
  }
  group = as.factor(group)
  nlibs = length(group)
  flibs <- rep(flibs, length.out=2)
  
  ##### Preparing the dataset
  
  if(verbose) message("Preparing dataset.\n") 
  
  # Sample parameters from dataset
  dataset <- getDataset(counts = dataset, 
                        drop.extreme.dispersion = drop.extreme.dispersion,
                        drop.low.lambda = drop.low.lambda)
  
  # Create new DGEList for storing dataset information
  dat <- new("DGEList", list(dataset = dataset, 
                             nTags = nTags, 
                             lib.size = lib.size, 
                             flibs = flibs, 
                             nlibs = nlibs, 
                             group = group, 
                             design = model.matrix(~group), 
                             pDiff= pDiff, 
                             pUp = pUp, 
                             foldDiff = foldDiff)) 
  
  ##### Simulate differentially expressed genes
  
  # If lib.size not already specified, sample them randomly from lib sizes sampled from real dataset, with some wiggle room
  if(is.null(lib.size)){
    set.seed(seed)
    dat$lib.size <- runif(nlibs, 
                          min = flibs[1]*median(dat$dataset$dataset.lib.size), 
                          max = flibs[2]*median(dat$dataset$dataset.lib.size))
  }
  
  # If nTags set to NULL, take them from real dataset
  if(is.null(nTags)){
    dat$nTags <- dat$dataset$dataset.nTags
  }
  
  # Sampling Lambda and Dispersion parameters to be used in simulated dataset
  # Adds Lambda and Dispersion matrices to dat object
  if(verbose) message("Sampling.\n")  
  dat <- sample.fun(dat, seed=seed) 
  
  # Alter Lambda's to incorporate differential expression
  if(verbose) message("Calculating differential expression.\n")	
  dat <- diff.fun(dat, seed=seed) 
  
  # Simulate data using these parameters
  if(verbose) message("Simulating data.\n")	
  dat <- sim.fun(dat, seed=seed)
  
  # Add dropouts
  if(add.dropout){
    if(verbose) message("Adding dropouts.\n")
    drop.method <- match.arg(drop.method)
    dat <- dropout.fun(dat, pDropout, drop.method, drop.lambda, seed=seed)
  }
  
  dat
}