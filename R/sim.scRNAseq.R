#' Simulating scRNA-seq expression data
#' Extract parameters from real dataset to simulate scRNA-seq expression data with differentially expressed genes; options for imposing dropouts on the data.
#' @param dataset Numeric matrix of read counts from which to extract parameters.
#' @param nlibs Desired number of replicates per group.
#' @param nTags Desired number of genes. If not provided, simulated dataset will contain as many genes as the real dataset.
#' @param design Model matrix (without an intercept) that you would like to simulate from
#' @param beta Set of coefficients for the model matrix (must have same number of columns as mod)
#' @param lib.size Numeric vector giving the total count (sequence depth) for each library. If not provided, library sizes are extracted from the real dataset.
#' @param flibs Some wiggle room for the library sizes extracted from real
#'   dataset.
#' @param drop.low.lambda Logical, whether to drop low lambdas.
#' @param drop.extreme.dispersion Proportion of extreme dispersions to drop.
#' @param add.dropout Logical, whether to impose dropouts.
#' @param pDropout Numeric vector of dropout probabilities to consider. Will only be used if \code{add.dropout=TRUE}.
#' @param drop.method Dropout method, with two options: zero-inflation (default), or low-rate Poisson. Will only be used if \code{add.dropout=TRUE}.
#' @param drop.lambda Mean of low-rate Poisson for the Poisson-dropout scenario. Will only be used if \code{add.dropout=TRUE} and \code{drop.method="poisson"}.
#' @param verbose Logical, whether to print simulation progress to screen.
#' @param seed Optional seed, if it desired to replicate the simulation.
#' @export
#' @examples
#' # Simulate toy dataset with 3 reps per group and 20 genes
#' # data(prostatedata)
#' #simData <- sim.scRNAseq(dataset=prostatedata, nreps=3, nTags=20, drop.method=c("zero"))
#' NULL
#' @return DGElist object containing all simulation settings and results, including the count matrix and dropout-imposed count matrix.

sim.scRNAseq <- function(dataset, nlibs=NULL, nTags=NULL,
                         design=NULL, beta=NULL,
                         drop.low.lambda=TRUE, drop.extreme.dispersion=0.1,
                         lib.size=NULL, flibs=c(0.7,1.3),
                         #add.dropout=TRUE, pDropout=c(0.1,0.2,0.5,0.8), drop.method=c("zero", "poisson"), drop.lambda=1,
                         verbose=TRUE, seed=NULL){

  require(edgeR)

  # Make sure that model matrix doesn't have intercept
  if( !is.null(design) ){
    if( any(apply(design,2,function(x){all(x==1)})) ){
      stop("Model matrix should not have an intercept.\n")
    }
    # Make sure betas have same columns as model matrix
    if( !is.null(beta) & (ncol(design) != ncol(beta)) ){
      stop("Beta coefficients must have same number of columns as model matrix.\n")
    }
  }

  ##### Preparing the dataset

  if(verbose) cat("Extracting parameters from dataset.\n")

  # Sample parameters from dataset
  dataset.params <- getDataset(counts = dataset,
                        drop.extreme.dispersion = drop.extreme.dispersion,
                        drop.low.lambda = drop.low.lambda)

  # Set nTags and nlibs when neither design or beta provided
  if(is.null(design) & is.null(beta)){
    if(is.null(nTags)){
      nTags <- dataset.params$nTags
    }
    if(is.null(nlibs)){
      nlibs <- dataset.params$nlibs
    }

    # Set nTags and nlibs when only design provided
  } else if(!is.null(design) & is.null(beta)){
    if(is.null(nlibs)){
      nlibs = dim(design)[1]
    } else if(nlibs!=dim(design)[1]){
      message("nlibs must be same size as design...\n setting nlibs to be same size of design.")
      nlibs = dim(design)[1]
    }
    if(is.null(nTags)){
      nTags <- dataset.params$nTags
    }

    # Set nTags and nlibs when only beta provided
  } else if(is.null(design) & !is.null(beta)){
    if(is.null(nTags)){
      nTags = dim(beta)[1]
    } else if(!is.null(nTags) & nTags!=dim(beta)[1]){
      message("nTags must be same size as beta...\n setting nTags to be same size of beta.")
      nTags = dim(beta)[1]
    }
    if(is.null(nlibs)){
      nlibs <- dataset.params$nlibs
    }

    # Set nTags and nlibs when both design and beta provided
  } else if(!is.null(design) & !is.null(beta)){
    if(is.null(nlibs)){
      nlibs = dim(design)[1]
    } else if(nlibs!=dim(design)[1]){
      message("nlibs must be same size as design...\n setting nlibs to be same size of design.")
      nlibs = dim(design)[1]
    }
    if(is.null(nTags)){
      nTags = dim(beta)[1]
    } else if(nTags!=dim(beta)[1]){
      message("nTags must be same size as beta...\n setting nTags to be same size of beta.")
      nTags = dim(beta)[1]
    }
  }

  # If lib.size not already specified, sample them randomly from lib sizes sampled from real dataset, with some wiggle room
  if(is.null(lib.size)){
    set.seed(seed)
    flibs <- rep(flibs, length.out=2)
    lib.size <- runif(nlibs, min = flibs[1]*median(dataset.params$lib.size), max = flibs[2]*median(dataset.params$lib.size))
  }

  # Create new DGEList for storing dataset information
  dat <- new("DGEList", list(dataset.params = dataset.params,
                             nTags = nTags,
                             nlibs = nlibs,
                             lib.size = lib.size,
                             design = design,
                             beta = beta))

  ##### Resampling parameters and getting the count matrix

  # Sampling Lambda and Dispersion parameters to be used in simulated dataset
  # Adds Lambda and Dispersion matrices to dat object
  if(verbose) cat("Re-sampling parameters.\n")
  dat <- sample.fun(dat, seed=seed)

  dat <- sim.fun(dat, seed=seed)

  # Add dropouts
  # if(add.dropout){
  #   if(verbose) message("Adding dropouts.\n")
  #   drop.method <- match.arg(drop.method)
  #   dat <- dropout.fun(dat, pDropout, drop.method, drop.lambda, seed=seed)
  #   colnames(dat$counts.dropouts) <- samplenames
  # }

  dat
}
