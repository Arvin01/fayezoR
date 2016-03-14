#' Apply ZINB on dataset.
#' This is a general function which applies ZINB on a given dataset.
#' @param dataset Numeric count matrix of read counts.
#' @param filename Filename of the dataset; ZINB results on that dataset will be named accordingly
#' @param group Vector or factor giving the experimental group/condition for each sample/library in the pilot dataset. If not provided, groups are assumed equally divided.
#' @export
#' @examples
#'
#' @return Applies ZINB to given dataset, and returns results table

apply.zinb <- function(dataset, group=NULL){

  require(pscl)

  G <- dim(dataset)[1]
  N <- dim(dataset)[2]
  libsizes <- colSums(dataset)

  if (is.null(group)){
    group = as.factor(rep(c("A", "B"), each=N/2))
  }
  group = as.factor(group)

  data.df <- as.data.frame(cbind(t(dataset), group))

  # Gets p-values from a call of zeroinfl()
  # Input is a gene vector y
  zeroinf.p <- function(y){
    fit <- NA
    try(fit <- zeroinfl(y ~ group | 1, data=data.df, dist="negbin", offset=log(libsizes)))
    if (class(fit)=="zeroinfl"){
      return(coef(summary(fit))$count[2,])
    }
    else{
      return(NA)
    }
  }

  # Apply ZINB on each gene
  ZINBres <- apply(data.df[,1:G], 2, zeroinf.p)

  # Remove ZINBres which are NAs
  lengths <- sapply(ZINBres, length)
  ind.NA <- which(lengths!=4)
  ZINBres.noNA <- ZINBres[-ind.NA]
  names <- names(ZINBres.noNA)

  # Getting results in matrix form
  ind.keep <- which(lengths==4)
  ZINBres.df <- as.data.frame(matrix(unlist(ZINBres.noNA), ncol=4, byrow=TRUE), row.names=ind.keep)
  names(ZINBres.df) <- c("Beta", "SE", "ZValue", "PValue")

  # Putting in adjusted p-value
  ZINBres.df$FDR <- p.adjust(ZINBres.df$PValue, method="BH")
  ZINBres.df
}
