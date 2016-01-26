#' Ratliff scRNA-seq prostate data
#'
#' Data from an unpublished scRNA-seq experiment on human prostate cancer cell lines, consisting of 209 cells from a Knockdown group (in which a particular gene was knocked out) and 190 cells from a Control group. Knockdown sample names begin with "S" and Control sample names begin with "C". Genes with average counts < 5 have been filtered out. 
#'
#' @name prostatedata
#' @docType data
#' @usage load(prostatedata)
#' @format Numeric count matrix.
#' @examples
#' 
"prostatedata"