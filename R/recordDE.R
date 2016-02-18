#' Record the number of differentially expressed genes found by each of the edgeR results
#' @param reps Numerical vector of replicates per group considered.
#' @param depths Numerical vector of depths considered (e.g. c(0.5, 1, 1.5)*1e6)
#' @param nsims Number of simulations that were run per reps and depths
#' @export
#' @examples
#'
#' @return Doesn't return anything, just writes results to file as "recordDE.txt". Table contains the number of DE genes for each reps x depth x sim dataset on which edgeR was applied.

recordDE <- function(reps, depths, nsims){
  recordDE.df <- data.frame(reps=0, depth=0, filename=NA, numDE=0)
  for (rep in as.character(reps)){
    for (depth in sprintf("%.1f", depths/1e6)){
      for (sim in 1:nsims){
        # load edgeR output (edgeRres) which was saved as RData
        filename <- paste("N", rep, "_D", depth, "_sim", sim, ".RData", sep="")
        load(filename)
        # get the number of differentially expressed genes at FDR 0.05
        numDE <- sum(edgeRres$FDR <= 0.05, na.rm=TRUE)
        # record the information
        recordDE.df <- rbind(recordDE.df, c(rep, depth, filename, numDE))
      }
    }
  }
  # Write the recorded DE information to a table
  recordDE.df <- recordDE.df[-1,]
  write.table(recordDE.df, file="recordDE.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}
