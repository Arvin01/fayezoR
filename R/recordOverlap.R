#' Records overlaps with gold standard DE gene list, and gets accuracy measures
#' @param reps Numerical vector of replicates per group considered.
#' @param depths Numerical vector of depths considered (e.g. c(0.5, 1, 1.5)*1e6)
#' @param nsimss Number of simulations that were run per reps and depths
#' @export
#' @examples
#'
#' @return

# recordDE.overlap as "recordDE.TP.FP.varyFDR.txt"
# agg.overlap as agg.DE.TPR.FPR.byRepDepthFDR.txt"

recordOverlap <- function(reps, depths, nsims){
  recordDE.overap <- data.frame(reps=NA, depth=NA, filename=NA, FDR=NA, numDE=NA, TP=NA, FP=NA)

  # Get the "gold standard" from original "full" dataset
  goldstandard.data <- read.csv(paste("edgeR.rep", max(reps), ".depth", sprintf("%.1f", max(depths*1e-6)), "M.sim1.csv", sep=""))
  goldstandard.DElist <- goldstandard.data[ goldstandard.data$FDR <= 0.001, 1]
  numPos <- length(goldstandard.DElist)
  numNeg <- dim(goldstandard.data)[1]-numPos


}

recordOverlap <- function(reps, depths, nsims){
  recordDE.overlap <- data.frame(reps=NA, depth=NA, software=NA, filepath=NA, FDR=NA, numDE=NA, TP=NA, FP=NA)
  # Get the "gold standard" from highest depth, highest reps at FDR 0.001 to compute power and ROC curve
  goldstandard.data <- read.csv(paste("edgeR.rep", max(reps), ".depth", sprintf("%.1f", max(depths*1e-6)), "M.sim1.csv", sep=""))
  goldstandard.DElist <- goldstandard.data[ goldstandard.data$FDR <= 0.001, 1]
  numPos <- length(goldstandard.DElist)
  numNeg <- dim(goldstandard.data)[1]-numPos

  for (rep in as.character(reps)){
    for (depth in sprintf("%.1f", depths/1e6)){
      for (sim in 1:nsims){
        filepath <- paste("edgeR",  ".rep", rep, ".depth", depth, "M.sim", sim, ".csv", sep="")
        try(edgeR.res <- read.csv(filepath, header=TRUE) )
        if (dim(edgeR.res)[1] == 0 ) {
          next;
        }
        #vary FDR level
        for (fdr.var in seq(from=0,to=1,by=0.01)){
          numDE <- sum(edgeR.res$FDR <= fdr.var, na.rm=TRUE)
          # pos = genes that were both in gold list and found DE in edgeR (true pos)
          TP <- length(intersect(edgeR.res[ edgeR.res$FDR <= fdr.var, 1], goldstandard.DElist))
          # neg = genes that were found DE in edgeR but not on gold list (false pos)
          FP <- length(setdiff(edgeR.res[ edgeR.res$FDR <= fdr.var, 1], goldstandard.DElist))
          recordDE.overlap <- rbind(recordDE.overlap, c(rep, depth, "edgeR", filepath, fdr.var, numDE, TP, FP))
          #cat ("recordOverlap()", "\t", filepath, "\t", "FDR level", fdr.var, "\t", "numDE", numDE, "\n")
        }
      }
    }
  }
  recordDE.overlap <- recordDE.overlap[-1,]

  # Saving things to plot: reps, depth, numDE, FDR level, TP, FP
  to.plot <- data.frame(reps=as.numeric(recordDE.overlap$reps), depth=as.numeric(recordDE.overlap$depth), FDR=as.numeric(recordDE.overlap$FDR), numDE=as.numeric(recordDE.overlap$numDE), TPR=as.numeric(recordDE.overlap$TP)/numPos, FPR=as.numeric(recordDE.overlap$FP)/numNeg)

  # Aggregating things to plot by taking mean over numreps, numreads, and FDR
  agg.overlap <- aggregate(to.plot, by=list(to.plot$reps, to.plot$depth, to.plot$FDR), FUN=mean)[,-(1:3)]

  write.table(agg.overlap, "agg.DE.TPR.FPR.byRepDepthFDR.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(recordDE.overlap, file="recordDE.TP.FP.varyFDR.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}
