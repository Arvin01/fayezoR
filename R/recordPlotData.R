#' Record the number of differentially expressed genes found by each of the edgeR results
#' @param reps Numerical vector of replicates per group considered.
#' @param depths Numerical vector of depths considered, in millions
#' @param nsims Number of simulations that were run per reps and depths
#' @export
#' @examples
#'
#' @return Writes plotting data to file

recordPlotData <- function(reps, depths, nsims){

  setwd("/scratch/radon/z/zheng64/ExpDesignApp/Results")

  ########## RECORD DE DATA

  recordDE.df <- data.frame(reps=0, depth=0, numDE=0)
  for (rep in as.character(reps)){
    for (depth in as.character(depths)){
      for (sim in 1:nsims){
        # load edgeR output (edgeRres) which was saved as RData
        filename <- paste("N", rep, "_D", depth, "_sim", sim, ".RData", sep="")
        load(filename)
        # get the number of differentially expressed genes at FDR 0.05
        edgeRres$FDR <- p.adjust(edgeRres$PValue, method="BH")
        numDE <- sum(edgeRres$FDR <= 0.05, na.rm=TRUE)
        # record the information
        recordDE.df <- rbind(recordDE.df, c(rep, depth, numDE))
      }
    }
  }
  # Write the recorded DE information to a table
  recordDE.df <- recordDE.df[-1,]

  # Aggregating taking mean over repsh/depth
  agg.DE <- aggregate(sapply(recordDE.df, as.numeric), by=list(recordDE.df$reps, recordDE.df$depth), FUN=mean)[,-(1:2)]

  write.table(recordDE.df, file="../Plots/recordDE.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(agg.DE, file="../Plots/aggDE.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

  ########## RECORD OVERLAP DATA

  recordDE.overlap <- data.frame(reps=NA, depth=NA, FDR=NA, numDE=NA, TPR=NA, FPR=NA)

  # Get goldstandard DE list as those found by edgeR
  # load("fullData.RData")
  # edgeRres$FDR <- p.adjust(edgeRres$PValue, method="BH")
  # goldstandard.DElist <- which(edgeRres$FDR <= 0.05)
  # numPos <- length(goldstandard.DElist)
  # numNeg <- dim(edgeRres)[1]-numPos

  # Alternatively, get goldstandard DE list as the true simulated DE genes
  load("../Datasets/fullData.RData")
  goldstandard.DElist <- 1:1000
  numPos <- length(goldstandard.DElist)
  numNeg <- fullData$nTags-numPos

  for (rep in as.character(reps)){
    for (depth in as.character(depths)){
      for (sim in 1:nsims){

        # load edgeR output (edgeRres) which was saved as RData
        filename <- paste("N", rep, "_D", depth, "_sim", sim, ".RData", sep="")
        load(filename)
        edgeRres$FDR <- p.adjust(edgeRres$PValue, method="BH")

        #vary FDR level
        for (fdr.var in seq(from=0,to=1,by=0.01)){
          edgeRres$FDR <- p.adjust(edgeRres$PValue, method="BH")
          numDE <- sum(edgeRres$FDR <= fdr.var, na.rm=TRUE)
          # TP = genes that were both in gold list and found DE in edgeR
          TP <- length(intersect(which(edgeRres$FDR <= fdr.var), goldstandard.DElist))
          # FP = genes that were found DE in edgeR but not on gold list
          FP <- length(setdiff(which(edgeRres$FDR <= fdr.var), goldstandard.DElist))
          recordDE.overlap <- rbind(recordDE.overlap, c(rep, depth, fdr.var, numDE, TP/numPos, FP/numNeg))
        }
      }
    }
  }
  recordDE.overlap <- recordDE.overlap[-1,]

  # Aggregating things to plot by taking mean over numreps, numreads, and FDR
  agg.overlap <- aggregate(sapply(recordDE.overlap, as.numeric), by=list(recordDE.overlap$reps, recordDE.overlap$depth, recordDE.overlap$FDR), FUN=mean)[,-(1:3)]

  write.table(agg.overlap, "../Plots/agg.overlap.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  write.table(recordDE.overlap, file="../Plots/recordDE.overlap.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}
