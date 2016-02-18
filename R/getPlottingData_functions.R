

### ---------------------------------------------------------------------
### RECORDDE()
### Gets numDE results for every output file
### ---------------------------------------------------------------------

### Writes to file
# recordDE.df as "recordDE.FDR.05.txt"

recordDE <- function(reps, depths, nSim){
  recordDE.df <- data.frame(reps=0, depth=0, software=NA, filepath=NA, numDE=0)
  for (rep in as.character(reps)){
    for (depth in sprintf("%.1f", depths/1e6)){
      for (sim in 1:nSim){
        # filepath to edgeR output
        filepath <- paste("edgeR.rep", rep, ".depth", depth, "M.sim", sim, ".csv", sep="")    
        # try reading in the filename, give error message if didn't work
        try(edgeR.res <- read.csv(filepath, header=TRUE))
        if (dim(edgeR.res)[1] == 0){
          next;
        }
        # get the number of differentially expressed genes at FDR 0.05
        numDE <- sum(edgeR.res$FDR <= 0.05, na.rm=TRUE)
        # record the information
        recordDE.df <- rbind(recordDE.df, c(rep, depth, "edgeR", filepath, numDE))
        #cat ("recordDE()", "\t", filepath, "\t", numDE, "\n")
      }
    }
  }
  # Write the recorded DE information to a table
  recordDE.df <- recordDE.df[-1,]
  write.table(recordDE.df, file="recordDE.FDR.05.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

### ---------------------------------------------------------------------
### RECORDOVERLAP()
### Records overlaps with gold standard DE gene list, and gets accuracy measures
### ---------------------------------------------------------------------

### Writes to file
# recordDE.overlap as "recordDE.TP.FP.varyFDR.txt"
# agg.overlap as agg.DE.TPR.FPR.byRepDepthFDR.txt"

recordOverlap <- function(reps, depths, nSim){
  recordDE.overlap <- data.frame(reps=NA, depth=NA, software=NA, filepath=NA, FDR=NA, numDE=NA, TP=NA, FP=NA)
  # Get the "gold standard" from highest depth, highest reps at FDR 0.001 to compute power and ROC curve
  goldstandard.data <- read.csv(paste("edgeR.rep", max(reps), ".depth", sprintf("%.1f", max(depths*1e-6)), "M.sim1.csv", sep=""))
  goldstandard.DElist <- goldstandard.data[ goldstandard.data$FDR <= 0.001, 1]
  numPos <- length(goldstandard.DElist)
  numNeg <- dim(goldstandard.data)[1]-numPos
  
  for (rep in as.character(reps)){
    for (depth in sprintf("%.1f", depths/1e6)){
      for (sim in 1:nSim){
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

### ---------------------------------------------------------------------
### RECORDCPM()
### Gets logCPM results for top/med/low 100 genes across depth/rep combos
### ---------------------------------------------------------------------

### Writes to file
# record_logCPM as "recordlogCPM.txt"
# record_logCPM_results as "agg.recordlogCPM.txt"

recordCPM <- function(reps, depths, nSim){
  # Getting top/middle/low 100 logCPM and logFC genes from gold standard dataset
  goldstandard <- read.csv(paste("edgeR.rep", max(reps), ".depth", sprintf("%.1f",max(depths*1e-6)), "M.sim1.csv", sep=""))
  logCPM_top100 <- goldstandard[ order( goldstandard[,"logCPM"],decreasing=TRUE )[1:100], 2]
  logCPM_medium100 <- goldstandard[ order( goldstandard[,"logCPM"],decreasing=TRUE )[5001:5100], 2]
  logCPM_low100 <- goldstandard[ order( goldstandard[,"logCPM"],decreasing=TRUE )[12001:12100], 2]
  
  nrows = length(depths)*length(reps)*nSim*length(logCPM_top100)*3
  
  # For each rep/depth combo, get the logCPM of the top/middle/low 100 genes
  record_logCPM <- data.frame( gene=rep(NA,nrows), reps=rep(NA,nrows), depth=rep(NA,nrows), sim=rep(NA,nrows), logCPM=rep(NA,nrows) )
  i <- 1
  for (rep in as.character(reps)){
    for (depth in sprintf("%.1f", depths/1e6)){
      for (sim in 1:nSim ){
        file_path <- paste("edgeR", ".rep", rep, ".depth", depth, "M.sim", sim, ".csv", sep="")
        edgeR_res <- read.csv(file_path, header=TRUE)
        logCPM1 <- edgeR_res[ match(logCPM_top100,edgeR_res$geneID), "logCPM"]
        logCPM2 <- edgeR_res[ match(logCPM_medium100,edgeR_res$geneID), "logCPM"]
        logCPM3 <- edgeR_res[ match(logCPM_low100,edgeR_res$geneID), "logCPM"]
        record_logCPM[i:(i+299),] <-  cbind(c(logCPM_top100, logCPM_medium100, logCPM_low100), rep, depth, sim, c(logCPM1,logCPM2,logCPM3))
        i <- i+300
	#cat("recordlogCPM", "\t", "Reps", rep, "\t", "Depth", depth, "\t", "sim", sim, "\n")
      }
    }  
  }
  write.table(record_logCPM, "recordlogCPM.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  
  nrows = length(depths)*length(reps)*length(logCPM_top100)*3
  record_logCPM_results <- data.frame( gene=rep(NA, nrows), reps=rep(NA,nrows), depth=rep(NA,nrows), logCPM_std=rep(NA,nrows), logCPM_mean=rep(NA,nrows), logCPM_mean_diff=rep(NA,nrows) )
  i <- 1
  for (x in reps){
    for (y in sprintf("%.1f", depths/1e06)){
      #cat("recordlogCPM_results", "Reps", x, "\t", "Depths", y, "\n")
      for(g in c(logCPM_top100, logCPM_medium100, logCPM_low100)){
        #cat( "Reps", x, "\t", "Depths", y, "\t", "i", i, "\n" )
        
        # Get the logCPMs for this gene in each reps/depth combination 
        # There should be one per simulation
        logCPM <- as.numeric( subset(record_logCPM, gene==g & reps==x & depth==y, select=logCPM)[,1] )
        
        standard_logCPM <- goldstandard[ goldstandard$geneID==g, "logCPM" ]
        # Report the SD, mean of logCPM for this gene/reps/depth combination
        record_logCPM_results[i,] <- c( g, x, y, sd(logCPM), mean(logCPM), mean(abs(logCPM-standard_logCPM)))
        i<-i+1
      }
    }
  }
  write.table(record_logCPM_results, "agg.recordlogCPM.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)  
}
