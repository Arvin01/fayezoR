#' Custom multiplot function for ggplot2
#' @param reps Numerical vector of replicates per gorup to plot
#' @param depths Numerical vector of depths to plot, in millions
#' @param cols Number of columns to plot (default)
#' @param width, height The width and height of the graphics region in inches; default values are 10
#' @param plotROC Logical (default TRUE), whether to make ROC plots
#' @param plotDE Logical (default TRUE), whether to make numDE plots
#' @param plotPower Logical (default TRUE), whether to make power plots
#' @export
#' @examples
#'
#' @return Makes plots of ROC curves, number of detected DE features, and power, all saved as PDFs

makePlots <- function(reps, depths, plotROC=TRUE, plotDE=TRUE, plotPower=TRUE, cols=2, width=10, height=10){

  require(ggplot2)

  #recordDE: reps, depth, numDE
  recordDE <- as.data.frame(read.table("recordDE.txt", header=T))

  #aggDE: reps, depth, numDE (aggregated)
  aggDE <- as.data.frame(read.table("aggDE.txt", header=T))

  #recordDE.overlap: reps, depth, FDR, numDE, TPR, FPR
  recordDE.overlap <- as.data.frame(read.table("recordDE.overlap.txt", header=T))

  #agg.overlap: reps, depth, FDR, numDE, TPR, FPR (aggregated)
  agg.overlap <- as.data.frame(read.table("agg.overlap.txt", header=T))

  ##### ROC curves, one for each level of replication; lines for each depth
  if (plotROC){
    message("Plotting ROC curves...")
    i <- 1
    plotlist <- list()
    for (rep in reps){
      p <- ggplot(subset(agg.overlap, reps==rep & depth %in% depths),
                  aes(x=FPR, y=TPR, color=factor(depth), group=factor(depth)))
      p <- p + geom_line() +
        scale_x_continuous("False Positive Rate", breaks=seq(0,1,0.2)) +
        scale_y_continuous("True Positive Rate", breaks=seq(0,1,0.2)) +
        scale_colour_brewer("# Reads (M)", breaks=depths, palette="Set1") +
        ggtitle(paste("ROC Curve for", rep, "Replicates", sep=" "))
      assign(paste("roc", rep, sep=""), p)
      plotlist[[i]] <- p
      i <- i+1
    }

    pdf("ROC_byReps.pdf", width=width, height=height)
    multiplot(plotlist=plotlist, cols=cols)
    dev.off()
  }

  ##### Number of DE genes found vs. depth/rep
  if(plotDE){
    message("Plotting number of DE genes...")
    p <- ggplot(recordDE,
                aes(x=depth, y=numDE, color=factor(reps), group=factor(reps)))
    p <- p + stat_smooth(method="loess", span=0.5) +
      scale_x_continuous("Number of Reads (M)") +
      scale_y_continuous("Number of DE genes (FDR 0.05)")  +
      scale_colour_brewer("# Reps", breaks=reps, palette="Set1") +
      ggtitle("# DE genes vs. Reps vs. # Reads")

    #pdf("numDE_RepsDepths.pdf")
    ggsave("numDE_RepsDepths.pdf", p)
    #dev.off()
  }

  ##### Power vs. depth/rep
  if(plotPower){
    message("Plotting power...")
    # Table for power calculation
    # Using FDR 0.05, plot the power curves (TP)
    power_plot <- subset(recordDE.overlap, FDR==0.05)

    p <- ggplot(power_plot,
                aes(x=depth, y=TPR, color=factor(reps), group=factor(reps)))
    p <- p  + stat_smooth(method="loess", span=0.5) +
      scale_x_continuous("# Reads", breaks=depths) +
      scale_y_continuous("Power at FDR 0.05") +
      scale_colour_brewer("# Reps", breaks=reps, palette="Set1" ) +
      ggtitle("Power at FDR 0.05")

    #pdf("Power_RepsDepth.pdf")
    ggsave("Power_RepsDepth.pdf", p)
    #dev.off()
  }
}
