#' plotSNPData
#'
#' @description
#' Plots the raw SNP array data using karyoploteR
#'
#' @details
#' Creates a plot with the LRR and BAF values along the genome
#'
#' @usage plotSNPData(snps, main="Raw Data", chromosomes="canonical", zoom=NULL, lrr.min=-4, lrr.max=2, total.height=1, bottom=0, margin=0.05, points.cex=0.3, labels.cex=1.5, main.cex=2, axis.cex=1.2, chr.cex=1.5)
#'
#' @param snps The SNP array data
#' @param main (defaults to "Raw Data")
#' @param chromosomes (defaults to "canonical")
#' @param zoom (defaults to NULL)
#' @param lrr.min (defaults to -4)
#' @param lrr.max (defaults to 2)
#' @param total.height (defaults to 1)
#' @param bottom (defaults to 0)
#' @param margin (defaults to 0.05)
#' @param points.cex (defaults to 0.3)
#' @param labels.cex (defaults to 1.5)
#' @param main.cex (defaults to 2)
#' @param axis.cex (defaults to 1.2)
#' @param chr.cex (defaults to 1.5)
#'
#'
#' @return
#' Invisibly returns the karyoplot object representing the plot. With it
#' it is possible to add other elements to the plot using standrad karyoploteR
#' functions
#'
#' @examples
#'
#'
#'
#' @export plotSNPData
#'
#' @import karyoploteR
#'


plotSNPData <- function(snps, main="Raw Data", chromosomes="canonical", zoom=NULL, lrr.min=-4, lrr.max=2, total.height=1, bottom=0, margin=0.05, points.cex=0.3, labels.cex=1.5, main.cex=2, axis.cex=1.2, chr.cex=1.5) {
  pp <- getDefaultPlotParams(plot.type=4)
  pp$bottommargin <- 20
  pp$data1inmargin <- 3
  kp <- plotKaryotype(chromosomes=chromosomes, plot.type=4, ideogram.plotter = NULL, labels.plotter = NULL, zoom=zoom, plot.params = pp)
  kpAddCytobandsAsLine(kp)
  kpAddChromosomeNames(kp, srt=45, cex=chr.cex)
  kpAddMainTitle(kp, main = main, cex=main.cex)

  #Plot LRR
  lrr.r0 <- bottom
  lrr.r1 <- (total.height-margin)/2+lrr.r0

  below.min <- snps$lrr < lrr.min
  above.max <- snps$lrr > lrr.max

  kpAddLabels(kp, r0=lrr.r0, r1=lrr.r1, labels = "LRR", srt=90, label.margin = 0.03, pos = 3, cex=labels.cex)
  kpAxis(kp, r0=lrr.r0, r1=lrr.r1, ymin=lrr.min, ymax=lrr.max, tick.pos = c(ceiling(lrr.min):floor(lrr.max)), cex=axis.cex)
  #plot the "in range" points
  kpPoints(kp, data=snps[!(below.min | above.max)], y=snps$lrr[!(below.min | above.max)], ymin=lrr.min, ymax=lrr.max,
           r0=lrr.r0, r1=lrr.r1, col="#333333", pch=16, cex=points.cex)

  #and plot the out of range points in red
  #Note: the rnorm is just a jitter in the y axis so they are more visible
  if(any(above.max)) {
    kpPoints(kp, data=snps[above.max], y=lrr.max - rnorm(length(which(above.max)), 0.02, 0.01), ymin=lrr.min, ymax=lrr.max,
             r0=lrr.r0, r1=lrr.r1, pch=16, cex=points.cex, col="red")
  }
  if(any(below.min)) {
    kpPoints(kp, data=snps[below.min], y=lrr.min + rnorm(length(which(below.min)), 0.02, 0.01) , ymin=lrr.min, ymax=lrr.max,
             r0=lrr.r0, r1=lrr.r1, pch=16, cex=points.cex, col="red")
  }

  #Plot BAF
  baf.r0 <- lrr.r1 + margin
  baf.r1 <- (total.height-margin)/2 + baf.r0

  kpAddLabels(kp, r0=baf.r0, r1=baf.r1, labels = "BAF", srt=90, label.margin = 0.03, pos = 3, cex=labels.cex)
  kpAxis(kp, r0=baf.r0, r1=baf.r1, ymin=0, ymax=1, cex=axis.cex)
  kpPoints(kp, data=snps, y=snps$baf, ymin=0, ymax=1, r0=baf.r0, r1=baf.r1, col="#333333", pch=16, cex=points.cex)

  invisible(kp)

}
