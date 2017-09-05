#' plotFocals
#'
#' @description
#' Plots the SNP array data with the identified focal gains and losses on top.
#'
#' @details
#' Creates a plot with the LRR and BAF values along the genome and adds
#' rectangles showing the identified focal gains and deletions. Optionally,
#' the SNPs in each of the focal events may be highlighted.
#'
#' @usage plotFocals(find.focal.results,
#' plot.gain=TRUE, plot.gain.points=FALSE, gain.col="#FFBD07AA", gain.border.col="#FFBD07AA",
#' plot.del=TRUE, plot.del.points=FALSE, del.col="#00A6EDAA", del.border.col="#00A6EDAA",
#' plot.mean.lrr=FALSE, mean.lrr.col="#B0FFA5AA", mean.lrr.border.col=NA, mean.lrr.line.col="#FF7C30",
#' main="Focal Gains and Losses", chromosomes="canonical", zoom=NULL,
#' lrr.min=-4, lrr.max=2, total.height=1, bottom=0, margin=0.05,
#' focal.points.cex=0.6, points.cex=0.3, labels.cex=1.5, main.cex=2, axis.cex=1.2, chr.cex=1.5)
#'
#'
#' @param find.focal.results The results from findFocals
#' @param plot.gain (defaults to TRUE)
#' @param plot.gain.points (defaults to FALSE)
#' @param gain.col (defaults to "#FFBD07AA")
#' @param gain.border.col (defaults to "#FFBD07AA")
#' @param plot.del (defaults to TRUE)
#' @param plot.del.points (defaults to FALSE)
#' @param del.col (defaults to "#00A6EDAA")
#' @param del.border.col (defaults to "#00A6EDAA")
#' @param plot.mean.lrr (defaults to FALSE)
#' @param mean.lrr.col (defaults to "#B0FFA5AA")
#' @param mean.lrr.border.col (defaults to NA)
#' @param mean.lrr.line.col (defaults to "#FF7C30")
#' @param main (defaults to "Focal Gains and Losses")
#' @param chromosomes (defaults to "canonical")
#' @param zoom (defaults to NULL)
#' @param lrr.min (defaults to -4)
#' @param lrr.max (defaults to 2)
#' @param total.height (defaults to 1)
#' @param bottom (defaults to 0)
#' @param margin (defaults to 0.05)
#' @param focal.points.cex (defaults to 0.6)
#' @param points.cex (defaults to 0.3)
#' @param labels.cex (defaults to 1.5)
#' @param main.cex (defaults to 2)
#' @param axis.cex (defaults to 1.2)
#' @param chr.cex (defaults to 1.5)
#'
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
#' @export plotFocals
#'

plotFocals <- function(find.focal.results,
                       plot.gain=TRUE, plot.gain.points=FALSE, gain.col="#FFBD07AA", gain.border.col="#FFBD07AA",
                       plot.del=TRUE, plot.del.points=FALSE, del.col="#00A6EDAA", del.border.col="#00A6EDAA",
                       plot.mean.lrr=FALSE, mean.lrr.col="#B0FFA5AA", mean.lrr.border.col=NA, mean.lrr.line.col="#FF7C30",
                       main="Focal Gains and Losses", chromosomes="canonical", zoom=NULL,
                       lrr.min=-4, lrr.max=2, total.height=1, bottom=0, margin=0.05,
                       focal.points.cex=0.6, points.cex=0.3, labels.cex=1.5, main.cex=2, axis.cex=1.2, chr.cex=1.5) {


  kp <- plotSNPData(snps = find.focal.results$snps, main = main, chromosomes = chromosomes, zoom=zoom,
                    lrr.min=lrr.min, lrr.max=lrr.max, total.height = total.height, bottom = bottom, margin = margin,
                    points.cex = points.cex, labels.cex = labels.cex, main.cex = main.cex, axis.cex = axis.cex, chr.cex = chr.cex)

  lrr.r0 <- bottom
  lrr.r1 <- (total.height-margin)/2+lrr.r0
  baf.r0 <- lrr.r1 + margin
  baf.r1 <- (total.height-margin)/2 + baf.r0


  if(plot.gain.points==TRUE) {
    gain.snps <- find.focal.results$snps[find.focal.results$snps.in.gains]
    if(length(gain.snps)>0) {
      kpPoints(kp, data=gain.snps, y=gain.snps$lrr, col=gain.col,  r0=lrr.r0, r1=lrr.r1, ymin=lrr.min, ymax=lrr.max, cex=focal.points.cex, pch=16)
      kpPoints(kp, data=gain.snps, y=gain.snps$baf, col=gain.col,  r0=baf.r0, r1=baf.r1, ymin=0, ymax=1, cex=focal.points.cex, pch=16)
    }
  }
  if(plot.del.points==TRUE) {
    del.snps <- find.focal.results$snps[find.focal.results$snps.in.deletions]
    if(length(del.snps)>0) {
      kpPoints(kp, data=del.snps, y=del.snps$lrr, col=del.col,  r0=lrr.r0, r1=lrr.r1, ymin=lrr.min, ymax=lrr.max, cex=focal.points.cex, pch=16)
      kpPoints(kp, data=del.snps, y=del.snps$baf, col=del.col,  r0=baf.r0, r1=baf.r1, ymin=0, ymax=1, cex=focal.points.cex, pch=16)
    }
  }

  if(plot.mean.lrr==TRUE) {
    kpPlotRibbon(kp,  data=find.focal.results$snps,
                 y0=find.focal.results$snps$mean.lrr - find.focal.results$params$fdel.sd.lim*find.focal.results$lrr.sd,
                 y1=find.focal.results$snps$mean.lrr + find.focal.results$params$fgain.sd.lim*find.focal.results$lrr.sd,
                 col=mean.lrr.col, border=mean.lrr.border.col,
                 r0=lrr.r0, r1=lrr.r1, ymin=lrr.min, ymax=lrr.max)
    kpLines(kp,  data=find.focal.results$snps, y=find.focal.results$snps$mean.lrr, col=mean.lrr.line.col,
            r0=lrr.r0, r1=lrr.r1, ymin=lrr.min, ymax=lrr.max)
  }

  if(plot.gain==TRUE) {
    kpRect(kp, data=find.focal.results$gains, r0=bottom, r1=bottom+total.height, y0=0, y1=1, border=gain.border.col, col=gain.col)
  }
  if(plot.del==TRUE) {
    kpRect(kp, data=find.focal.results$deletions, r0=bottom, r1=bottom+total.height, y0=0, y1=1, border=del.border.col, col=del.col)
  }


  invisible(kp)
}

