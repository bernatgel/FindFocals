#' findFocals
#'
#' @description
#' Find clusters of extrem LRR values and determine if they might be classified
#' as a focal gain or deletion.
#'
#' @details
#' Given a GRanges with SNP data as the one created by loadSNPData, it
#' processess the LRR values to detect clusters of outliers, that is,
#' clusters of SNP a number of standard deviations above or below the
#' mean LRR. Once detected, it will try to extend the detected clusters until
#' the addition of a new SNP fails to comply with the required thresholds.
#'
#'
#' @usage findFocals(snps, sample.name, ignore.chr.Y=NULL, min.snps.in.cluster=5, max.extension=6,
#' fdel.sd.lim=4, fdel.min.num.below.thr=6, fdel.min.pct.below.thr=0.3, fdel.max.pct.above.mean=0.1,
#' fgain.sd.lim=3.5, fgain.min.num.above.thr=6, fgain.min.pct.above.thr=0.4, fgain.max.pct.below.mean=0.2,
#' cbs.min.segment.length=2e6, cbs.alpha=0.01, verbose=TRUE)
#'
#'
#' @param snps The SNPs data
#' @param sample.name The name of the sample
#' @param ignore.chr.Y (defaults to NULL)
#' @param min.snps.in.cluster (defaults to 5)
#' @param max.extension (defaults to 6)
#' @param fdel.sd.lim (defaults to 4)
#' @param fdel.min.num.below.thr (defaults to 6)
#' @param fdel.min.pct.below.thr (defaults to 0.3)
#' @param fdel.max.pct.above.mean (defaults to 0.1)
#' @param fgain.sd.lim (defaults to 3.5)
#' @param fgain.min.num.above.thr (defaults to 6)
#' @param fgain.min.pct.above.thr (defaults to 0.4)
#' @param fgain.max.pct.below.mean (defaults to 0.2)
#' @param cbs.min.segment.length (defaults to 2e6)
#' @param cbs.alpha (defaults to 0.01)
#' @param verbose (defaults to TRUE)
#'
#'
#' @return
#' A "FindFocalsResults" object with the results of the analysis. The object is
#' a list with the results for gains and deletions in different slots.
#'
#' @examples
#'
#'
#'
#' @export findFocals
#' @importFrom DNAcopy CNA
#' @importFrom DNAcopy smooth.CNA
#' @importFrom DNAcopy segment
#' @importFrom zoo rollapply
#' @importFrom stats median rnorm sd
#' @importFrom utils read.table
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @import GenomicRanges
#' @import regioneR
#'
#'
#'

findFocals <- function(snps, sample.name, ignore.chr.Y=NULL, min.snps.in.cluster=5, max.extension=6,
                       fdel.sd.lim=4, fdel.min.num.below.thr=6, fdel.min.pct.below.thr=0.3, fdel.max.pct.above.mean=0.1,
                       fgain.sd.lim=3.5, fgain.min.num.above.thr=6, fgain.min.pct.above.thr=0.4, fgain.max.pct.below.mean=0.2,
                       cbs.min.segment.length=2e6, cbs.alpha=0.01, verbose=TRUE) {

  if(verbose==TRUE) message("Segmenting LRR using CBS...")
  CNA.object <- DNAcopy::CNA(genomdat = snps$lrr, chrom = as.character(seqnames(snps)), maploc = seq_len(length(snps)), data.type="logratio", sampleid = sample.name)
  CNA.object <- DNAcopy::smooth.CNA(x=CNA.object, smooth.region=10, outlier.SD.scale=2.5, smooth.SD.scale=2, trim=0.1)
  lrr.segments <- DNAcopy::segment(CNA.object, alpha = cbs.alpha, nperm = 1000, p.method = "hybrid", min.width=5, kmax=25, nmin=200,
                          eta=0.05, trim = 0.025, undo.splits = "sdundo", undo.SD = 2, verbose=ifelse(verbose==TRUE, 1, 0))

  segments <- lrr.segments$output
  segments$chrom <- as.character(segments$chrom)

  segments <- cbind(segments[,-1], start=start(snps[segments$loc.start]), end=end(snps[segments$loc.end]), stringsAsFactors=FALSE)
  lrr.segs <- toGRanges(segments[,c("chrom", "start", "end", "seg.mean", "num.mark")])
  seqlevelsStyle(lrr.segs) <- "UCSC" #If snps object was created with loadSNPData, it should already be UCSC, but just in case ensure that it is

  #Filter out the segments smaller than a threshold
  lrr.segs <- lrr.segs[width(lrr.segs)>cbs.min.segment.length]


  if(verbose==TRUE) message("Computing the segmented LRR mean...")
  #Compute the lrr.mean for each point as:
  #   1 - if it overlaps a segment, the segment mean
  #   2 - If it does not, the mean of the 2 flanking segments

  mean.lrr <- rep(NA, length(snps))

  #Part 1: mean.lrr is seg.mean for the snps overlapping a segment
  for(i in c(1:length(lrr.segs))) {
    s <- lrr.segs[i]
    mean.lrr[overlapsAny(snps, s)] <- s$seg.mean
  }

  #Part 2: create a gaps granges and set mean.lrr to the mean of the two surrounding segments to it
  seg.gaps <- gaps(x=lrr.segs)
  for(i in c(1:length(seg.gaps))) {
    gg <- seg.gaps[i]
    lrr.mean <- mean(c(subsetByOverlaps(lrr.segs, extendRegions(gg, extend.start = 1))$seg.mean,
                       subsetByOverlaps(lrr.segs, extendRegions(gg, extend.end = 1))$seg.mean))
    mean.lrr[overlapsAny(snps, gg)] <- lrr.mean
  }

  #Some snps might still without value because they are "out of the chromosome".
  #Do it "by hand" finding the last segment of that chromosome
  chrs <- unique(as.character(seqnames(snps[which(is.na(mean.lrr))])))
  for(chr in chrs) {
    chr.segs <- lrr.segs[as.character(seqnames(lrr.segs))==chr]
    lrr.mean <- chr.segs[length(chr.segs)]$seg.mean
    mean.lrr[which(is.na(mean.lrr) & as.character(seqnames(snps))==chr)] <- lrr.mean
  }

  #and add the mean.lrr to the snps object
  snps$mean.lrr <- mean.lrr

  if(verbose==TRUE) message("Computing the global LRR standard deviation...")
  #Compute the global sd
  #get the median of the rolling sd's of lrr to get the "expected spread of lrr data"
  #NOTE: Leave the Y chromosome out, since it can affect the computation in female genomes
  in.chrY <- as.character(seqnames(snps))=="chrY"
  lrr.sd <- median(zoo::rollapply(data = snps[!in.chrY]$lrr, FUN=sd, width=301))

  #and compute the Y sd separately (to filter the chr out is necessary)
  y.lrr.sd <- median(zoo::rollapply(data = snps[in.chrY]$lrr, FUN=sd, width=301))

  if(verbose==TRUE) message("Detecting focal deletions...")
  #Now, detect the clusters below the expected mean

  #For speed, extract lrr and mean.lrr from snps only once
  snps.lrr <- snps$lrr
  snps.mean.lrr <- snps$mean.lrr



  #Get the points below the threshold (mean - a number of sd's)
  #and for every one of them, study if they form a cluster with most of the points very below the mean
  below.thr.idx <- which(snps$lrr < snps$mean.lrr - fdel.sd.lim * lrr.sd)

  #for each point below the index, try to expand to create a cluster.
  #The expansion is done one below.thr snps at a time, and iteratively changing
  #from left extension to right extension
  focal.del <- GRanges()
  in.focal.del <- c()
  for(i in c(1:length(below.thr.idx))) {
    #check if it has any "close" neighbours
    has.neighbours <- FALSE
    if(i>1) {
      if(below.thr.idx[i-1]>below.thr.idx[i]-max.extension) {
        has.neighbours <- TRUE
      }
    }
    if(i<length(below.thr.idx) & below.thr.idx[i+1] < below.thr.idx[i]+max.extension) {
      has.neighbours <- TRUE
    }
    if(has.neighbours) {
      #Try to "expand the region". Extend it so at least fdel.min.pct.below.thr part if below the threshold and we extend it a max of max.extension points not below the threshold
      k.left <- below.thr.idx[i]
      k.right <- below.thr.idx[i]
      done <- FALSE
      while(!done) {
        done <- TRUE
        #First try to the extend to the left
        for(m.left in c(1:max.extension)) {
          if((k.left  - m.left) > 1) {
            if(snps.lrr[k.left-m.left] < snps.mean.lrr[k.left  - m.left] - fdel.sd.lim * lrr.sd) {
              possible.cluster <- (k.left-m.left):k.right
              num.below.thr <- length(which(snps.lrr[possible.cluster] < snps.mean.lrr[possible.cluster]-fdel.sd.lim*lrr.sd))
              num.above.mean <- length(which(snps.lrr[possible.cluster] > snps.mean.lrr[possible.cluster]))
              if(num.below.thr > length(possible.cluster)*fdel.min.pct.below.thr &&
                 num.above.mean < length(possible.cluster)*fdel.max.pct.above.mean) {
                done <- FALSE
                k.left <- k.left - m.left
              }
            }
          }
        }
        #And then try to the extend to the right
        for(m.right in c(1:max.extension)) {
          if((k.right + m.right) < length(snps)) { #WARNING: - TODO: Check we do not cross the chromosome (or arm?) boundaries!
            if(snps.lrr[k.right+m.right] < snps.mean.lrr[k.right+m.right] - fdel.sd.lim*lrr.sd) {
              possible.cluster <- k.left:(k.right+m.right)
              num.below.thr <- length(which(snps.lrr[possible.cluster] < snps.mean.lrr[possible.cluster]-fdel.sd.lim*lrr.sd))
              num.above.mean <- length(which(snps.lrr[possible.cluster] > snps.mean.lrr[possible.cluster]))
              if(num.below.thr > length(possible.cluster)*fdel.min.pct.below.thr &&
                 num.above.mean < length(possible.cluster)*fdel.max.pct.above.mean) {
                done <- FALSE
                k.right <- k.right - m.right
              }
            }
          }
        }
      }
      if(k.right - k.left > min.snps.in.cluster) { #If there are at least 5 elements in the extended cluster
        num.below.thr <- length(which(snps.lrr[k.left:k.right] < snps.mean.lrr[k.left:k.right]-fdel.sd.lim*lrr.sd))
        if(num.below.thr >= fdel.min.num.below.thr) {
          new.focal.del <- snps[k.left]
          end(new.focal.del) <- end(snps[k.right])
          focal.del <- c(focal.del, new.focal.del)
          in.focal.del <- c(in.focal.del, k.left:k.right)
        }
      }
    }
  }

  #And the clusters above the expected mean
  #get the indexes of the snps above a threshold (a number of sd above the mean)
  above.thr.idx <- which(snps$lrr > snps$mean.lrr + fgain.sd.lim*lrr.sd)

  focal.gains <- GRanges()
  in.focal.gain <- c()
  for(i in c(1:length(above.thr.idx))) {
    #check if it has any "close" neighbours
    has.neighbours <- FALSE
    if(i>1) {
      if(above.thr.idx[i-1]>above.thr.idx[i]-max.extension) {
        has.neighbours <- TRUE
      }
    }
    if(i<length(above.thr.idx) & above.thr.idx[i+1] < above.thr.idx[i]+max.extension) {
      has.neighbours <- TRUE
    }
    if(has.neighbours) {
      #Try to "expand the region". Extend it so at least fgain.min.pct.above.thr part if below the threshold and we extend it a max of max.extension points not below the threshold
      k.left <- above.thr.idx[i]
      k.right <- above.thr.idx[i]
      done <- FALSE
      while(!done) {
        done <- TRUE
        #First try to the extend to the left
        for(m.left in c(1:max.extension)) {
          if((k.left  - m.left) > 1) {
            if(snps.lrr[k.left-m.left] > snps.mean.lrr[k.left  - m.left]+fgain.sd.lim*lrr.sd) {
              possible.cluster <- (k.left-m.left):k.right
              num.above.thr <- length(which(snps.lrr[possible.cluster] > snps.mean.lrr[possible.cluster]+fgain.sd.lim*lrr.sd))
              num.below.mean <- length(which(snps.lrr[possible.cluster] < snps.mean.lrr[possible.cluster]))
              if(num.above.thr > length(possible.cluster)*fgain.min.pct.above.thr &&
                 num.below.mean < length(possible.cluster)*fgain.max.pct.below.mean) {
                done <- FALSE
                k.left <- k.left - m.left
              }
            }
          }
        }
        #And then try to the extend to the right
        for(m.right in c(1:max.extension)) {
          if((k.right + m.right) < length(snps)) {
            if(snps.lrr[k.right+m.right] > snps.mean.lrr[k.right+m.right]+fgain.sd.lim*lrr.sd) {
              possible.cluster <- k.left:(k.right+m.right)
              num.above.thr <- length(which(snps.lrr[possible.cluster] > snps.mean.lrr[possible.cluster]+fgain.sd.lim*lrr.sd))
              num.below.mean <- length(which(snps.lrr[possible.cluster] < snps.mean.lrr[possible.cluster]))
              if(num.above.thr > length(possible.cluster)*fgain.min.pct.above.thr &&
                 num.below.mean < length(possible.cluster)*fgain.max.pct.below.mean) {
                done <- FALSE
                k.right <- k.right - m.right
              }
            }
          }
        }
      }
      if(k.right - k.left > min.snps.in.cluster) { #If at least 5 elements in the extended cluster
        num.above.thr <- length(which(snps.lrr[k.left:k.right] > snps.mean.lrr[k.left:k.right]+fgain.sd.lim*lrr.sd))
        if(num.above.thr >= fgain.min.num.above.thr) {
          new.focal.gain <- snps[k.left]
          end(new.focal.gain) <- end(snps[k.right])
          focal.gains <- c(focal.gains, new.focal.gain)
          in.focal.gain <- c(in.focal.gain, k.left:k.right)
        }
      }
    }
  }


  #If the Y chromsome sd is higher than the rest, it usually means it's a female genome and Y data is only noise. Ignore it.
  if(is.null(ignore.chr.Y)) ignore.chr.Y <- ifelse(y.lrr.sd > 2*lrr.sd, TRUE, FALSE)
  if(ignore.chr.Y==TRUE ) {
    if(verbose==TRUE) message("Removing focal gains and deletions in chromosome Y...")
    focal.del <- focal.del[as.character(seqnames(focal.del))!="chrY"]
    focal.gains <- focal.gains[as.character(seqnames(focal.gains))!="chrY"]
    in.focal.del <- in.focal.del[as.character(seqnames(snps[in.focal.del]))!="chrY"]
    in.focal.gain <- in.focal.gain[as.character(seqnames(snps[in.focal.gain]))!="chrY"]
  }

  ffr <- list(snps=snps, gains=focal.gains, deletions=focal.del,
              snps.in.gains=in.focal.gain, snps.in.deletions=in.focal.del,
              lrr.sd=lrr.sd, y.lrr.sd=y.lrr.sd,
              params=list(fdel.sd.lim=fdel.sd.lim, fgain.sd.lim=fgain.sd.lim))
  class(ffr) <- "FindFocalsResults"
  return(ffr)
}
