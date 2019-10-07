# scores.R
#
# Functions to identify the globally-altered arms and computes various scores.
#
# Author: Yann Christinat
# Date: 19.9.2019


#' Get all globally-altered chromosome arms.
#'
#' @details By default uses the sum of all alterations and set the arm as globally altered if \>80\% of the arm
#' is altered. Does not account for alteration type and copy number.
#' Will run the function \code{trim_to_coverage} on the segments.
#'
#' @param segments A \code{GRanges} object containing the segments.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on each chromosome arm.
#' @param threshold The minimum percentage of the arm to be considered as globally altered. Defaults
#' to 80\%.
#'
#' @return A list of globally-altered chromosome arms with the percentage of arm altered.
#' @export
#'
#' @import GenomicRanges
#'
#' @examples
#' arms <- armlevel_alt(segs.chas_example, oncoscan_na33.cov, 0.9)
armlevel_alt <- function(segments, kit.coverage, threshold = 0.8){
  is.cn_segment(segments, raise_error = TRUE)

  vals <- sapply(levels(seqnames(kit.coverage)), function(arm){
    if(!(arm %in% seqnames(segments))){
      return(0)
    }
    # Get the number of bases covered by the kit
    armlen <- sum(width(kit.coverage[seqnames(kit.coverage) == arm]))

    # Run trim_to_coverage on arm segments
    segs <- trim_to_coverage(segments[seqnames(segments) == arm], kit.coverage)
    if(length(segs)==0){
      return(0)
    }


    # Set the copy number to NA to make sure overlapping segments will be merged, then merge
    segs$cn <- NA
    segs.flat <- merge_segments(segs, 1/1000)

    # Get the number of bases altered
    segslen <- sum(width(segs.flat))

    return(segslen/armlen)
  })

  return(vals[vals >= threshold])
}


#' Compute the number of Large-scale State Transitions (LSTs).
#'
#' @details Procedure based on the paper from Popova et al, Can. Res. 2012 (PMID: 22933060). First segments
#' smaller than 3Mb are removed, then segments are smoothed with respect to copy number at a distance of 3Mb.
#' The number of LSTs is the number of breakpoints (breakpoints closer than 3Mb are merged) that have a segment
#' larger or equal to 10Mb on each side. This score was linked to BRCA1/2-deficient tumors.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy number and copy number types.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on each chromosome arm.
#'
#' @return An integer representing the number of LSTs.
#' @export
#'
#' @examples
#' score_lst(segs.chas_example, oncoscan_na33.cov)
score_lst <- function(segments, kit.coverage){
  #Prune all segments < 3Mb
  segs.min3mb <- prune_by_size(segments, 3000)
  #Merge segments by CN at a distance of 3Mb
  segs.merged <- merge_segments(segs.min3mb, 3000)

  #Compute the number of Large-scale State Transitions for each arm
  vals <- lapply(unique(seqnames(segs.merged)), function(arm){
    arm.cov <- kit.coverage[seqnames(kit.coverage) == arm]
    segs <- segs.merged[seqnames(segs.merged) == arm]

    # Get breakpoints and start/end of arm
    starts <- start(segs)
    ends <- end(segs)
    breakpoints <- sort(unique(c(start(arm.cov), starts, ends+1, end(arm.cov)+1)))

    if(length(breakpoints) < 3){
      # Since we add the start/end of arm, there has to be at least 3 elements in 'breakpoints'
      # to have LSTs
      return(0)
    }

    # Smooth breakpoints at 3Mb
    for(i in 2:length(breakpoints)){
      if(breakpoints[i]-breakpoints[i-1] < 3*10^6){ # Merge breakpoints
        midpoint <- round((breakpoints[i]+breakpoints[i-1])/2)
        breakpoints[i-1] <- NA
        breakpoints[i] <- midpoint
      }
    }
    bp.smooth <- breakpoints[!is.na(breakpoints)]

    if(length(bp.smooth) < 3){
      # Since we add the start/end of arm, there has to be at least 3 elements in 'breakpoints'
      # to have LSTs
      return(0)
    }

    # Get the number of LSTs
    n.lst <- 0
    for(i in 2:(length(bp.smooth)-1)){
      n.lst <- ifelse(bp.smooth[i]-bp.smooth[i-1] >= 10*10^6 &
                        bp.smooth[i+1]-bp.smooth[i] >= 10*10^6,
                      n.lst+1, n.lst)
    }
    return(n.lst)
  })
  return(sum(unlist(vals)))
}


#' Compute the number HR deficiency-associated LOH regions.
#'
#' @details Procedure based on the paper from Abkevich et al., Br J Cancer 2012 (PMID: 23047548). All LOH
#' segments larger than 15Mb but excluding chromosome with a global LOH alteration (to compute with the
#' \code{armlevel_alt} function on LOH segments only). This score was linked to BRCA1/2-deficient tumors.
#' Note that the function will merge overlapping or neighbor LOH segments (at a distance of 1bp).
#'
#' @param segments A \code{GRanges} object containing the segments, their copy number and copy number types.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on each chromosome arm.
#' @param armlevel.loh A list of arms with global/arm-level LOH alteration.
#'
#' @return An integer representing the number of HRD-LOH regions.
#' @export
#'
#' @examples
#' armlevel.loh <- armlevel_alt(segs.chas_example[segs.chas_example$cn.type == cntype.loh],
#'                              kit.coverage = oncoscan_na33.cov)
#' score_loh(segs.chas_example, oncoscan_na33.cov, names(armlevel.loh))
score_loh <- function(segments, kit.coverage, armlevel.loh){
  is.cn_segment(segments, raise_error = TRUE)

  if(length(armlevel.loh) > 0){
    # Check on each chromosome if both arms are within the parameter 'armlevel.loh'. Also account
    # for the fact that a chromosome may have only one arm covered by the kit.
    chromnames <- unique(unlist(strsplit(as.character(seqnames(kit.coverage)), '[pq]')))
    arms.to_ban <- lapply(chromnames, function(chrom){
      arms <- paste0(chrom, c('p','q'))
      arms.covered <- intersect(arms, seqnames(kit.coverage))
      arms.found <- unique(intersect(arms.covered, armlevel.loh))
      if(length(arms.found) == length(arms.covered)){
        return(as.character(arms.covered))
      }

    })
    arms.to_ban <- do.call("c", arms.to_ban)
  }
  else {
    arms.to_ban <- c()
  }
  # Get LOH segments larger than 15Mb and not in 'arms.to_ban'
  segs.loh <- segments[segments$cn.subtype %in% c(cntype.loh, cntype.hetloss) &
                             !(seqnames(segments) %in% arms.to_ban)]
  segs.merged <- reduce(segs.loh)
  segs.loh_15m <- segs.loh[width(segs.merged) > 15*10^6]
  return(length(segs.loh_15m))
}


#' Compute the number of large tandem duplication (TDplus).
#'
#' @details Procedure based on the paper from Popova et al., Cancer Res 2016 (PMID: 26787835). The TDplus
#' score is defined as the number of regions larger than 1Mb but smaller or equal to 10Mb with a gain of one
#' or two copies (\code{cntype.gain} in the field \code{cn.subtype}). This score was linked to CDK12-deficient
#' tumors. They also identified as second category of tandem duplication whose size is smaller or equal than
#' 1Mb and around 300Kb but could not link it to a phenotype. Note that due to its resolution the Oncoscan assay
#' will most likely miss this second category. Nonetheless it is reported by the function.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy number and copy number types.
#'
#' @return A list of integer containing the TDplus score (\code{"TDplus"}) and the small TD score (\code{"TD"}).
#' @export
#'
#' @examples
#' score_td(segs.chas_example)
score_td <- function(segments){
  is.cn_segment(segments, raise_error = TRUE)

  segs.gain <- segments[segments$cn.subtype == cntype.gain]
  segs.width <- width(segs.gain)
  segs.tdplus <- segs.gain[segs.width > 1*10^6 & segs.width <= 10*10^6]
  segs.td <- segs.gain[segs.width <= 1*10^6]

  return(list(TDplus=length(segs.tdplus), TD=length(segs.td)))
}


