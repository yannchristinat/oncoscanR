#!/usr/bin/Rscript
# oncoscan-workflow.R
#
# Standard workflow for the detection of arm-level alterations.
#
# Author: Yann Christinat
# Date: 1.10.2019

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(oncoscanR))

#' Retrieve arguments from command line.
#' Expects that the first argument is the ChAS file and the second the gender (F or M).
args <- commandArgs(TRUE)
if(length(args) != 2){
  stop("The first argument has to be the ChAS file and the second the gender (F or M).")
}

chas.fn <- args[1]
gender <- args[2]
if(!(gender %in% c('M', 'F'))){
  stop("The gender (second argument) has to be F or M.")
}

#' Remove the 21p arm from the Oncoscan coverage as it is only partly covered and we don't
#' want to return results on this arm.
oncoscan.cov <- oncoscan_na33.cov[seqnames(oncoscan_na33.cov) != '21p']

#' Load the ChAS file and assign subtypes.
segments <- load_chas(chas.fn, oncoscan.cov)
segments$cn.subtype <- get_cn_subtype(segments, gender)

#' Clean the segments: resctricted to Oncoscan coverage, LOH not overlapping with copy loss
#' segments, smooth&merge segments within 300kb and prune segments smaller than 300kb.
segs.clean <- trim_to_coverage(segments, oncoscan.cov) %>%
  adjust_loh() %>%
  merge_segments() %>%
  prune_by_size()

#' Split segments by type: Loss, LOH, gain or amplification and get the arm-level alterations.
#' Note that the segments with copy gains include all amplified segments.
armlevel.loss <- segs.clean[segs.clean$cn.type == cntype.loss] %>%
  armlevel_alt(kit.coverage = oncoscan.cov)
armlevel.loh <- segs.clean[segs.clean$cn.type == cntype.loh] %>%
  armlevel_alt(kit.coverage = oncoscan.cov)
armlevel.gain <- segs.clean[segs.clean$cn.type == cntype.gain] %>%
  armlevel_alt(kit.coverage = oncoscan.cov)
armlevel.amp <- segs.clean[segs.clean$cn.subtype %in% c(cntype.strongamp, cntype.weakamp)] %>%
  armlevel_alt(kit.coverage = oncoscan.cov)

#' Remove amplified segments from armlevel.gain
armlevel.gain <- armlevel.gain[!(names(armlevel.gain) %in% names(armlevel.amp))]

#' Get the number of LST, LOH, TDplus and TD
n.lst <- score_lst(segs.clean, oncoscan.cov)
n.loh <- score_loh(segs.clean, oncoscan.cov, names(armlevel.loh))
n.td <- score_td(segs.clean)

#' Get the alterations into a single list and print it in a JSON format.
armlevel_alt.list <- list(AMP=names(armlevel.amp),
                          LOSS=names(armlevel.loss),
                          LOH=names(armlevel.loh),
                          GAIN=names(armlevel.gain))
scores.list <- c(list(LST=n.lst, LOH=n.loh), n.td)
print(toJSON(list(armlevel=armlevel_alt.list,
                  scores=scores.list,
                  gender=gender,
                  file=chas.fn),
             auto_unbox=TRUE, pretty=TRUE))

