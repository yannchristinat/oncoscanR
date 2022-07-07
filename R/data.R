# objects_doc.R Documentation of all R objects included in the package
# Author: Yann Christinat
# Date: 30.6.2022


#' GenomicRanges object of the chromosomal arms coverage for the oncoscan assay
#' (based on file extdata/Oncoscan.na33.r2.cov.processed.bed).
#'
#' @source \code{oncoscan_na33.cov <- get_oncoscan_coverage_from_bed(
#'     system.file('extdata', 'Oncoscan.na33.r2.cov.processed.bed',
#'                 package = 'oncoscanR'))}
#'
#' @format A \code{GRanges} object containing the region covered on each
#' chromosome arm.
"oncoscan_na33.cov"


#' Expected segments from loading the ChAS file 'chas_example.txt'.
#'
#' @source
#' segs.filename <- system.file('extdata', 'chas_example.txt', 
#' package = 'oncoscanR')
#' mykit.cov <- get_oncoscan_coverage_from_probes()
#' segs.chas_example <- load_chas(segs.filename, kit.coverage = mykit.cov)
#'
#' @format A \code{GRanges} object containing the segments, their copy number
#' (field \code{cn}) and their copy number types (field \code{cn.type}).
"segs.chas_example"

#' Accepted types of CN for the segments
#'   - 'Gain': 1-2 extra copies
#'   - 'Weak amplification': 3-7 extra copies
#'   - 'Strong amplification': 8 or more extra copies
#'   - 'Heterozygote loss': Loss of one copy out of two
#'   - 'Homozygote loss': Loss of all copies
#'   - 'LOH': copy-neutral loss of one parental allele
#' @source
#' cntypes <- list(LOH='LOH', Gain='Gain', Loss='Loss')
#'
#' @noRd
"cntypes"

