# data.R Documentation of all R objects included in the package
# Author: Yann Christinat
# Date: 14.06.2022

#' Arm coverage of the Oncoscan na33.r1 assay
#'
#' @source Processing of the file OncoScan.na33.r1.annot.csv (obtained from
#' https://www.affymetrix.com/analysis/downloads/na33/genotyping/) by the
#' function \code{oncoscan_na33.cov <- get_oncoscan_coverage_from_probes()}.
#'
#' @format A \code{GRanges} object containing the regions covered on each
#' chromosome arm. Arms not covered by the kit are absent from the object.
"oncoscan_na33.cov"


#' Expected segments from loading the ChAS file 'chas_example.txt'.
#'
#' @source
#' segs.filename <- system.file('extdata', 'chas_example.txt', package = 'oncoscanR')
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
#'   @source
#'   cntypes <- list(LOH='LOH', Gain='Gain', Loss='Loss')
#'
#'   @noRd
"cntypes"

