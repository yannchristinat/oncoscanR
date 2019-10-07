# objects_doc.R
#
# Documentation of all R objects included in the package
#
# Author: Yann Christinat
# Date: 23.9.2019
#
# The following commands were used to integrate the data into the package:
# use_data(oncoscan_na33.cov, internal=FALSE)
# use_data(segs.chas_example, internal=FALSE)
# use_data(cntype.loh, cntype.gain, cntype.loss, cntype.hetloss, cntype.homloss, cntype.amp, cntype.weakamp,
#          cntype.strongamp)

#' Arm coverage of the Oncoscan na33.r1 assay
#'
#' @source Processing of the file OncoScan.na33.r1.annot.csv (obtained on Affymetrix website) by the function
#' \code{oncoscan_na33.cov <- get_oncoscan_coverage_from_probes()}.
#'
#' @format A \code{GRanges} object containing the regions covered on each chromosome arm. Arms not covered by
#' the kit are absent from the object.
#'
#' @examples
#' \dontrun{
#'  oncoscan_na33.cov
#' }
"oncoscan_na33.cov"


#' copy number type "LOH" (copy neutral loss of heterozygosity).
#'
#' @source \code{cntype.loh <- "LOH"}
#' @format character
"cntype.loh"


#' copy number type "Gain" (one or more additional copies).
#' In some contextes can be defined as the gain of 1-2 extra copies but not more.
#'
#' @source \code{cntype.gain <- "Gain"}
#' @format character
"cntype.gain"


#' copy number type "Loss" (loss of one or more copies; depends on gender for the sexual chromosomes).
#'
#' @source \code{cntype.loss <- "Loss"}
#' @format character
"cntype.loss"


#' copy number type "Heterozygous loss" (loss of one of the two copies).
#'
#' @source \code{cntype.hetloss <- "Heterozygous loss"}
#' @format character
"cntype.hetloss"


#' copy number type "Homozygous loss" (loss of the two copies or only one on sexual chromosomes if the
#' subject is male).
#'
#' @source \code{cntype.homloss <- "Homozygous loss"}
#' @format character
"cntype.homloss"


#' copy number type "Amplification" (3 additional copies or more).
#'
#' @source \code{cntype.amp <- "Amplification"}
#' @format character
"cntype.amp"


#' copy number type "Weak amplification" (3 to 7 additional copies).
#'
#' @source \code{cntype.weakamp <- "Weak amplification"}
#' @format character
"cntype.weakamp"


#' copy number type "Strong amplification" (8 additional copies or more).
#'
#' @source \code{cntype.strongamp <- "Strong amplification"}
#' @format character
"cntype.strongamp"


#' Expected segments from loading the ChAS file 'chas_example.txt'.
#'
#' @source
#' segs.filename <- system.file("extdata", "chas_example.txt", package = "oncoscanR")
#' mykit.cov <- get_oncoscan_coverage_from_probes()
#' segs.chas_example <- load_chas(segs.filename, kit.coverage = mykit.cov)
#'
#' @format A \code{GRanges} object containing the segments, their copy number (field \code{cn}) and
#' their copy number types (field \code{cn.type}).
"segs.chas_example"

