# objects_doc.R Documentation of all R objects included in the package Author:
# Yann Christinat Date: 23.9.2019 The following commands were used to integrate
# the data into the package: use_data(oncoscan_na33.cov, internal=FALSE)
# use_data(segs.chas_example, internal=FALSE) use_data(cntype.loh, cntype.gain,
# cntype.loss, cntype.hetloss, cntype.homloss, cntype.amp, cntype.weakamp,
# cntype.strongamp)

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

# Set the objects as global variables otherwise R CMD CHECK complains
#utils::globalVariables(c("cntype.gain", "cntype.amp", "cntype.strongamp",
#    "cntype.weakamp", "cntype.loh", "cntype.loss", "cntype.homloss",
#    "cntype.hetloss", "oncoscan_na33.cov", "segs.chas_example"))
