# workflow_oncoscan.R Functions to run the complete workflow from input files
# to scores and arm-level alterations.  Author: Yann Christinat Date:
# 23.06.2020

#' Run the standard workflow for Oncoscan ChAS files.
#'
#' @details Identifies the globally altered arms (\>=90\% of arm altered), computes the HRD and
#' TD+ scores. The amplification is defined as a CN subtype \code{cntype.weakamp} or
#' \code{cntype.strongamp}. An arm is gained if of CN type \code{cntype.gain} unless the arm is
#' amplified.
#'
#' @param chas.fn Path to the text-export ChAS file
#' @param gender Gender of the sample (M or F)
#'
#' @return A list of lists with the following elements:
#' \code{armlevel = list(AMP= list of arms, GAIN= list of arms, LOSS= list of arms, LOH= list of arms),
#' scores = list(LST= number, LOH= number, TDplus= number, TD= number),
#' gender = gender as given by the parameter,
#' file = path of the ChAS file as given by the parameter)}
#'
#' @export
#'
#' @import magrittr
#'
#' @examples
#' segs.filename <- system.file('extdata', 'chas_example.txt', package = 'oncoscanR')
#' workflow_oncoscan.run(segs.filename, 'M')
workflow_oncoscan.run <- function(chas.fn, gender) {
    if (!(gender %in% c("M", "F"))) {
        stop("The gender (second argument) has to be F or M.")
    }

    # Remove the 21p arm from the Oncoscan coverage as it is only partly
    # covered and we don't want to return results on this arm.
    oncoscan.cov <- oncoscanR::oncoscan_na33.cov[seqnames(oncoscanR::oncoscan_na33.cov) !=
        "21p"]

    # Load the ChAS file and assign subtypes.
    segments <- load_chas(chas.fn, oncoscan.cov)
    segments$cn.subtype <- get_cn_subtype(segments, gender)

    # Clean the segments: resctricted to Oncoscan coverage, LOH not overlapping
    # with copy loss segments, smooth&merge segments within 300kb and prune
    # segments smaller than 300kb.
    segs.clean <- trim_to_coverage(segments, oncoscan.cov) %>%
        adjust_loh() %>%
        merge_segments() %>%
        prune_by_size()

    # Split segments by type: Loss, LOH, gain or amplification and get the
    # arm-level alterations.  Note that the segments with copy gains include
    # all amplified segments.
    armlevel.loss <- segs.clean[segs.clean$cn.type == oncoscanR::cntype.loss] %>%
        armlevel_alt(kit.coverage = oncoscan.cov)
    armlevel.loh <- segs.clean[segs.clean$cn.type == oncoscanR::cntype.loh] %>%
        armlevel_alt(kit.coverage = oncoscan.cov)
    armlevel.gain <- segs.clean[segs.clean$cn.type == oncoscanR::cntype.gain] %>%
        armlevel_alt(kit.coverage = oncoscan.cov)
    armlevel.amp <- segs.clean[segs.clean$cn.subtype %in% c(oncoscanR::cntype.strongamp, oncoscanR::cntype.weakamp)] %>%
        armlevel_alt(kit.coverage = oncoscan.cov)

    # Remove amplified segments from armlevel.gain
    armlevel.gain <- armlevel.gain[!(names(armlevel.gain) %in% names(armlevel.amp))]

    # Get the number of nLST and TDplus
    wgd <- score_estwgd(segs.clean, oncoscan.cov)  # Get the avg CN, including 21p
    hrd <- score_nlst(segs.clean, wgd["WGD"], oncoscan.cov)

    n.td <- score_td(segs.clean)

    mbalt <- score_mbalt(segs.clean, oncoscan.cov, loh.rm=TRUE)

    hrd.label <- hrd["HRD"]
    if(mbalt['sample']/mbalt['kit'] < 0.01) hrd.label <- paste(hrd["HRD"], "(no tumor?)")

    # Get the alterations into a single list and print it in a JSON format.
    armlevel_alt.list <- list(AMP = sort(names(armlevel.amp)), LOSS = sort(names(armlevel.loss)),
        LOH = sort(names(armlevel.loh)), GAIN = sort(names(armlevel.gain)))
    scores.list <- list(HRD = paste0(hrd.label, ", nLST=", hrd["nLST"]), TDplus = n.td$TDplus,
        avgCN = substr(as.character(wgd["avgCN"]), 1, 4))

    return(list(armlevel = armlevel_alt.list, scores = scores.list, gender = gender,
        file = basename(chas.fn)))
}
