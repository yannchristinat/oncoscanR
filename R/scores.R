# scores.R Functions to identify the globally-altered arms and computes various
# scores.  Author: Yann Christinat Date: 19.9.2019


#' Get all globally-altered chromosome arms.
#'
#' @details By default uses the sum of all alterations and set the arm as
#' globally altered if \>80\% of the arm is altered. Does not account for
#' alteration type and copy number.
#' Will run the function \code{trim_to_coverage} on the segments.
#'
#' @param segments A \code{GRanges} object containing the segments.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm.
#' @param threshold The minimum percentage of the arm to be considered as
#' globally altered. Defaults to 80\%.
#'
#' @return A list of globally-altered chromosome arms with the percentage of arm
#'  altered.
#' @export
#'
#' @import GenomicRanges
#'
#' @examples
#' arms <- armlevel_alt(segs.chas_example, oncoscan_na33.cov, 0.9)
armlevel_alt <- function(segments, kit.coverage, threshold = 0.9) {
    is_cn_segment(segments, raise_error = TRUE)

    vals <- vapply(levels(seqnames(kit.coverage)), function(arm) {
        if (!(arm %in% as.vector(seqnames(segments)))) {
            return(0)
        }
        # Get the number of bases covered by the kit
        armlen <- sum(width(kit.coverage[as.vector(seqnames(kit.coverage)) == arm]))

        # Run trim_to_coverage on arm segments
        segs <- trim_to_coverage(segments[as.vector(seqnames(segments)) == arm],
            kit.coverage)
        if (length(segs) == 0) {
            return(0)
        }


        # Set the copy number to NA to make sure overlapping segments will be
        # merged, then merge
        segs$cn <- NA
        segs.flat <- merge_segments(segs, 1/1000)

        # Get the number of bases altered
        segslen <- sum(width(segs.flat))

        return(segslen/armlen)
    }, c(0), USE.NAMES = TRUE)
    return(vals[vals >= threshold])
}


#' Compute the number of Large-scale State Transitions (LSTs).
#'
#' @details Procedure based on the paper from Popova et al, Can. Res. 2012
#' (PMID: 22933060). First segments smaller than 3Mb are removed, then segments
#' are smoothed with respect to copy number at a distance of 3Mb.
#' The number of LSTs is the number of breakpoints (breakpoints closer than 3Mb
#' are merged) that have a segment larger or equal to 10Mb on each side. This
#' score was linked to BRCA1/2-deficient tumors.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm.
#'
#' @return An integer representing the number of LSTs.
#' @export
#'
#' @examples
#' score_lst(segs.chas_example, oncoscan_na33.cov)
score_lst <- function(segments, kit.coverage) {
    # Prune all segments < 3Mb
    segs.min3mb <- prune_by_size(segments, 3000)
    if (length(segs.min3mb) == 0) {
        return(0)
    }

    # Merge segments by CN at a distance of 3Mb
    segs.merged <- merge_segments(segs.min3mb, 3000)

    # Compute the number of Large-scale State Transitions for each arm
    vals <- lapply(unique(seqnames(segs.merged)), function(arm) {
        arm.cov <- kit.coverage[as.vector(seqnames(kit.coverage)) == arm]
        segs <- segs.merged[as.vector(seqnames(segs.merged)) == arm]

        # Get breakpoints and start/end of arm
        starts <- start(segs)
        ends <- end(segs)
        breakpoints <- sort(unique(c(start(arm.cov), starts, ends + 1, end(arm.cov) +
            1)))

        if (length(breakpoints) < 3) {
            # Since we add the start/end of arm, there has to be at least 3
            # elements in 'breakpoints' to have LSTs
            return(0)
        }

        # Smooth breakpoints at 3Mb
        for (i in 2:length(breakpoints)) {
            if (breakpoints[i] - breakpoints[i - 1] < 3 * 10^6) {
                # Merge breakpoints
                midpoint <- round((breakpoints[i] + breakpoints[i - 1])/2)
                breakpoints[i - 1] <- NA
                breakpoints[i] <- midpoint
            }
        }
        bp.smooth <- breakpoints[!is.na(breakpoints)]

        if (length(bp.smooth) < 3) {
            # Since we add the start/end of arm, there has to be at least 3
            # elements in 'breakpoints' to have LSTs
            return(0)
        }

        # Get the number of LSTs
        n.lst <- 0
        for (i in 2:(length(bp.smooth) - 1)) {
            n.lst <- ifelse(bp.smooth[i] - bp.smooth[i - 1] >= 10 * 10^6 & bp.smooth[i +
                1] - bp.smooth[i] >= 10 * 10^6, n.lst + 1, n.lst)
        }
        return(n.lst)
    })
    return(sum(unlist(vals)))
}


#' Compute the number HR deficiency-associated LOH regions.
#'
#' @details Procedure based on the paper from Abkevich et al., Br J Cancer 2012
#' (PMID: 23047548). All LOH segments larger than 15Mb but excluding chromosome
#' with a global LOH alteration (to compute with the \code{armlevel_alt}
#' function on LOH segments only). This score was linked to BRCA1/2-deficient
#' tumors.
#' Note that the function will merge overlapping or neighbor LOH segments (at a
#' distance of 1bp).
#'
#' Of note the cn.subtype has to be defined first via the function
#' \code{get_cn_subtype}.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm.
#' @param arms.loh A list of arms with global/arm-level LOH alteration.
#' @param arms.hetloss A list of arms with global/arm-level heterozygous
#' losses.
#'
#' @return An integer representing the number of HRD-LOH regions.
#' @export
#'
#' @examples
#' segs.chas_example$cn.subtype <- get_cn_subtype(segs.chas_example, 'F')
#' armlevel.loh <- armlevel_alt(segs.chas_example[segs.chas_example$cn.type == cntype.loh],
#'                              kit.coverage = oncoscan_na33.cov)
#' armlevel.hetloss <- armlevel_alt(segs.chas_example[segs.chas_example$cn.subtype == cntype.hetloss],
#'                              kit.coverage = oncoscan_na33.cov)
#' score_loh(segs.chas_example, oncoscan_na33.cov, names(armlevel.loh), names(armlevel.hetloss))
score_loh <- function(segments, kit.coverage, arms.loh, arms.hetloss) {
    is_cn_segment(segments, raise_error = TRUE)

    if (length(segments) == 0) {
        return(0)
    }

    if (is.null(segments$cn.subtype)) {
        stop("Segments are missing the field 'cn.subtype'.")
    }

    # Test if any arms in arms.loh or arms.loss is not in kit.coverage
    if(length(c(arms.loh, arms.hetloss)) > 0){
        unknown.arms <- setdiff(c(arms.loh, arms.hetloss), seqnames(kit.coverage))
        if (length(unknown.arms)>0) {
            msg <- paste('Some arms were not found in the kit:',
                         paste(as.character(unknown.arms), collapse = TRUE))
            stop(msg)
        }
    }

    arms.to_ban <- c()
    # Get chromosomes with LOH on both arms
    if (length(arms.loh) > 0) {
        # Check on each chromosome if both arms are within the parameter
        # 'arms.loh'. Also account for the fact that a chromosome may have
        # only one arm covered by the kit.
        chromnames <- unique(unlist(strsplit(as.character(seqnames(kit.coverage)),
            "[pq]")))
        arms_found.to_ban <- lapply(chromnames, function(chrom) {
            arms <- paste0(chrom, c("p", "q"))

            arms.covered <- intersect(arms, as.vector(seqnames(kit.coverage)))

            arms.found <- unique(intersect(arms.covered, arms.loh))
            if (length(arms.found) == length(arms.covered)) {
                return(as.character(arms.covered))
            }
        })
        arms.to_ban <- append(arms.to_ban, do.call("c", arms_found.to_ban))
    }
    # Get chromosomes with het loss on both arms
    if (length(arms.hetloss) > 0) {
        # Check on each chromosome if both arms are within the parameter
        # 'arms.loh'. Also account for the fact that a chromosome may have
        # only one arm covered by the kit.
        chromnames <- unique(unlist(strsplit(as.character(seqnames(kit.coverage)),
            "[pq]")))
        arms_found.to_ban <- lapply(chromnames, function(chrom) {
            arms <- paste0(chrom, c("p", "q"))

            arms.covered <- intersect(arms, as.vector(seqnames(kit.coverage)))

            arms.found <- unique(intersect(arms.covered, arms.hetloss))
            if (length(arms.found) == length(arms.covered)) {
                return(as.character(arms.covered))
            }
        })
        arms.to_ban <- append(arms.to_ban, do.call("c", arms_found.to_ban))
    }
    # Get LOH segments larger than 15Mb and not in 'arms.to_ban'
    segs.loh <- segments[segments$cn.subtype %in% c(oncoscanR::cntype.loh, oncoscanR::cntype.hetloss) &
        !(as.vector(seqnames(segments)) %in% arms.to_ban)]
    segs.merged <- merge_segments(segs.loh, 1/1000)
    segs.loh_15m <- segs.merged[width(segs.merged) > 15 * 10^6]
    return(length(segs.loh_15m))
}


#' Compute the number of large tandem duplication (TDplus).
#'
#' @details Procedure based on the paper from Popova et al., Cancer Res 2016
#' (PMID: 26787835). The TDplus score is defined as the number of regions larger
#'  than 1Mb but smaller or equal to 10Mb with a gain of one or two copies
#'  (\code{cntype.gain} in the field \code{cn.subtype}). This score was linked
#'  to CDK12-deficient tumors. They also identified as second category of tandem
#'   duplication whose size is smaller or equal than 1Mb and around 300Kb but
#'   could not link it to a phenotype. Note that due to its resolution the
#'   Oncoscan assay will most likely miss this second category. Nonetheless it
#'   is reported by the function.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#'
#' @return A list of integer containing the TDplus score (\code{'TDplus'}) and
#' the small TD score (\code{'TD'}).
#' @export
#'
#' @examples
#' score_td(segs.chas_example)
score_td <- function(segments) {
    is_cn_segment(segments, raise_error = TRUE)

    if (length(segments) == 0) {
        return(list(TDplus = 0, TD = 0))
    }

    if (is.null(segments$cn.type)) {
        stop("Segments are missing the field 'cn.subtype'.")
    }

    segs.gain <- segments[segments$cn.subtype == oncoscanR::cntype.gain]
    segs.width <- width(segs.gain)
    segs.tdplus <- segs.gain[segs.width > 1 * 10^6 & segs.width <= 10 * 10^6]
    segs.td <- segs.gain[segs.width <= 1 * 10^6]

    return(list(TDplus = length(segs.tdplus), TD = length(segs.td)))
}



#' Compute the average copy number variation across the genome.
#'
#' @details Compute the weighted average (by segment length) of the copy number
#' variation. LOH segments and sexual chromosomes are excluded. Copy number
#' variation is rounded to the next level (1.67 -> 1 but 2.33 -> 3).
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm.
#'
#' @return A decimal value
#' @export
#'
#' @examples
#' score_avgcn(segs.chas_example, oncoscan_na33.cov)
score_avgcn <- function(segments, kit.coverage) {
    autosomes <- as.vector(seqnames(kit.coverage)[!(as.vector(seqnames(kit.coverage)) %in%
        c("Xp", "Xq", "Yp", "Yq"))])

    # Total number of bases (in Mbp) covered by the kit
    kitwidth.noXY <- sum(width(kit.coverage[as.vector(seqnames(kit.coverage)) %in%
        autosomes])/10^6)

    # Select segments on autosomes and with non-copy neutral alterations
    segs <- segments[as.vector(seqnames(segments)) %in% autosomes & segments$cn.type %in%
        c(oncoscanR::cntype.gain, oncoscanR::cntype.loss)]
    if (length(segs) == 0) {
        return(2)
    }

    segs.w <- width(segs)/10^6  # Width of each segment in Mbp

    # Round the copy number in case of subclones
    cn <- segs$cn
    cn[cn < 2] <- floor(cn[cn < 2])
    cn[cn > 2] <- ceiling(cn[cn > 2])

    avgcn <- (sum(segs.w * cn) + 2 * (kitwidth.noXY - sum(segs.w)))/kitwidth.noXY
    return(avgcn)
}


#' Estimates the number of whole-genome doubling events (WGD).
#'
#' @details Based on the publication from Carter et al. (Nature Biotechnology
#' 2012; PubMed ID: 22544022).
#' On a pan-cancer cohort, they observed that tumors that underwent one
#' whole-genome doubling event had a ploidy (average copy number) between 2.2
#' and 3.4. This function relies on the function \code{score_avgcn} to compute
#' the ploidy.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm.
#'
#' @return A named list with two values: WGD (whole-genome doubling events) and
#' avgCN (the average copy number). WGD values are 0 for no WGD event, 1 for one
#'  WGD event, 2 for several WGD events.
#' @export
#'
#' @examples
#' score_estwgd(segs.chas_example, oncoscan_na33.cov)
score_estwgd <- function(segments, kit.coverage) {
    # Get the average copy number
    avgcn <- score_avgcn(segments, kit.coverage)

    # Estimate the number of WGD events
    wgd.est <- ifelse(avgcn > 3.4, 2, ifelse(avgcn > 2.2, 1, 0))

    return(c(WGD = wgd.est, avgCN = avgcn))
}


#' Compute the number of LSTs, normalized by the number of WGD events.
#'
#' @details Compute the number of LSTs in non-LOH segments via the
#' \code{score_lst} function and subtract the extra noise induced by WGD events:
#'  nLST = LST - 7*W/2 where W is the number of WGD events.
#' A sample is HRD positive (deficient in HR pathway) if nLST is greater or equal
#' to the threshold (15 by default).
#' This score was linked to BRCA1/2-deficient tumors.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param n.wgd Number of whole-genome doubling events (0 if diploid).
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm.
#' @param threshold A number above which the test is returned positive (>=).
#'
#' @return A named list with the number of nLSTs and the corresponding label
#' ('Positive', 'Negative').
#' @export
#'
#' @examples
#' w <- score_estwgd(segs.chas_example, oncoscan_na33.cov)
#' score_nlst(segs.chas_example, w['WGD'], oncoscan_na33.cov)
score_nlst <- function(segments, n.wgd, kit.coverage, threshold=15) {
    lst.noLOH <- score_lst(segments[segments$cn.type != oncoscanR::cntype.loh], kit.coverage)

    nlst <- max(0, lst.noLOH - 3.5 * as.numeric(n.wgd))
    label <- ifelse(nlst >= threshold, "Positive", "Negative")
    return(c(nLST = nlst, HRD = label))
}


#' Computes the total number of Mbp altered.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm.
#' @param loh.rm A boolean (TRUE by default) to indicate whether LOH segments
#' should be excluded.
#'
#' @return A named list representing the Mbp altered in the sample and the total
#'  Mbp of the kit.
#' @export
#'
#' @examples
#' score_mbalt(segs.chas_example, oncoscan_na33.cov)
#' score_mbalt(segs.chas_example, oncoscan_na33.cov, FALSE)
score_mbalt <- function(segments, kit.coverage, loh.rm = TRUE) {
    # Compute the number Mbp present in the whole kit
    mb.kit <- round(sum(width(kit.coverage)/10^6))

    # Exclude (or not) the LOH segments
    segs <- segments
    if (loh.rm) {
        segs <- segments[segments$cn.type != oncoscanR::cntype.loh]
    }

    if (length(segs) == 0) {
        return(c(sample = 0, kit = mb.kit))
    } else {
        segs.w <- width(reduce(segs))/10^6  # Width of each segment in Mbp
        mb.alt <- round(sum(segs.w))

        return(c(sample = mb.alt, kit = mb.kit))
    }
}


#' Compute the genomic LOH score.
#'
#' @details The percentage genomic LOH score is computed as described in the
#' FoundationFocus CDx BRCA LOH assay; i.e. the percentage of bases covered by
#' the Oncoscan that display a loss of heterozygosity independently of the
#' number of copies, excluding chromosomal arms that have a global LOH (>=90% of
#' arm length).
#' To compute with the \code{armlevel_alt} function on LOH segments only).
#' This score was linked to BRCA1/2-deficient tumors.
#'
#' Of note the cn.subtype has to be defined first via the function
#' \code{get_cn_subtype}.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param arms.loh A list of arms with global/arm-level LOH alteration.
#' @param arms.hetloss A list of arms with global/arm-level heterozygous
#' loss.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm.
#'
#' @return An integer representing the percentage of LOH bases.
#' @export
#'
#' @examples
#' segs.chas_example$cn.subtype <- get_cn_subtype(segs.chas_example, 'F')
#' armlevel.loh <- armlevel_alt(segs.chas_example[segs.chas_example$cn.type == cntype.loh],
#'                              kit.coverage = oncoscan_na33.cov)
#' armlevel.hetloss <- armlevel_alt(segs.chas_example[segs.chas_example$cn.subtype == cntype.hetloss],
#'                              kit.coverage = oncoscan_na33.cov)
#' score_gloh(segs.chas_example, names(armlevel.loh), names(armlevel.hetloss), oncoscan_na33.cov)
score_gloh <- function(segments, arms.loh, arms.hetloss, kit.coverage) {
    # Test if any arms in arms.loh or arms.loss is not in kit.coverage
    if(length(c(arms.loh, arms.hetloss)) > 0){
        unknown.arms <- setdiff(c(arms.loh, arms.hetloss), seqnames(kit.coverage))
        if (length(unknown.arms)>0) {
          msg <- paste('Some arms were not found in the kit:',
                       paste(as.character(unknown.arms), collapse = TRUE))
          stop(msg)
        }
    }

    arms2discard <- c(arms.loh, arms.hetloss)
    sel <- segments$cn.subtype %in% c(oncoscanR::cntype.hetloss, oncoscanR::cntype.loh) &
        !(as.vector(seqnames(segments)) %in% arms2discard)
    width.loh <- IRanges::width(segments[sel])
    return(sum(width.loh)/sum(width(kit.coverage)))
}
