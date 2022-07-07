# input.R Functions to handle loading/cleaning of segments and reference files.
# Author: Yann Christinat Date: 19.9.2019

#' Load a ChAS text export file.
#'
#' @details The ChAS file is expected to have the following column names:
#' 'CN State' (number or empty), 'Type' (expected value: 'Gain', 'Loss' or
#' 'LOH') and 'Full Location' (in the format 'chr:start-end').
#'
#' The segments are attributed to each chromosome arm and split if necessary.
#'
#' @param filename Path to the ChAS file.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm by the kit.
#'
#' @return A \code{GRanges} object containing the segments, their copy number
#' (field \code{cn}), their copy
#' number types (field \code{cntype}). \code{cntype} contains either 'Gain',
#' 'Loss' or 'LOH'.
#' If the file contains twice the same segment or does not respect the format
#' specifications, then an error is
#' raised. NB. The chromosome name is in the format '1' and not 'chr1' and will
#' be transformed if needed.
#'
#' @export
#'
#' @import readr
#' @import GenomicRanges
#' @import IRanges
#' @importFrom methods is
#'
#' @examples
#' segs.filename <- system.file('extdata', 'chas_example.txt',
#'   package = 'oncoscanR')
#' segs.chas_example <- load_chas(segs.filename, oncoscan_na33.cov)
load_chas <- function(filename, kit.coverage) {
    # Entry check on arguments
    stopifnot(is(kit.coverage, "GRanges"))

    # Reads in the ChAS file
    cols <- cols_only(`CN State` = col_number(), Type = col_character(),
                      `Full Location` = col_character())
    oncoscan_table <- read_tsv(filename, comment = "#", col_names = TRUE,
                               col_types = cols)

    if (dim(oncoscan_table)[2] != 3) {
        stop("Parsing ChAS file failed.")
    }
    if (dim(oncoscan_table)[1] == 0) {
        warning("No segments loaded!")
        return(GRanges())
    }

    # Process the data into GRanges segments
    segs <- process_chas(oncoscan_table, kit.coverage)

    # Test for duplicated entries
    errmsg <- paste("The file", filename, "contains duplicated entries.")
    for (arm in unique(seqnames(segs))) {
        arm_segs <- sort(segs[seqnames(segs) == arm])
        if (length(arm_segs) > 1) {
            for (i in 2:length(arm_segs)) {
                segA <- arm_segs[i - 1]
                segB <- arm_segs[i]
                if (IRanges::start(segA) == IRanges::start(segB) &
                    IRanges::end(segA) == IRanges::end(segB) &
                    segA$cn.type == segB$cn.type) {
                    if (segA$cn.type == cntypes$LOH) {
                        stop(errmsg)
                    } else if (segA$cn == segB$cn) {
                        stop(errmsg)
                    }
                }
            }
        }
    }

    if (length(segs) == 0) {
        warning("No segments loaded!")
    }

    return(segs)
}


#' Process ChAS table.
#'
#' @details Used in the load_chas function.
#'
#' @param oncoscan_table A tibble with the following column names:
#' 'CN State' (number or empty), 'Type' (expected value: 'Gain', 'Loss' or
#' 'LOH') and 'Full Location' (in the format 'chr:start-end').
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm by the kit.
#'
#' @return A \code{GRanges} object containing the segments, their copy number
#' (field \code{cn}), their copy
#' number types (field \code{cntype}). \code{cntype} contains either 'Gain',
#' 'Loss' or 'LOH'.
#'
#' @noRd
process_chas <- function(oncoscan_table, kit.coverage){
    # Allocate the place for the GRanges segments. At most we will end up with
    # twice as much segments as present in the raw data.
    segments_list <- vector(mode = "list", length = 2 * dim(oncoscan_table)[1])

    # Parse throuh all lines of the data
    counter <- 0
    for (i in seq(dim(oncoscan_table)[1])) {
        counter <- counter + 1
        # Start: Extract chr no, start, end, copy number

        # Full location from oncoscan file
        loc_cord <- oncoscan_table$`Full Location`[i]
        loc_cord_list <- strsplit(loc_cord, split = ":", fixed = TRUE)[[1]]

        # seg chr no based on oncoscan file
        seg_chr <- loc_cord_list[1]
        seg_chr <- gsub("chr", "", seg_chr)  #remove 'chr' if present

        # seg coordinates
        coord_list <- strsplit(loc_cord_list[2], split = "-", fixed = TRUE)[[1]]
        seg_start <- as.numeric(coord_list[1])
        seg_end <- as.numeric(coord_list[2])

        seg_cn <- oncoscan_table$`CN State`[i]
        seg_cntype <- oncoscan_table$Type[i]

        # Test if copy number type is correct
        if (!(seg_cntype %in% cntypes)) {
            msg <- paste("The column \"Type\" should contain only the",
                         "following values:", cntypes, "whereas",
                         seg_cntype, "was found!")
            stop(msg)
        }

        # Get the arms
        parm <- kit.coverage[seqnames(kit.coverage) == paste0(seg_chr, "p")]
        qarm <- kit.coverage[seqnames(kit.coverage) == paste0(seg_chr, "q")]

        kit.arms <- levels(seqnames(kit.coverage))
        if (length(parm) > 0 || length(qarm) > 0) {
            if (length(parm) == 0 || seg_start > end(parm)) {
                # Then the segment is only in the q arm
                sn <- factor(paste0(seg_chr, "q"), levels = kit.arms)
                rg <- IRanges(start = seg_start, end = seg_end)
                seg <- GRanges(seqnames=sn, ranges=rg,
                               cn=seg_cn, cn.type=seg_cntype)
                segments_list[[counter]] <- seg
            } else if (length(qarm) == 0 || seg_end < start(qarm)) {
                # Then the segment is only in the p arm
                sn <- factor(paste0(seg_chr, "p"), levels = kit.arms)
                rg <- IRanges(start = seg_start, end = seg_end)
                seg <- GRanges(seqnames=sn, ranges=rg,
                               cn=seg_cn, cn.type=seg_cntype)
                segments_list[[counter]] <- seg
            } else {
                # Create a segment for each arm
                sn <- factor(paste0(seg_chr, c("p", "q")), levels = kit.arms)
                rg <- IRanges(start = rep(seg_start, 2), end = rep(seg_end, 2))
                seg <- GRanges(seqnames = sn, ranges = rg)

                # Get the segments overlapping with the arms
                o <- IRanges::findOverlapPairs(c(parm, qarm), seg)
                new_segs <- pintersect(o)
                new_segs$cn <- rep(seg_cn, length(new_segs))
                new_segs$cn.type <- rep(seg_cntype, length(new_segs))
                new_segs$hit <- NULL

                segments_list[[counter]] <- new_segs
            }
        }
    }

    segs <- do.call("c", unlist(segments_list))
    return(segs)
}


#' Load the oncoscan coverage BED file into a GenomicRanges object.
#'
#' @details Expects the following columns from the BED file (no header):
#'   1. Name of the chromosomal arm (e.g. "1p")
#'   2. Start position of the arm
#'   3. End position of the arm
#'
#' @param filename Path to the coverage BED file.
#'
#' @return A \code{GRanges} object containing the regions covered on each
#' chromosome arm.
#'
#' @export
#'
#' @import GenomicRanges
#' @import IRanges
#' @importFrom utils read.table
#'
#' @examples
#' oncoscan_na33.cov <- get_oncoscan_coverage_from_bed(
#'        system.file('extdata', 'OncoScan.na33.r2.cov.processed.bed',
#'        package = 'oncoscanR'))
get_oncoscan_coverage_from_bed <- function(filename) {
    # Read the annotation file
    dat <- read.table(filename, sep = '\t', header = FALSE)
    colnames(dat) <- c('Arm', 'Start', 'Stop')

    cov <- GRanges(seqnames = factor(dat$Arm),
                   ranges = IRanges(start = dat$Start, end = dat$Stop))
    return(cov)
}


#' Return all segments of type LOH, independently of the copy number.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#'
#' @return A \code{GRanges} object containing the selected segments, their copy
#' number and copy number types.
#' @export
#'
#' @examples
#' segs.loh <- get_loh_segments(segs.chas_example)
get_loh_segments <- function(segments) {
    is_cn_segment(segments, raise_error = TRUE)

    return(segments[segments$cn.type == cntypes$LOH])
}

#' Return all segments with loss of 1 or 2 copies.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#'
#' @return A \code{GRanges} object containing the selected segments, their copy
#' number and copy number types.
#' @export
#'
#' @examples
#' segs.loh <- get_loh_segments(segs.chas_example)
get_loss_segments <- function(segments) {
    is_cn_segment(segments, raise_error = TRUE)

    return(segments[segments$cn.type == cntypes$Loss])
}

#' Return all segments with heterozygous loss.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#'
#' @return A \code{GRanges} object containing the selected segments, their copy
#' number and copy number types.
#' @export
#'
#' @examples
#' segs.hetloss <- get_hetloss_segments(segs.chas_example)
get_hetloss_segments <- function(segments) {
    is_cn_segment(segments, raise_error = TRUE)

    return(segments[segments$cn.type == cntypes$Loss & segments$cn >= 1])
}

#' Return all segments with homozygous loss.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#'
#' @return A \code{GRanges} object containing the selected segments, their copy
#' number and copy number types.
#' @export
#'
#' @examples
#' get_homloss_segments <- get_homloss_segments(segs.chas_example)
get_homloss_segments <- function(segments) {
    is_cn_segment(segments, raise_error = TRUE)

    return(segments[segments$cn.type == cntypes$Loss & segments$cn < 1])
}

#' Return all segments with gain of copies.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#'
#' @return A \code{GRanges} object containing the selected segments, their copy
#' number and copy number types.
#' @export
#'
#' @examples
#' segs.gain <- get_gain_segments(segs.chas_example)
get_gain_segments <- function(segments) {
    is_cn_segment(segments, raise_error = TRUE)

    return(segments[segments$cn.type == cntypes$Gain])
}

#' Return all segments with an amplification (5 or more copies).
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#'
#' @return A \code{GRanges} object containing the selected segments, their copy
#' number and copy number types.
#' @export
#'
#' @examples
#' segs.amp <- get_amp_segments(segs.chas_example)
get_amp_segments <- function(segments) {
    is_cn_segment(segments, raise_error = TRUE)

    return(segments[segments$cn.type == cntypes$Gain & segments$cn >= 5])
}


#' Trim segments with respect to the kit's coverage.
#'
#' @details All segments that are not entirely contained within the kit coverage
#'  will be trimmed to the coverage's limits.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param kit.coverage A \code{GRanges} object containing the regions covered on
#'  each chromosome arm.
#'
#' @return A \code{GRanges} object containing the cleaned segments, their copy
#' number and copy number types.
#' @export
#'
#' @import GenomicRanges
#' @import IRanges
#'
#' @examples
#' segs.trimmed <- trim_to_coverage(segs.chas_example, oncoscan_na33.cov)
trim_to_coverage <- function(segments, kit.coverage) {
    is_cn_segment(segments, raise_error = TRUE)

    if (length(segments) == 0) {
        return(segments)
    }

    # Apply on each arm...
    segs.clean <- lapply(unique(seqnames(segments)), function(arm) {
        # Trim segments wrt coverage
        arm.cov <- kit.coverage[seqnames(kit.coverage) == arm]
        o <- IRanges::findOverlapPairs(segments[seqnames(segments) == arm],
                                       arm.cov)
        new_segs <- sort(pintersect(o))
        if (length(new_segs) == 0) {
            return(new_segs)
        }
        new_segs$hit <- NULL
        return(new_segs)
    })
    return(do.call("c", unlist(segs.clean)))
}


#' Merge segments with respect to the kit resolution and the copy number.
#'
#' @details If two segments are at a distance smaller than the resolution, then
#' the segments are merged if the
#' share the same \code{cn} value. Note that the function does not look at the
#' copy number type or subtype but
#' only at the actual copy number to decide whether segments can be merged.
#'
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param kit.resolution Number >0 indicating the minimum segment size
#' detectable by the technique (in kilobases).
#' Defaults to the Oncoscan assay resolution outside of cancer genes: 300Kb.
#'
#' @return A \code{GRanges} object containing the cleaned segments, their copy
#' number and copy number types.
#' @export
#'
#' @import GenomicRanges
#'
#' @examples
#' segs.merged <- merge_segments(segs.chas_example)
#' segs.merged_50k <- merge_segments(segs.chas_example, 50)
merge_segments <- function(segments, kit.resolution = 300) {
    is_cn_segment(segments, raise_error = TRUE)

    if (kit.resolution < 1 / 1000) {
        stop("Kit resolution has to be greater than zero.")
    }

    if (length(segments) == 0) {
        return(segments)
    }

    # Go through each arm...
    segs.merged <- lapply(unique(seqnames(segments)), function(arm) {
        # Go through each copy number
        armsegs <- segments[seqnames(segments) == arm]
        armsegs.merged <- lapply(unique(armsegs$cn), function(cn) {
            segs <- NULL
            if (is.na(cn)) {
                segs <- armsegs[is.na(armsegs$cn)]
            } else {
                segs <- armsegs[!is.na(armsegs$cn) & armsegs$cn == cn]
            }
            new_segs <-
                reduce(segs, min.gapwidth = kit.resolution * 1000)
            # Re-set the copy number, cn type and subtype
            new_segs$cn <- cn
            new_segs$cn.type <- segs[1]$cn.type
            if (!is.null(segs$cn.subtype)) {
                new_segs$cn.subtype <- segs[1]$cn.subtype
            }
            return(new_segs)
        })
        do.call("c", unlist(armsegs.merged))
    })
    return(do.call("c", unlist(segs.merged)))
}


#' Trim LOH segments with respect to loss segments.
#'
#' @details LOH segments completely contained within (or equal to)  a copy loss
#' segment are deleted.
#' LOH segments partially overlapping (on one end only) with a copy loss segment
#'  are trimmed to remove the overlap or split into several segments.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#'
#' @return A \code{GRanges} object containing the cleaned segments, their copy
#' number and copy number types.
#' @export
#'
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#'
#' @examples
#' segs.adj <- adjust_loh(segs.chas_example)
adjust_loh <- function(segments) {
    is_cn_segment(segments, raise_error = TRUE)
    if (length(segments) == 0) { return(segments) }

    # Apply on each arm
    loh.adj <- lapply(unique(seqnames(segments)), function(arm) {
        segs.loh <- segments[segments$cn.type == cntypes$LOH &
                                 seqnames(segments) == arm]

        if (length(segs.loh) == 0) { return(GRanges()) }
        segs.loss <- segments[segments$cn.type == cntypes$Loss &
                                  seqnames(segments) == arm]

        pos.all <- sort(unique(c(start(segs.loh), end(segs.loh),
                                 start(segs.loss), end(segs.loss))))
        dt <- data.frame(
            row.names = pos.all,
            loss = factor(rep('nd', length(pos.all)),
                          levels = c('nd', 'start', 'end')),
            loh = factor(rep('nd', length(pos.all)),
                         levels = c('nd', 'start', 'end'))
        )
        dt[as.character(start(segs.loh)), 'loh'] <- 'start'
        dt[as.character(end(segs.loh)), 'loh'] <- 'end'
        dt[as.character(start(segs.loss)), 'loss'] <- 'start'
        dt[as.character(end(segs.loss)), 'loss'] <- 'end'

        loh.toadd <- getLOHtoadd(dt)

        if (length(loh.toadd) == 0) {
            return(GRanges())
        }
        else {
            return(GRanges(seqnames = rep(arm, length(loh.toadd)),
                           ranges = loh.toadd
            ))
        }
    })

    new.loh <- do.call("c", unlist(loh.adj))
    new.loh$cn <- as.numeric(NA)
    new.loh$cn.type <- as.character(cntypes$LOH)
    if (!is.null(segments$cn.subtype)) {
        new.loh$cn.subtype <- as.character(cntypes$LOH)
    }

    return(c(segments[segments$cn.type != cntypes$LOH], new.loh))
}


#' Given a list of segments, trims or split the LOH segments if they overlap
#' with a loss.
#'
#' @details Used in the adjust_loh function
#'
#' @param dt A data.frame with the position as row names and two columns: 'loh'
#' and 'loss' with values 'end' or 'start', indicating the starts and ends of
#' the segments.
#'
#' @return a list of IRanges segments
#'
#' @noRd
getLOHtoadd <- function(dt){
    in.loh <- FALSE
    in.loss <- FALSE
    lohseg.start <- NULL

    loh.toadd <- IRanges()
    for (i in seq(dim(dt)[1])) {
        # Update status whether we are in a loss
        in.loss <- ifelse(dt[i, 'loss'] == 'start', TRUE,
                          ifelse(dt[i, 'loss'] == 'end', FALSE, in.loss))

        # Update status whether we are in a LOH
        in.loh <- ifelse(dt[i, 'loh'] == 'start', TRUE,
                         ifelse(dt[i, 'loh'] == 'end', FALSE, in.loh))

        if (in.loh & !in.loss) {
            # Start of a LOH segment
            lohseg.start <- ifelse(dt[i, 'loss'] == 'end',
                                   as.numeric(rownames(dt)[i]) + 1,
                                   as.numeric(rownames(dt)[i]))
        }
        else if ((dt[i, 'loh'] == 'nd' && dt[i, 'loss'] == 'start' && in.loh) ||
                 (dt[i, 'loh'] == 'end' && dt[i, 'loss'] == 'start') ||
                 (dt[i, 'loh'] == 'end' && dt[i, 'loss'] == 'nd' && !in.loss)) {
            # End of a LOH segment
            lohseg.end <- ifelse(dt[i, 'loss'] == 'start',
                                 as.numeric(rownames(dt)[i]) - 1,
                                 as.numeric(rownames(dt)[i]))

            # Add segment
            loh.toadd <- append(loh.toadd,
                                IRanges(start = lohseg.start, end = lohseg.end))
            }
    }

    return(loh.toadd)
}

#' Remove segments smaller than the kit resolution.
#'
#' @param segments A \code{GRanges} object containing the segments, their copy
#' number and copy number types.
#' @param threshold Number indicating the minimum segment size to be kept (in
#' kilobases).
#' Defaults to the Oncoscan assay resolution outside of cancer genes: 300Kb.
#'
#' @return A \code{GRanges} object containing the cleaned segments, their copy
#' number and copy number types.
#' @export
#'
#' @import IRanges
#'
#' @examples
#' segs.300k <- prune_by_size(segs.chas_example)
#' segs.50k <- prune_by_size(segs.chas_example, 50)
prune_by_size <- function(segments, threshold = 300) {
    is_cn_segment(segments, raise_error = TRUE)

    if (threshold < 0) {
        stop("Threshold has to be greater than or equal to zero.")
    }
    if (length(segments) == 0) {
        return(segments)
    }

    return(segments[width(segments) >= threshold * 1000])
}
