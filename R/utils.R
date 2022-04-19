# utils.R Functions that are too general to fit in any other R file.
# Author: Yann Christinat Date: 24.9.2019


#' Compare two CNV segments and test if equal
#'
#' @param x A \code{GRanges} object one segment, with optional fields \code{cn},
#'  \code{cn.type} and \code{cn.subtype}
#' @param y A \code{GRanges} object one segment, with optional fields \code{cn},
#'  \code{cn.type} and \code{cn.subtype}
#'
#' @return A boolean. TRUE if the segments' chromosome/arm, start/end and CN
#' tags match.
#'
#' @export
#'
#' @importFrom GenomicRanges seqnames GRanges
#' @importFrom IRanges start end
#'
#' @examples
#' library(GenomicRanges)
#' library(IRanges)
#' s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
#'               cn = 5, cn.type = cntype.gain, cn.subtype = cntype.weakamp)
#' s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
#'               cn = 6, cn.type = cntype.gain, cn.subtype = cntype.weakamp)
#' same_segments(s1,s2)
same_segments <- function(x, y) {
    # Test the class of the parameters
    if (!is(x, "GRanges") || !is(y, "GRanges")) {
        stop("Function 'compare_segments' can only handle GRanges objects.")
    }

    # Ensures that the GRanges contain only 1 segment each
    if (length(x) != length(y) || length(x) != 1) {
        msg <- paste("Function 'compare_segments' can only handle GRanges with only one segment. Found:",
                     length(x), length(y))
        stop(msg)
    }

    # Replace any NA by NULL in cn, cn.type and cn.subtype
    if (!is.null(x$cn) && is.na(x$cn)) {
        x$cn <- NULL
    }
    if (!is.null(x$cn.type) && is.na(x$cn.type)) {
        x$cn.type <- NULL
    }
    if (!is.null(x$cn.subtype) && is.na(x$cn.subtype)) {
        x$cn.subtype <- NULL
    }

    if (!is.null(y$cn) && is.na(y$cn)) {
        y$cn <- NULL
    }
    if (!is.null(y$cn.type) && is.na(y$cn.type)) {
        y$cn.type <- NULL
    }
    if (!is.null(y$cn.subtype) && is.na(y$cn.subtype)) {
        y$cn.subtype <- NULL
    }


    # Compare segment position and arm
    if (as.character(seqnames(x)) != as.character(seqnames(y)) ||
        start(x) != start(y) || end(x) != end(y)) {
        return(FALSE)
    } else {
        # Same segment: Test cn (if set)
        if ((is.null(x$cn) && !is.null(y$cn)) ||
            (!is.null(x$cn) && is.null(y$cn))) {
            # Only one cn is null
            return(FALSE)
        } else if (!is.null(x$cn) && !is.null(y$cn)) {
            # cn is set on both segments
            if (x$cn != y$cn) {
                return(FALSE)
            }
        }

        # Test cn.type (if set)
        if ((is.null(x$cn.type) && !is.null(y$cn.type)) ||
            (!is.null(x$cn.type) && is.null(y$cn.type))) {
            # Only one cn.type is null
            return(FALSE)
        } else if (!is.null(x$cn.type) & !is.null(y$cn.type)) {
            # cn is set on both segments
            if (x$cn.type != y$cn.type) {
                return(FALSE)
            }
        }

        # Test cn.subtype (if set)
        if ((is.null(x$cn.subtype) && !is.null(y$cn.subtype)) ||
            (!is.null(x$cn.subtype) && is.null(y$cn.subtype))) {
            # Only one cn.subtype is null
            return(FALSE)
        } else if (!is.null(x$cn.subtype) && !is.null(y$cn.subtype)) {
            # cn is set on both segments
            if (x$cn.subtype != y$cn.subtype) {
                return(FALSE)
            }
        }
    }
    return(TRUE)
}


#' Test if the object is a CNV segment
#'
#' @param obj Object to test
#' @param raise_error Boolean. If TRUE then raises an error if the test fails.
#'
#' @return Boolean. TRUE if the object is of class \code{GRanges} with fields
#' \code{cn} and \code{cn.type}.
#' @export
#'
#' @examples
#' library(GenomicRanges)
#' library(IRanges)
#' s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
#'               cn = 6, cn.type = cntype.gain, cn.subtype = cntype.weakamp)
#' is_cn_segment(s1)
#' s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
#'               cn = 6)
#' is_cn_segment(s2, raise_error = FALSE)
is_cn_segment <- function(obj, raise_error = TRUE) {
    test <- is(obj, "GRanges") && (length(obj) == 0 || (!is.null(obj$cn) &&
        !is.null(obj$cn.type)))
    if (!test && raise_error) {
        stop("The parameter has to be a CNV segment (GRanges object with 'cn' and 'cn.type' columns).")
    }
    return(test)
}


#' Small function to test if two sets have the same CN segments
#'
#' @param grA A \code{GRanges} object with fields \code{cn} and \code{cn.type}
#' @param grB A \code{GRanges} object with fields \code{cn} and \code{cn.type}
#'
#' @return TRUE if the two sets contain exactly the same CN segments
#' @export
#'
#' @examples
#' library(GenomicRanges)
#' library(IRanges)
#' s1 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
#'               cn = 5, cn.type = cntype.gain, cn.subtype = cntype.weakamp)
#' s2 <- GRanges(seqnames = '1p', ranges = IRanges(start = 1, end = 10),
#'               cn = 6, cn.type = cntype.gain, cn.subtype = cntype.weakamp)
#' same_segmentsets(s1,s2)
same_segmentsets <- function(grA, grB) {
    if (length(grA) != length(grB)) {
        return(FALSE)
    }

    found <- 0
    for (i in seq_along(grA)) {
        segA <- grA[i]
        for (j in seq_along(grB)) {
            segB <- grB[j]
            if (same_segments(segA, segB)) {
                found <- found + 1
            }
        }
    }
    return(found == length(grA))
}
