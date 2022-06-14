#' Compare two CNV segments and test if equal
#'
#' @param x A \code{GRanges} object one segment, with optional fields \code{cn},
#'  \code{cn.type} and \code{cn.subtype}
#' @param y A \code{GRanges} object one segment, with optional fields \code{cn},
#'  \code{cn.type} and \code{cn.subtype}
#'
#' @return A boolean. TRUE if the segments' chromosome/arm, start/end and CN
#' tags match.
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

    # Replace any NA by NULL in cn, cntype
    if (!is.null(x$cn) && is.na(x$cn)) {
        x$cn <- NULL
    }
    if (!is.null(x$cn.type) && is.na(x$cn.type)) {
        x$cn.type <- NULL
    }

    if (!is.null(y$cn) && is.na(y$cn)) {
        y$cn <- NULL
    }
    if (!is.null(y$cn.type) && is.na(y$cn.type)) {
        y$cn.type <- NULL
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
    }
    return(TRUE)
}


#' Small function to test if two sets have the same CN segments
#'
#' @param grA A \code{GRanges} object with fields \code{cn} and \code{cn.type}
#' @param grB A \code{GRanges} object with fields \code{cn} and \code{cn.type}
#'
#' @return TRUE if the two sets contain exactly the same CN segments
#'
#' @noRd
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


#' Loads the results of the armlevel_alt function and compare it to the results
#' from a text file
armlevel.test <- function(dat.armlevel, armlevel.fn){
    # Get ground truth
    dat.true <- read.csv(armlevel.fn, header = TRUE, row.names = 1,
                         stringsAsFactors = FALSE, colClasses = 'character')

    #Check AMP, GAIN, LOSS and LOH
    amp.test <- identical(sort(dat.armlevel[['AMP']]),
                          sort(unlist(strsplit(dat.true['AMPL','Arms'], ','))))
    gain.test <- identical(sort(dat.armlevel[['GAIN']]),
                           sort(unlist(strsplit(dat.true['GAIN','Arms'], ','))))
    loh.test <- identical(sort(dat.armlevel[['LOH']]),
                          sort(unlist(strsplit(dat.true['LOH','Arms'], ','))))
    loss.test <- identical(sort(dat.armlevel[['LOSS']]),
                           sort(unlist(strsplit(dat.true['LOSS','Arms'], ','))))
    return(c(amp.test, gain.test, loh.test, loss.test))
}
