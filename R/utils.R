# utils.R Functions that are too general to fit in any other R file.
# Author: Yann Christinat Date: 24.9.2019


#' Test if the object is a CNV segment (GRange object with metadata: cn and cntype)
#'
#' @param obj Object to test
#' @param raise_error Boolean. If TRUE then raises an error if the test fails.
#'
#' @return Boolean. TRUE if the object is of class \code{GRanges} with fields
#' \code{cn} and \code{cn.type}.
#'
#' @noRd
#'
#' @importFrom methods is
is_cn_segment <- function(obj, raise_error = TRUE) {
    test <- is(obj, "GRanges") && (length(obj) == 0 || (!is.null(obj$cn) &&
        !is.null(obj$cn.type)))
    if (!test && raise_error) {
        stop("The parameter has to be a CNV segment (GRanges object with 'cn' and 'cntype' columns).")
    }
    return(test)
}


