#' Select missing cells by classification
#'
#' Return a logical mask for the remaining missing values in selected
#' feature-condition states. Use it to restrict downstream replacements without
#' selecting observed values or synthetic rescue values.
#'
#' @param x An \code{imputefinder_result} from
#'   \code{\link{classify_missingness}}.
#' @param states A character vector containing exact state names:
#'   \code{"MNAR"}, \code{"MAR"}, or \code{"absent"}. The default selects all
#'   three. An empty character vector selects no cells.
#'
#' @details
#' A cell is selected only when it is currently missing and its retained
#' feature-condition block has a selected state. Dropped features are excluded.
#' Absence remains separate from MAR and MNAR; selecting it does not assign a
#' missingness mechanism or an imputation method.
#'
#' For experiment results, the recorded \code{assay_index} identifies the assay
#' that was classified. Older single-assay results also work. Re-run
#' \code{classify_missingness()} for older multi-assay results without this
#' record. Keep the result's assay order unchanged.
#'
#' Create masks before downstream processing. Align any replacement matrix to
#' the mask's feature and sample names before assigning values. The function
#' does not normalize, impute, change the result, or restore rescued cells.
#'
#' @return A base logical matrix with the dimensions, names, and order of
#'   \code{x$data}. Experiment results use the selected assay. The matrix
#'   contains no missing logical values. An empty retained result has zero rows.
#'
#' @seealso \code{\link{classify_missingness}}
#' @export
#' @examples
#' x <- rbind(
#'     on_off = c(NA, NA, NA, NA, 15, 16, 15, 17),
#'     mar = c(14, 15, NA, 16, 14, NA, 15, 16),
#'     low_anchor = c(8, NA, NA, NA, 14, 15, 16, 17)
#' )
#' colnames(x) <- paste0("s", seq_len(ncol(x)))
#' fit <- classify_missingness(x, rep(c("A", "B"), each = 4),
#'     cutoffs = c(A = 12, B = 12), rescue = "none")
#' missingness_mask(fit, "MNAR")
#' missingness_mask(fit, "MAR")
#' missingness_mask(fit, "absent")
#' stopifnot(identical(missingness_mask(fit), is.na(fit$data)))
missingness_mask <- function(x, states = c("MNAR", "MAR", "absent")) {
    if (!inherits(x, "imputefinder_result")) {
        stop("`x` must be an imputefinder_result.", call. = FALSE)
    }
    valid <- is.character(states) && !is.object(states) &&
        is.null(dim(states)) && !anyNA(states) &&
        all(states %in% c("MNAR", "MAR", "absent"))
    if (!valid) {
        stop("`states` must contain only MNAR, MAR, or absent.", call. = FALSE)
    }

    data <- .missingness_mask_data(x)
    groups <- .validate_result_groups(data, x$groups_by_sample)
    classified <- x$classifications
    fields <- c("feature", "condition", "state")
    if (!is.data.frame(classified) || !all(fields %in% names(classified)) ||
        !all(vapply(classified[fields], is.character, logical(1L))) ||
        anyNA(classified[fields])) {
        stop("`x` must contain non-missing feature-condition classifications.",
            call. = FALSE)
    }
    .validate_classification_states(classified)

    mask <- is.na(data)
    for (condition in unique(unname(groups))) {
        block <- classified[classified$condition == condition &
            classified$feature %in% rownames(data), , drop = FALSE]
        index <- match(rownames(data), block$feature)
        if (nrow(block) != nrow(data) || anyDuplicated(block$feature) ||
            anyNA(index)) {
            stop("`x` must have one classification per retained ",
                "feature and condition.", call. = FALSE)
        }
        mask[!block$state[index] %in% states, which(groups == condition)] <- FALSE
    }
    mask
}

.missingness_mask_data <- function(x) {
    data <- x$data
    if (methods::is(data, "SummarizedExperiment")) {
        count <- length(SummarizedExperiment::assays(data, withDimnames = FALSE))
        index <- x$assay_index
        if (is.null(index) && count == 1L) {
            index <- 1L
        }
        valid <- is.numeric(index) && !is.object(index) &&
            is.null(dim(index)) && length(index) == 1L &&
            is.finite(index) && index >= 1 && index <= count &&
            index == floor(index)
        if (!valid) {
            stop("`x` needs a valid recorded assay. Re-run ",
                "classify_missingness() on the original input.", call. = FALSE)
        }
        data <- SummarizedExperiment::assay(data, index)
    }
    if (!is.matrix(data) || !is.numeric(data)) {
        stop("The classified data must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(data) > 0L) {
        .validate_axis_names(rownames(data), nrow(data), "feature")
    }
    .validate_axis_names(colnames(data), ncol(data), "sample")
    data
}
