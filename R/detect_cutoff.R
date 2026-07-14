.resolve_manual_cutoffs <- function(statistics, cutoffs, conditions = NULL) {
    .validate_cutoff_statistics(statistics)
    if (is.null(conditions)) {
        conditions <- sort(unique(statistics$condition), method = "radix")
    } else {
        conditions <- .normalise_cutoff_conditions(conditions)
    }
    unknown_statistics <- setdiff(unique(statistics$condition), conditions)
    if (length(unknown_statistics) > 0L) {
        stop(
            "Cutoff conditions must include every statistics condition.",
            call. = FALSE
        )
    }
    manual <- .normalise_manual_cutoffs(cutoffs, conditions)

    resolved <- stats::setNames(rep(NA_real_, length(conditions)), conditions)
    diagnostics <- stats::setNames(
        vector("list", length(conditions)),
        conditions
    )

    for (condition_index in seq_along(conditions)) {
        condition <- conditions[[condition_index]]
        condition_statistics <- statistics[
            statistics$condition == condition,
            ,
            drop = FALSE
        ]
        has_missing <- any(condition_statistics$missing_count > 0L)

        if (!has_missing) {
            diagnostics[[condition_index]] <- .not_needed_cutoff_diagnostic()
        } else if (condition %in% names(manual)) {
            cutoff <- unname(manual[[condition]])
            resolved[[condition]] <- cutoff
            diagnostics[[condition_index]] <- .manual_cutoff_diagnostic(
                condition,
                cutoff,
                condition_statistics$mean_intensity
            )
        }
    }

    list(cutoffs = resolved, diagnostics = diagnostics)
}

.normalise_cutoff_conditions <- function(conditions) {
    valid <- is.character(conditions) &&
        is.null(dim(conditions)) &&
        !anyNA(conditions) &&
        all(nzchar(conditions)) &&
        !anyDuplicated(conditions)
    if (!valid) {
        stop("Cutoff conditions must be unique, non-empty labels.", call. = FALSE)
    }

    sort(conditions, method = "radix")
}

.normalise_manual_cutoffs <- function(cutoffs, conditions) {
    if (is.null(cutoffs)) {
        return(stats::setNames(numeric(), character()))
    }

    valid_shape <- is.numeric(cutoffs) &&
        is.null(dim(cutoffs)) &&
        length(cutoffs) > 0L &&
        !is.null(names(cutoffs))
    if (!valid_shape) {
        stop(
            "`cutoffs` must be NULL or a non-empty named numeric vector.",
            call. = FALSE
        )
    }

    cutoff_names <- names(cutoffs)
    valid_names <- length(cutoff_names) == length(cutoffs) &&
        !anyNA(cutoff_names) &&
        all(nzchar(cutoff_names)) &&
        !anyDuplicated(cutoff_names)
    if (!valid_names) {
        stop(
            "Manual cutoff names must be unique, non-empty condition labels.",
            call. = FALSE
        )
    }
    if (any(!is.finite(cutoffs))) {
        stop(
            "Manual cutoffs must contain only finite values.",
            call. = FALSE
        )
    }

    unknown <- sort(setdiff(cutoff_names, conditions), method = "radix")
    if (length(unknown) > 0L) {
        stop(
            sprintf(
                "Manual cutoff names must match condition labels exactly: %s.",
                paste(sprintf("`%s`", unknown), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    stats::setNames(as.numeric(cutoffs), cutoff_names)
}

.manual_cutoff_diagnostic <- function(condition, cutoff, mean_intensity) {
    observed_range <- range(mean_intensity)
    warnings <- character()
    if (cutoff < observed_range[[1L]] || cutoff > observed_range[[2L]]) {
        warnings <- sprintf(
            paste0(
                "Manual cutoff %s is outside the observed feature-mean ",
                "range [%s, %s] for condition `%s`."
            ),
            format(cutoff, trim = TRUE),
            format(observed_range[[1L]], trim = TRUE),
            format(observed_range[[2L]], trim = TRUE),
            condition
        )
    }

    list(
        source = "manual",
        method = "manual",
        method_version = NA_character_,
        quality = list(
            observed_mean_range = stats::setNames(
                observed_range,
                c("minimum", "maximum")
            )
        ),
        warnings = warnings
    )
}

.not_needed_cutoff_diagnostic <- function() {
    list(
        source = "not_needed",
        method = NA_character_,
        method_version = NA_character_,
        quality = list(),
        warnings = character()
    )
}

.validate_cutoff_statistics <- function(statistics) {
    required <- c("condition", "missing_count", "mean_intensity")
    if (!is.data.frame(statistics) || !all(required %in% names(statistics))) {
        stop(
            "`statistics` must contain condition, missing_count, and mean_intensity.",
            call. = FALSE
        )
    }
    if (anyNA(statistics$condition) || any(!nzchar(statistics$condition)) ||
        any(!is.finite(statistics$mean_intensity))) {
        stop("`statistics` contains invalid cutoff inputs.", call. = FALSE)
    }

    invisible(statistics)
}
