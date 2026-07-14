.assign_condition_states <- function(statistics, cutoffs) {
    .validate_state_statistics(statistics)
    conditions <- sort(unique(statistics$condition), method = "radix")
    cutoffs <- .normalise_resolved_cutoffs(cutoffs, conditions)

    for (condition in conditions) {
        condition_rows <- statistics$condition == condition
        has_missing <- any(statistics$missing_count[condition_rows] > 0L)
        if (!has_missing) {
            cutoffs[[condition]] <- NA_real_
        } else if (is.na(cutoffs[[condition]])) {
            stop(
                sprintf(
                    "Condition `%s` has incomplete features but no resolved cutoff.",
                    condition
                ),
                call. = FALSE
            )
        }
    }

    row_cutoffs <- unname(cutoffs[statistics$condition])
    state <- rep(NA_character_, nrow(statistics))
    complete <- statistics$missing_count == 0L
    state[complete] <- "complete"

    incomplete <- !complete
    below_cutoff <- incomplete &
        statistics$mean_intensity < row_cutoffs
    state[below_cutoff] <- "MNAR"

    mar_side <- incomplete & !below_cutoff
    strict_majority <-
        statistics$observed_count > statistics$sample_count / 2
    state[mar_side & strict_majority] <- "MAR"
    state[mar_side & !strict_majority] <- "insufficient"

    data.frame(
        feature = statistics$feature,
        condition = statistics$condition,
        sample_count = as.integer(statistics$sample_count),
        observed_count = as.integer(statistics$observed_count),
        missing_count = as.integer(statistics$missing_count),
        missing_fraction = as.numeric(statistics$missing_fraction),
        mean_intensity = as.numeric(statistics$mean_intensity),
        cutoff = as.numeric(row_cutoffs),
        state = state,
        seeded = as.logical(statistics$seeded),
        retained = rep(NA, nrow(statistics)),
        drop_reason = rep(NA_character_, nrow(statistics)),
        stringsAsFactors = FALSE
    )
}

.normalise_resolved_cutoffs <- function(cutoffs, conditions) {
    valid <- is.numeric(cutoffs) &&
        is.null(dim(cutoffs)) &&
        length(cutoffs) == length(conditions) &&
        !is.null(names(cutoffs)) &&
        !anyNA(names(cutoffs)) &&
        all(nzchar(names(cutoffs))) &&
        !anyDuplicated(names(cutoffs)) &&
        setequal(names(cutoffs), conditions) &&
        !any(is.nan(cutoffs)) &&
        !any(is.infinite(cutoffs))
    if (!valid) {
        stop(
            "Resolved cutoffs must contain one finite-or-NA value per condition.",
            call. = FALSE
        )
    }

    stats::setNames(as.numeric(cutoffs[conditions]), conditions)
}

.validate_state_statistics <- function(statistics) {
    required <- c(
        "feature",
        "condition",
        "sample_count",
        "observed_count",
        "missing_count",
        "missing_fraction",
        "mean_intensity",
        "seeded"
    )
    if (!is.data.frame(statistics) || !all(required %in% names(statistics))) {
        stop("`statistics` does not satisfy the state-input schema.", call. = FALSE)
    }
    valid_counts <- statistics$sample_count > 0L &
        statistics$observed_count >= 0L &
        statistics$missing_count >= 0L &
        statistics$observed_count + statistics$missing_count ==
            statistics$sample_count
    valid <- !anyNA(statistics[, required]) &&
        all(valid_counts) &&
        all(is.finite(statistics$mean_intensity))
    if (!valid) {
        stop("`statistics` contains invalid state inputs.", call. = FALSE)
    }

    invisible(statistics)
}
