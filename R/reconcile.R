.assign_condition_states <- function(statistics, cutoffs) {
    .validate_state_statistics(statistics)
    row_cutoffs <- .classification_row_cutoffs(statistics, cutoffs)
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

.classification_row_cutoffs <- function(statistics, cutoffs) {
    conditions <- sort(unique(statistics$condition), method = "radix")
    cutoffs <- .normalise_resolved_cutoffs(cutoffs, conditions)
    for (condition in conditions) {
        rows <- statistics$condition == condition
        has_missing <- any(statistics$missing_count[rows] > 0L)
        if (!has_missing) {
            cutoffs[[condition]] <- NA_real_
        } else if (is.na(cutoffs[[condition]])) {
            stop(
                "Condition `",
                condition,
                "` has incomplete features but no resolved cutoff.",
                call. = FALSE
            )
        }
    }

    unname(cutoffs[statistics$condition])
}

.reconcile_condition_states <- function(classifications, feature_status) {
    .validate_reconciliation_inputs(classifications, feature_status)

    features <- feature_status$feature[is.na(feature_status$retained)]
    conditions <- sort(unique(classifications$condition), method = "radix")
    state_matrix <- matrix(
        NA_character_,
        nrow = length(features),
        ncol = length(conditions),
        dimnames = list(features, conditions)
    )
    state_matrix[cbind(
        match(classifications$feature, features),
        match(classifications$condition, conditions)
    )] <- classifications$state

    decision <- .reconciled_feature_decision(state_matrix, conditions)
    retained <- decision$retained
    drop_reason <- decision$drop_reason

    status_rows <- match(features, feature_status$feature)
    feature_status$retained[status_rows] <- retained
    feature_status$drop_reason[status_rows] <- drop_reason
    classification_features <- match(classifications$feature, features)
    classifications$retained <- retained[classification_features]
    classifications$drop_reason <- drop_reason[classification_features]

    if (anyNA(classifications$retained)) {
        stop(
            "Reconciliation produced an incomplete feature audit.",
            call. = FALSE
        )
    }

    list(
        classifications = classifications,
        feature_status = feature_status
    )
}

.reconciled_feature_decision <- function(state_matrix, conditions) {
    insufficient <- state_matrix == "insufficient"
    has_insufficient <- rowSums(insufficient) > 0L
    all_mnar <- rowSums(state_matrix == "MNAR") == ncol(state_matrix)
    retained <- !has_insufficient & !all_mnar
    drop_reason <- rep(NA_character_, nrow(state_matrix))
    drop_reason[all_mnar] <- "MNAR_all_conditions"
    drop_reason[has_insufficient] <- apply(
        insufficient[has_insufficient, , drop = FALSE],
        1L,
        function(flag) {
            paste0(
                "insufficient:",
                paste(conditions[flag], collapse = ",")
            )
        }
    )

    list(retained = retained, drop_reason = drop_reason)
}

.retained_condition_groups <- function(classifications, feature_order) {
    .validate_reconciled_classifications(classifications, feature_order)
    conditions <- sort(unique(classifications$condition), method = "radix")
    retained_features <- feature_order[
        feature_order %in%
            unique(classifications$feature[classifications$retained])
    ]

    groups <- lapply(
        conditions,
        function(condition) {
            condition_rows <- classifications$condition == condition &
                classifications$retained
            condition_states <- stats::setNames(
                classifications$state[condition_rows],
                classifications$feature[condition_rows]
            )
            states <- unname(condition_states[retained_features])

            list(
                MNAR = retained_features[states == "MNAR"],
                MAR = retained_features[states == "MAR"],
                complete = retained_features[states == "complete"],
                MAR_or_complete = retained_features[
                    states %in% c("MAR", "complete")
                ]
            )
        }
    )

    stats::setNames(groups, conditions)
}

.validate_reconciliation_inputs <- function(classifications, feature_status) {
    .validate_classification_states(classifications)
    required_classification <- c(
        "feature",
        "condition",
        "state",
        "retained",
        "drop_reason"
    )
    required_status <- c("feature", "retained", "drop_reason")
    valid_status <- .valid_reconciliation_status(
        feature_status,
        required_status
    )
    valid_classifications <- .valid_reconciliation_classifications(
        classifications,
        required_classification
    )
    if (!valid_status || !valid_classifications) {
        stop(
            "Reconciliation inputs do not satisfy the audit schemas.",
            call. = FALSE
        )
    }

    .validate_reconciliation_coverage(classifications, feature_status)
    invisible(classifications)
}

.validate_reconciliation_coverage <- function(
    classifications,
    feature_status
) {
    all_missing <- !is.na(feature_status$drop_reason) &
        feature_status$drop_reason == "all_missing"
    valid_initial_status <-
        all(!feature_status$retained[all_missing]) &&
        all(
            is.na(feature_status$retained[!all_missing]) &
                is.na(feature_status$drop_reason[!all_missing])
        )
    classified_features <- unique(classifications$feature)
    expected_features <- feature_status$feature[!all_missing]
    if (!isTRUE(valid_initial_status) ||
        !setequal(classified_features, expected_features) ||
        any(classifications$feature %in% feature_status$feature[all_missing])) {
        stop(
            "Reconciliation must preserve all-missing precedence and ",
            "classify every survivor.",
            call. = FALSE
        )
    }

    conditions <- sort(unique(classifications$condition), method = "radix")
    coverage <- table(
        factor(classifications$feature, levels = expected_features),
        factor(classifications$condition, levels = conditions)
    )
    if (length(conditions) == 0L || any(coverage != 1L)) {
        stop(
            "Every surviving feature must have exactly one state per ",
            "condition.",
            call. = FALSE
        )
    }

    invisible(classifications)
}

.valid_reconciliation_status <- function(feature_status, required) {
    is.data.frame(feature_status) &&
        identical(names(feature_status), required) &&
        is.character(feature_status$feature) &&
        is.logical(feature_status$retained) &&
        is.character(feature_status$drop_reason) &&
        !anyNA(feature_status$feature) &&
        all(nzchar(feature_status$feature)) &&
        !anyDuplicated(feature_status$feature)
}

.valid_reconciliation_classifications <- function(classifications, required) {
    all(required %in% names(classifications)) &&
        is.character(classifications$feature) &&
        is.character(classifications$condition) &&
        is.logical(classifications$retained) &&
        is.character(classifications$drop_reason) &&
        !anyNA(classifications$feature) &&
        !anyNA(classifications$condition) &&
        !anyNA(classifications$state) &&
        all(nzchar(classifications$feature)) &&
        all(nzchar(classifications$condition)) &&
        !anyDuplicated(
            classifications[, c("feature", "condition"), drop = FALSE]
        )
}

.validate_reconciled_classifications <- function(
    classifications,
    feature_order
) {
    .validate_classification_states(classifications)
    valid_order <- is.character(feature_order) &&
        is.null(dim(feature_order)) &&
        !anyNA(feature_order) &&
        all(nzchar(feature_order)) &&
        !anyDuplicated(feature_order)
    required <- c("feature", "condition", "state", "retained", "drop_reason")
    valid_classifications <- all(required %in% names(classifications)) &&
        is.logical(classifications$retained) &&
        !anyNA(classifications$retained) &&
        all(classifications$feature %in% feature_order)
    if (!valid_order || !valid_classifications) {
        stop(
            "Reconciled classifications do not satisfy the group schema.",
            call. = FALSE
        )
    }

    retained <- unique(classifications$feature[classifications$retained])
    conditions <- unique(classifications$condition)
    retained_rows <- classifications[classifications$retained, , drop = FALSE]
    coverage <- table(
        factor(retained_rows$feature, levels = retained),
        factor(retained_rows$condition, levels = conditions)
    )
    valid_states <- retained_rows$state %in% c("complete", "MAR", "MNAR")
    if (any(coverage != 1L) || any(!valid_states)) {
        stop(
            "Every retained feature must have one eligible state per ",
            "condition.",
            call. = FALSE
        )
    }

    invisible(classifications)
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
            "Resolved cutoffs must contain one finite-or-NA value per ",
            "condition.",
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
        stop(
            "`statistics` does not satisfy the state-input schema.",
            call. = FALSE
        )
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
