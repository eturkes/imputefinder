.missingness_states <- function() {
    c("complete", "MNAR", "MAR", "insufficient")
}

.empty_classifications <- function() {
    data.frame(
        feature = character(),
        condition = character(),
        sample_count = integer(),
        observed_count = integer(),
        missing_count = integer(),
        missing_fraction = numeric(),
        mean_intensity = numeric(),
        cutoff = numeric(),
        state = character(),
        seeded = logical(),
        retained = logical(),
        drop_reason = character(),
        stringsAsFactors = FALSE
    )
}

.empty_feature_status <- function() {
    data.frame(
        feature = character(),
        retained = logical(),
        drop_reason = character(),
        stringsAsFactors = FALSE
    )
}

.empty_seed_log <- function() {
    data.frame(
        feature = character(),
        condition = character(),
        sample = character(),
        old_value = numeric(),
        inserted_value = numeric(),
        seed = integer(),
        stringsAsFactors = FALSE
    )
}

.empty_condition_groups <- function(conditions) {
    groups <- lapply(
        conditions,
        function(condition) {
            list(
                MNAR = character(),
                MAR = character(),
                complete = character(),
                MAR_or_complete = character()
            )
        }
    )
    stats::setNames(groups, conditions)
}

.new_imputefinder_result <- function(
    data,
    groups_by_sample,
    call,
    classifications = .empty_classifications(),
    groups = NULL,
    feature_status = .empty_feature_status(),
    cutoffs = NULL,
    cutoff_diagnostics = NULL,
    profiles = NULL,
    seed_log = .empty_seed_log()
) {
    groups_by_sample <- .validate_result_groups(data, groups_by_sample)
    conditions <- sort(unique(unname(groups_by_sample)), method = "radix")

    if (is.null(groups)) {
        groups <- .empty_condition_groups(conditions)
    }
    if (is.null(cutoffs)) {
        cutoffs <- stats::setNames(rep(NA_real_, length(conditions)), conditions)
    }
    if (is.null(cutoff_diagnostics)) {
        cutoff_diagnostics <- stats::setNames(
            vector("list", length(conditions)),
            conditions
        )
    }
    if (is.null(profiles)) {
        profiles <- stats::setNames(
            vector("list", length(conditions)),
            conditions
        )
    }

    .validate_classification_states(classifications)

    structure(
        list(
            data = data,
            classifications = classifications,
            groups = groups,
            feature_status = feature_status,
            cutoffs = cutoffs,
            cutoff_diagnostics = cutoff_diagnostics,
            profiles = profiles,
            seed_log = seed_log,
            groups_by_sample = groups_by_sample,
            call = call
        ),
        class = "imputefinder_result"
    )
}

.validate_result_groups <- function(data, groups_by_sample) {
    sample_names <- colnames(data)
    valid_groups <- is.atomic(groups_by_sample) &&
        is.null(dim(groups_by_sample)) &&
        length(groups_by_sample) == length(sample_names) &&
        !anyNA(groups_by_sample)
    if (!valid_groups || !identical(names(groups_by_sample), sample_names)) {
        stop(
            "`groups_by_sample` names must match output sample order.",
            call. = FALSE
        )
    }

    labels <- as.character(groups_by_sample)
    if (anyNA(labels) || any(!nzchar(labels))) {
        stop("`groups_by_sample` must contain non-empty labels.", call. = FALSE)
    }

    stats::setNames(labels, sample_names)
}

.validate_classification_states <- function(classifications) {
    if (!is.data.frame(classifications) ||
        !"state" %in% names(classifications)) {
        stop("`classifications` must contain a state column.", call. = FALSE)
    }

    states <- classifications$state[!is.na(classifications$state)]
    if (!is.character(states) ||
        any(!states %in% .missingness_states())) {
        stop(
            "Classification states must use the stable state vocabulary.",
            call. = FALSE
        )
    }

    invisible(classifications)
}
