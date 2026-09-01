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
    groups <- if (is.null(groups)) {
        .empty_condition_groups(conditions)
    } else {
        groups
    }
    cutoffs <- .default_result_cutoffs(cutoffs, conditions)
    cutoff_diagnostics <- .default_condition_records(
        cutoff_diagnostics,
        conditions
    )
    profiles <- .default_condition_records(profiles, conditions)
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

.default_result_cutoffs <- function(cutoffs, conditions) {
    if (is.null(cutoffs)) {
        stats::setNames(rep(NA_real_, length(conditions)), conditions)
    } else {
        cutoffs
    }
}

.default_condition_records <- function(records, conditions) {
    if (is.null(records)) {
        stats::setNames(vector("list", length(conditions)), conditions)
    } else {
        records
    }
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

#' Present an ImputeFinder Result
#'
#' Print a compact overview or calculate structured feature, condition-state,
#' drop-reason, cutoff, and seed-insertion counts for an ImputeFinder result.
#'
#' @param x An \code{imputefinder_result}.
#' @param object An \code{imputefinder_result}.
#' @param ... Unused.
#'
#' @return \code{print.imputefinder_result()} returns \code{x} invisibly.
#'   \code{summary.imputefinder_result()} returns a list of class
#'   \code{summary.imputefinder_result}. Its count tables include all classified
#'   feature-condition blocks and the retained subset separately.
#'   \code{print.summary.imputefinder_result()} returns \code{x} invisibly.
#'
#' @seealso \code{\link{classify_missingness}},
#'   \code{\link{plot_missingness}}, and \code{vignette("imputefinder")}.
#' @name imputefinder_result
#' @examples
#' x <- rbind(
#'     on_off = c(NA, NA, NA, NA, 15, 16, 15, 17),
#'     mar = c(14, 15, NA, 16, 14, NA, 15, 16),
#'     low_anchor = c(8, NA, NA, NA, 14, 15, 16, 17)
#' )
#' colnames(x) <- paste0("sample", seq_len(ncol(x)))
#' group <- rep(c("A", "B"), each = 4)
#' fit <- classify_missingness(x, group, cutoffs = c(A = 12, B = 12))
#'
#' fit
#' summary(fit)
NULL

#' @rdname imputefinder_result
#' @export
print.imputefinder_result <- function(x, ...) {
    result_summary <- summary(x)
    .print_result_overview(result_summary)

    for (index in seq_len(nrow(result_summary$retained_states))) {
        states <- result_summary$retained_states[index, , drop = FALSE]
        cutoff <- result_summary$cutoffs[index, , drop = FALSE]
        cat(sprintf(
            paste0(
                "%s: retained MNAR=%d, MAR=%d, complete=%d; ",
                "cutoff=%s\n"
            ),
            states$condition,
            states$MNAR,
            states$MAR,
            states$complete,
            .format_result_cutoff(cutoff$cutoff, cutoff$source)
        ))
    }

    invisible(x)
}

#' @rdname imputefinder_result
#' @export
summary.imputefinder_result <- function(object, ...) {
    conditions <- .validate_result_for_presentation(object)
    feature_count <- nrow(object$feature_status)
    retained_count <- sum(object$feature_status$retained)

    structure(
        list(
            call = object$call,
            features = c(
                total = as.integer(feature_count),
                retained = as.integer(retained_count),
                dropped = as.integer(feature_count - retained_count)
            ),
            samples = c(
                total = as.integer(length(object$groups_by_sample)),
                conditions = as.integer(length(conditions))
            ),
            seed_insertions = as.integer(nrow(object$seed_log)),
            states = .result_state_counts(
                object$classifications,
                conditions
            ),
            retained_states = .result_state_counts(
                object$classifications,
                conditions,
                retained_only = TRUE
            ),
            drops = .result_drop_counts(object$feature_status),
            cutoffs = .result_cutoff_summary(object, conditions)
        ),
        class = "summary.imputefinder_result"
    )
}

#' @rdname imputefinder_result
#' @export
print.summary.imputefinder_result <- function(x, ...) {
    .print_result_overview(x, summary = TRUE)
    cat("States (all classified features):\n")
    for (index in seq_len(nrow(x$states))) {
        states <- x$states[index, , drop = FALSE]
        cat(sprintf(
            "%s: complete=%d, MNAR=%d, MAR=%d, insufficient=%d\n",
            states$condition,
            states$complete,
            states$MNAR,
            states$MAR,
            states$insufficient
        ))
    }

    if (nrow(x$drops) == 0L) {
        cat("Dropped: none\n")
    } else {
        drop_text <- sprintf("%s=%d", x$drops$reason, x$drops$count)
        cat("Dropped: ", paste(drop_text, collapse = "; "), "\n", sep = "")
    }

    cutoff_text <- vapply(
        seq_len(nrow(x$cutoffs)),
        function(index) {
            sprintf(
                "%s=%s",
                x$cutoffs$condition[[index]],
                .format_result_cutoff(
                    x$cutoffs$cutoff[[index]],
                    x$cutoffs$source[[index]]
                )
            )
        },
        character(1L)
    )
    cat("Cutoffs: ", paste(cutoff_text, collapse = "; "), "\n", sep = "")

    invisible(x)
}

.validate_result_for_presentation <- function(result) {
    conditions <- names(result$groups)
    required <- c(
        "classifications", "feature_status", "cutoffs",
        "cutoff_diagnostics", "seed_log", "groups_by_sample", "call"
    )
    valid <- inherits(result, "imputefinder_result") &&
        is.list(result) &&
        all(required %in% names(result)) &&
        is.character(conditions) &&
        length(conditions) > 0L &&
        identical(names(result$cutoffs), conditions) &&
        identical(names(result$cutoff_diagnostics), conditions) &&
        is.data.frame(result$classifications) &&
        is.data.frame(result$feature_status) &&
        is.data.frame(result$seed_log)
    if (!valid) {
        stop("Result does not satisfy the presentation schema.", call. = FALSE)
    }
    .validate_classification_states(result$classifications)

    conditions
}

.result_state_counts <- function(
    classifications,
    conditions,
    retained_only = FALSE
) {
    if (retained_only) {
        classifications <- classifications[
            classifications$retained,
            ,
            drop = FALSE
        ]
    }
    states <- .missingness_states()
    counts <- table(
        factor(classifications$condition, levels = conditions),
        factor(classifications$state, levels = states)
    )
    columns <- lapply(
        seq_along(states),
        function(index) as.integer(counts[, index])
    )
    names(columns) <- states

    data.frame(
        c(list(condition = conditions), columns),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

.result_drop_counts <- function(feature_status) {
    reasons <- feature_status$drop_reason[!feature_status$retained]
    if (length(reasons) == 0L) {
        return(data.frame(
            reason = character(),
            count = integer(),
            stringsAsFactors = FALSE
        ))
    }

    unique_reasons <- unique(reasons)
    insufficient <- sort(
        unique_reasons[startsWith(unique_reasons, "insufficient:")],
        method = "radix"
    )
    fixed <- c("all_missing", "MNAR_all_conditions")
    other <- sort(
        setdiff(unique_reasons, c(fixed, insufficient)),
        method = "radix"
    )
    ordered_reasons <- c(
        intersect("all_missing", unique_reasons),
        insufficient,
        intersect("MNAR_all_conditions", unique_reasons),
        other
    )
    counts <- table(factor(reasons, levels = ordered_reasons))

    data.frame(
        reason = ordered_reasons,
        count = as.integer(counts),
        stringsAsFactors = FALSE
    )
}

.result_cutoff_summary <- function(result, conditions) {
    diagnostics <- result$cutoff_diagnostics[conditions]
    sources <- vapply(
        diagnostics,
        function(diagnostic) {
            valid <- is.list(diagnostic) &&
                is.character(diagnostic$source) &&
                length(diagnostic$source) == 1L &&
                !is.na(diagnostic$source) &&
                nzchar(diagnostic$source)
            if (!valid) {
                stop(
                    "Result does not satisfy the presentation schema.",
                    call. = FALSE
                )
            }
            diagnostic$source
        },
        character(1L)
    )

    data.frame(
        condition = conditions,
        cutoff = as.numeric(result$cutoffs[conditions]),
        source = unname(sources),
        stringsAsFactors = FALSE
    )
}

.print_result_overview <- function(result_summary, summary = FALSE) {
    label <- if (summary) {
        "<summary.imputefinder_result>"
    } else {
        "<imputefinder_result>"
    }
    cat(label, "\n", sep = "")
    cat(sprintf(
        "Features: %d/%d retained; %d dropped\n",
        result_summary$features[["retained"]],
        result_summary$features[["total"]],
        result_summary$features[["dropped"]]
    ))
    condition_label <- if (result_summary$samples[["conditions"]] == 1L) {
        "condition"
    } else {
        "conditions"
    }
    cat(sprintf(
        "Samples: %d across %d %s; seed insertions: %d\n",
        result_summary$samples[["total"]],
        result_summary$samples[["conditions"]],
        condition_label,
        result_summary$seed_insertions
    ))
}

.format_result_cutoff <- function(cutoff, source) {
    if (is.na(cutoff) && identical(source, "not_needed")) {
        return("not needed")
    }
    if (is.na(cutoff)) {
        return(sprintf("unavailable (%s)", source))
    }

    formatted <- format(cutoff, digits = 7L, trim = TRUE)
    sprintf("%s (%s)", formatted, source)
}
