.BASE_FIT_COMPATIBILITY_SCHEMA <- "base_fit_compatibility_v1"

.classic_compatibility_tolerance <- function() {
    tolerance <- sqrt(.Machine$double.eps)
    c(absolute = tolerance, relative = tolerance)
}

.validate_base_fit_result <- function(base_fit) {
    tryCatch(
        .validate_classic_result_shape(base_fit),
        error = function(condition) {
            .abort_sidecar(
                "`base_fit` must be a complete imputefinder result.",
                "imputefinder_base_fit_schema_error",
                condition_class = unname(class(condition)),
                condition_message = conditionMessage(condition)
            )
        }
    )
    invisible(base_fit)
}

.base_fit_source <- function(diagnostic) {
    if (is.list(diagnostic) && "source" %in% names(diagnostic)) {
        .sidecar_scalar_character(diagnostic$source)
    } else {
        NA_character_
    }
}

.valid_base_fit_cutoff <- function(source, cutoff) {
    valid_value <- is.numeric(cutoff) &&
        length(cutoff) == 1L &&
        !is.nan(cutoff) &&
        !is.infinite(cutoff)
    if (!valid_value) {
        return(FALSE)
    }
    if (identical(source, "not_needed")) {
        is.na(cutoff)
    } else {
        is.finite(cutoff)
    }
}

.base_fit_cutoff_policy <- function(base_fit) {
    .validate_base_fit_result(base_fit)
    conditions <- names(base_fit$cutoffs)
    sources <- vapply(
        base_fit$cutoff_diagnostics,
        .base_fit_source,
        character(1L)
    )
    allowed <- c("manual", "automatic", "not_needed")
    for (condition in conditions) {
        source <- unname(sources[[condition]])
        cutoff <- unname(base_fit$cutoffs[[condition]])
        if (!source %in% allowed ||
            !.valid_base_fit_cutoff(source, cutoff)) {
            .abort_sidecar(
                "`base_fit` has an invalid cutoff-source policy.",
                "imputefinder_base_fit_policy_error",
                condition = condition,
                source = source
            )
        }
    }
    manual_names <- conditions[sources[conditions] == "manual"]
    manual <- stats::setNames(
        as.numeric(base_fit$cutoffs[manual_names]),
        manual_names
    )

    list(sources = sources[conditions], manual = manual)
}

.without_finite_doubles <- function(value) {
    if (is.double(value)) {
        value[is.finite(value)] <- 0
        return(value)
    }
    if (is.list(value)) {
        for (index in seq_along(value)) {
            value[[index]] <- .without_finite_doubles(value[[index]])
        }
    }
    value
}

.classic_matrix_fingerprint <- function(data) {
    if (!is.matrix(data) || !is.numeric(data)) {
        return(NULL)
    }
    tryCatch(
        .matrix_fingerprint(data),
        error = function(condition) NULL
    )
}

.classic_data_identity <- function(data) {
    if (is.matrix(data)) {
        return(list(
            object = data,
            numeric_fingerprint = .classic_matrix_fingerprint(data)
        ))
    }
    assays <- SummarizedExperiment::assays(data, withDimnames = TRUE)
    fingerprints <- lapply(assays, .classic_matrix_fingerprint)

    list(object = data, numeric_fingerprints = fingerprints)
}

.classic_exact_components <- function(fit) {
    list(
        schema = list(class = class(fit), fields = names(fit)),
        data = .classic_data_identity(fit$data),
        classifications = .without_finite_doubles(fit$classifications),
        groups = fit$groups,
        feature_status = fit$feature_status,
        seed_log = fit$seed_log,
        groups_by_sample = fit$groups_by_sample,
        call = fit$call
    )
}

.classic_canonical_components <- function(fit) {
    list(
        cutoffs = .without_finite_doubles(fit$cutoffs),
        cutoff_diagnostics = .without_finite_doubles(
            fit$cutoff_diagnostics
        ),
        profiles = .without_finite_doubles(fit$profiles)
    )
}

.classic_tolerance_components <- function(fit) {
    list(
        classifications = fit$classifications,
        cutoffs = fit$cutoffs,
        cutoff_diagnostics = fit$cutoff_diagnostics,
        profiles = fit$profiles
    )
}

.different_components <- function(candidate, reference) {
    components <- union(names(candidate), names(reference))
    different <- vapply(
        components,
        function(component) {
            !identical(candidate[[component]], reference[[component]])
        },
        logical(1L)
    )
    components[different]
}

.double_vectors_equal <- function(candidate, reference, tolerance) {
    if (!is.double(candidate) || !is.double(reference) ||
        !identical(attributes(candidate), attributes(reference)) ||
        length(candidate) != length(reference)) {
        return(FALSE)
    }
    same_special <- identical(is.na(candidate), is.na(reference)) &&
        identical(is.nan(candidate), is.nan(reference)) &&
        identical(is.infinite(candidate), is.infinite(reference))
    if (!same_special) {
        return(FALSE)
    }
    infinite <- is.infinite(candidate)
    if (any(infinite) &&
        !identical(candidate[infinite], reference[infinite])) {
        return(FALSE)
    }
    finite <- is.finite(candidate)
    bound <- tolerance[["absolute"]] +
        tolerance[["relative"]] * pmax(
            abs(candidate[finite]),
            abs(reference[finite])
        )
    all(abs(candidate[finite] - reference[finite]) <= bound)
}

.numeric_trees_equal <- function(candidate, reference, tolerance) {
    if (is.double(candidate) || is.double(reference)) {
        return(.double_vectors_equal(candidate, reference, tolerance))
    }
    if (is.list(candidate) || is.list(reference)) {
        if (!is.list(candidate) || !is.list(reference) ||
            length(candidate) != length(reference)) {
            return(FALSE)
        }
        equal <- mapply(
            .numeric_trees_equal,
            candidate,
            reference,
            MoreArgs = list(tolerance = tolerance),
            SIMPLIFY = TRUE,
            USE.NAMES = FALSE
        )
        return(all(equal))
    }
    TRUE
}

.different_numeric_components <- function(candidate, reference, tolerance) {
    components <- union(names(candidate), names(reference))
    different <- vapply(
        components,
        function(component) {
            candidate_value <- candidate[[component]]
            reference_value <- reference[[component]]
            same_structure <- identical(
                .without_finite_doubles(candidate_value),
                .without_finite_doubles(reference_value)
            )
            if (!same_structure) {
                return(FALSE)
            }
            !.numeric_trees_equal(
                candidate_value,
                reference_value,
                tolerance
            )
        },
        logical(1L)
    )
    components[different]
}

.new_base_fit_report <- function(
    exact = character(),
    canonical = character(),
    tolerance = character()
) {
    compatible <- !length(c(exact, canonical, tolerance))
    list(
        schema = .BASE_FIT_COMPATIBILITY_SCHEMA,
        compatible = compatible,
        exact = unname(exact),
        canonical = unname(canonical),
        tolerance = unname(tolerance),
        numeric_tolerance = .classic_compatibility_tolerance()
    )
}

.compare_classic_fits <- function(candidate, reference) {
    .validate_base_fit_result(candidate)
    .validate_base_fit_result(reference)
    numeric_tolerance <- .classic_compatibility_tolerance()

    .new_base_fit_report(
        exact = .different_components(
            .classic_exact_components(candidate),
            .classic_exact_components(reference)
        ),
        canonical = .different_components(
            .classic_canonical_components(candidate),
            .classic_canonical_components(reference)
        ),
        tolerance = .different_numeric_components(
            .classic_tolerance_components(candidate),
            .classic_tolerance_components(reference),
            numeric_tolerance
        )
    )
}

.base_fit_differences <- function(report) {
    list(
        exact = report$exact,
        canonical = report$canonical,
        tolerance = report$tolerance
    )
}

.abort_base_fit_mismatch <- function(report, recomputed = NULL) {
    .abort_sidecar(
        "`base_fit` is incompatible with the recomputed classic path.",
        "imputefinder_base_fit_mismatch_error",
        report = report,
        differences = .base_fit_differences(report),
        recomputed = recomputed
    )
}

.base_fit_recomputation_cutoffs <- function(policy) {
    if (length(policy$manual) == 0L) NULL else policy$manual
}

.accept_compatible_base_fit <- function(
    base_fit,
    prepared,
    original,
    seed
) {
    policy <- .base_fit_cutoff_policy(base_fit)
    recomputed <- .attempt_input_first_classic(
        prepared = prepared,
        original = original,
        call = base_fit$call,
        cutoffs = .base_fit_recomputation_cutoffs(policy),
        seed = seed
    )
    if (inherits(recomputed, "imputefinder_classic_failure")) {
        report <- .new_base_fit_report(exact = "recomputation")
        .abort_base_fit_mismatch(report, recomputed)
    }
    report <- .compare_classic_fits(base_fit, recomputed)
    if (!report$compatible) {
        .abort_base_fit_mismatch(report)
    }

    base_fit
}
