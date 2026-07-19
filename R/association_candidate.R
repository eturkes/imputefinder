.ASSOCIATION_CANDIDATE_ARTIFACT_SCHEMA <- "association_candidate_artifact_v1"
.ASSOCIATION_QUASIBINOMIAL_CANDIDATE <- "a_fraction_quasibinomial"
.ASSOCIATION_CANDIDATES <- c(
    "a_fraction_ols_hc3_cr2",
    "a_fraction_freedman_lane",
    .ASSOCIATION_QUASIBINOMIAL_CANDIDATE
)
.ASSOCIATION_CANDIDATE_ARTIFACT_FIELDS <- c(
    "schema", "protocol", "candidate", "input_sha256", "response",
    "hypotheses", "support", "outcomes", "multiplicity", "diagnostics"
)
.ASSOCIATION_CANDIDATE_HYPOTHESIS_FIELDS <- c(
    "hypothesis", "stratum", "coefficient", "label", "term_id", "term",
    "kind", "components", "component_encodings", "role", "eligible",
    "estimable", "family_member"
)
.ASSOCIATION_AVAILABLE_OUTCOME_FIELDS <- c(
    "status", "quantity", "candidate", "stratum", "hypothesis", "effect",
    "standard_error", "conf_low", "conf_high", "statistic", "reference_df",
    "raw_p", "adjusted_p", "flag", "log_odds",
    "log_odds_standard_error", "log_odds_conf_low", "log_odds_conf_high",
    "permutation_count", "diagnostics"
)
.ASSOCIATION_OUTCOME_DIAGNOSTIC_FIELDS <- c(
    "variance_method", "rank", "residual_df", "leverage_max",
    "satterthwaite_df", "dispersion", "converged", "boundary",
    "allowable_transformations", "evaluated_transformations",
    "exceedance_count", "permutation_mode"
)
.ASSOCIATION_MULTIPLICITY_FIELDS <- c(
    "stratum", "method", "level", "family_size", "available_count",
    "unavailable_count"
)
.ASSOCIATION_CANDIDATE_DIAGNOSTIC_FIELDS <- c(
    "input_sha256", "strata", "seed_manifest", "warnings"
)
.ASSOCIATION_STRATUM_DIAGNOSTIC_FIELDS <- c(
    "stratum", "acquisition", "sample_count", "rank", "coefficient_count",
    "residual_df", "testable_count", "unavailable_count"
)
.ASSOCIATION_SEED_MANIFEST_FIELDS <- c(
    "protocol_id", "candidate_id", "acquisition", "hypothesis_id",
    "draw_id", "seed", "seed_nonce", "map_retry", "map_sha256"
)
.ASSOCIATION_UNAVAILABLE_CODES <- c(
    "association_no_observable_features",
    "association_nonestimable",
    "association_low_independent_support",
    "association_low_block_support",
    "association_low_interaction_support",
    "association_incompatible_unit_positions",
    "association_low_residual_df",
    "association_low_reference_df",
    "association_low_permutation_resolution",
    "association_incompatible_permutation",
    "association_degenerate_response",
    "association_singular_covariance",
    "association_quasibinomial_boundary",
    "association_quasibinomial_nonconvergence",
    "association_quasibinomial_paired_scope",
    "association_numerical_failure",
    "association_no_testable_hypotheses"
)
.ASSOCIATION_SUPPORT_CODES <- c(
    "association_low_independent_support",
    "association_low_block_support",
    "association_low_interaction_support",
    "association_incompatible_unit_positions"
)
.ASSOCIATION_ROBUST_UNAVAILABLE_CODES <- c(
    "association_nonestimable",
    .ASSOCIATION_SUPPORT_CODES,
    "association_low_residual_df",
    "association_low_reference_df",
    "association_degenerate_response",
    "association_singular_covariance",
    "association_numerical_failure"
)
.ASSOCIATION_ROBUST_FIT_UNAVAILABLE_CODES <- c(
    "association_low_reference_df",
    "association_singular_covariance",
    "association_numerical_failure"
)
.ASSOCIATION_FREEDMAN_LANE_UNAVAILABLE_CODES <- c(
    .ASSOCIATION_ROBUST_UNAVAILABLE_CODES,
    "association_low_permutation_resolution",
    "association_incompatible_permutation"
)
.ASSOCIATION_FREEDMAN_LANE_FIT_UNAVAILABLE_CODES <- c(
    .ASSOCIATION_ROBUST_FIT_UNAVAILABLE_CODES,
    "association_low_permutation_resolution",
    "association_incompatible_permutation"
)
.ASSOCIATION_QUASIBINOMIAL_UNAVAILABLE_CODES <- c(
    "association_nonestimable",
    .ASSOCIATION_SUPPORT_CODES,
    "association_low_residual_df",
    "association_degenerate_response",
    "association_singular_covariance",
    "association_quasibinomial_boundary",
    "association_quasibinomial_nonconvergence",
    "association_quasibinomial_paired_scope",
    "association_numerical_failure"
)
.ASSOCIATION_QUASIBINOMIAL_FIT_UNAVAILABLE_CODES <- c(
    "association_singular_covariance",
    "association_quasibinomial_boundary",
    "association_quasibinomial_nonconvergence",
    "association_quasibinomial_paired_scope",
    "association_numerical_failure"
)

.association_matrix_tolerance <- function(x) {
    if (!is.matrix(x) || !is.numeric(x) || !length(x) ||
        any(!is.finite(x))) {
        return(NA_real_)
    }
    max(dim(x)) * .Machine$double.eps * max(1, max(abs(x)))
}

.association_clean_psd <- function(x) {
    if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x)) {
        return(NULL)
    }
    x <- (x + t(x)) / 2
    tolerance <- .association_matrix_tolerance(x)
    if (!is.finite(tolerance)) {
        return(NULL)
    }
    decomposition <- tryCatch(
        eigen(x, symmetric = TRUE),
        error = function(error) NULL
    )
    if (is.null(decomposition) ||
        any(!is.finite(decomposition$values)) ||
        any(decomposition$values < -tolerance)) {
        return(NULL)
    }
    values <- decomposition$values
    selected <- abs(values) <= tolerance
    if (!any(selected)) {
        return(list(matrix = x, tolerance = tolerance))
    }
    vectors <- decomposition$vectors[, selected, drop = FALSE]
    correction <- sweep(vectors, 2L, values[selected], `*`) %*% t(vectors)
    cleaned <- x - correction
    cleaned <- (cleaned + t(cleaned)) / 2
    dimnames(cleaned) <- dimnames(x)
    list(matrix = cleaned, tolerance = tolerance)
}

.association_candidate_failure <- function(code) {
    list(ok = FALSE, code = code)
}

.abort_association_candidate_artifact <- function(message) {
    .abort_association(
        message,
        "imputefinder_association_artifact_error"
    )
}

.association_scalar_character <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.association_sha256 <- function(x) {
    .association_scalar_character(x) && grepl("^[0-9a-f]{64}$", x)
}

.association_double_scalar <- function(x, finite = TRUE) {
    is.double(x) && length(x) == 1L &&
        if (finite) is.finite(x) else is.na(x)
}

.association_integer_scalar <- function(x, minimum = NULL, na = FALSE) {
    valid <- is.integer(x) && length(x) == 1L
    if (!valid) {
        return(FALSE)
    }
    if (na) {
        return(is.na(x))
    }
    !is.na(x) && (is.null(minimum) || x >= minimum)
}

.association_numeric_close <- function(x, y) {
    is.finite(x) && is.finite(y) &&
        abs(x - y) <= sqrt(.Machine$double.eps) *
            max(1, abs(x), abs(y))
}

.empty_association_seed_manifest <- function() {
    data.frame(
        protocol_id = character(),
        candidate_id = character(),
        acquisition = character(),
        hypothesis_id = character(),
        draw_id = integer(),
        seed = integer(),
        seed_nonce = integer(),
        map_retry = integer(),
        map_sha256 = character(),
        stringsAsFactors = FALSE
    )[.ASSOCIATION_SEED_MANIFEST_FIELDS]
}

.association_unavailable_specification <- function(code) {
    specifications <- list(
        association_nonestimable = c(
            "The encoded coefficient is not estimable in the rebuilt design.",
            "an estimable encoded coefficient"
        ),
        association_low_independent_support = c(
            "The encoded coefficient has too little independent-unit support.",
            "at least four independent units on each required side"
        ),
        association_low_block_support = c(
            "The encoded coefficient has too little complete-block support.",
            "at least six complete blocks on the required sides"
        ),
        association_low_interaction_support = c(
            "The encoded interaction has too little joint-cell support.",
            "the frozen interaction cell and complete-block support"
        ),
        association_incompatible_unit_positions = c(
            "The declared units do not admit the required restricted positions.",
            "compatible declared-unit positions"
        ),
        association_low_residual_df = c(
            "The full design leaves fewer than three residual degrees of freedom.",
            "at least three residual degrees of freedom"
        ),
        association_low_reference_df = c(
            "The CR2 scalar reference distribution has fewer than five degrees of freedom.",
            "CR2 Satterthwaite degrees of freedom of at least five"
        ),
        association_low_permutation_resolution = c(
            "The constrained-null design admits fewer than 20 transformations.",
            "at least 20 distinct constrained-null transformations"
        ),
        association_incompatible_permutation = c(
            "The declared design does not admit the frozen constrained-null maps.",
            "a design with compatible constrained-null permutation maps"
        ),
        association_degenerate_response = c(
            "The detection-fraction response is constant.",
            "a nonconstant detection-fraction response"
        ),
        association_singular_covariance = c(
            "The scalar association covariance is nonpositive or singular.",
            "a positive scalar association covariance"
        ),
        association_quasibinomial_boundary = c(
            "The grouped counts or fitted probabilities reach the frozen quasibinomial boundary.",
            "strictly interior counts and fitted probabilities"
        ),
        association_quasibinomial_nonconvergence = c(
            "The frozen quasibinomial fit did not converge.",
            "convergence within 100 IRLS iterations"
        ),
        association_quasibinomial_paired_scope = c(
            "The empirical quasibinomial candidate is not licensed for blocked units.",
            "independent biological units"
        ),
        association_numerical_failure = c(
            "The frozen association calculation failed its numerical checks.",
            "a finite, numerically stable association fit"
        )
    )
    specification <- specifications[[code]]
    if (is.null(specification)) {
        .abort_association_candidate_artifact(
            "No unavailable-result specification exists for this code."
        )
    }
    list(message = specification[[1L]], requires = specification[[2L]])
}

.new_association_candidate_unavailable <- function(hypothesis, code) {
    specification <- .association_unavailable_specification(code)
    .new_unavailable(
        quantity = hypothesis,
        code = code,
        message = specification$message,
        requires = specification$requires
    )
}

.association_candidate_multiplicity <- function(
    response,
    hypotheses,
    available
) {
    strata <- unique(response$stratum)
    rows <- lapply(strata, function(stratum) {
        selected <- hypotheses$stratum == stratum
        count <- as.integer(sum(available[selected]))
        data.frame(
            stratum = stratum,
            method = "Holm",
            level = 0.05,
            family_size = count,
            available_count = count,
            unavailable_count = as.integer(sum(!available[selected])),
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    })
    output <- do.call(rbind, unname(rows))
    row.names(output) <- NULL
    output[.ASSOCIATION_MULTIPLICITY_FIELDS]
}

.association_candidate_stratum_diagnostics <- function(
    preparation,
    hypotheses,
    available
) {
    rows <- lapply(preparation$strata, function(stratum) {
        selected <- hypotheses$stratum == stratum$stratum
        rank <- stratum$core$rank$rank
        data.frame(
            stratum = stratum$stratum,
            acquisition = stratum$acquisition,
            sample_count = as.integer(length(stratum$samples)),
            rank = rank,
            coefficient_count = as.integer(ncol(stratum$core$model$matrix)),
            residual_df = as.integer(length(stratum$samples) - rank),
            testable_count = as.integer(sum(available[selected])),
            unavailable_count = as.integer(sum(!available[selected])),
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    })
    output <- do.call(rbind, unname(rows))
    row.names(output) <- NULL
    output[.ASSOCIATION_STRATUM_DIAGNOSTIC_FIELDS]
}

.valid_association_candidate_hypotheses <- function(hypotheses, response) {
    valid <- is.data.frame(hypotheses) && nrow(hypotheses) > 0L &&
        identical(names(hypotheses), .ASSOCIATION_CANDIDATE_HYPOTHESIS_FIELDS) &&
        is.character(hypotheses$hypothesis) &&
        all(grepl("^a_[0-9a-f]{64}$", hypotheses$hypothesis)) &&
        !anyDuplicated(hypotheses$hypothesis) &&
        all(vapply(hypotheses[c(
            "stratum", "coefficient", "label", "term_id", "term", "kind",
            "role"
        )], is.character, logical(1L))) &&
        !anyNA(hypotheses[c(
            "stratum", "coefficient", "label", "term_id", "term", "kind",
            "role"
        )]) &&
        all(vapply(hypotheses[c(
            "stratum", "coefficient", "label", "term_id", "term", "kind",
            "role"
        )], function(x) all(nzchar(x)), logical(1L))) &&
        is.list(hypotheses$components) &&
        inherits(hypotheses$components, "AsIs") &&
        is.list(hypotheses$component_encodings) &&
        inherits(hypotheses$component_encodings, "AsIs") &&
        all(hypotheses$kind %in% c("main", "interaction")) &&
        all(hypotheses$role %in% c("condition", "nuisance", "interaction")) &&
        all((hypotheses$kind == "main") ==
            (hypotheses$role %in% c("condition", "nuisance"))) &&
        all(lengths(hypotheses$components)[hypotheses$kind == "main"] == 1L) &&
        all(lengths(hypotheses$components)[hypotheses$kind == "interaction"] >= 2L) &&
        all(grepl("^coef_[0-9]{4,}$", hypotheses$coefficient)) &&
        all(grepl("^term_[0-9]{4,}$", hypotheses$term_id)) &&
        is.logical(hypotheses$eligible) && all(hypotheses$eligible) &&
        is.logical(hypotheses$estimable) && !anyNA(hypotheses$estimable) &&
        is.logical(hypotheses$family_member) &&
        !anyNA(hypotheses$family_member) &&
        all(hypotheses$stratum %in% response$stratum) &&
        !anyDuplicated(hypotheses[c("stratum", "coefficient")])
    if (!valid) {
        return(FALSE)
    }
    component_valid <- all(vapply(seq_len(nrow(hypotheses)), function(index) {
        components <- hypotheses$components[[index]]
        encodings <- hypotheses$component_encodings[[index]]
        is.character(components) && length(components) > 0L &&
            !anyNA(components) && all(nzchar(components)) &&
            !anyDuplicated(components) &&
            identical(components, sort(components, method = "radix")) &&
            is.character(encodings) && identical(names(encodings), components) &&
            all(encodings %in% c("numeric", "treatment"))
    }, logical(1L)))
    expected_ids <- vapply(seq_len(nrow(hypotheses)), function(index) {
        .association_hypothesis_id(
            hypotheses$stratum[[index]],
            hypotheses$coefficient[[index]],
            hypotheses$label[[index]],
            hypotheses$term_id[[index]]
        )
    }, character(1L))
    strata <- unique(response$stratum)
    canonical <- do.call(order, list(
        match(hypotheses$stratum, strata),
        hypotheses$coefficient,
        method = "radix"
    ))
    component_valid && identical(hypotheses$hypothesis, expected_ids) &&
        identical(canonical, seq_len(nrow(hypotheses)))
}

.valid_association_candidate_support <- function(support, hypotheses) {
    valid <- is.data.frame(support) &&
        identical(names(support), .ASSOCIATION_SUPPORT_FIELDS) &&
        identical(support$hypothesis, hypotheses$hypothesis) &&
        all(vapply(support[c(
            "hypothesis", "design", "side_definition"
        )], is.character, logical(1L))) &&
        !anyNA(support[c("hypothesis", "design", "side_definition")]) &&
        all(nzchar(support$side_definition)) &&
        all(support$design %in% c("independent", "blocked")) &&
        all(vapply(support[c(
            "side_count_low", "side_count_high", "complete_block_count",
            "cell_count_min"
        )], is.integer, logical(1L))) &&
        all(support$side_count_low >= 0L) &&
        all(support$side_count_high >= 0L) &&
        all(is.na(support$complete_block_count) |
            support$complete_block_count >= 0L) &&
        all(is.na(support$cell_count_min) | support$cell_count_min >= 0L) &&
        all(is.na(support$complete_block_count[
            support$design == "independent"
        ])) &&
        all(!is.na(support$complete_block_count[
            support$design == "blocked"
        ])) &&
        all(is.na(support$cell_count_min[hypotheses$kind == "main"])) &&
        all(!is.na(support$cell_count_min[
            hypotheses$kind == "interaction"
        ])) &&
        is.list(support$numeric_components) &&
        inherits(support$numeric_components, "AsIs") &&
        is.logical(support$eligible) && !anyNA(support$eligible) &&
        is.character(support$code) &&
        identical(is.na(support$code), support$eligible) &&
        all(is.na(support$code) | support$code %in% .ASSOCIATION_SUPPORT_CODES)
    if (!valid) {
        return(FALSE)
    }
    numeric_valid <- all(vapply(seq_len(nrow(support)), function(index) {
        numeric <- support$numeric_components[[index]]
        expected <- names(hypotheses$component_encodings[[index]])[
            hypotheses$component_encodings[[index]] == "numeric"
        ]
        is.data.frame(numeric) && identical(
            names(numeric),
            c(
                "component", "median", "minimum", "maximum",
                "extrapolates_0_1"
            )
        ) && is.character(numeric$component) &&
            is.double(numeric$median) && is.double(numeric$minimum) &&
            is.double(numeric$maximum) &&
            is.logical(numeric$extrapolates_0_1) && !anyNA(numeric) &&
            all(nzchar(numeric$component)) &&
            !anyDuplicated(numeric$component) &&
            identical(numeric$component, expected) &&
            all(is.finite(numeric$median)) &&
            all(is.finite(numeric$minimum)) &&
            all(is.finite(numeric$maximum)) &&
            all(numeric$minimum <= numeric$median) &&
            all(numeric$median <= numeric$maximum) &&
            identical(
                numeric$extrapolates_0_1,
                numeric$minimum > 0 | numeric$maximum < 1
            )
    }, logical(1L)))
    expected_eligible <- vapply(seq_len(nrow(support)), function(index) {
        if (identical(
            support$code[[index]],
            "association_incompatible_unit_positions"
        )) {
            return(FALSE)
        }
        blocked <- identical(support$design[[index]], "blocked")
        interaction <- identical(hypotheses$kind[[index]], "interaction")
        if (interaction) {
            if (blocked) {
                support$cell_count_min[[index]] >= 6L &&
                    support$complete_block_count[[index]] >= 6L
            } else {
                support$cell_count_min[[index]] >= 4L
            }
        } else if (blocked) {
            support$side_count_low[[index]] >= 6L &&
                support$side_count_high[[index]] >= 6L &&
                support$complete_block_count[[index]] >= 6L
        } else {
            support$side_count_low[[index]] >= 4L &&
                support$side_count_high[[index]] >= 4L
        }
    }, logical(1L))
    expected_code <- vapply(seq_len(nrow(support)), function(index) {
        if (expected_eligible[[index]]) {
            return(NA_character_)
        }
        if (identical(
            support$code[[index]],
            "association_incompatible_unit_positions"
        )) {
            return("association_incompatible_unit_positions")
        }
        if (identical(hypotheses$kind[[index]], "interaction")) {
            "association_low_interaction_support"
        } else if (identical(support$design[[index]], "blocked")) {
            "association_low_block_support"
        } else {
            "association_low_independent_support"
        }
    }, character(1L))
    numeric_valid && identical(support$eligible, expected_eligible) &&
        identical(support$code, expected_code)
}

.valid_association_robust_diagnostics <- function(
    diagnostics,
    outcome,
    design,
    candidate
) {
    valid <- is.list(diagnostics) &&
        identical(names(diagnostics), .ASSOCIATION_OUTCOME_DIAGNOSTIC_FIELDS) &&
        .association_scalar_character(diagnostics$variance_method) &&
        diagnostics$variance_method %in% c("HC3", "CR2") &&
        .association_integer_scalar(diagnostics$rank, 1L) &&
        .association_integer_scalar(diagnostics$residual_df, 3L) &&
        .association_double_scalar(diagnostics$leverage_max) &&
        diagnostics$leverage_max >= 0 && diagnostics$leverage_max <= 1 &&
        .association_double_scalar(diagnostics$dispersion, FALSE) &&
        is.logical(diagnostics$converged) &&
        length(diagnostics$converged) == 1L && is.na(diagnostics$converged) &&
        is.logical(diagnostics$boundary) &&
        length(diagnostics$boundary) == 1L && is.na(diagnostics$boundary) &&
        .association_scalar_character(diagnostics$permutation_mode)
    if (!valid) {
        return(FALSE)
    }
    reference_valid <- if (identical(design, "independent")) {
        identical(diagnostics$variance_method, "HC3") &&
            .association_double_scalar(diagnostics$satterthwaite_df, FALSE) &&
            identical(outcome$reference_df, as.double(diagnostics$residual_df))
    } else {
        identical(diagnostics$variance_method, "CR2") &&
            .association_double_scalar(diagnostics$satterthwaite_df) &&
            diagnostics$satterthwaite_df >= 5 &&
            .association_numeric_close(
                outcome$reference_df,
                diagnostics$satterthwaite_df
            )
    }
    if (!reference_valid) {
        return(FALSE)
    }
    freedman_lane <- identical(
        candidate,
        .ASSOCIATION_FREEDMAN_LANE_CANDIDATE
    )
    if (!freedman_lane) {
        return(
            .association_double_scalar(
                diagnostics$allowable_transformations,
                FALSE
            ) && .association_integer_scalar(
                diagnostics$evaluated_transformations,
                na = TRUE
            ) && .association_integer_scalar(
                diagnostics$exceedance_count,
                na = TRUE
            ) && identical(diagnostics$permutation_mode, "none") &&
                .association_integer_scalar(
                    outcome$permutation_count,
                    na = TRUE
                )
        )
    }
    allowable <- diagnostics$allowable_transformations
    valid <- .association_double_scalar(allowable) && allowable >= 20 &&
        allowable == floor(allowable) &&
        .association_integer_scalar(
            diagnostics$evaluated_transformations,
            1L
        ) && .association_integer_scalar(
            diagnostics$exceedance_count,
            0L
        ) && .association_integer_scalar(outcome$permutation_count, 20L)
    if (!valid) {
        return(FALSE)
    }
    if (allowable <= .ASSOCIATION_EXACT_PERMUTATION_MAX) {
        valid <- identical(diagnostics$permutation_mode, "exact") &&
            diagnostics$evaluated_transformations == allowable &&
            outcome$permutation_count ==
                diagnostics$evaluated_transformations &&
            diagnostics$exceedance_count >= 1L &&
            diagnostics$exceedance_count <=
                diagnostics$evaluated_transformations
        expected_p <- diagnostics$exceedance_count /
            diagnostics$evaluated_transformations
    } else {
        valid <- identical(
            diagnostics$permutation_mode,
            "monte_carlo"
        ) && diagnostics$evaluated_transformations ==
            .ASSOCIATION_MONTE_CARLO_DRAWS &&
            outcome$permutation_count == 10000L &&
            diagnostics$exceedance_count <=
                .ASSOCIATION_MONTE_CARLO_DRAWS
        expected_p <- (1 + diagnostics$exceedance_count) / 10000
    }
    valid && identical(outcome$raw_p, as.double(expected_p))
}

.valid_association_quasibinomial_diagnostics <- function(
    diagnostics,
    outcome,
    design
) {
    is.list(diagnostics) &&
        identical(names(diagnostics), .ASSOCIATION_OUTCOME_DIAGNOSTIC_FIELDS) &&
        identical(design, "independent") &&
        identical(diagnostics$variance_method, "quasibinomial") &&
        .association_integer_scalar(diagnostics$rank, 1L) &&
        .association_integer_scalar(diagnostics$residual_df, 3L) &&
        .association_double_scalar(diagnostics$leverage_max, FALSE) &&
        .association_double_scalar(diagnostics$satterthwaite_df, FALSE) &&
        .association_double_scalar(diagnostics$dispersion) &&
        diagnostics$dispersion > 0 &&
        identical(diagnostics$converged, TRUE) &&
        identical(diagnostics$boundary, FALSE) &&
        .association_double_scalar(
            diagnostics$allowable_transformations,
            FALSE
        ) &&
        .association_integer_scalar(
            diagnostics$evaluated_transformations,
            na = TRUE
        ) &&
        .association_integer_scalar(
            diagnostics$exceedance_count,
            na = TRUE
        ) &&
        identical(diagnostics$permutation_mode, "none") &&
        .association_integer_scalar(outcome$permutation_count, na = TRUE) &&
        identical(
            outcome$reference_df,
            as.double(diagnostics$residual_df)
        )
}

.valid_association_available_outcome <- function(
    outcome,
    key,
    candidate,
    hypothesis,
    support
) {
    valid <- is.list(outcome) &&
        identical(class(outcome), "imputefinder_association") &&
        identical(names(outcome), .ASSOCIATION_AVAILABLE_OUTCOME_FIELDS) &&
        identical(outcome$status, "available") &&
        identical(outcome$quantity, key) &&
        identical(outcome$candidate, candidate) &&
        identical(outcome$stratum, hypothesis$stratum) &&
        identical(outcome$hypothesis, key) &&
        all(vapply(outcome[c(
            "effect", "standard_error", "conf_low", "conf_high", "statistic",
            "reference_df", "raw_p", "adjusted_p"
        )], .association_double_scalar, logical(1L))) &&
        outcome$standard_error > 0 && outcome$conf_low <= outcome$effect &&
        outcome$conf_high >= outcome$effect && outcome$reference_df >= 3 &&
        outcome$raw_p >= 0 && outcome$raw_p <= 1 &&
        outcome$adjusted_p >= 0 && outcome$adjusted_p <= 1 &&
        is.logical(outcome$flag) && length(outcome$flag) == 1L &&
        !is.na(outcome$flag) &&
        identical(outcome$flag, outcome$adjusted_p <= 0.05)
    robust <- identical(candidate, "a_fraction_ols_hc3_cr2")
    freedman_lane <- identical(
        candidate,
        .ASSOCIATION_FREEDMAN_LANE_CANDIDATE
    )
    quasibinomial <- identical(
        candidate,
        .ASSOCIATION_QUASIBINOMIAL_CANDIDATE
    )
    if (!valid || (!robust && !freedman_lane && !quasibinomial)) {
        return(FALSE)
    }
    secondary <- unlist(outcome[c(
        "log_odds", "log_odds_standard_error", "log_odds_conf_low",
        "log_odds_conf_high"
    )], use.names = FALSE)
    if (!is.double(secondary)) {
        return(FALSE)
    }
    critical <- stats::qt(0.975, outcome$reference_df)
    statistic_valid <- if (freedman_lane) {
        outcome$statistic >= 0 && .association_numeric_close(
            outcome$statistic,
            (outcome$effect / outcome$standard_error)^2
        )
    } else if (quasibinomial) {
        .association_numeric_close(
            outcome$statistic,
            outcome$log_odds / outcome$log_odds_standard_error
        ) && .association_numeric_close(
            outcome$raw_p,
            2 * stats::pt(-abs(outcome$statistic), outcome$reference_df)
        )
    } else {
        .association_numeric_close(
            outcome$statistic,
            outcome$effect / outcome$standard_error
        ) && .association_numeric_close(
            outcome$raw_p,
            2 * stats::pt(-abs(outcome$statistic), outcome$reference_df)
        )
    }
    secondary_valid <- if (quasibinomial) {
        all(is.finite(secondary)) && outcome$log_odds_standard_error > 0 &&
            .association_numeric_close(
                outcome$log_odds_conf_low,
                outcome$log_odds - critical *
                    outcome$log_odds_standard_error
            ) && .association_numeric_close(
                outcome$log_odds_conf_high,
                outcome$log_odds + critical *
                    outcome$log_odds_standard_error
            )
    } else {
        all(is.na(secondary))
    }
    diagnostics_valid <- if (quasibinomial) {
        .valid_association_quasibinomial_diagnostics(
            outcome$diagnostics,
            outcome,
            support$design
        )
    } else {
        .valid_association_robust_diagnostics(
            outcome$diagnostics,
            outcome,
            support$design,
            candidate
        )
    }
    secondary_valid && statistic_valid &&
        .association_numeric_close(
            outcome$conf_low,
            outcome$effect - critical * outcome$standard_error
        ) &&
        .association_numeric_close(
            outcome$conf_high,
            outcome$effect + critical * outcome$standard_error
        ) && diagnostics_valid
}

.valid_association_candidate_unavailable <- function(outcome, key, candidate) {
    allowed <- if (identical(candidate, "a_fraction_ols_hc3_cr2")) {
        .ASSOCIATION_ROBUST_UNAVAILABLE_CODES
    } else if (identical(
        candidate,
        .ASSOCIATION_FREEDMAN_LANE_CANDIDATE
    )) {
        .ASSOCIATION_FREEDMAN_LANE_UNAVAILABLE_CODES
    } else if (identical(
        candidate,
        .ASSOCIATION_QUASIBINOMIAL_CANDIDATE
    )) {
        .ASSOCIATION_QUASIBINOMIAL_UNAVAILABLE_CODES
    } else {
        character()
    }
    valid <- is.list(outcome) &&
        identical(class(outcome), "imputefinder_unavailable") &&
        identical(names(outcome), .UNAVAILABLE_FIELDS) &&
        identical(outcome$status, "unavailable") &&
        identical(outcome$quantity, key) &&
        .association_scalar_character(outcome$code) &&
        outcome$code %in% allowed &&
        .association_scalar_character(outcome$message) &&
        is.character(outcome$requires) && length(outcome$requires) > 0L &&
        !anyNA(outcome$requires) && all(nzchar(outcome$requires)) &&
        !anyDuplicated(outcome$requires)
    valid && identical(
        outcome,
        .new_association_candidate_unavailable(key, outcome$code)
    )
}

.valid_association_candidate_multiplicity <- function(
    multiplicity,
    response,
    hypotheses,
    available
) {
    strata <- unique(response$stratum)
    valid <- is.data.frame(multiplicity) &&
        identical(names(multiplicity), .ASSOCIATION_MULTIPLICITY_FIELDS) &&
        identical(multiplicity$stratum, strata) && !anyNA(multiplicity) &&
        all(multiplicity$method == "Holm") &&
        all(multiplicity$level == 0.05) &&
        all(vapply(multiplicity[c(
            "family_size", "available_count", "unavailable_count"
        )], is.integer, logical(1L))) &&
        all(multiplicity$family_size >= 0L) &&
        all(multiplicity$available_count >= 0L) &&
        all(multiplicity$unavailable_count >= 0L)
    if (!valid) {
        return(FALSE)
    }
    all(vapply(seq_along(strata), function(index) {
        selected <- hypotheses$stratum == strata[[index]]
        count <- as.integer(sum(available[selected]))
        multiplicity$family_size[[index]] == count &&
            multiplicity$available_count[[index]] == count &&
            multiplicity$unavailable_count[[index]] ==
                as.integer(sum(!available[selected]))
    }, logical(1L)))
}

.valid_association_freedman_lane_seed_manifest <- function(seed, artifact) {
    monte_carlo <- vapply(artifact$outcomes, function(outcome) {
        identical(class(outcome), "imputefinder_association") &&
            identical(outcome$diagnostics$permutation_mode, "monte_carlo")
    }, logical(1L))
    monte_carlo_ids <- artifact$hypotheses$hypothesis[monte_carlo]
    if (identical(seed, .empty_association_seed_manifest())) {
        return(!length(monte_carlo_ids))
    }
    valid <- is.data.frame(seed) && nrow(seed) > 0L &&
        identical(names(seed), .ASSOCIATION_SEED_MANIFEST_FIELDS) &&
        all(vapply(seed[c(
            "protocol_id", "candidate_id", "acquisition", "hypothesis_id",
            "map_sha256"
        )], is.character, logical(1L))) &&
        all(vapply(seed[c(
            "draw_id", "seed", "seed_nonce", "map_retry"
        )], is.integer, logical(1L))) && !anyNA(seed) &&
        all(seed$protocol_id == .ASSOCIATION_PROTOCOL_ID) &&
        all(seed$candidate_id == .ASSOCIATION_FREEDMAN_LANE_CANDIDATE) &&
        all(seed$hypothesis_id %in% monte_carlo_ids) &&
        all(seed$draw_id > 0L) && all(seed$seed > 0L) &&
        all(seed$seed_nonce >= 0L) && all(seed$map_retry >= 0L) &&
        all(seed$map_retry <= .ASSOCIATION_MAP_RETRY_LIMIT) &&
        all(grepl("^[0-9a-f]{64}$", seed$map_sha256)) &&
        !anyDuplicated(seed$seed)
    if (!valid) {
        return(FALSE)
    }
    ordered <- do.call(order, c(
        unname(seed[c(
            "candidate_id", "acquisition", "hypothesis_id", "draw_id"
        )]),
        list(method = "radix")
    ))
    if (!identical(ordered, seq_len(nrow(seed)))) {
        return(FALSE)
    }
    instances <- unique(seed[c(
        "candidate_id", "acquisition", "hypothesis_id"
    )])
    row.names(instances) <- NULL
    if (!identical(
        sort(instances$hypothesis_id, method = "radix"),
        sort(monte_carlo_ids, method = "radix")
    )) {
        return(FALSE)
    }
    all(vapply(seq_len(nrow(instances)), function(index) {
        instance <- instances[index, , drop = FALSE]
        selected <- seed$candidate_id == instance$candidate_id &
            seed$acquisition == instance$acquisition &
            seed$hypothesis_id == instance$hypothesis_id
        rows <- seed[selected, , drop = FALSE]
        hypothesis_index <- match(
            instance$hypothesis_id,
            artifact$hypotheses$hypothesis
        )
        stratum <- artifact$hypotheses$stratum[[hypothesis_index]]
        acquisition <- unique(
            artifact$response$acquisition[
                artifact$response$stratum == stratum
            ]
        )
        identical(rows$draw_id, seq_len(.ASSOCIATION_MONTE_CARLO_DRAWS)) &&
            !anyDuplicated(rows$map_sha256) &&
            identical(instance$acquisition, acquisition)
    }, logical(1L)))
}

.valid_association_candidate_diagnostics <- function(
    diagnostics,
    artifact,
    available
) {
    if (!is.list(diagnostics)) {
        return(FALSE)
    }
    strata <- unique(artifact$response$stratum)
    records <- diagnostics$strata
    seed <- diagnostics$seed_manifest
    seed_valid <- if (identical(
        artifact$candidate,
        "a_fraction_ols_hc3_cr2"
    ) || identical(
        artifact$candidate,
        .ASSOCIATION_QUASIBINOMIAL_CANDIDATE
    )) {
        identical(seed, .empty_association_seed_manifest())
    } else if (identical(
        artifact$candidate,
        .ASSOCIATION_FREEDMAN_LANE_CANDIDATE
    )) {
        .valid_association_freedman_lane_seed_manifest(seed, artifact)
    } else {
        FALSE
    }
    valid <- is.list(diagnostics) &&
        identical(names(diagnostics), .ASSOCIATION_CANDIDATE_DIAGNOSTIC_FIELDS) &&
        identical(diagnostics$input_sha256, artifact$input_sha256) &&
        is.data.frame(records) &&
        identical(names(records), .ASSOCIATION_STRATUM_DIAGNOSTIC_FIELDS) &&
        identical(records$stratum, strata) &&
        all(vapply(records[c("stratum", "acquisition")],
            is.character, logical(1L))) &&
        all(vapply(records[c(
            "sample_count", "rank", "coefficient_count", "residual_df",
            "testable_count", "unavailable_count"
        )], is.integer, logical(1L))) &&
        !anyNA(records) && all(records$sample_count > 0L) &&
        all(records$rank > 0L) && all(records$coefficient_count > 0L) &&
        all(records$residual_df == records$sample_count - records$rank) &&
        all(records$testable_count >= 0L) &&
        all(records$unavailable_count >= 0L) &&
        seed_valid &&
        is.character(diagnostics$warnings) &&
        !anyNA(diagnostics$warnings) && all(nzchar(diagnostics$warnings)) &&
        !anyDuplicated(diagnostics$warnings) &&
        identical(
            diagnostics$warnings,
            sort(diagnostics$warnings, method = "radix")
        )
    if (!valid) {
        return(FALSE)
    }
    expected_warnings <- sort(unique(unname(vapply(
        artifact$outcomes[!available],
        `[[`,
        character(1L),
        "code"
    ))), method = "radix")
    if (!identical(diagnostics$warnings, expected_warnings)) {
        return(FALSE)
    }
    all(vapply(seq_along(strata), function(index) {
        stratum <- strata[[index]]
        response_selected <- artifact$response$stratum == stratum
        hypothesis_selected <- artifact$hypotheses$stratum == stratum
        available_selected <- which(hypothesis_selected & available)
        outcome_diagnostics <- lapply(
            artifact$outcomes[available_selected],
            `[[`,
            "diagnostics"
        )
        outcome_join <- !length(outcome_diagnostics) ||
            all(vapply(outcome_diagnostics, function(outcome) {
                identical(outcome$rank, records$rank[[index]]) &&
                    identical(
                        outcome$residual_df,
                        records$residual_df[[index]]
                    )
            }, logical(1L)))
        identical(
            records$acquisition[[index]],
            unique(artifact$response$acquisition[response_selected])
        ) && records$sample_count[[index]] == sum(response_selected) &&
            records$coefficient_count[[index]] >= sum(hypothesis_selected) &&
            records$testable_count[[index]] ==
                sum(available[hypothesis_selected]) &&
            records$unavailable_count[[index]] ==
                sum(!available[hypothesis_selected]) && outcome_join
    }, logical(1L)))
}

.valid_association_freedman_lane_provenance <- function(
    artifact,
    preparation,
    available
) {
    records <- tryCatch(
        .association_freedman_lane_records(preparation),
        error = function(error) NULL
    )
    if (is.null(records) || length(records) != nrow(artifact$hypotheses)) {
        return(FALSE)
    }
    seeds <- tryCatch(
        .association_freedman_lane_seeds(
            records,
            preparation$input_sha256
        ),
        error = function(error) NULL
    )
    if (is.null(seeds)) {
        return(FALSE)
    }
    seed_manifests <- vector("list", length(records))
    for (index in seq_along(records)) {
        record <- records[[index]]
        outcome <- artifact$outcomes[[index]]
        if (!is.null(record$code)) {
            if (available[[index]] || !identical(
                outcome$code,
                record$code
            )) {
                return(FALSE)
            }
            next
        }
        generated <- NULL
        if (identical(record$permutations$mode, "monte_carlo")) {
            generated <- tryCatch(
                .association_freedman_lane_generate_maps(record, seeds),
                error = function(error) NULL
            )
            if (is.null(generated)) {
                return(FALSE)
            }
            if (!generated$ok) {
                if (available[[index]] || !identical(
                    outcome$code,
                    "association_numerical_failure"
                )) {
                    return(FALSE)
                }
                next
            }
            record$maps <- generated$maps
        }
        evaluated <- .association_freedman_lane_evaluate(record)
        if (!evaluated$ok) {
            if (available[[index]] || !identical(
                outcome$code,
                "association_numerical_failure"
            )) {
                return(FALSE)
            }
            next
        }
        if (!available[[index]]) {
            return(FALSE)
        }
        diagnostics <- outcome$diagnostics
        observed <- evaluated$result
        valid <- all(vapply(c(
            "effect", "standard_error", "conf_low", "conf_high",
            "reference_df"
        ), function(field) {
            .association_numeric_close(outcome[[field]], observed[[field]])
        }, logical(1L))) &&
            identical(outcome$statistic, evaluated$statistic) &&
            identical(outcome$raw_p, evaluated$raw_p) &&
            identical(
                outcome$permutation_count,
                evaluated$permutation_count
            ) && identical(
                diagnostics$variance_method,
                observed$diagnostics$variance_method
            ) && identical(diagnostics$rank, observed$diagnostics$rank) &&
            identical(
                diagnostics$residual_df,
                observed$diagnostics$residual_df
            ) && .association_numeric_close(
                diagnostics$leverage_max,
                observed$diagnostics$leverage_max
            ) && identical(
                is.na(diagnostics$satterthwaite_df),
                is.na(observed$diagnostics$satterthwaite_df)
            ) && (is.na(diagnostics$satterthwaite_df) ||
                .association_numeric_close(
                    diagnostics$satterthwaite_df,
                    observed$diagnostics$satterthwaite_df
                )) && identical(
                diagnostics$allowable_transformations,
                evaluated$allowable_transformations
            ) && identical(
                diagnostics$evaluated_transformations,
                evaluated$evaluated_transformations
            ) && identical(
                diagnostics$exceedance_count,
                evaluated$exceedance_count
            ) && identical(
                diagnostics$permutation_mode,
                evaluated$permutation_mode
            )
        if (!valid) {
            return(FALSE)
        }
        if (!is.null(generated)) {
            seed_manifests[[index]] <- generated$manifest
        }
    }
    expected_seed <- .association_bind_freedman_lane_seed_manifests(
        seed_manifests
    )
    identical(artifact$diagnostics$seed_manifest, expected_seed)
}

.valid_association_quasibinomial_provenance <- function(
    artifact,
    preparation,
    available
) {
    records <- tryCatch(
        .association_quasibinomial_records(preparation),
        error = function(error) NULL
    )
    if (is.null(records) || length(records) != nrow(artifact$hypotheses)) {
        return(FALSE)
    }
    comparison_fields <- setdiff(
        .ASSOCIATION_AVAILABLE_OUTCOME_FIELDS,
        c("adjusted_p", "flag")
    )
    all(vapply(seq_along(records), function(index) {
        record <- records[[index]]
        outcome <- artifact$outcomes[[index]]
        if (!is.null(record$code)) {
            return(!available[[index]] && identical(
                outcome$code,
                record$code
            ))
        }
        if (!available[[index]] || is.null(record$result) ||
            !isTRUE(record$result$ok)) {
            return(FALSE)
        }
        expected <- .new_association_quasibinomial_outcome(
            record$result,
            record$hypothesis
        )
        identical(
            outcome[comparison_fields],
            expected[comparison_fields]
        )
    }, logical(1L)))
}

.validate_association_candidate_artifact <- function(
    artifact,
    preparation
) {
    if (missing(preparation)) {
        .abort_association_candidate_artifact(
            "Association candidate validation requires its preparation."
        )
    }
    .validate_association_preparation(preparation)
    valid <- is.list(artifact) &&
        identical(
            class(artifact),
            "imputefinder_association_candidate_artifact"
        ) &&
        identical(names(artifact), .ASSOCIATION_CANDIDATE_ARTIFACT_FIELDS) &&
        identical(artifact$schema, .ASSOCIATION_CANDIDATE_ARTIFACT_SCHEMA) &&
        identical(artifact$protocol, .association_protocol()) &&
        .association_scalar_character(artifact$candidate) &&
        artifact$candidate %in% .ASSOCIATION_CANDIDATES &&
        .association_sha256(artifact$input_sha256) &&
        .valid_association_response(artifact$response)
    if (!valid || !.valid_association_candidate_hypotheses(
        artifact$hypotheses,
        artifact$response
    ) || !.valid_association_candidate_support(
        artifact$support,
        artifact$hypotheses
    )) {
        .abort_association_candidate_artifact(
            "Stored association candidate header or joins are malformed."
        )
    }
    expected_hypotheses <- artifact$hypotheses[
        .ASSOCIATION_HYPOTHESIS_FIELDS
    ]
    provenance_valid <- identical(
        artifact$protocol,
        preparation$protocol
    ) && identical(
        artifact$input_sha256,
        preparation$input_sha256
    ) && identical(
        artifact$response,
        preparation$response
    ) && identical(
        expected_hypotheses,
        preparation$hypotheses
    ) && identical(artifact$support, preparation$support)
    if (!provenance_valid) {
        .abort_association_candidate_artifact(
            "Stored association candidate provenance is detached."
        )
    }
    outcomes <- artifact$outcomes
    if (!is.list(outcomes) ||
        !identical(names(outcomes), artifact$hypotheses$hypothesis)) {
        .abort_association_candidate_artifact(
            "Stored association candidate outcome keys are malformed."
        )
    }
    available <- unname(vapply(outcomes, function(outcome) {
        identical(class(outcome), "imputefinder_association")
    }, logical(1L)))
    outcome_valid <- vapply(seq_along(outcomes), function(index) {
        outcome <- outcomes[[index]]
        key <- names(outcomes)[[index]]
        if (available[[index]]) {
            .valid_association_available_outcome(
                outcome,
                key,
                artifact$candidate,
                artifact$hypotheses[index, , drop = FALSE],
                artifact$support[index, , drop = FALSE]
            )
        } else {
            .valid_association_candidate_unavailable(
                outcome,
                key,
                artifact$candidate
            )
        }
    }, logical(1L))
    if (!all(outcome_valid) ||
        !identical(artifact$hypotheses$family_member, available) ||
        any(available & (!artifact$hypotheses$estimable |
            !artifact$support$eligible))) {
        .abort_association_candidate_artifact(
            "Stored association candidate outcomes are malformed."
        )
    }
    if (identical(
        artifact$candidate,
        .ASSOCIATION_FREEDMAN_LANE_CANDIDATE
    ) && !.valid_association_freedman_lane_provenance(
        artifact,
        preparation,
        available
    )) {
        .abort_association_candidate_artifact(
            "Stored Freedman-Lane provenance or disposition is malformed."
        )
    }
    if (identical(
        artifact$candidate,
        .ASSOCIATION_QUASIBINOMIAL_CANDIDATE
    ) && !.valid_association_quasibinomial_provenance(
        artifact,
        preparation,
        available
    )) {
        .abort_association_candidate_artifact(
            "Stored quasibinomial provenance or disposition is malformed."
        )
    }
    quasibinomial <- identical(
        artifact$candidate,
        .ASSOCIATION_QUASIBINOMIAL_CANDIDATE
    )
    for (index in which(!available)) {
        stratum <- preparation$strata[[
            artifact$hypotheses$stratum[[index]]
        ]]
        residual_df <- length(stratum$samples) - stratum$core$rank$rank
        degenerate <- length(unique(unname(stratum$response))) == 1L
        response_rows <- artifact$response$stratum ==
            artifact$hypotheses$stratum[[index]]
        boundary <- any(
            artifact$response$detected_count[response_rows] <= 0L |
                artifact$response$detected_count[response_rows] >=
                    artifact$response$globally_observable_count[response_rows]
        )
        expected <- if (!artifact$hypotheses$estimable[[index]]) {
            "association_nonestimable"
        } else if (!artifact$support$eligible[[index]]) {
            artifact$support$code[[index]]
        } else if (residual_df < 3L) {
            "association_low_residual_df"
        } else if (quasibinomial && identical(
            artifact$support$design[[index]],
            "blocked"
        )) {
            "association_quasibinomial_paired_scope"
        } else if (quasibinomial && boundary) {
            "association_quasibinomial_boundary"
        } else if (degenerate) {
            "association_degenerate_response"
        } else {
            NULL
        }
        if (!is.null(expected) &&
            !identical(outcomes[[index]]$code, expected)) {
            .abort_association_candidate_artifact(
                "Stored association unavailable precedence is malformed."
            )
        }
        fit_codes <- if (quasibinomial) {
            .ASSOCIATION_QUASIBINOMIAL_FIT_UNAVAILABLE_CODES
        } else if (identical(
            artifact$candidate,
            .ASSOCIATION_FREEDMAN_LANE_CANDIDATE
        )) {
            .ASSOCIATION_FREEDMAN_LANE_FIT_UNAVAILABLE_CODES
        } else {
            .ASSOCIATION_ROBUST_FIT_UNAVAILABLE_CODES
        }
        if (is.null(expected) &&
            !outcomes[[index]]$code %in% fit_codes) {
            .abort_association_candidate_artifact(
                "Stored association fit-stage disposition is malformed."
            )
        }
        if (is.null(expected) &&
            identical(outcomes[[index]]$code, "association_low_reference_df") &&
            !identical(artifact$support$design[[index]], "blocked")) {
            .abort_association_candidate_artifact(
                "Stored association reference-df scope is malformed."
            )
        }
    }
    for (stratum in unique(artifact$hypotheses$stratum)) {
        selected <- artifact$hypotheses$stratum == stratum & available
        if (any(selected)) {
            raw <- vapply(outcomes[selected], `[[`, numeric(1L), "raw_p")
            adjusted <- vapply(
                outcomes[selected],
                `[[`,
                numeric(1L),
                "adjusted_p"
            )
            if (!identical(adjusted, stats::p.adjust(raw, method = "holm"))) {
                .abort_association_candidate_artifact(
                    "Stored association Holm arithmetic is malformed."
                )
            }
        }
    }
    diagnostic_valid <- .valid_association_candidate_multiplicity(
        artifact$multiplicity,
        artifact$response,
        artifact$hypotheses,
        available
    ) && .valid_association_candidate_diagnostics(
        artifact$diagnostics,
        artifact,
        available
    )
    if (!diagnostic_valid) {
        .abort_association_candidate_artifact(
            "Stored association multiplicity or diagnostics are malformed."
        )
    }
    expected_strata <- do.call(rbind, lapply(
        preparation$strata,
        function(stratum) {
            rank <- stratum$core$rank$rank
            data.frame(
                stratum = stratum$stratum,
                acquisition = stratum$acquisition,
                sample_count = as.integer(length(stratum$samples)),
                rank = rank,
                coefficient_count = as.integer(ncol(
                    stratum$core$model$matrix
                )),
                residual_df = as.integer(length(stratum$samples) - rank),
                stringsAsFactors = FALSE
            )
        }
    ))
    row.names(expected_strata) <- NULL
    provenance_valid <- provenance_valid && identical(
        artifact$diagnostics$strata[c(
            "stratum", "acquisition", "sample_count", "rank",
            "coefficient_count", "residual_df"
        )],
        expected_strata
    )
    if (!provenance_valid) {
        .abort_association_candidate_artifact(
            "Stored association candidate provenance is detached."
        )
    }
    invisible(artifact)
}
