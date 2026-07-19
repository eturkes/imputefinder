.association_quasibinomial_counts <- function(stratum, response) {
    selected <- match(stratum$samples, response$sample)
    if (anyNA(selected)) {
        return(NULL)
    }
    rows <- response[selected, , drop = FALSE]
    totals <- unique(rows$globally_observable_count)
    valid <- length(totals) == 1L && is.integer(totals) && totals > 0L &&
        is.integer(rows$detected_count) &&
        identical(rows$sample, stratum$samples) &&
        identical(
            unname(rows$detection_fraction),
            unname(stratum$response)
        )
    if (!valid) {
        return(NULL)
    }
    list(
        detected = rows$detected_count,
        total = totals[[1L]]
    )
}

.association_quasibinomial_stats_glm_fit <- function(...) {
    stats::glm.fit(...)
}

.association_quasibinomial_glm_fit <- function(z, detected, total) {
    suppressWarnings(.association_quasibinomial_stats_glm_fit(
        x = z,
        y = cbind(detected, total - detected),
        family = stats::quasibinomial(),
        start = NULL,
        etastart = NULL,
        control = stats::glm.control(
            epsilon = 1e-10,
            maxit = 100L,
            trace = FALSE
        )
    ))
}

.association_quasibinomial_weighted_svd <- function(weighted_design) {
    .design_svd(weighted_design)
}

.association_quasibinomial_fit <- function(stratum, response) {
    counts <- .association_quasibinomial_counts(stratum, response)
    if (is.null(counts)) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    if (any(counts$detected <= 0L | counts$detected >= counts$total)) {
        return(.association_candidate_failure(
            "association_quasibinomial_boundary"
        ))
    }
    x <- stratum$core$model$matrix
    decomposition <- tryCatch(
        .design_svd(x),
        error = function(error) NULL
    )
    rank <- stratum$core$rank$rank
    n <- nrow(x)
    residual_df <- as.integer(n - rank)
    if (is.null(decomposition) || decomposition$rank != rank || rank < 1L ||
        any(!is.finite(x))) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    u <- decomposition$decomposition$u[, seq_len(rank), drop = FALSE]
    basis <- decomposition$decomposition$v[
        ,
        seq_len(rank),
        drop = FALSE
    ]
    singular_values <- decomposition$decomposition$d[seq_len(rank)]
    z <- sweep(u, 2L, singular_values, `*`)
    if (any(!is.finite(z))) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    fit <- tryCatch(
        .association_quasibinomial_glm_fit(
            z,
            counts$detected,
            counts$total
        ),
        error = function(error) NULL
    )
    if (is.null(fit) || !is.logical(fit$converged) ||
        length(fit$converged) != 1L || is.na(fit$converged) ||
        !is.logical(fit$boundary) || length(fit$boundary) != 1L ||
        is.na(fit$boundary)) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    if (!fit$converged) {
        return(.association_candidate_failure(
            "association_quasibinomial_nonconvergence"
        ))
    }
    if (fit$boundary) {
        return(.association_candidate_failure(
            "association_quasibinomial_boundary"
        ))
    }
    coefficients <- as.double(fit$coefficients)
    mu <- as.double(fit$fitted.values)
    valid_fit <- length(coefficients) == rank && length(mu) == n &&
        all(is.finite(coefficients)) && all(is.finite(mu))
    if (!valid_fit) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    probability_floor <- sqrt(.Machine$double.eps)
    if (any(mu < probability_floor | mu > 1 - probability_floor)) {
        return(.association_candidate_failure(
            "association_quasibinomial_boundary"
        ))
    }
    weights <- as.double(counts$total * mu * (1 - mu))
    if (any(!is.finite(weights)) || any(weights <= 0)) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    weighted_design <- sweep(z, 1L, sqrt(weights), `*`)
    weighted <- tryCatch(
        .association_quasibinomial_weighted_svd(weighted_design),
        error = function(error) NULL
    )
    if (is.null(weighted) || weighted$rank != rank) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    weighted_values <- weighted$decomposition$d[seq_len(rank)]
    weighted_basis <- weighted$decomposition$v[
        ,
        seq_len(rank),
        drop = FALSE
    ]
    inverse_root <- sweep(weighted_basis, 2L, weighted_values, `/`)
    information_inverse <- tcrossprod(inverse_root)
    dispersion <- sum(
        (counts$detected - counts$total * mu)^2 / weights
    ) / residual_df
    if (!is.finite(dispersion) || dispersion <= 0 ||
        any(!is.finite(information_inverse))) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    covariance <- .association_clean_psd(dispersion * information_inverse)
    if (is.null(covariance)) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    list(
        ok = TRUE,
        rank = rank,
        residual_df = residual_df,
        basis = basis,
        z = z,
        coefficients = coefficients,
        covariance = covariance$matrix,
        covariance_tolerance = covariance$tolerance,
        fitted_probability = mu,
        dispersion = as.double(dispersion),
        converged = TRUE,
        boundary = FALSE
    )
}

.association_quasibinomial_contrast <- function(fit, hypothesis, stratum) {
    columns <- colnames(stratum$core$model$matrix)
    column <- match(hypothesis$coefficient, columns)
    if (!is.list(fit) || !isTRUE(fit$ok) || is.na(column)) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    axis <- numeric(length(columns))
    axis[[column]] <- 1
    contrast <- as.vector(crossprod(fit$basis, axis))
    observed_column <- stratum$core$model$matrix[, column]
    z0 <- fit$z - tcrossprod(observed_column, contrast)
    z1 <- sweep(z0, 2L, contrast, `+`)
    mu0 <- stats::plogis(as.vector(z0 %*% fit$coefficients))
    mu1 <- stats::plogis(as.vector(z1 %*% fit$coefficients))
    effect <- mean(mu1 - mu0)
    gradient <- colMeans(
        sweep(z1, 1L, mu1 * (1 - mu1), `*`) -
            sweep(z0, 1L, mu0 * (1 - mu0), `*`)
    )
    log_odds <- as.double(crossprod(contrast, fit$coefficients))
    effect_variance <- as.double(crossprod(
        gradient,
        fit$covariance %*% gradient
    ))
    log_odds_variance <- as.double(crossprod(
        contrast,
        fit$covariance %*% contrast
    ))
    values <- c(
        contrast, z0, z1, mu0, mu1, effect, gradient, log_odds,
        effect_variance, log_odds_variance
    )
    if (any(!is.finite(values))) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    if (effect_variance <= fit$covariance_tolerance ||
        log_odds_variance <= fit$covariance_tolerance) {
        return(.association_candidate_failure(
            "association_singular_covariance"
        ))
    }
    standard_error <- sqrt(effect_variance)
    log_odds_standard_error <- sqrt(log_odds_variance)
    reference_df <- as.double(fit$residual_df)
    statistic <- log_odds / log_odds_standard_error
    raw_p <- 2 * stats::pt(-abs(statistic), reference_df)
    critical <- stats::qt(0.975, reference_df)
    output <- c(
        standard_error, log_odds_standard_error, reference_df, statistic,
        raw_p, critical,
        effect - critical * standard_error,
        effect + critical * standard_error,
        log_odds - critical * log_odds_standard_error,
        log_odds + critical * log_odds_standard_error
    )
    if (any(!is.finite(output))) {
        return(.association_candidate_failure("association_numerical_failure"))
    }
    list(
        ok = TRUE,
        effect = as.double(effect),
        standard_error = as.double(standard_error),
        conf_low = as.double(effect - critical * standard_error),
        conf_high = as.double(effect + critical * standard_error),
        statistic = as.double(statistic),
        reference_df = reference_df,
        raw_p = as.double(raw_p),
        log_odds = log_odds,
        log_odds_standard_error = as.double(log_odds_standard_error),
        log_odds_conf_low = as.double(
            log_odds - critical * log_odds_standard_error
        ),
        log_odds_conf_high = as.double(
            log_odds + critical * log_odds_standard_error
        ),
        diagnostics = list(
            variance_method = "quasibinomial",
            rank = fit$rank,
            residual_df = fit$residual_df,
            leverage_max = NA_real_,
            satterthwaite_df = NA_real_,
            dispersion = fit$dispersion,
            converged = fit$converged,
            boundary = fit$boundary,
            allowable_transformations = NA_real_,
            evaluated_transformations = NA_integer_,
            exceedance_count = NA_integer_,
            permutation_mode = "none"
        )
    )
}

.association_quasibinomial_records <- function(preparation) {
    hypotheses <- preparation$hypotheses
    records <- vector("list", nrow(hypotheses))
    names(records) <- hypotheses$hypothesis
    for (stratum in preparation$strata) {
        selected <- which(hypotheses$stratum == stratum$stratum)
        if (!length(selected)) {
            next
        }
        support <- preparation$support[selected, , drop = FALSE]
        residual_df <- as.integer(
            length(stratum$samples) - stratum$core$rank$rank
        )
        blocked <- identical(stratum$core$units$resampling_role, "block")
        counts <- .association_quasibinomial_counts(
            stratum,
            preparation$response
        )
        boundary <- is.null(counts) || any(
            counts$detected <= 0L | counts$detected >= counts$total
        )
        degenerate <- length(unique(unname(stratum$response))) == 1L
        fit <- NULL
        for (position in seq_along(selected)) {
            index <- selected[[position]]
            hypothesis <- hypotheses[index, , drop = FALSE]
            code <- if (!hypothesis$estimable[[1L]]) {
                "association_nonestimable"
            } else if (!support$eligible[[position]]) {
                support$code[[position]]
            } else if (residual_df < 3L) {
                "association_low_residual_df"
            } else if (blocked) {
                "association_quasibinomial_paired_scope"
            } else if (boundary) {
                "association_quasibinomial_boundary"
            } else if (degenerate) {
                "association_degenerate_response"
            } else {
                NULL
            }
            record <- list(
                stratum = stratum,
                hypothesis = hypothesis,
                code = code,
                fit = NULL,
                result = NULL
            )
            if (!is.null(code)) {
                records[[index]] <- record
                next
            }
            if (is.null(fit)) {
                fit <- .association_quasibinomial_fit(
                    stratum,
                    preparation$response
                )
            }
            result <- if (fit$ok) {
                .association_quasibinomial_contrast(
                    fit,
                    hypothesis,
                    stratum
                )
            } else {
                fit
            }
            if (!result$ok) {
                record$code <- result$code
            }
            record$fit <- fit
            record$result <- result
            records[[index]] <- record
        }
    }
    records
}

.new_association_quasibinomial_outcome <- function(result, hypothesis) {
    structure(
        list(
            status = "available",
            quantity = hypothesis$hypothesis,
            candidate = .ASSOCIATION_QUASIBINOMIAL_CANDIDATE,
            stratum = hypothesis$stratum,
            hypothesis = hypothesis$hypothesis,
            effect = result$effect,
            standard_error = result$standard_error,
            conf_low = result$conf_low,
            conf_high = result$conf_high,
            statistic = result$statistic,
            reference_df = result$reference_df,
            raw_p = result$raw_p,
            adjusted_p = NA_real_,
            flag = NA,
            log_odds = result$log_odds,
            log_odds_standard_error = result$log_odds_standard_error,
            log_odds_conf_low = result$log_odds_conf_low,
            log_odds_conf_high = result$log_odds_conf_high,
            permutation_count = NA_integer_,
            diagnostics = result$diagnostics
        ),
        class = "imputefinder_association"
    )
}

.run_association_quasibinomial <- function(preparation) {
    .validate_association_preparation(preparation)
    hypotheses <- preparation$hypotheses
    records <- .association_quasibinomial_records(preparation)
    outcomes <- lapply(records, function(record) {
        if (is.null(record$code)) {
            .new_association_quasibinomial_outcome(
                record$result,
                record$hypothesis
            )
        } else {
            .new_association_candidate_unavailable(
                record$hypothesis$hypothesis,
                record$code
            )
        }
    })
    names(outcomes) <- hypotheses$hypothesis
    available <- unname(vapply(outcomes, function(outcome) {
        identical(class(outcome), "imputefinder_association")
    }, logical(1L)))
    for (stratum in unique(hypotheses$stratum)) {
        selected <- which(hypotheses$stratum == stratum & available)
        if (!length(selected)) {
            next
        }
        raw <- vapply(outcomes[selected], `[[`, numeric(1L), "raw_p")
        adjusted <- stats::p.adjust(raw, method = "holm")
        for (position in seq_along(selected)) {
            index <- selected[[position]]
            outcomes[[index]]$adjusted_p <- unname(adjusted[[position]])
            outcomes[[index]]$flag <- outcomes[[index]]$adjusted_p <= 0.05
        }
    }
    hypotheses$family_member <- available
    multiplicity <- .association_candidate_multiplicity(
        preparation$response,
        hypotheses,
        available
    )
    warnings <- sort(unique(unname(vapply(
        outcomes[!available],
        `[[`,
        character(1L),
        "code"
    ))), method = "radix")
    diagnostics <- list(
        input_sha256 = preparation$input_sha256,
        strata = .association_candidate_stratum_diagnostics(
            preparation,
            hypotheses,
            available
        ),
        seed_manifest = .empty_association_seed_manifest(),
        warnings = warnings
    )
    artifact <- structure(
        list(
            schema = .ASSOCIATION_CANDIDATE_ARTIFACT_SCHEMA,
            protocol = preparation$protocol,
            candidate = .ASSOCIATION_QUASIBINOMIAL_CANDIDATE,
            input_sha256 = preparation$input_sha256,
            response = preparation$response,
            hypotheses = hypotheses,
            support = preparation$support,
            outcomes = outcomes,
            multiplicity = multiplicity,
            diagnostics = diagnostics
        ),
        class = "imputefinder_association_candidate_artifact"
    )
    .validate_association_candidate_artifact(artifact, preparation)
    artifact
}
