.ASSOCIATION_ROBUST_CANDIDATE <- "a_fraction_ols_hc3_cr2"

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

.association_cr2_adjustment <- function(block) {
    if (!is.matrix(block) || !is.numeric(block) ||
        nrow(block) != ncol(block)) {
        return(NULL)
    }
    block <- (block + t(block)) / 2
    tolerance <- .association_matrix_tolerance(block)
    if (!is.finite(tolerance)) {
        return(NULL)
    }
    decomposition <- tryCatch(
        eigen(block, symmetric = TRUE),
        error = function(error) NULL
    )
    if (is.null(decomposition) ||
        any(!is.finite(decomposition$values)) ||
        any(decomposition$values < -tolerance)) {
        return(NULL)
    }
    keep <- decomposition$values > tolerance
    if (!any(keep)) {
        return(matrix(0, nrow(block), ncol(block)))
    }
    vectors <- decomposition$vectors[, keep, drop = FALSE]
    inverse_root <- sweep(
        vectors,
        2L,
        sqrt(decomposition$values[keep]),
        `/`
    )
    inverse_root %*% t(vectors)
}

.association_robust_failure <- function(code) {
    list(ok = FALSE, code = code)
}

.association_cr2_reference_df <- function(projection) {
    if (!is.matrix(projection) || !is.numeric(projection) ||
        !length(projection) || any(!is.finite(projection))) {
        return(NA_real_)
    }
    scale <- max(abs(projection))
    if (!is.finite(scale) || scale <= 0) {
        return(NA_real_)
    }
    projection <- projection / scale
    norms <- colSums(projection^2)
    denominator <- sum(crossprod(projection)^2)
    if (!is.finite(denominator) || denominator <= 0) {
        return(NA_real_)
    }
    as.double(sum(norms)^2 / denominator)
}

.association_residual_projector <- function(left_basis) {
    if (!is.matrix(left_basis) || !is.numeric(left_basis) ||
        nrow(left_basis) < ncol(left_basis) ||
        any(!is.finite(left_basis))) {
        return(NULL)
    }
    n <- nrow(left_basis)
    rank <- ncol(left_basis)
    if (rank == n) {
        return(matrix(0, nrow = n, ncol = n))
    }
    decomposition <- tryCatch(
        qr(left_basis, LAPACK = TRUE),
        error = function(error) NULL
    )
    if (is.null(decomposition)) {
        return(NULL)
    }
    complete <- tryCatch(
        qr.Q(decomposition, complete = TRUE),
        error = function(error) NULL
    )
    if (is.null(complete) || !identical(dim(complete), c(n, n)) ||
        any(!is.finite(complete))) {
        return(NULL)
    }
    complement <- complete[, seq.int(rank + 1L, n), drop = FALSE]
    output <- tcrossprod(complement)
    output <- (output + t(output)) / 2
    if (any(!is.finite(output))) NULL else output
}

.association_robust_algebra <- function(stratum) {
    x <- stratum$core$model$matrix
    y <- unname(stratum$response)
    decomposition <- .design_svd(x)
    rank <- decomposition$rank
    n <- nrow(x)
    residual_df <- as.integer(n - rank)
    if (rank < 1L || rank != stratum$core$rank$rank ||
        any(!is.finite(x)) || any(!is.finite(y))) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    u <- decomposition$decomposition$u[, seq_len(rank), drop = FALSE]
    v <- decomposition$decomposition$v[, seq_len(rank), drop = FALSE]
    singular_values <- decomposition$decomposition$d[seq_len(rank)]
    inverse_singular_values <- 1 / singular_values
    z <- sweep(
        u,
        2L,
        singular_values,
        `*`
    )
    if (any(!is.finite(z)) || any(!is.finite(inverse_singular_values))) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    h_matrix <- tcrossprod(u)
    h_tolerance <- .association_matrix_tolerance(h_matrix)
    leverage <- diag(h_matrix)
    invalid_leverage <- !is.finite(h_tolerance) ||
        any(!is.finite(leverage)) ||
        any(leverage < -h_tolerance) ||
        any(leverage > 1 + h_tolerance)
    if (invalid_leverage) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    leverage[leverage < 0] <- 0
    leverage[leverage > 1] <- 1
    blocked <- identical(stratum$core$units$resampling_role, "block")
    if (!blocked && any(1 - leverage <= h_tolerance)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    residual_maker <- if (blocked) {
        .association_residual_projector(u)
    } else {
        diag(n) - h_matrix
    }
    if (is.null(residual_maker)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    residual <- as.vector(residual_maker %*% y)
    coefficients <- inverse_singular_values * as.vector(crossprod(u, y))
    if (any(!is.finite(residual)) || any(!is.finite(coefficients))) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    adjustments <- NULL
    if (!blocked) {
        scaled_residual <- residual / (1 - leverage)
        meat <- crossprod(u * scaled_residual, u * scaled_residual)
        method <- "HC3"
    } else {
        sample_units <- stratum$core$units$sample
        selected <- match(rownames(x), sample_units$sample)
        if (anyNA(selected)) {
            return(.association_robust_failure("association_numerical_failure"))
        }
        units <- sample_units$unit[selected]
        clusters <- unique(units)
        meat <- matrix(0, nrow = rank, ncol = rank)
        adjustments <- vector("list", length(clusters))
        names(adjustments) <- clusters
        for (index in seq_along(clusters)) {
            rows <- which(units == clusters[[index]])
            adjustment <- .association_cr2_adjustment(
                residual_maker[rows, rows, drop = FALSE]
            )
            if (is.null(adjustment) || any(!is.finite(adjustment))) {
                return(.association_robust_failure(
                    "association_numerical_failure"
                ))
            }
            score <- crossprod(
                u[rows, , drop = FALSE],
                adjustment %*% residual[rows]
            )
            meat <- meat + tcrossprod(as.vector(score))
            adjustments[[index]] <- list(
                rows = rows,
                matrix = adjustment
            )
        }
        method <- "CR2"
    }
    covariance <- sweep(meat, 1L, inverse_singular_values, `*`)
    covariance <- sweep(covariance, 2L, inverse_singular_values, `*`)
    covariance <- .association_clean_psd(covariance)
    if (is.null(covariance)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    list(
        ok = TRUE,
        method = method,
        rank = rank,
        residual_df = residual_df,
        leverage_max = as.double(max(leverage)),
        leverage = leverage,
        leverage_tolerance = h_tolerance,
        basis = v,
        left_basis = u,
        inverse_singular_values = inverse_singular_values,
        z = z,
        coefficients = coefficients,
        covariance = covariance$matrix,
        covariance_tolerance = covariance$tolerance,
        residual_maker = residual_maker,
        adjustments = adjustments
    )
}

.association_robust_refit <- function(fit, response) {
    valid <- is.list(fit) && isTRUE(fit$ok) &&
        is.double(response) &&
        length(response) == nrow(fit$left_basis) &&
        all(is.finite(response))
    if (!valid) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    residual <- as.vector(fit$residual_maker %*% response)
    coefficients <- fit$inverse_singular_values *
        as.vector(crossprod(fit$left_basis, response))
    if (any(!is.finite(residual)) || any(!is.finite(coefficients))) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    if (identical(fit$method, "HC3")) {
        if (!is.double(fit$leverage) ||
            length(fit$leverage) != length(response) ||
            !is.finite(fit$leverage_tolerance) ||
            any(1 - fit$leverage <= fit$leverage_tolerance)) {
            return(.association_robust_failure(
                "association_numerical_failure"
            ))
        }
        scaled_residual <- residual / (1 - fit$leverage)
        meat <- crossprod(
            fit$left_basis * scaled_residual,
            fit$left_basis * scaled_residual
        )
    } else if (identical(fit$method, "CR2")) {
        meat <- matrix(0, nrow = fit$rank, ncol = fit$rank)
        for (adjustment in fit$adjustments) {
            rows <- adjustment$rows
            score <- crossprod(
                fit$left_basis[rows, , drop = FALSE],
                adjustment$matrix %*% residual[rows]
            )
            meat <- meat + tcrossprod(as.vector(score))
        }
    } else {
        return(.association_robust_failure("association_numerical_failure"))
    }
    covariance <- sweep(
        meat,
        1L,
        fit$inverse_singular_values,
        `*`
    )
    covariance <- sweep(
        covariance,
        2L,
        fit$inverse_singular_values,
        `*`
    )
    covariance <- .association_clean_psd(covariance)
    if (is.null(covariance)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    output <- fit
    output$coefficients <- coefficients
    output$covariance <- covariance$matrix
    output$covariance_tolerance <- covariance$tolerance
    output
}

.association_robust_contrast <- function(fit, coefficient, columns) {
    column <- match(coefficient, columns)
    if (is.na(column)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    axis <- numeric(length(columns))
    axis[[column]] <- 1
    contrast <- as.vector(crossprod(fit$basis, axis))
    effect <- as.double(crossprod(contrast, fit$coefficients))
    variance <- as.double(crossprod(
        contrast,
        fit$covariance %*% contrast
    ))
    tolerance <- fit$covariance_tolerance
    if (!is.finite(effect) || !is.finite(variance) ||
        !is.finite(tolerance)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    if (variance <= tolerance) {
        return(.association_robust_failure("association_singular_covariance"))
    }
    standard_error <- sqrt(variance)
    if (identical(fit$method, "HC3")) {
        reference_df <- as.double(fit$residual_df)
        satterthwaite_df <- NA_real_
    } else {
        projections <- lapply(fit$adjustments, function(adjustment) {
            rows <- adjustment$rows
            as.vector(
                fit$residual_maker[, rows, drop = FALSE] %*%
                    adjustment$matrix %*%
                    fit$left_basis[rows, , drop = FALSE] %*%
                    (fit$inverse_singular_values * contrast)
            )
        })
        projection <- do.call(cbind, unname(projections))
        reference_df <- .association_cr2_reference_df(projection)
        # The registered global floor is a literal comparison on the computed
        # double; it has no tolerance band.
        if (!is.finite(reference_df)) {
            return(.association_robust_failure(
                "association_numerical_failure"
            ))
        }
        if (reference_df < 5) {
            return(.association_robust_failure(
                "association_low_reference_df"
            ))
        }
        reference_df <- as.double(reference_df)
        satterthwaite_df <- reference_df
    }
    statistic <- effect / standard_error
    raw_p <- 2 * stats::pt(-abs(statistic), reference_df)
    critical <- stats::qt(0.975, reference_df)
    values <- c(
        standard_error, reference_df, statistic, raw_p, critical,
        effect - critical * standard_error,
        effect + critical * standard_error
    )
    if (any(!is.finite(values))) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    list(
        ok = TRUE,
        effect = effect,
        standard_error = as.double(standard_error),
        conf_low = as.double(effect - critical * standard_error),
        conf_high = as.double(effect + critical * standard_error),
        statistic = as.double(statistic),
        reference_df = reference_df,
        raw_p = as.double(raw_p),
        diagnostics = list(
            variance_method = fit$method,
            rank = fit$rank,
            residual_df = fit$residual_df,
            leverage_max = fit$leverage_max,
            satterthwaite_df = satterthwaite_df,
            dispersion = NA_real_,
            converged = NA,
            boundary = NA,
            allowable_transformations = NA_real_,
            evaluated_transformations = NA_integer_,
            exceedance_count = NA_integer_,
            permutation_mode = "none"
        )
    )
}

.new_association_robust_outcome <- function(result, hypothesis) {
    structure(
        list(
            status = "available",
            quantity = hypothesis$hypothesis,
            candidate = .ASSOCIATION_ROBUST_CANDIDATE,
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
            log_odds = NA_real_,
            log_odds_standard_error = NA_real_,
            log_odds_conf_low = NA_real_,
            log_odds_conf_high = NA_real_,
            permutation_count = NA_integer_,
            diagnostics = result$diagnostics
        ),
        class = "imputefinder_association"
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

.run_association_ols_hc3_cr2 <- function(preparation) {
    .validate_association_preparation(preparation)
    hypotheses <- preparation$hypotheses
    outcomes <- vector("list", nrow(hypotheses))
    names(outcomes) <- hypotheses$hypothesis
    for (stratum in preparation$strata) {
        selected <- which(hypotheses$stratum == stratum$stratum)
        if (!length(selected)) {
            next
        }
        support <- preparation$support[selected, , drop = FALSE]
        residual_df <- as.integer(
            length(stratum$samples) - stratum$core$rank$rank
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
            } else if (degenerate) {
                "association_degenerate_response"
            } else {
                NULL
            }
            if (!is.null(code)) {
                outcomes[[index]] <- .new_association_candidate_unavailable(
                    hypothesis$hypothesis,
                    code
                )
                next
            }
            if (is.null(fit)) {
                fit <- .association_robust_algebra(stratum)
            }
            result <- if (fit$ok) {
                .association_robust_contrast(
                    fit,
                    hypothesis$coefficient,
                    colnames(stratum$core$model$matrix)
                )
            } else {
                fit
            }
            outcomes[[index]] <- if (result$ok) {
                .new_association_robust_outcome(result, hypothesis)
            } else {
                .new_association_candidate_unavailable(
                    hypothesis$hypothesis,
                    result$code
                )
            }
        }
    }
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
            candidate = .ASSOCIATION_ROBUST_CANDIDATE,
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
