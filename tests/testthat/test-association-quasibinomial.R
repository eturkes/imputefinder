association_quasi_detection_matrix <- function(counts, total) {
    samples <- sprintf("sample_%02d", seq_along(counts))
    output <- matrix(
        NA_real_,
        nrow = total,
        ncol = length(counts),
        dimnames = list(sprintf("feature_%03d", seq_len(total)), samples)
    )
    for (index in seq_along(counts)) {
        if (counts[[index]] > 0L) {
            positions <- (
                seq_len(counts[[index]]) - 1L + index - 1L
            ) %% total + 1L
            output[positions, index] <- positions
        }
    }
    output
}

association_quasi_preparation <- function(counts, total, metadata, ...) {
    data <- association_quasi_detection_matrix(counts, total)
    row.names(metadata) <- colnames(data)
    imputefinder:::.new_association_preparation(
        data,
        missingness_design(metadata, ...)
    )
}

association_quasi_outcome <- function(artifact, label) {
    selected <- match(label, artifact$hypotheses$label)
    artifact$outcomes[[selected]]
}

association_quasi_direct_contrast <- function(preparation, index) {
    hypothesis <- preparation$hypotheses[index, , drop = FALSE]
    stratum <- preparation$strata[[hypothesis$stratum]]
    selected <- match(stratum$samples, preparation$response$sample)
    counts <- preparation$response$detected_count[selected]
    total <- unique(preparation$response$globally_observable_count[selected])
    x <- stratum$core$model$matrix
    fit <- suppressWarnings(stats::glm.fit(
        x = x,
        y = cbind(counts, total - counts),
        family = stats::quasibinomial(),
        start = NULL,
        etastart = NULL,
        control = stats::glm.control(
            epsilon = 1e-10,
            maxit = 100L,
            trace = FALSE
        )
    ))
    stopifnot(fit$converged, ncol(x) == stratum$core$rank$rank)
    probability <- fit$fitted.values
    weights <- total * probability * (1 - probability)
    dispersion <- sum((counts - total * probability)^2 / weights) /
        (nrow(x) - ncol(x))
    covariance <- dispersion * solve(crossprod(x * sqrt(weights)))
    column <- match(hypothesis$coefficient, colnames(x))
    x0 <- x
    x1 <- x
    x0[, column] <- 0
    x1[, column] <- 1
    delta <- function(coefficients) {
        mean(
            stats::plogis(as.vector(x1 %*% coefficients)) -
                stats::plogis(as.vector(x0 %*% coefficients))
        )
    }
    step <- 1e-6 * pmax(1, abs(fit$coefficients))
    gradient <- vapply(seq_along(step), function(position) {
        offset <- numeric(length(step))
        offset[[position]] <- step[[position]]
        (delta(fit$coefficients + offset) -
            delta(fit$coefficients - offset)) / (2 * step[[position]])
    }, numeric(1L))
    list(
        effect = delta(fit$coefficients),
        standard_error = sqrt(as.double(crossprod(
            gradient,
            covariance %*% gradient
        ))),
        log_odds = as.double(fit$coefficients[[column]]),
        log_odds_standard_error = sqrt(covariance[column, column]),
        dispersion = as.double(dispersion)
    )
}

test_that("quasibinomial GLM invocation locks every declared control", {
    z <- cbind(intercept = 1, condition = c(0, 0, 1, 1))
    detected <- c(2L, 3L, 5L, 6L)
    captured <- testthat::with_mocked_bindings(
        imputefinder:::.association_quasibinomial_glm_fit(
            z,
            detected,
            8L
        ),
        .association_quasibinomial_stats_glm_fit = function(...) list(...),
        .package = "imputefinder"
    )
    expect_identical(
        names(captured),
        c("x", "y", "family", "start", "etastart", "control")
    )
    expect_identical(captured$x, z)
    expect_identical(captured$y, cbind(detected, 8L - detected))
    expect_identical(captured$family$family, "quasibinomial")
    expect_null(captured$start)
    expect_null(captured$etastart)
    expect_identical(
        captured$control,
        stats::glm.control(epsilon = 1e-10, maxit = 100L, trace = FALSE)
    )
})

test_that("quasibinomial reproduces the frozen grouped-count projection", {
    counts <- c(2L, 4L, 5L, 7L, 8L, 10L, 6L, 9L, 11L, 12L, 14L, 16L)
    metadata <- data.frame(
        condition = rep(c("A", "B"), each = 6L),
        stringsAsFactors = FALSE
    )
    preparation <- association_quasi_preparation(
        counts,
        20L,
        metadata,
        condition = "condition"
    )
    artifact <- imputefinder:::.run_association_quasibinomial(preparation)
    outcome <- association_quasi_outcome(artifact, "condition[B]")

    expect_identical(artifact$candidate, "a_fraction_quasibinomial")
    expect_s3_class(outcome, "imputefinder_association")
    expect_equal(outcome$effect, 0.26666666666666411, tolerance = 1e-13)
    expect_equal(
        outcome$standard_error,
        0.093230748766176541,
        tolerance = 1e-13
    )
    expect_equal(outcome$conf_low, 0.058935613140974785, tolerance = 1e-13)
    expect_equal(outcome$conf_high, 0.47439772019235343, tolerance = 1e-13)
    expect_equal(outcome$log_odds, 1.1155618469818709, tolerance = 1e-13)
    expect_equal(
        outcome$log_odds_standard_error,
        0.41055811687699634,
        tolerance = 1e-13
    )
    expect_equal(
        outcome$log_odds_conf_low,
        0.20078135576991374,
        tolerance = 1e-13
    )
    expect_equal(
        outcome$log_odds_conf_high,
        2.0303423381938281,
        tolerance = 1e-13
    )
    expect_equal(outcome$statistic, 2.7171837582158789, tolerance = 1e-13)
    expect_identical(outcome$reference_df, 10.0)
    expect_equal(outcome$raw_p, 0.021664804947531889, tolerance = 1e-13)
    expect_identical(outcome$adjusted_p, outcome$raw_p)
    expect_true(outcome$flag)
    expect_identical(outcome$permutation_count, NA_integer_)
    expect_identical(outcome$diagnostics$variance_method, "quasibinomial")
    expect_equal(
        outcome$diagnostics$dispersion,
        2.2895927601809909,
        tolerance = 1e-13
    )
    expect_true(outcome$diagnostics$converged)
    expect_false(outcome$diagnostics$boundary)
    expect_identical(
        artifact$diagnostics$seed_manifest,
        imputefinder:::.empty_association_seed_manifest()
    )
    expect_invisible(
        imputefinder:::.validate_association_candidate_artifact(
            artifact,
            preparation
        )
    )
})

test_that("numeric effects toggle only their encoded column from zero to one", {
    counts <- c(3L, 5L, 6L, 8L, 9L, 11L, 12L, 14L,
        5L, 7L, 10L, 12L, 13L, 15L, 17L, 19L)
    metadata <- data.frame(
        condition = rep(c("A", "B"), each = 8L),
        dose = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9,
            0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85),
        stringsAsFactors = FALSE
    )
    preparation <- association_quasi_preparation(
        counts,
        24L,
        metadata,
        condition = "condition",
        nuisance = "dose"
    )
    artifact <- imputefinder:::.run_association_quasibinomial(preparation)
    dose_index <- match("dose", artifact$hypotheses$label)
    outcome <- artifact$outcomes[[dose_index]]
    expect_s3_class(outcome, "imputefinder_association")

    expected <- association_quasi_direct_contrast(preparation, dose_index)
    expect_equal(outcome$effect, expected$effect, tolerance = 1e-12)
    expect_equal(
        outcome$standard_error,
        expected$standard_error,
        tolerance = 1e-8
    )
    expect_equal(outcome$log_odds, expected$log_odds, tolerance = 1e-12)
    expect_equal(
        outcome$log_odds_standard_error,
        expected$log_odds_standard_error,
        tolerance = 1e-12
    )
    expect_equal(
        outcome$diagnostics$dispersion,
        expected$dispersion,
        tolerance = 1e-12
    )
    expect_false(isTRUE(all.equal(outcome$effect, outcome$log_odds)))

    available <- vapply(artifact$outcomes, inherits, logical(1L),
        "imputefinder_association")
    raw <- vapply(artifact$outcomes[available], `[[`, numeric(1L), "raw_p")
    adjusted <- vapply(
        artifact$outcomes[available],
        `[[`,
        numeric(1L),
        "adjusted_p"
    )
    expect_identical(adjusted, stats::p.adjust(raw, method = "holm"))
})

test_that("interaction delta inference matches a raw-design oracle", {
    metadata <- expand.grid(
        replicate = seq_len(4L),
        dose = c(0.2, 0.8),
        condition = c("A", "B"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    preparation <- association_quasi_preparation(
        c(3L, 5L, 7L, 9L, 6L, 8L, 10L, 12L,
            4L, 7L, 9L, 11L, 8L, 10L, 13L, 15L),
        18L,
        metadata,
        condition = "condition",
        nuisance = "dose",
        interactions = list(c("condition", "dose"))
    )
    artifact <- imputefinder:::.run_association_quasibinomial(preparation)
    index <- which(artifact$hypotheses$kind == "interaction")
    expect_length(index, 1L)
    outcome <- artifact$outcomes[[index]]
    expect_s3_class(outcome, "imputefinder_association")
    expected <- association_quasi_direct_contrast(preparation, index)
    expect_equal(outcome$effect, expected$effect, tolerance = 1e-12)
    expect_equal(
        outcome$standard_error,
        expected$standard_error,
        tolerance = 1e-8
    )
    expect_equal(outcome$log_odds, expected$log_odds, tolerance = 1e-12)
    expect_equal(
        outcome$log_odds_standard_error,
        expected$log_odds_standard_error,
        tolerance = 1e-12
    )
})

test_that("quasibinomial preserves strata and retained-rank coordinates", {
    acquisition_counts <- c(
        2L, 4L, 5L, 7L, 3L, 6L, 8L, 9L,
        4L, 6L, 7L, 10L, 5L, 8L, 11L, 13L
    )
    acquisition_metadata <- data.frame(
        acquisition = rep(c("DDA", "DIA"), each = 8L),
        condition = rep(rep(c("A", "B"), each = 4L), 2L),
        stringsAsFactors = FALSE
    )
    acquisition_preparation <- association_quasi_preparation(
        acquisition_counts,
        16L,
        acquisition_metadata,
        condition = "condition",
        acquisition = "acquisition"
    )
    acquisition_artifact <- imputefinder:::.run_association_quasibinomial(
        acquisition_preparation
    )
    expect_identical(acquisition_artifact$multiplicity$family_size, c(1L, 1L))
    expect_true(all(vapply(
        acquisition_artifact$outcomes,
        inherits,
        logical(1L),
        "imputefinder_association"
    )))

    aliased_metadata <- data.frame(
        condition = rep(c("A", "B"), each = 6L),
        nuisance_a = rep(c("x", "y"), 6L),
        nuisance_b = rep(c("x", "y"), 6L),
        stringsAsFactors = FALSE
    )
    aliased_preparation <- association_quasi_preparation(
        c(2L, 4L, 5L, 7L, 8L, 10L, 6L, 9L, 11L, 12L, 14L, 15L),
        18L,
        aliased_metadata,
        condition = "condition",
        nuisance = c("nuisance_a", "nuisance_b")
    )
    expect_false(aliased_preparation$strata[[1L]]$core$rank$full_column_rank)
    aliased_artifact <- imputefinder:::.run_association_quasibinomial(
        aliased_preparation
    )
    condition <- association_quasi_outcome(
        aliased_artifact,
        "condition[B]"
    )
    expect_s3_class(condition, "imputefinder_association")
    unavailable <- aliased_artifact$outcomes[
        !aliased_artifact$hypotheses$estimable
    ]
    expect_identical(
        unname(vapply(unavailable, `[[`, character(1L), "code")),
        rep("association_nonestimable", 2L)
    )
})

test_that("one failed contrast leaves the remaining Holm family available", {
    counts <- c(3L, 5L, 6L, 8L, 9L, 11L, 12L, 14L,
        5L, 7L, 10L, 12L, 13L, 15L, 17L, 19L)
    metadata <- data.frame(
        condition = rep(c("A", "B"), each = 8L),
        dose = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9,
            0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85),
        stringsAsFactors = FALSE
    )
    preparation <- association_quasi_preparation(
        counts,
        24L,
        metadata,
        condition = "condition",
        nuisance = "dose"
    )
    original <- imputefinder:::.association_quasibinomial_contrast
    artifact <- testthat::with_mocked_bindings(
        imputefinder:::.run_association_quasibinomial(preparation),
        .association_quasibinomial_contrast = function(
            fit,
            hypothesis,
            stratum
        ) {
            if (identical(hypothesis$label, "dose")) {
                return(imputefinder:::.association_candidate_failure(
                    "association_singular_covariance"
                ))
            }
            original(fit, hypothesis, stratum)
        },
        .package = "imputefinder"
    )
    condition <- association_quasi_outcome(artifact, "condition[B]")
    dose <- association_quasi_outcome(artifact, "dose")
    expect_s3_class(condition, "imputefinder_association")
    expect_identical(condition$adjusted_p, condition$raw_p)
    expect_identical(dose$code, "association_singular_covariance")
    expect_identical(artifact$multiplicity$family_size, 1L)
    expect_identical(
        artifact$diagnostics$warnings,
        "association_singular_covariance"
    )
    expect_identical(
        artifact$diagnostics$seed_manifest,
        imputefinder:::.empty_association_seed_manifest()
    )
})

test_that("quasibinomial applies scope, boundary, and fit dispositions", {
    independent <- data.frame(
        condition = rep(c("A", "B"), each = 4L),
        stringsAsFactors = FALSE
    )
    boundary_preparation <- association_quasi_preparation(
        c(0L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
        12L,
        independent,
        condition = "condition"
    )
    boundary_artifact <- imputefinder:::.run_association_quasibinomial(
        boundary_preparation
    )
    boundary <- association_quasi_outcome(
        boundary_artifact,
        "condition[B]"
    )
    expect_s3_class(boundary, "imputefinder_unavailable")
    expect_identical(boundary$code, "association_quasibinomial_boundary")
    expect_identical(
        boundary$message,
        "The grouped counts or fitted probabilities reach the frozen quasibinomial boundary."
    )
    expect_identical(
        boundary$requires,
        "strictly interior counts and fitted probabilities"
    )
    forged_text <- boundary_artifact
    forged_text$outcomes[[1L]]$message <- "Arbitrary replacement."
    expect_error(
        imputefinder:::.validate_association_candidate_artifact(
            forged_text,
            boundary_preparation
        ),
        class = "imputefinder_association_artifact_error"
    )

    upper_boundary_preparation <- association_quasi_preparation(
        c(12L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
        12L,
        independent,
        condition = "condition"
    )
    upper_boundary <- association_quasi_outcome(
        imputefinder:::.run_association_quasibinomial(
            upper_boundary_preparation
        ),
        "condition[B]"
    )
    expect_identical(
        upper_boundary$code,
        "association_quasibinomial_boundary"
    )

    paired_metadata <- data.frame(
        condition = rep(c("A", "B"), 6L),
        subject = rep(sprintf("p%02d", seq_len(6L)), each = 2L),
        stringsAsFactors = FALSE
    )
    paired_preparation <- association_quasi_preparation(
        c(3L, 6L, 4L, 7L, 5L, 8L, 6L, 9L, 7L, 10L, 8L, 11L),
        14L,
        paired_metadata,
        condition = "condition",
        block = "subject"
    )
    paired <- association_quasi_outcome(
        imputefinder:::.run_association_quasibinomial(paired_preparation),
        "condition[B]"
    )
    expect_s3_class(paired, "imputefinder_unavailable")
    expect_identical(paired$code, "association_quasibinomial_paired_scope")

    degenerate_preparation <- association_quasi_preparation(
        rep(5L, 8L),
        10L,
        independent,
        condition = "condition"
    )
    degenerate <- association_quasi_outcome(
        imputefinder:::.run_association_quasibinomial(
            degenerate_preparation
        ),
        "condition[B]"
    )
    expect_identical(degenerate$code, "association_degenerate_response")

    zero_dispersion_preparation <- association_quasi_preparation(
        c(rep(2L, 4L), rep(8L, 4L)),
        10L,
        independent,
        condition = "condition"
    )
    zero_dispersion <- association_quasi_outcome(
        imputefinder:::.run_association_quasibinomial(
            zero_dispersion_preparation
        ),
        "condition[B]"
    )
    expect_identical(zero_dispersion$code, "association_singular_covariance")

    unclamped_preparation <- association_quasi_preparation(
        c(3L, 6L, 4L, 4L, 1L, 1L, 8L, 16L),
        20L,
        independent,
        condition = "condition"
    )
    unclamped <- association_quasi_outcome(
        imputefinder:::.run_association_quasibinomial(
            unclamped_preparation
        ),
        "condition[B]"
    )
    expect_s3_class(unclamped, "imputefinder_association")
    expect_lt(unclamped$conf_low, 0)
})

test_that("quasibinomial fit failures retain their exact disposition", {
    counts <- c(2L, 4L, 5L, 7L, 6L, 8L, 9L, 11L)
    metadata <- data.frame(
        condition = rep(c("A", "B"), each = 4L),
        stringsAsFactors = FALSE
    )
    preparation <- association_quasi_preparation(
        counts,
        14L,
        metadata,
        condition = "condition"
    )
    original_glm <- imputefinder:::.association_quasibinomial_glm_fit
    nonconverged <- testthat::with_mocked_bindings(
        imputefinder:::.run_association_quasibinomial(preparation),
        .association_quasibinomial_glm_fit = function(...) {
            fit <- original_glm(...)
            fit$converged <- FALSE
            fit
        },
        .package = "imputefinder"
    )
    expect_identical(
        association_quasi_outcome(nonconverged, "condition[B]")$code,
        "association_quasibinomial_nonconvergence"
    )

    malformed <- testthat::with_mocked_bindings(
        imputefinder:::.run_association_quasibinomial(preparation),
        .association_quasibinomial_glm_fit = function(...) {
            fit <- original_glm(...)
            fit$coefficients[[1L]] <- NA_real_
            fit
        },
        .package = "imputefinder"
    )
    expect_identical(
        association_quasi_outcome(malformed, "condition[B]")$code,
        "association_numerical_failure"
    )

    fitted_boundary <- testthat::with_mocked_bindings(
        imputefinder:::.run_association_quasibinomial(preparation),
        .association_quasibinomial_glm_fit = function(...) {
            fit <- original_glm(...)
            fit$fitted.values[[1L]] <- 0
            fit
        },
        .package = "imputefinder"
    )
    expect_identical(
        association_quasi_outcome(fitted_boundary, "condition[B]")$code,
        "association_quasibinomial_boundary"
    )

    probability_floor <- sqrt(.Machine$double.eps)
    stratum <- preparation$strata[[1L]]
    inclusive_boundary <- testthat::with_mocked_bindings(
        imputefinder:::.association_quasibinomial_fit(
            stratum,
            preparation$response
        ),
        .association_quasibinomial_glm_fit = function(...) {
            fit <- original_glm(...)
            fit$fitted.values[[1L]] <- probability_floor
            fit$fitted.values[[2L]] <- 1 - probability_floor
            fit
        },
        .package = "imputefinder"
    )
    expect_true(inclusive_boundary$ok)

    exact_phi_preparation <- association_quasi_preparation(
        c(2L, 4L, 6L, 8L, 3L, 5L, 7L, 9L),
        16L,
        metadata,
        condition = "condition"
    )
    exact_phi_stratum <- exact_phi_preparation$strata[[1L]]
    exact_zero_phi <- testthat::with_mocked_bindings(
        imputefinder:::.association_quasibinomial_fit(
            exact_phi_stratum,
            exact_phi_preparation$response
        ),
        .association_quasibinomial_glm_fit = function(z, detected, total) {
            fit <- original_glm(z, detected, total)
            fit$fitted.values <- detected / total
            fit
        },
        .package = "imputefinder"
    )
    expect_identical(exact_zero_phi$code, "association_numerical_failure")

    original_weighted_svd <-
        imputefinder:::.association_quasibinomial_weighted_svd
    rank_loss <- testthat::with_mocked_bindings(
        imputefinder:::.run_association_quasibinomial(preparation),
        .association_quasibinomial_weighted_svd = function(x) {
            decomposition <- original_weighted_svd(x)
            decomposition$rank <- decomposition$rank - 1L
            decomposition
        },
        .package = "imputefinder"
    )
    expect_identical(
        association_quasi_outcome(rank_loss, "condition[B]")$code,
        "association_numerical_failure"
    )
})

test_that("quasibinomial artifact validation replays fit and delta arithmetic", {
    counts <- c(2L, 4L, 5L, 7L, 8L, 10L, 6L, 9L, 11L, 12L, 14L, 16L)
    metadata <- data.frame(
        condition = rep(c("A", "B"), each = 6L),
        stringsAsFactors = FALSE
    )
    preparation <- association_quasi_preparation(
        counts,
        20L,
        metadata,
        condition = "condition"
    )
    artifact <- imputefinder:::.run_association_quasibinomial(preparation)

    shifted_effect <- artifact
    shifted <- shifted_effect$outcomes[[1L]]
    shifted$effect <- shifted$effect + 0.01
    shifted$conf_low <- shifted$conf_low + 0.01
    shifted$conf_high <- shifted$conf_high + 0.01
    shifted_effect$outcomes[[1L]] <- shifted

    shifted_log_odds <- artifact
    shifted <- shifted_log_odds$outcomes[[1L]]
    shifted$log_odds <- shifted$log_odds + 0.01
    shifted$log_odds_conf_low <- shifted$log_odds_conf_low + 0.01
    shifted$log_odds_conf_high <- shifted$log_odds_conf_high + 0.01
    shifted$statistic <- shifted$log_odds /
        shifted$log_odds_standard_error
    shifted$raw_p <- 2 * stats::pt(
        -abs(shifted$statistic),
        shifted$reference_df
    )
    shifted$adjusted_p <- shifted$raw_p
    shifted$flag <- shifted$adjusted_p <= 0.05
    shifted_log_odds$outcomes[[1L]] <- shifted

    changed_dispersion <- artifact
    changed_dispersion$outcomes[[1L]]$diagnostics$dispersion <-
        changed_dispersion$outcomes[[1L]]$diagnostics$dispersion + 0.01

    seeded <- artifact
    seeded$diagnostics$seed_manifest <- data.frame(
        protocol_id = "m13_a_association_protocol_v4",
        candidate_id = "a_fraction_quasibinomial",
        acquisition = "undeclared",
        hypothesis_id = artifact$hypotheses$hypothesis[[1L]],
        draw_id = 1L,
        seed = 1L,
        seed_nonce = 0L,
        map_retry = 0L,
        map_sha256 = paste(rep("0", 64L), collapse = ""),
        stringsAsFactors = FALSE
    )

    forged_failure <- artifact
    forged_failure$outcomes[[1L]] <-
        imputefinder:::.new_association_candidate_unavailable(
            artifact$hypotheses$hypothesis[[1L]],
            "association_quasibinomial_boundary"
        )
    forged_failure$hypotheses$family_member[[1L]] <- FALSE
    forged_failure$multiplicity$family_size[[1L]] <- 0L
    forged_failure$multiplicity$available_count[[1L]] <- 0L
    forged_failure$multiplicity$unavailable_count[[1L]] <- 1L
    forged_failure$diagnostics$strata$testable_count[[1L]] <- 0L
    forged_failure$diagnostics$strata$unavailable_count[[1L]] <- 1L
    forged_failure$diagnostics$warnings <-
        "association_quasibinomial_boundary"

    for (corrupted in list(
        shifted_effect,
        shifted_log_odds,
        changed_dispersion,
        seeded,
        forged_failure
    )) {
        expect_error(
            imputefinder:::.validate_association_candidate_artifact(
                corrupted,
                preparation
            ),
            class = "imputefinder_association_artifact_error"
        )
    }
})
