automatic_cutoff_statistics <- function(
    n = 800L,
    right = 12,
    left = right - 0.6,
    height = 0.72,
    mar_rate = 0.05,
    center = 13,
    scale = 2.2,
    score_offset = 0L,
    missing = NULL,
    condition = "A"
) {
    index <- seq_len(n)
    mean_intensity <- stats::qnorm(
        (index - 0.5) / n,
        mean = center,
        sd = scale
    )
    if (is.null(missing)) {
        score <- (((index + score_offset) * 104729) %% 1000003) / 1000003
        ramp <- pmin(pmax((right - mean_intensity) / (right - left), 0), 1)
        missing_probability <- 1 - (1 - mar_rate) * (1 - height * ramp)
        missing <- score < missing_probability
    }

    data.frame(
        feature = sprintf("feature_%05d", index),
        condition = rep(condition, n),
        sample_count = rep(8L, n),
        observed_count = ifelse(missing, 7L, 8L),
        missing_count = ifelse(missing, 1L, 0L),
        missing_fraction = ifelse(missing, 1 / 8, 0),
        mean_intensity = mean_intensity,
        seeded = rep(FALSE, n),
        stringsAsFactors = FALSE
    )
}

automatic_cutoff_profile <- function(...) {
    statistics <- automatic_cutoff_statistics(...)
    imputefinder:::.condition_missingness_profile(statistics)
}

automatic_cutoff_matrix_fixture <- function() {
    a <- automatic_cutoff_statistics(
        n = 1200L,
        right = 11,
        left = 10.2,
        height = 0.68,
        center = 12.5,
        scale = 2.4
    )
    b <- automatic_cutoff_statistics(
        n = 1200L,
        right = 14,
        left = 13.2,
        height = 0.68,
        center = 12.5,
        scale = 2.4,
        score_offset = 37L,
        condition = "B"
    )
    condition_block <- function(statistics) {
        block <- matrix(
            rep(statistics$mean_intensity, 4L),
            nrow = nrow(statistics)
        )
        missing <- statistics$missing_count > 0L
        block[cbind(which(missing), rep(4L, sum(missing)))] <- NA_real_
        block
    }
    x <- cbind(condition_block(a), condition_block(b))
    rownames(x) <- a$feature
    colnames(x) <- paste0("sample_", seq_len(ncol(x)))

    list(
        x = x,
        group = rep(c("A", "B"), each = 4L),
        expected = c(A = 10.698082335694377, B = 13.815875319932115)
    )
}

test_that("automatic cutoff records the selected algorithm and quality", {
    profile <- automatic_cutoff_profile()
    decision <- imputefinder:::.detect_automatic_cutoff(profile)

    expect_identical(decision$status, "ok")
    expect_equal(decision$cutoff, 11.7719792534904, tolerance = 1e-12)
    expect_identical(decision$method, "derivative_boundary")
    expect_identical(decision$method_version, "1")
    expect_identical(decision$reason, "")
    expect_identical(decision$warnings, character())

    quality <- decision$quality
    expect_identical(quality$evidence$feature_count, 800L)
    expect_identical(
        quality$evidence$class_counts,
        profile$metadata$class_counts
    )
    expect_true(quality$trend$descending)
    expect_lte(quality$trend$p_value, 0.001)
    expect_gte(quality$support$density_supported_grid_points, 25L)
    expect_gte(quality$derivative$credible_lobe_count, 1L)
    expect_gte(quality$selected_transition$drop, 0.05)
    expect_gte(quality$selected_transition$peak_to_noise, 1.5)
})

test_that("automatic cutoff passes frozen stress profiles", {
    fixtures <- list(
        heavy_mar = list(
            profile = automatic_cutoff_profile(
                n = 4000L,
                left = 11,
                height = 0.68,
                mar_rate = 0.25
            ),
            expected = 11.660167998616307,
            target = 12,
            tolerance = 1
        ),
        near_floor = list(
            profile = automatic_cutoff_profile(n = 96L),
            expected = 11.76302764656138,
            target = 12,
            tolerance = 1
        )
    )

    for (fixture in fixtures) {
        decision <- imputefinder:::.detect_automatic_cutoff(fixture$profile)
        observed_range <- fixture$profile$metadata$observed_mean_range

        expect_identical(decision$status, "ok")
        expect_equal(decision$cutoff, fixture$expected, tolerance = 1e-12)
        expect_lte(abs(decision$cutoff - fixture$target), fixture$tolerance)
        expect_true(is.finite(decision$cutoff))
        expect_gt(decision$cutoff, observed_range[["minimum"]])
        expect_lt(decision$cutoff, observed_range[["maximum"]])
    }
})

test_that("automatic cutoff enforces total and per-class evidence floors", {
    below_total <- imputefinder:::.detect_automatic_cutoff(
        automatic_cutoff_profile(n = 79L)
    )
    expect_identical(below_total$status, "unidentifiable")
    expect_true(is.na(below_total$cutoff))
    expect_match(below_total$reason, "at least 80 feature blocks", fixed = TRUE)
    expect_identical(below_total$quality$evidence$feature_count, 79L)

    missing <- seq_len(80L) <= 11L
    below_class <- imputefinder:::.detect_automatic_cutoff(
        automatic_cutoff_profile(n = 80L, missing = missing)
    )
    expect_identical(below_class$status, "unidentifiable")
    expect_match(
        below_class$reason,
        "at least 12 missing and 12 complete feature blocks",
        fixed = TRUE
    )
    expect_identical(
        below_class$quality$evidence$class_counts,
        c(missing = 11L, complete = 69L)
    )

    complete_only <- imputefinder:::.detect_automatic_cutoff(
        automatic_cutoff_profile(n = 100L, missing = rep(FALSE, 100L))
    )
    expect_identical(complete_only$status, "unidentifiable")
    expect_true(is.na(complete_only$cutoff))
    expect_identical(
        complete_only$quality$evidence$class_counts,
        c(missing = 0L, complete = 100L)
    )

    missing_only <- imputefinder:::.detect_automatic_cutoff(
        automatic_cutoff_profile(n = 100L, missing = rep(TRUE, 100L))
    )
    expect_identical(missing_only$status, "unidentifiable")
    expect_true(is.na(missing_only$cutoff))
    expect_identical(
        missing_only$quality$evidence$class_counts,
        c(missing = 100L, complete = 0L)
    )
})

test_that("flat and unsupported profiles fail without endpoint fallbacks", {
    flat_profile <- automatic_cutoff_profile(
        n = 2000L,
        missing = seq_len(2000L) %% 8L == 0L
    )
    flat <- imputefinder:::.detect_automatic_cutoff(flat_profile)

    expect_identical(flat$status, "unidentifiable")
    expect_true(is.na(flat$cutoff))
    expect_match(
        flat$reason,
        "No statistically credible descending missingness trend",
        fixed = TRUE
    )
    expect_false(flat$quality$trend$credible)

    unsupported_profile <- automatic_cutoff_profile()
    unsupported_profile$grid$supported[] <- FALSE
    unsupported <- imputefinder:::.detect_automatic_cutoff(unsupported_profile)

    expect_identical(unsupported$status, "unidentifiable")
    expect_true(is.na(unsupported$cutoff))
    expect_match(
        unsupported$reason,
        "no supported observed-intensity grid",
        fixed = TRUE
    )
})

test_that("automatic cutoff decision is exact under evidence row order", {
    statistics <- automatic_cutoff_statistics()
    baseline <- imputefinder:::.detect_automatic_cutoff(
        imputefinder:::.condition_missingness_profile(statistics)
    )
    permutation <- c(seq(2L, nrow(statistics), by = 2L),
        seq(1L, nrow(statistics), by = 2L))
    permuted <- imputefinder:::.detect_automatic_cutoff(
        imputefinder:::.condition_missingness_profile(
            statistics[permutation, , drop = FALSE]
        )
    )

    expect_identical(permuted, baseline)
})

test_that("public automatic cutoffs are condition-specific and overridable", {
    fixture <- automatic_cutoff_matrix_fixture()
    original <- fixture$x
    result <- classify_missingness(fixture$x, fixture$group)

    expect_identical(fixture$x, original)
    expect_equal(result$cutoffs, fixture$expected, tolerance = 1e-12)
    expect_identical(
        vapply(result$cutoff_diagnostics, `[[`, character(1L), "source"),
        c(A = "automatic", B = "automatic")
    )
    expect_identical(result$cutoff_diagnostics$A$method, "derivative_boundary")
    expect_identical(result$cutoff_diagnostics$A$method_version, "1")
    expect_true(all(
        result$classifications$cutoff[
            result$classifications$condition == "A"
        ] == result$cutoffs[["A"]]
    ))

    partial <- classify_missingness(
        fixture$x,
        fixture$group,
        cutoffs = c(A = 10.5)
    )
    expect_identical(partial$cutoffs[["A"]], 10.5)
    expect_equal(
        partial$cutoffs[["B"]],
        fixture$expected[["B"]],
        tolerance = 1e-12
    )
    expect_identical(partial$cutoff_diagnostics$A$source, "manual")
    expect_identical(partial$cutoff_diagnostics$B$source, "automatic")
})

test_that("public automatic decisions are exact under input order", {
    fixture <- automatic_cutoff_matrix_fixture()
    baseline <- classify_missingness(fixture$x, fixture$group)
    feature_order <- rev(seq_len(nrow(fixture$x)))
    sample_order <- rev(seq_len(ncol(fixture$x)))
    permuted <- classify_missingness(
        fixture$x[feature_order, sample_order, drop = FALSE],
        fixture$group[sample_order]
    )

    expect_identical(permuted$cutoffs, baseline$cutoffs)
    expect_identical(
        permuted$cutoff_diagnostics,
        baseline$cutoff_diagnostics
    )
    for (condition in names(baseline$profiles)) {
        expect_identical(
            permuted$profiles[[condition]]$grid,
            baseline$profiles[[condition]]$grid
        )
        expect_identical(
            permuted$profiles[[condition]]$metadata,
            baseline$profiles[[condition]]$metadata
        )
    }

    baseline_state <- stats::setNames(
        baseline$classifications$state,
        paste(
            baseline$classifications$feature,
            baseline$classifications$condition,
            sep = "\r"
        )
    )
    permuted_state <- stats::setNames(
        permuted$classifications$state,
        paste(
            permuted$classifications$feature,
            permuted$classifications$condition,
            sep = "\r"
        )
    )
    expect_identical(
        unname(permuted_state[names(baseline_state)]),
        unname(baseline_state)
    )
})

test_that("automatic failure carries its condition and diagnostic profile", {
    fixture <- normative_fixture()
    error <- tryCatch(
        classify_missingness(
            fixture$x,
            fixture$group,
            cutoffs = c(A = 12)
        ),
        imputefinder_cutoff_error = identity
    )

    expect_s3_class(error, "imputefinder_cutoff_error")
    expect_s3_class(error, "imputefinder_cutoff_unidentifiable")
    expect_identical(error$condition, "B")
    expect_match(error$reason, "at least 80 feature blocks", fixed = TRUE)
    expect_match(
        conditionMessage(error),
        "supply a manual cutoff for `B`",
        fixed = TRUE
    )
    expect_identical(error$diagnostic$source, "automatic")
    expect_identical(error$diagnostic$method, "derivative_boundary")
    expect_identical(error$diagnostic$method_version, "1")
    expect_identical(error$diagnostic$profile, error$profile$metadata)
    expect_identical(error$profile$metadata$feature_count, 6L)
})
