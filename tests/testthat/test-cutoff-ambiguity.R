ambiguous_cutoff_fixture <- function() {
    n <- 4000L
    intensity <- seq(6, 22, length.out = n)
    probability <- 0.1 + 0.2 * stats::plogis((10 - intensity) / 0.2) +
        0.5 * stats::plogis((16 - intensity) / 0.2)
    quantiles <- (seq_len(n) * ((sqrt(5) - 1) / 2)) %% 1
    x <- matrix(
        rep(intensity, 8L),
        ncol = 8L,
        dimnames = list(paste0("f", seq_len(n)), paste0("s", 1:8))
    )
    x[quantiles < probability, 1L] <- NA_real_
    list(x = x, group = rep(c("A", "B"), each = 4L))
}

test_that("competing credible cutoff transitions expose the selection policy", {
    fixture <- ambiguous_cutoff_fixture()
    expect_warning(
        result <- classify_missingness(fixture$x, fixture$group),
        class = "imputefinder_cutoff_ambiguity_warning"
    )
    diagnostic <- result$cutoff_diagnostics$A
    expect_match(
        diagnostic$warnings,
        "multiple credible downward transitions",
        fixed = TRUE
    )
    expect_match(diagnostic$warnings, "lowest-intensity", fixed = TRUE)
    expect_identical(
        diagnostic$quality$derivative$credible_lobe_count,
        2L
    )
    expect_lt(diagnostic$quality$selected_transition$peak_intensity, 12)
    expect_lt(result$cutoffs[["A"]], 12)
    expect_match(
        plot_missingness(result, "A")$labels$caption,
        "multiple credible downward transitions",
        fixed = TRUE
    )

    manual <- classify_missingness(
        fixture$x, fixture$group, cutoffs = c(A = 15)
    )
    expect_identical(manual$cutoffs[["A"]], 15)
    expect_identical(manual$cutoff_diagnostics$A$source, "manual")
    expect_identical(manual$cutoff_diagnostics$B$source, "not_needed")
})

test_that("strict ambiguity errors preserve alternatives and allow manual override", {
    fixture <- ambiguous_cutoff_fixture()
    error <- tryCatch(
        classify_missingness(fixture$x, fixture$group, ambiguity = "error"),
        imputefinder_cutoff_ambiguous = identity
    )
    expect_s3_class(error, "imputefinder_cutoff_error")
    expect_identical(error$condition, "A")
    expect_identical(error$diagnostic$ambiguity_policy, "error")
    alternatives <- error$diagnostic$quality$candidate_transitions
    expect_identical(nrow(alternatives), 2L)
    expect_true(all(is.finite(alternatives$cutoff)))
    expect_true(all(diff(alternatives$peak_intensity) > 0))
    expect_identical(sum(alternatives$selected), 1L)
    expect_type(error$profile$raw$mean_intensity, "double")

    expect_silent(lowest <- classify_missingness(fixture$x, fixture$group,
        ambiguity = "lowest"))
    expect_identical(lowest$cutoff_diagnostics$A$ambiguity_policy, "lowest")
    expect_identical(lowest$cutoffs[["A"]], alternatives$cutoff[[1]])
    expect_silent(manual <- classify_missingness(fixture$x, fixture$group,
        ambiguity = "error", cutoffs = c(A = 15)))
    expect_identical(manual$cutoffs[["A"]], 15)
})
