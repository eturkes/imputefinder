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
    result <- classify_missingness(fixture$x, fixture$group)
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
