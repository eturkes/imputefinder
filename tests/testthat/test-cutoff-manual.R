manual_cutoff_statistics <- function() {
    fixture <- normative_fixture()
    post_seed_statistics(fixture$x, fixture$group)$statistics
}

test_that("manual cutoffs are named, finite, and optionally partial", {
    statistics <- manual_cutoff_statistics()

    partial <- imputefinder:::.resolve_manual_cutoffs(
        statistics,
        c(B = 12)
    )
    expect_identical(partial$cutoffs, c(A = NA_real_, B = 12))
    expect_null(partial$diagnostics$A)
    expect_identical(partial$diagnostics$B$source, "manual")
    expect_identical(partial$diagnostics$B$warnings, character())

    complete <- imputefinder:::.resolve_manual_cutoffs(
        statistics,
        c(B = 13, A = 11)
    )
    expect_identical(complete$cutoffs, c(A = 11, B = 13))
    expect_identical(names(complete$diagnostics), c("A", "B"))
})

test_that("invalid manual cutoff specifications fail explicitly", {
    statistics <- manual_cutoff_statistics()
    invalid_shape <- list(
        c(12, 12),
        c(A = "12"),
        numeric(),
        matrix(c(A = 12), nrow = 1L)
    )
    for (cutoffs in invalid_shape) {
        expect_error(
            imputefinder:::.resolve_manual_cutoffs(statistics, cutoffs),
            "`cutoffs` must be NULL or a non-empty named numeric vector.",
            fixed = TRUE
        )
    }

    for (cutoffs in list(c(A = NA_real_), c(A = NaN), c(A = Inf))) {
        expect_error(
            imputefinder:::.resolve_manual_cutoffs(statistics, cutoffs),
            "Manual cutoffs must contain only finite values.",
            fixed = TRUE
        )
    }

    expect_error(
        imputefinder:::.resolve_manual_cutoffs(statistics, c(A = 12, A = 13)),
        "Manual cutoff names must be unique, non-empty condition labels.",
        fixed = TRUE
    )
    expect_error(
        imputefinder:::.resolve_manual_cutoffs(statistics, c(other = 12)),
        "Manual cutoff names must match condition labels exactly: `other`.",
        fixed = TRUE
    )
})

test_that("out-of-range manual cutoffs are retained with diagnostics", {
    statistics <- manual_cutoff_statistics()

    resolved <- expect_silent(
        imputefinder:::.resolve_manual_cutoffs(
            statistics,
            c(A = -100, B = 100)
        )
    )

    expect_identical(resolved$cutoffs, c(A = -100, B = 100))
    expect_match(
        resolved$diagnostics$A$warnings,
        "outside the observed feature-mean range",
        fixed = TRUE
    )
    expect_match(
        resolved$diagnostics$B$warnings,
        "outside the observed feature-mean range",
        fixed = TRUE
    )
})

test_that("a condition without missing values needs no cutoff", {
    x <- rbind(
        low = c(1, 2, 5, NA),
        high = c(3, 4, 6, 7)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    groups <- stats::setNames(rep(c("A", "B"), each = 2L), colnames(x))
    statistics <- imputefinder:::.feature_condition_statistics(
        x,
        groups,
        imputefinder:::.empty_seed_log()
    )

    resolved <- imputefinder:::.resolve_manual_cutoffs(
        statistics,
        c(A = -100, B = 5.5)
    )
    classified <- imputefinder:::.assign_condition_states(
        statistics,
        resolved$cutoffs
    )

    expect_identical(resolved$cutoffs, c(A = NA_real_, B = 5.5))
    expect_identical(resolved$diagnostics$A$source, "not_needed")
    expect_identical(resolved$diagnostics$A$warnings, character())
    expect_true(all(classified$state[classified$condition == "A"] == "complete"))
    expect_true(all(is.na(classified$cutoff[classified$condition == "A"])))
})
