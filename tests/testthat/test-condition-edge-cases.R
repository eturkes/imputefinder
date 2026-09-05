test_that("cancellation cannot make classifications depend on sample order", {
    x <- rbind(
        target = c(-1e20, 1e20, 1, NA, rep(16, 4)),
        low = rep(-1, 8),
        high = rep(16, 8)
    )
    colnames(x) <- paste0("s", 1:8)
    group <- stats::setNames(rep(c("A", "B"), each = 4L), colnames(x))
    original <- classify_missingness(x, group, cutoffs = c(A = 0.2))
    sample_order <- c(1L, 3L, 2L, 4:8)
    reordered <- classify_missingness(
        x[, sample_order, drop = FALSE], group, cutoffs = c(A = 0.2)
    )

    expect_identical(reordered$classifications, original$classifications)
    expect_identical(reordered$feature_status, original$feature_status)
    expect_identical(reordered$profiles, original$profiles)
    expect_identical(
        reordered$data[, colnames(x), drop = FALSE], original$data
    )
    target <- original$classifications$feature == "target" &
        original$classifications$condition == "A"
    expect_identical(original$classifications$state[target], "MAR")
    expect_equal(original$classifications$mean_intensity[target], 1 / 3)
})

test_that("singleton rescue records why a previously absent block is complete", {
    x <- rbind(on_off = c(NA, 15), low = c(8, 8))
    colnames(x) <- c("s1", "s2")
    result <- classify_missingness(x, c("A", "B"), cutoffs = c(A = 12))
    target <- result$classifications$feature == "on_off" &
        result$classifications$condition == "A"

    expect_identical(result$classifications$state[target], "complete")
    expect_true(result$classifications$seeded[target])
    expect_identical(result$classifications$missing_count[target], 0L)
    expect_identical(result$seed_log$feature, "on_off")
    expect_identical(result$seed_log$inserted_value, 8)
    expect_identical(result$cutoff_diagnostics$A$source, "not_needed")
})
