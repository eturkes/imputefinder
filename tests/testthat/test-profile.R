test_that("post-seed statistics follow the long result contract", {
    fixture <- normative_fixture()
    calculated <- post_seed_statistics(fixture$x, fixture$group)
    statistics <- calculated$statistics
    rescued <- calculated$rescued

    expect_identical(
        names(statistics),
        c(
            "feature",
            "condition",
            "sample_count",
            "observed_count",
            "missing_count",
            "missing_fraction",
            "mean_intensity",
            "seeded"
        )
    )
    expect_identical(
        vapply(statistics, typeof, character(1L)),
        c(
            feature = "character",
            condition = "character",
            sample_count = "integer",
            observed_count = "integer",
            missing_count = "integer",
            missing_fraction = "double",
            mean_intensity = "double",
            seeded = "logical"
        )
    )
    expect_identical(
        statistics$feature,
        rep(rownames(rescued$data), each = 2L)
    )
    expect_identical(
        statistics$condition,
        rep(c("A", "B"), times = nrow(rescued$data))
    )
    expect_true(all(is.finite(statistics$mean_intensity)))
})

test_that("post-seed statistics audit counts, arithmetic means, and seeds", {
    fixture <- normative_fixture()
    statistics <- post_seed_statistics(fixture$x, fixture$group)$statistics

    on_off_a <- statistics[
        statistics$feature == "on_off" & statistics$condition == "A",
        ,
        drop = FALSE
    ]
    expect_identical(on_off_a$sample_count, 4L)
    expect_identical(on_off_a$observed_count, 1L)
    expect_identical(on_off_a$missing_count, 3L)
    expect_identical(on_off_a$missing_fraction, 0.75)
    expect_identical(on_off_a$mean_intensity, 8)
    expect_true(on_off_a$seeded)

    mar_both_a <- statistics[
        statistics$feature == "mar_both" & statistics$condition == "A",
        ,
        drop = FALSE
    ]
    expect_identical(mar_both_a$observed_count, 3L)
    expect_identical(mar_both_a$missing_count, 1L)
    expect_equal(mar_both_a$mean_intensity, mean(c(14, 15, 16)))
    expect_false(mar_both_a$seeded)

    expect_false("globally_absent" %in% statistics$feature)
})

test_that("statistics reject a surviving all-missing condition block", {
    x <- rbind(
        invalid = c(NA, NA, 5, 6),
        valid = c(1, 2, 3, 4)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    groups <- stats::setNames(rep(c("A", "B"), each = 2L), colnames(x))

    expect_error(
        imputefinder:::.feature_condition_statistics(
            x,
            groups,
            imputefinder:::.empty_seed_log()
        ),
        "Every surviving feature-condition block must have a finite mean",
        fixed = TRUE
    )
})
