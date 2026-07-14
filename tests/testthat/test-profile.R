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

test_that("profiles implement the count-weighted density formula", {
    fixture <- normative_fixture()
    statistics <- post_seed_statistics(fixture$x, fixture$group)$statistics
    profiles <- missingness_profiles(statistics)

    expect_identical(names(profiles), c("A", "B"))
    expect_identical(names(profiles$A), c("raw", "grid", "metadata"))
    expect_identical(
        names(profiles$A$grid),
        c(
            "intensity",
            "missing_density",
            "complete_density",
            "weighted_missing_density",
            "weighted_complete_density",
            "missing_proportion",
            "supported"
        )
    )
    expect_identical(
        profiles$A$metadata$class_counts,
        c(missing = 5L, complete = 1L)
    )

    profile <- profiles$A
    grid <- profile$grid
    expected <- with(
        grid,
        profile$metadata$class_counts[["missing"]] * missing_density /
            (
                profile$metadata$class_counts[["missing"]] *
                    missing_density +
                    profile$metadata$class_counts[["complete"]] *
                    complete_density
            )
    )
    expect_true(all(grid$supported))
    expect_equal(grid$missing_proportion, expected, tolerance = 1e-14)

    missing_means <- profile$raw$mean_intensity[profile$raw$has_missing]
    expected_density <- stats::density(
        missing_means,
        bw = profile$metadata$bandwidth,
        kernel = profile$metadata$kernel,
        n = profile$metadata$grid_points,
        from = min(grid$intensity),
        to = max(grid$intensity)
    )
    expect_equal(
        grid$missing_density,
        pmax(expected_density$y, 0),
        tolerance = 1e-14
    )
})

test_that("repeated means produce finite count-weighted profiles", {
    statistics <- data.frame(
        feature = c("missing_1", "complete", "missing_2"),
        condition = rep("A", 3L),
        sample_count = rep(2L, 3L),
        observed_count = c(1L, 2L, 1L),
        missing_count = c(1L, 0L, 1L),
        missing_fraction = c(0.5, 0, 0.5),
        mean_intensity = rep(10, 3L),
        seeded = c(FALSE, FALSE, TRUE),
        stringsAsFactors = FALSE
    )

    profile <- missingness_profiles(statistics)$A

    expect_true(is.finite(profile$metadata$bandwidth))
    expect_gt(profile$metadata$bandwidth, 0)
    expect_true(all(profile$grid$supported))
    expect_true(all(is.finite(profile$grid$missing_density)))
    expect_true(all(is.finite(profile$grid$complete_density)))
    expect_equal(
        profile$grid$missing_proportion,
        rep(2 / 3, profile$metadata$grid_points),
        tolerance = 1e-12
    )
})

test_that("one-class profile proportions are explicit", {
    statistics <- data.frame(
        feature = rep(c("one", "two"), each = 2L),
        condition = rep(c("complete_only", "missing_only"), times = 2L),
        sample_count = rep(2L, 4L),
        observed_count = c(2L, 1L, 2L, 1L),
        missing_count = c(0L, 1L, 0L, 1L),
        missing_fraction = c(0, 0.5, 0, 0.5),
        mean_intensity = c(1, 1, 2, 2),
        seeded = rep(FALSE, 4L),
        stringsAsFactors = FALSE
    )

    profiles <- missingness_profiles(statistics)
    complete_grid <- profiles$complete_only$grid
    missing_grid <- profiles$missing_only$grid

    expect_identical(
        profiles$complete_only$metadata$class_counts,
        c(missing = 0L, complete = 2L)
    )
    expect_identical(
        profiles$missing_only$metadata$class_counts,
        c(missing = 2L, complete = 0L)
    )
    expect_identical(
        profiles$complete_only$metadata$profile_type,
        "complete_only"
    )
    expect_identical(
        profiles$missing_only$metadata$profile_type,
        "missing_only"
    )
    expect_true(all(complete_grid$supported))
    expect_true(all(missing_grid$supported))
    expect_true(all(complete_grid$missing_proportion == 0))
    expect_true(all(missing_grid$missing_proportion == 1))
})

test_that("public results store raw evidence and profile grids", {
    result <- classify_normative_fixture()

    expect_identical(names(result$profiles), c("A", "B"))
    expect_identical(
        result$profiles$A$raw$feature,
        result$classifications$feature[
            result$classifications$condition == "A"
        ]
    )
    expect_identical(
        result$profiles$A$raw$has_missing,
        result$profiles$A$raw$missing_count > 0L
    )
    expect_true(
        result$profiles$A$raw$seeded[
            result$profiles$A$raw$feature == "on_off"
        ]
    )
    expect_identical(result$profiles$A$metadata$grid_points, 512L)
    expect_equal(
        result$profiles$A$metadata$observed_mean_range,
        c(minimum = 8, maximum = 20)
    )
    expect_identical(
        result$cutoff_diagnostics$A$profile,
        result$profiles$A$metadata
    )
})

test_that("profile grids are invariant to statistics row order", {
    fixture <- normative_fixture()
    statistics <- post_seed_statistics(fixture$x, fixture$group)$statistics
    baseline <- missingness_profiles(statistics)
    permuted <- missingness_profiles(statistics[rev(seq_len(nrow(statistics))), ])

    for (condition in names(baseline)) {
        expect_identical(
            baseline[[condition]]$grid,
            permuted[[condition]]$grid
        )
        expect_identical(
            baseline[[condition]]$metadata,
            permuted[[condition]]$metadata
        )
    }
})
