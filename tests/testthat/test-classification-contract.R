test_that("classification has no absent prototype plot-helper dependency", {
    implementation <- paste(deparse(body(classify_missingness)), collapse = "\n")

    expect_false(grepl("plot_detect_custom", implementation, fixed = TRUE))
})

test_that("a feature absent in one condition is rescued and audited", {
    result <- classify_normative_fixture()
    fixture <- normative_fixture()
    a_samples <- colnames(fixture$x)[fixture$group == "A"]
    seeded_a <- result$data["on_off", a_samples]
    log_row <- result$seed_log[
        result$seed_log$feature == "on_off" &
            result$seed_log$condition == "A",
        ,
        drop = FALSE
    ]

    expect_equal(sum(!is.na(seeded_a)), 1L)
    expect_equal(unname(seeded_a[!is.na(seeded_a)]), 8)
    expect_equal(nrow(log_row), 1L)
    expect_equal(log_row$inserted_value, 8)
    expect_identical(classification_row(result, "on_off", "A")$state, "MNAR")
    expect_false(any(result$seed_log$feature == "globally_absent"))
})

test_that("strict majority is derived from each condition size", {
    result <- classify_normative_fixture()

    expect_identical(classification_row(result, "mar_both", "A")$state, "MAR")
    expect_identical(classification_row(result, "mar_both", "B")$state, "MAR")
    expect_identical(
        classification_row(result, "sparse_mar", "A")$state,
        "insufficient"
    )
})

test_that("classifications remain condition-specific", {
    result <- classify_normative_fixture()

    expect_identical(classification_row(result, "on_off", "A")$state, "MNAR")
    expect_identical(classification_row(result, "on_off", "B")$state, "complete")
    expect_contains(result$groups$A$MNAR, "on_off")
    expect_contains(result$groups$B$complete, "on_off")
    expect_false("on_off" %in% result$groups$B$MNAR)
})

test_that("features classified MNAR in every condition are excluded", {
    result <- classify_normative_fixture()
    status <- result$feature_status[
        result$feature_status$feature == "all_mnar",
        ,
        drop = FALSE
    ]

    expect_false(status$retained)
    expect_identical(status$drop_reason, "MNAR_all_conditions")
    expect_false("all_mnar" %in% rownames(result$data))
})
