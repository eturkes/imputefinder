new_result_skeleton <- function(data, groups_by_sample) {
    imputefinder:::.new_imputefinder_result(
        data = data,
        groups_by_sample = groups_by_sample,
        call = quote(classify_missingness(x = x))
    )
}

test_that("the state vocabulary is stable and exhaustive", {
    expect_identical(
        imputefinder:::.missingness_states(),
        c("complete", "MNAR", "MAR", "insufficient")
    )
})

test_that("the result skeleton has typed auditable components", {
    fixture <- matrix_input_fixture()
    groups_by_sample <- setNames(fixture$group, colnames(fixture$x))

    result <- new_result_skeleton(fixture$x, groups_by_sample)

    expect_s3_class(result, "imputefinder_result")
    expect_identical(
        names(result),
        c(
            "data",
            "classifications",
            "groups",
            "feature_status",
            "cutoffs",
            "cutoff_diagnostics",
            "profiles",
            "seed_log",
            "groups_by_sample",
            "call"
        )
    )
    expect_identical(result$data, fixture$x)
    expect_identical(
        names(result$classifications),
        c(
            "feature",
            "condition",
            "sample_count",
            "observed_count",
            "missing_count",
            "missing_fraction",
            "mean_intensity",
            "cutoff",
            "state",
            "seeded",
            "retained",
            "drop_reason"
        )
    )
    expect_identical(
        vapply(result$classifications, typeof, character(1L)),
        c(
            feature = "character",
            condition = "character",
            sample_count = "integer",
            observed_count = "integer",
            missing_count = "integer",
            missing_fraction = "double",
            mean_intensity = "double",
            cutoff = "double",
            state = "character",
            seeded = "logical",
            retained = "logical",
            drop_reason = "character"
        )
    )
    expect_identical(
        names(result$feature_status),
        c("feature", "retained", "drop_reason")
    )
    expect_identical(
        names(result$seed_log),
        c(
            "feature",
            "condition",
            "sample",
            "old_value",
            "inserted_value",
            "seed"
        )
    )
    expect_identical(names(result$groups), c("A", "B"))
    expect_identical(
        names(result$groups$A),
        c("MNAR", "MAR", "complete", "MAR_or_complete")
    )
    expect_true(all(lengths(result$groups$A) == 0L))
    expect_identical(result$cutoffs, c(A = NA_real_, B = NA_real_))
    expect_identical(names(result$cutoff_diagnostics), c("A", "B"))
    expect_identical(names(result$profiles), c("A", "B"))
    expect_identical(result$groups_by_sample, groups_by_sample)
    expect_identical(result$call, quote(classify_missingness(x = x)))
})

test_that("the result constructor preserves a SummarizedExperiment output", {
    fixture <- summarized_experiment_input_fixture()
    prepared <- prepare_input(fixture$se, group_col = "condition")
    restored <- restore_output_data(prepared$data, prepared, fixture$se)

    result <- new_result_skeleton(restored, prepared$groups_by_sample)

    expect_s4_class(result$data, "SummarizedExperiment")
    expect_identical(result$data, fixture$se)
})

test_that("the result constructor rejects invalid state and group schemas", {
    fixture <- matrix_input_fixture()
    groups_by_sample <- setNames(fixture$group, colnames(fixture$x))
    classifications <- imputefinder:::.empty_classifications()
    classifications[1L, ] <- list(
        "protein_b",
        "A",
        1L,
        1L,
        0L,
        0,
        1,
        1,
        "unknown",
        FALSE,
        TRUE,
        NA_character_
    )

    expect_error(
        imputefinder:::.new_imputefinder_result(
            data = fixture$x,
            groups_by_sample = groups_by_sample,
            call = quote(classify_missingness(x = x)),
            classifications = classifications
        ),
        "Classification states must use the stable state vocabulary",
        fixed = TRUE
    )
    expect_error(
        new_result_skeleton(
            fixture$x,
            setNames(fixture$group, c("s2", "other", "s3"))
        ),
        "`groups_by_sample` names must match output sample order",
        fixed = TRUE
    )
})

test_that("the public classifier exposes the decided argument contract", {
    expect_identical(
        names(formals(classify_missingness)),
        c("x", "group", "group_col", "assay", "cutoffs", "seed")
    )
})
