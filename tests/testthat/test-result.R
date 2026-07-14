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

test_that("summary exposes exact feature, state, drop, and cutoff counts", {
    result <- classify_normative_fixture()

    result_summary <- summary(result)

    expect_s3_class(result_summary, "summary.imputefinder_result")
    expect_identical(
        names(result_summary),
        c(
            "call",
            "features",
            "samples",
            "seed_insertions",
            "states",
            "retained_states",
            "drops",
            "cutoffs"
        )
    )
    expect_identical(
        result_summary$features,
        c(total = 7L, retained = 4L, dropped = 3L)
    )
    expect_identical(
        result_summary$samples,
        c(total = 8L, conditions = 2L)
    )
    expect_identical(result_summary$seed_insertions, 1L)
    expect_identical(
        result_summary$states,
        data.frame(
            condition = c("A", "B"),
            complete = c(1L, 4L),
            MNAR = c(3L, 1L),
            MAR = c(1L, 1L),
            insufficient = c(1L, 0L),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(
        result_summary$retained_states,
        data.frame(
            condition = c("A", "B"),
            complete = c(1L, 3L),
            MNAR = c(2L, 0L),
            MAR = c(1L, 1L),
            insufficient = c(0L, 0L),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(
        result_summary$drops,
        data.frame(
            reason = c(
                "all_missing",
                "insufficient:A",
                "MNAR_all_conditions"
            ),
            count = rep(1L, 3L),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(
        result_summary$cutoffs,
        data.frame(
            condition = c("A", "B"),
            cutoff = c(12, 12),
            source = c("manual", "manual"),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(result_summary$call, result$call)
})

test_that("result and summary print methods are concise and return invisibly", {
    result <- classify_normative_fixture()
    result_visibility <- NULL
    result_output <- capture.output(
        result_visibility <- withVisible(print(result))
    )

    expect_false(result_visibility$visible)
    expect_identical(result_visibility$value, result)
    expect_identical(
        result_output,
        c(
            "<imputefinder_result>",
            "Features: 4/7 retained; 3 dropped",
            "Samples: 8 across 2 conditions; seed insertions: 1",
            paste0(
                "A: retained MNAR=2, MAR=1, complete=1; ",
                "cutoff=12 (manual)"
            ),
            paste0(
                "B: retained MNAR=0, MAR=1, complete=3; ",
                "cutoff=12 (manual)"
            )
        )
    )

    result_summary <- summary(result)
    summary_visibility <- NULL
    summary_output <- capture.output(
        summary_visibility <- withVisible(print(result_summary))
    )

    expect_false(summary_visibility$visible)
    expect_identical(summary_visibility$value, result_summary)
    expect_identical(
        summary_output,
        c(
            "<summary.imputefinder_result>",
            "Features: 4/7 retained; 3 dropped",
            "Samples: 8 across 2 conditions; seed insertions: 1",
            "States (all classified features):",
            "A: complete=1, MNAR=3, MAR=1, insufficient=1",
            "B: complete=4, MNAR=1, MAR=1, insufficient=0",
            paste0(
                "Dropped: all_missing=1; insufficient:A=1; ",
                "MNAR_all_conditions=1"
            ),
            "Cutoffs: A=12 (manual); B=12 (manual)"
        )
    )
})

test_that("presentation handles no drops and conditions without cutoffs", {
    x <- rbind(
        low = c(1, 1, 2, 2),
        high = c(10, 10, 20, 20)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    result <- classify_missingness(x, c("B", "B", "A", "A"))

    result_summary <- summary(result)
    output <- capture.output(print(result_summary))

    expect_identical(
        result_summary$drops,
        data.frame(
            reason = character(),
            count = integer(),
            stringsAsFactors = FALSE
        )
    )
    expect_true(any(output == "Dropped: none"))
    expect_true(any(output == "Cutoffs: A=not needed; B=not needed"))
})
