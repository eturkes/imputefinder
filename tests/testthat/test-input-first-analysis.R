input_first_design <- function(fixture, acquisition = FALSE) {
    sample_data <- data.frame(
        condition = fixture$group,
        row.names = colnames(fixture$x)
    )
    if (acquisition) {
        sample_data$acquisition_mode <- rep("DIA", nrow(sample_data))
    }
    sample_data <- sample_data[rev(rownames(sample_data)), , drop = FALSE]

    missingness_design(
        sample_data,
        condition = "condition",
        acquisition = if (acquisition) "acquisition_mode" else NULL
    )
}

test_that("input-first matrix construction freezes pre-rescue evidence", {
    fixture <- normative_fixture()
    design <- input_first_design(fixture, acquisition = TRUE)
    original_x <- fixture$x
    original_design <- design
    analysis_call <- quote(
        analyze_missingness(x, design, cutoffs = cutoffs, seed = 19L)
    )

    analysis <- imputefinder:::.analyze_missingness_input_first(
        x = fixture$x,
        design = design,
        cutoffs = fixture$cutoffs,
        seed = 19L,
        call = analysis_call
    )
    expected <- classify_missingness(
        fixture$x,
        group = stats::setNames(fixture$group, colnames(fixture$x)),
        cutoffs = fixture$cutoffs,
        seed = 19L
    )
    expected$call <- analysis_call

    expect_identical(
        names(formals(imputefinder:::.analyze_missingness_input_first)),
        c(
            "x",
            "design",
            "assay",
            "base_fit",
            "cutoffs",
            "seed",
            "modules",
            "call"
        )
    )
    expect_true("analyze_missingness" %in% getNamespaceExports("imputefinder"))
    expect_s3_class(analysis, "imputefinder_analysis")
    expect_identical(analysis$classic, expected)
    expect_identical(
        rownames(analysis$design$declared$sample_data),
        colnames(fixture$x)
    )
    expect_identical(
        analysis$input$acquisition,
        stats::setNames(rep("DIA", ncol(fixture$x)), colnames(fixture$x))
    )
    expect_identical(names(analysis$sentinel), "pre_rescue")
    expect_identical(
        names(analysis$sentinel$pre_rescue),
        c("schema", "feature_condition", "sample")
    )
    expect_identical(
        analysis$sentinel$pre_rescue$schema,
        "pre_rescue_evidence_v1"
    )
    expect_identical(
        names(analysis$sentinel$pre_rescue$feature_condition),
        c(
            "feature",
            "condition",
            "sample_count",
            "observed_count",
            "missing_count",
            "missing_fraction",
            "mean_intensity"
        )
    )
    expect_identical(
        names(analysis$sentinel$pre_rescue$sample),
        c(
            "sample",
            "condition",
            "feature_count",
            "observed_count",
            "missing_count",
            "missing_fraction",
            "mean_intensity"
        )
    )

    feature_condition <- analysis$sentinel$pre_rescue$feature_condition
    on_off_a <- feature_condition[
        feature_condition$feature == "on_off" &
            feature_condition$condition == "A",
        ,
        drop = FALSE
    ]
    globally_absent <- feature_condition[
        feature_condition$feature == "globally_absent",
        ,
        drop = FALSE
    ]
    expect_identical(on_off_a$observed_count, 0L)
    expect_identical(on_off_a$missing_count, 4L)
    expect_identical(on_off_a$missing_fraction, 1)
    expect_identical(on_off_a$mean_intensity, NA_real_)
    expect_identical(globally_absent$observed_count, c(0L, 0L))
    expect_true(all(is.na(globally_absent$mean_intensity)))

    sample <- analysis$sentinel$pre_rescue$sample
    expect_identical(sample$sample, colnames(fixture$x))
    expect_identical(sample$condition, fixture$group)
    expect_identical(sample$feature_count, rep(7L, 8L))
    expect_identical(analysis$provenance$seeds, list(classic_rescue = 19L))
    expect_identical(fixture$x, original_x)
    expect_identical(design, original_design)
    expect_invisible(imputefinder:::.validate_imputefinder_analysis(analysis))
})

test_that("pre-rescue evidence survives an automatic classic failure", {
    fixture <- normative_fixture()
    design <- input_first_design(fixture)
    original <- fixture$x
    analysis_call <- quote(analyze_missingness(x, design))

    analysis <- imputefinder:::.analyze_missingness_input_first(
        fixture$x,
        design,
        call = analysis_call
    )

    expect_s3_class(analysis$classic, "imputefinder_classic_failure")
    expect_identical(analysis$classic$stage, "cutoff")
    expect_true(
        "imputefinder_cutoff_unidentifiable" %in% analysis$classic$class
    )
    expect_match(
        analysis$classic$message,
        "Automatic cutoff is unidentifiable",
        fixed = TRUE
    )
    expect_identical(analysis$classic$call, analysis_call)
    expect_identical(
        analysis$provenance$failures,
        list(classic = analysis$classic)
    )
    expect_identical(
        imputefinder:::.unpack_original_mask(
            analysis$input$original_mask,
            analysis$input$dimensions,
            analysis$input$feature_names,
            analysis$input$sample_names
        ),
        is.na(fixture$x)
    )
    expect_identical(
        analysis$sentinel$pre_rescue,
        imputefinder:::.new_pre_rescue_evidence(
            fixture$x,
            stats::setNames(fixture$group, colnames(fixture$x))
        )
    )
    expect_false("seeded" %in% names(
        analysis$sentinel$pre_rescue$feature_condition
    ))
    expect_identical(fixture$x, original)

    round_trip <- unserialize(serialize(analysis, NULL, version = 3L))
    expect_identical(round_trip, analysis)
    expect_invisible(
        imputefinder:::.validate_imputefinder_analysis(round_trip)
    )
})

test_that("pre-rescue evidence survives a rescue-stage classic failure", {
    x <- matrix(
        c(NA, NA, 10, 11),
        nrow = 2L,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    design <- missingness_design(
        data.frame(
            condition = c("A", "B"),
            row.names = colnames(x)
        ),
        condition = "condition"
    )

    analysis <- imputefinder:::.analyze_missingness_input_first(
        x,
        design,
        cutoffs = c(A = 5, B = 5)
    )

    expect_s3_class(analysis$classic, "imputefinder_classic_failure")
    expect_identical(analysis$classic$stage, "rescue")
    expect_match(
        analysis$classic$message,
        "Condition `A` has no finite intensity",
        fixed = TRUE
    )
    expect_identical(
        analysis$sentinel$pre_rescue$sample$observed_count,
        c(0L, 2L)
    )
    expect_identical(
        is.na(analysis$sentinel$pre_rescue$sample$mean_intensity),
        c(TRUE, FALSE)
    )
    expect_invisible(imputefinder:::.validate_imputefinder_analysis(analysis))
})

test_that("input-first SummarizedExperiment routing is aligned and exact", {
    fixture <- normative_fixture()
    se <- rich_summarized_experiment(fixture$x, fixture$group)
    sample_data <- as.data.frame(SummarizedExperiment::colData(se))
    sample_data$acquisition_mode <- rep("DDA", nrow(sample_data))
    design <- missingness_design(
        sample_data[rev(rownames(sample_data)), , drop = FALSE],
        condition = "condition",
        acquisition = "acquisition_mode"
    )
    original_se <- se
    original_design <- design
    analysis_call <- quote(
        analyze_missingness(
            se,
            design,
            assay = "intensity",
            cutoffs = cutoffs,
            seed = 7L
        )
    )

    analysis <- imputefinder:::.analyze_missingness_input_first(
        se,
        design,
        assay = "intensity",
        cutoffs = fixture$cutoffs,
        seed = 7L,
        call = analysis_call
    )
    expected <- classify_missingness(
        se,
        group_col = "condition",
        assay = "intensity",
        cutoffs = fixture$cutoffs,
        seed = 7L
    )
    expected$call <- analysis_call

    expect_identical(analysis$classic, expected)
    expect_identical(analysis$input$representation, "SummarizedExperiment")
    expect_identical(analysis$input$assay, "intensity")
    expect_identical(
        analysis$input$fingerprint,
        imputefinder:::.matrix_fingerprint(fixture$x)
    )
    expect_identical(
        analysis$sentinel$pre_rescue,
        imputefinder:::.new_pre_rescue_evidence(
            fixture$x,
            stats::setNames(fixture$group, colnames(fixture$x))
        )
    )
    expect_identical(se, original_se)
    expect_identical(design, original_design)

    mismatched <- design
    mismatched$sample_data$condition[[1L]] <- "mismatched"
    error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            se,
            mismatched,
            assay = "intensity",
            cutoffs = fixture$cutoffs
        ),
        imputefinder_design_condition_error = identity
    )
    expect_s3_class(error, "imputefinder_design_condition_error")
    expect_s3_class(error, "imputefinder_design_error")
    expect_identical(error$mismatched_samples, "s8")
    expect_identical(se, original_se)
})

test_that("input-first construction rejects design and assay misalignment", {
    fixture <- normative_fixture()
    design <- input_first_design(fixture)

    wrong_samples <- design
    rownames(wrong_samples$sample_data)[[1L]] <- "unexpected"
    alignment_error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            fixture$x,
            wrong_samples,
            cutoffs = fixture$cutoffs
        ),
        imputefinder_design_alignment_error = identity
    )
    expect_s3_class(alignment_error, "imputefinder_design_alignment_error")
    expect_identical(alignment_error$missing_samples, "unexpected")
    expect_identical(alignment_error$unexpected_samples, "s8")

    se <- rich_summarized_experiment(fixture$x, fixture$group)
    names(SummarizedExperiment::colData(se))[
        names(SummarizedExperiment::colData(se)) == "condition"
    ] <- "other_condition"
    condition_error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            se,
            design,
            assay = "intensity",
            cutoffs = fixture$cutoffs
        ),
        imputefinder_design_condition_error = identity
    )
    expect_s3_class(condition_error, "imputefinder_design_condition_error")
    expect_identical(condition_error$condition_role, "condition")
    expect_identical(condition_error$mismatched_samples, character())

    expect_error(
        imputefinder:::.analyze_missingness_input_first(
            fixture$x,
            design,
            assay = "intensity",
            cutoffs = fixture$cutoffs
        ),
        "`assay` must be NULL for matrix input",
        fixed = TRUE
    )
})

test_that("pre-rescue schema validation rejects drift and corruption", {
    fixture <- normative_fixture()
    analysis <- imputefinder:::.analyze_missingness_input_first(
        fixture$x,
        input_first_design(fixture),
        cutoffs = fixture$cutoffs
    )

    wrong_schema <- analysis
    wrong_schema$sentinel$pre_rescue$schema <- "pre_rescue_evidence_v2"
    lifecycle_error <- tryCatch(
        imputefinder:::.validate_imputefinder_analysis(wrong_schema),
        imputefinder_pre_rescue_lifecycle_error = identity
    )
    expect_s3_class(
        lifecycle_error,
        "imputefinder_pre_rescue_lifecycle_error"
    )
    expect_identical(
        lifecycle_error$expected_schema,
        "pre_rescue_evidence_v1"
    )
    expect_identical(
        lifecycle_error$actual_schema,
        "pre_rescue_evidence_v2"
    )

    bad_count <- analysis
    bad_count$sentinel$pre_rescue$sample$observed_count[[1L]] <- 99L
    expect_error(
        imputefinder:::.validate_imputefinder_analysis(bad_count),
        class = "imputefinder_pre_rescue_schema_error"
    )

    bad_mean <- analysis
    row <- which(
        bad_mean$sentinel$pre_rescue$feature_condition$observed_count == 0L
    )[[1L]]
    bad_mean$sentinel$pre_rescue$feature_condition$mean_intensity[[row]] <- 0
    expect_error(
        imputefinder:::.validate_imputefinder_analysis(bad_mean),
        class = "imputefinder_pre_rescue_schema_error"
    )

    malformed_sentinel <- analysis
    malformed_sentinel$sentinel$extra <- list()
    expect_error(
        imputefinder:::.validate_imputefinder_analysis(malformed_sentinel),
        class = "imputefinder_analysis_schema_error"
    )
})

test_that("input-first execution preserves caller RNG and stable v1 output", {
    fixture <- normative_fixture()
    design <- input_first_design(fixture)
    baseline <- classify_normative_fixture()
    baseline_bytes <- serialize(baseline, NULL, version = 3L)

    with_preserved_random_state({
        RNGkind("Knuth-TAOCP-2002", "Box-Muller", "Rejection")
        set.seed(4812L)
        kind_before <- RNGkind()
        seed_before <- .Random.seed

        analysis <- imputefinder:::.analyze_missingness_input_first(
            fixture$x,
            design,
            cutoffs = fixture$cutoffs,
            seed = 31L
        )

        expect_s3_class(analysis$classic, "imputefinder_result")
        expect_identical(RNGkind(), kind_before)
        expect_identical(.Random.seed, seed_before)
    })

    after <- classify_normative_fixture()
    expect_identical(serialize(after, NULL, version = 3L), baseline_bytes)
    expect_identical(
        names(formals(classify_missingness)),
        c("x", "group", "group_col", "assay", "cutoffs", "seed")
    )
})
