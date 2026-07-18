public_analysis_fixture <- function() {
    fixture <- normative_fixture()
    design <- missingness_design(
        data.frame(
            condition = fixture$group,
            row.names = colnames(fixture$x)
        ),
        condition = "condition"
    )

    list(fixture = fixture, design = design)
}

test_that("the public analyzer exposes the frozen complete signature", {
    candidate <- public_analysis_fixture()
    original_x <- candidate$fixture$x
    original_design <- candidate$design

    analysis <- analyze_missingness(
        candidate$fixture$x,
        candidate$design,
        cutoffs = candidate$fixture$cutoffs,
        seed = 19L
    )

    expect_true("analyze_missingness" %in% getNamespaceExports("imputefinder"))
    expect_identical(
        names(formals(analyze_missingness)),
        c(
            "x",
            "design",
            "assay",
            "base_fit",
            "cutoffs",
            "seed",
            "modules"
        )
    )
    expect_identical(
        formals(analyze_missingness)$modules,
        quote(c("sentinel", "stability"))
    )
    expect_s3_class(analysis, "imputefinder_analysis")
    expect_s3_class(analysis$classic, "imputefinder_result")
    expect_identical(names(analysis$sentinel), "pre_rescue")
    expect_s3_class(analysis$stability, "imputefinder_unavailable")
    expect_identical(analysis$stability$status, "unavailable")
    expect_identical(analysis$stability$quantity, "stability")
    expect_identical(
        analysis$stability$code,
        "module_pending_validation"
    )
    expect_identical(
        analysis$stability$requires,
        "robustness_certificate_v1"
    )
    expect_identical(
        analysis$classic$call,
        analysis$provenance$call
    )
    expect_identical(candidate$fixture$x, original_x)
    expect_identical(candidate$design, original_design)
    expect_invisible(imputefinder:::.validate_imputefinder_analysis(analysis))

    printed <- capture.output(visibility <- withVisible(print(analysis)))
    expect_false(visibility$visible)
    expect_identical(visibility$value, analysis)
    expect_lte(length(printed), 6L)
    expect_match(printed[[1L]], "<imputefinder_analysis>", fixed = TRUE)
    expect_true(any(grepl("stability: unavailable", printed, fixed = TRUE)))
    expect_false(any(grepl("globally_absent", printed, fixed = TRUE)))
})

test_that("module selection routes selected, unavailable, and skipped slots", {
    candidate <- public_analysis_fixture()
    arguments <- list(
        x = candidate$fixture$x,
        design = candidate$design,
        cutoffs = candidate$fixture$cutoffs
    )

    sentinel <- do.call(
        analyze_missingness,
        c(arguments, list(modules = "sentinel"))
    )
    stability <- do.call(
        analyze_missingness,
        c(arguments, list(modules = "stability"))
    )
    skipped <- do.call(
        analyze_missingness,
        c(arguments, list(modules = character()))
    )

    expect_type(sentinel$sentinel, "list")
    expect_null(sentinel$stability)
    expect_null(stability$sentinel)
    expect_s3_class(stability$stability, "imputefinder_unavailable")
    expect_null(skipped$sentinel)
    expect_null(skipped$stability)
    expect_identical(
        sentinel$classic$feature_status,
        stability$classic$feature_status
    )
    expect_identical(
        stability$classic$feature_status,
        skipped$classic$feature_status
    )
})

test_that("module requests reject ambiguity through a typed boundary", {
    candidate <- public_analysis_fixture()
    invalid <- list(
        NULL,
        factor("sentinel"),
        c("sentinel", NA_character_),
        c("sentinel", "sentinel"),
        "unknown"
    )

    for (modules in invalid) {
        expect_error(
            analyze_missingness(
                candidate$fixture$x,
                candidate$design,
                cutoffs = candidate$fixture$cutoffs,
                modules = modules
            ),
            class = "imputefinder_analysis_module_error"
        )
    }

    error <- tryCatch(
        analyze_missingness(
            candidate$fixture$x,
            candidate$design,
            cutoffs = candidate$fixture$cutoffs,
            modules = "unknown"
        ),
        imputefinder_analysis_module_error = identity
    )
    expect_s3_class(error, "imputefinder_analysis_error")
    expect_identical(error$supported, c("sentinel", "stability"))
    expect_identical(error$problem, "unknown")
})

test_that("selected pre-rescue evidence survives public classic failure", {
    candidate <- public_analysis_fixture()

    analysis <- analyze_missingness(
        candidate$fixture$x,
        candidate$design,
        modules = c("stability", "sentinel")
    )

    expect_s3_class(analysis$classic, "imputefinder_classic_failure")
    expect_identical(analysis$classic$stage, "cutoff")
    expect_identical(names(analysis$sentinel), "pre_rescue")
    expect_s3_class(analysis$stability, "imputefinder_unavailable")
    expect_identical(
        analysis$provenance$failures,
        list(classic = analysis$classic)
    )
})

test_that("stored module status rejects corruption and invented output", {
    candidate <- public_analysis_fixture()
    analysis <- analyze_missingness(
        candidate$fixture$x,
        candidate$design,
        cutoffs = candidate$fixture$cutoffs
    )

    wrong_status <- analysis
    wrong_status$stability$code <- "available"
    expect_error(
        imputefinder:::.validate_imputefinder_analysis(wrong_status),
        class = "imputefinder_analysis_schema_error"
    )

    invented <- analysis
    invented$stability <- list(schema = "robustness_certificate_v1")
    expect_error(
        imputefinder:::.validate_imputefinder_analysis(invented),
        class = "imputefinder_analysis_schema_error"
    )
})

test_that("the public base-fit path preserves its compatible result exactly", {
    candidate <- public_analysis_fixture()
    fit <- classify_missingness(
        candidate$fixture$x,
        candidate$fixture$group,
        cutoffs = candidate$fixture$cutoffs,
        seed = 19L
    )
    original_fit <- fit

    analysis <- analyze_missingness(
        candidate$fixture$x,
        candidate$design,
        base_fit = fit,
        seed = 19L,
        modules = character()
    )

    expect_identical(analysis$classic, fit)
    expect_identical(fit, original_fit)
    expect_null(analysis$sentinel)
    expect_null(analysis$stability)
    expect_error(
        analyze_missingness(
            candidate$fixture$x,
            candidate$design,
            base_fit = fit,
            cutoffs = candidate$fixture$cutoffs,
            seed = 19L
        ),
        class = "imputefinder_base_fit_conflict_error"
    )
})
