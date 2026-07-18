base_fit_design <- function(fixture) {
    missingness_design(
        data.frame(
            condition = fixture$group,
            row.names = colnames(fixture$x)
        ),
        condition = "condition"
    )
}

base_fit_fixture <- function(seed = 19L) {
    fixture <- normative_fixture()
    groups <- stats::setNames(fixture$group, colnames(fixture$x))
    fit <- classify_missingness(
        fixture$x,
        group = groups,
        cutoffs = fixture$cutoffs,
        seed = seed
    )

    list(
        fixture = fixture,
        design = base_fit_design(fixture),
        fit = fit,
        seed = seed
    )
}

test_that("a recomputed compatible base fit is reused exactly", {
    candidate <- base_fit_fixture()
    original_x <- candidate$fixture$x
    original_design <- candidate$design
    original_fit <- candidate$fit
    analysis_call <- quote(
        analyze_missingness(x, design, base_fit = fit, seed = 19L)
    )

    analysis <- imputefinder:::.analyze_missingness_input_first(
        x = candidate$fixture$x,
        design = candidate$design,
        base_fit = candidate$fit,
        seed = candidate$seed,
        call = analysis_call
    )
    report <- imputefinder:::.compare_classic_fits(
        candidate$fit,
        analysis$classic
    )

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
    expect_identical(analysis$classic, candidate$fit)
    expect_identical(analysis$provenance$seeds, list(classic_rescue = 19L))
    expect_identical(
        analysis$sentinel$pre_rescue,
        imputefinder:::.new_pre_rescue_evidence(
            candidate$fixture$x,
            stats::setNames(
                candidate$fixture$group,
                colnames(candidate$fixture$x)
            )
        )
    )
    expect_identical(
        names(report),
        c(
            "schema",
            "compatible",
            "exact",
            "canonical",
            "tolerance",
            "numeric_tolerance"
        )
    )
    expect_identical(report$schema, "base_fit_compatibility_v1")
    expect_true(report$compatible)
    expect_identical(report$exact, character())
    expect_identical(report$canonical, character())
    expect_identical(report$tolerance, character())
    expect_identical(
        names(report$numeric_tolerance),
        c("absolute", "relative")
    )
    expect_identical(candidate$fixture$x, original_x)
    expect_identical(candidate$design, original_design)
    expect_identical(candidate$fit, original_fit)
    expect_invisible(imputefinder:::.validate_imputefinder_analysis(analysis))
})

test_that("base-fit recomputation preserves the caller RNG exactly", {
    candidate <- base_fit_fixture()

    with_preserved_random_state({
        RNGkind("Knuth-TAOCP-2002", "Box-Muller", "Rejection")
        set.seed(9127L)
        kind_before <- RNGkind()
        seed_before <- .Random.seed

        analysis <- imputefinder:::.analyze_missingness_input_first(
            candidate$fixture$x,
            candidate$design,
            base_fit = candidate$fit,
            seed = candidate$seed
        )

        expect_s3_class(analysis$classic, "imputefinder_result")
        expect_identical(RNGkind(), kind_before)
        expect_identical(.Random.seed, seed_before)
    })
})

test_that("base-fit cutoff policy reproduces manual and automatic sources", {
    manual <- base_fit_fixture()
    manual_policy <- imputefinder:::.base_fit_cutoff_policy(manual$fit)

    expect_identical(
        manual_policy,
        list(
            sources = c(A = "manual", B = "manual"),
            manual = c(A = 12, B = 12)
        )
    )

    automatic <- automatic_cutoff_matrix_fixture()
    design <- missingness_design(
        data.frame(
            condition = automatic$group,
            row.names = colnames(automatic$x)
        ),
        condition = "condition"
    )
    fit <- classify_missingness(
        automatic$x,
        group = automatic$group,
        seed = 7L
    )
    original_fit <- fit
    policy <- imputefinder:::.base_fit_cutoff_policy(fit)
    analysis <- imputefinder:::.analyze_missingness_input_first(
        automatic$x,
        design,
        base_fit = fit,
        seed = 7L
    )

    expect_identical(policy$sources, c(A = "automatic", B = "automatic"))
    expect_identical(policy$manual, stats::setNames(numeric(), character()))
    expect_identical(analysis$classic, fit)
    expect_identical(fit, original_fit)

    complete_x <- matrix(
        seq_len(16),
        nrow = 4L,
        dimnames = list(paste0("f", 1:4), paste0("s", 1:4))
    )
    complete_groups <- rep(c("A", "B"), each = 2L)
    complete_fit <- classify_missingness(
        complete_x,
        group = complete_groups
    )
    complete_policy <- imputefinder:::.base_fit_cutoff_policy(complete_fit)
    expect_identical(
        complete_policy$sources,
        c(A = "not_needed", B = "not_needed")
    )
    expect_identical(
        complete_policy$manual,
        stats::setNames(numeric(), character())
    )
})

test_that("base fit and explicit cutoffs are conflicting specifications", {
    candidate <- base_fit_fixture()
    original_fit <- candidate$fit

    error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            candidate$fixture$x,
            candidate$design,
            base_fit = candidate$fit,
            cutoffs = candidate$fixture$cutoffs,
            seed = candidate$seed
        ),
        imputefinder_base_fit_conflict_error = identity
    )

    expect_s3_class(error, "imputefinder_base_fit_conflict_error")
    expect_s3_class(error, "imputefinder_analysis_error")
    expect_identical(error$arguments, c("base_fit", "cutoffs"))
    expect_identical(candidate$fit, original_fit)
})

test_that("recomputation rejects all identity and seed drift", {
    candidate <- base_fit_fixture(seed = 1L)

    changed_x <- candidate$fixture$x
    changed_x["mar_both", "s1"] <- changed_x["mar_both", "s1"] + 0.5
    original_changed_x <- changed_x
    original_fit <- candidate$fit
    input_error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            changed_x,
            candidate$design,
            base_fit = candidate$fit,
            seed = 1L
        ),
        imputefinder_base_fit_mismatch_error = identity
    )
    expect_s3_class(input_error, "imputefinder_base_fit_mismatch_error")
    expect_false(input_error$report$compatible)
    expect_true("data" %in% input_error$report$exact)
    expect_true("classifications" %in% input_error$report$tolerance)
    expect_identical(changed_x, original_changed_x)
    expect_identical(candidate$fit, original_fit)

    changed_design <- candidate$design
    changed_design$sample_data$condition[1:2] <- "C"
    design_error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            candidate$fixture$x,
            changed_design,
            base_fit = candidate$fit,
            seed = 1L
        ),
        imputefinder_base_fit_mismatch_error = identity
    )
    expect_s3_class(design_error, "imputefinder_base_fit_mismatch_error")
    expect_s3_class(
        design_error$recomputed,
        "imputefinder_classic_failure"
    )
    expect_identical(design_error$report$exact, "recomputation")

    se <- rich_summarized_experiment(
        candidate$fixture$x,
        candidate$fixture$group
    )
    representation_error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            se,
            candidate$design,
            assay = "intensity",
            base_fit = candidate$fit,
            seed = 1L
        ),
        imputefinder_base_fit_mismatch_error = identity
    )
    expect_s3_class(
        representation_error,
        "imputefinder_base_fit_mismatch_error"
    )
    expect_true("data" %in% representation_error$report$exact)

    seed_error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            candidate$fixture$x,
            candidate$design,
            base_fit = candidate$fit,
            seed = 2L
        ),
        imputefinder_base_fit_mismatch_error = identity
    )
    expect_s3_class(seed_error, "imputefinder_base_fit_mismatch_error")
    expect_true(any(c("data", "seed_log") %in% seed_error$report$exact))

    zero_x <- matrix(
        c(0, 1, 2, 3),
        nrow = 1L,
        dimnames = list("f1", paste0("z", 1:4))
    )
    zero_groups <- rep(c("A", "B"), each = 2L)
    zero_design <- missingness_design(
        data.frame(
            condition = zero_groups,
            row.names = colnames(zero_x)
        ),
        condition = "condition"
    )
    zero_fit <- classify_missingness(zero_x, group = zero_groups)
    zero_x[[1L]] <- -0
    zero_error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            zero_x,
            zero_design,
            base_fit = zero_fit
        ),
        imputefinder_base_fit_mismatch_error = identity
    )
    expect_s3_class(zero_error, "imputefinder_base_fit_mismatch_error")
    expect_identical(zero_error$report$exact, "data")
})

test_that("malformed base fits fail through the specialized schema boundary", {
    candidate <- base_fit_fixture()
    partial <- candidate$fit
    partial$profiles <- NULL

    error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            candidate$fixture$x,
            candidate$design,
            base_fit = partial,
            seed = candidate$seed
        ),
        imputefinder_base_fit_schema_error = identity
    )

    expect_s3_class(error, "imputefinder_base_fit_schema_error")
    expect_s3_class(error, "imputefinder_analysis_error")
    expect_true("imputefinder_analysis_schema_error" %in% error$condition_class)
    expect_match(error$condition_message, "complete imputefinder result")
})

test_that("the comparator separates exact, canonical, and tolerance drift", {
    candidate <- base_fit_fixture()
    reference <- candidate$fit

    discrete <- reference
    discrete$classifications$state[[1L]] <- "MAR"
    discrete_report <- imputefinder:::.compare_classic_fits(
        discrete,
        reference
    )
    expect_false(discrete_report$compatible)
    expect_identical(discrete_report$exact, "classifications")

    canonical <- reference
    canonical$profiles$A$metadata$warnings <- "changed"
    canonical_report <- imputefinder:::.compare_classic_fits(
        canonical,
        reference
    )
    expect_false(canonical_report$compatible)
    expect_identical(canonical_report$canonical, "profiles")

    tolerance <- unname(
        imputefinder:::.classic_compatibility_tolerance()[["absolute"]]
    )
    close <- reference
    close$profiles$A$grid$intensity[[1L]] <-
        close$profiles$A$grid$intensity[[1L]] + tolerance / 4
    close_report <- imputefinder:::.compare_classic_fits(close, reference)
    expect_true(close_report$compatible)
    expect_identical(close_report$tolerance, character())

    distant <- reference
    distant$profiles$A$grid$intensity[[1L]] <-
        distant$profiles$A$grid$intensity[[1L]] + tolerance * 100
    distant_report <- imputefinder:::.compare_classic_fits(
        distant,
        reference
    )
    expect_false(distant_report$compatible)
    expect_true(any(startsWith(distant_report$tolerance, "profiles")))
})

test_that("scientific mismatches fail with the complete comparison report", {
    candidate <- base_fit_fixture()
    changed <- candidate$fit
    changed$groups$A$MNAR <- character()

    error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            candidate$fixture$x,
            candidate$design,
            base_fit = changed,
            seed = candidate$seed
        ),
        imputefinder_base_fit_mismatch_error = identity
    )

    expect_s3_class(error, "imputefinder_base_fit_mismatch_error")
    expect_identical(error$report$schema, "base_fit_compatibility_v1")
    expect_false(error$report$compatible)
    expect_identical(error$report$exact, "groups")
    expect_identical(
        error$differences,
        list(
            exact = "groups",
            canonical = character(),
            tolerance = character()
        )
    )
})

test_that("malformed cutoff policy fails before scientific comparison", {
    candidate <- base_fit_fixture()
    malformed <- candidate$fit
    malformed$cutoff_diagnostics$A$source <- "estimated"

    error <- tryCatch(
        imputefinder:::.analyze_missingness_input_first(
            candidate$fixture$x,
            candidate$design,
            base_fit = malformed,
            seed = candidate$seed
        ),
        imputefinder_base_fit_policy_error = identity
    )

    expect_s3_class(error, "imputefinder_base_fit_policy_error")
    expect_identical(error$condition, "A")
    expect_identical(error$source, "estimated")
})

test_that("fit-only input artifacts are portable unavailable records", {
    candidate <- base_fit_fixture()
    original_fit <- candidate$fit
    artifacts <- c(
        "dropped_row_cells",
        "original_mask",
        "pre_rescue_evidence"
    )

    unavailable <- lapply(
        artifacts,
        function(artifact) {
            imputefinder:::.fit_only_artifact(candidate$fit, artifact)
        }
    )

    for (index in seq_along(artifacts)) {
        value <- unavailable[[index]]
        expect_s3_class(value, "imputefinder_unavailable")
        expect_identical(
            names(value),
            c("status", "quantity", "code", "message", "requires")
        )
        expect_identical(value$status, "unavailable")
        expect_identical(value$quantity, artifacts[[index]])
        expect_identical(value$code, "original_input_required")
        expect_identical(value$requires, c("x", "missingness_design"))
        expect_match(value$message, "original `x`", fixed = TRUE)
        expect_identical(
            unserialize(serialize(value, NULL, version = 3L)),
            value
        )
        expect_invisible(imputefinder:::.validate_unavailable(value))
    }
    expect_identical(candidate$fit, original_fit)
})
