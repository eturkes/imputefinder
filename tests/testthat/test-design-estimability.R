design_core_crossed_fixture <- function(perfect = FALSE) {
    sample <- paste0("s", seq_len(8L))
    condition <- rep(c("A", "B"), each = 4L)
    batch <- if (perfect) {
        rep(c("x", "y"), each = 4L)
    } else {
        rep(c("x", "y"), times = 4L)
    }
    metadata <- data.frame(
        condition = condition,
        batch = batch,
        row.names = sample,
        stringsAsFactors = FALSE
    )

    missingness_design(
        metadata,
        condition = "condition",
        nuisance = "batch"
    )
}

design_core_paired_fixture <- function() {
    metadata <- expand.grid(
        technical = c("t1", "t2"),
        subject = paste0("p", seq_len(3L)),
        condition = c("A", "B"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    rownames(metadata) <- sprintf("sample_%02d", seq_len(nrow(metadata)))

    missingness_design(
        metadata,
        condition = "condition",
        nuisance = "technical",
        block = "subject"
    )
}

test_that("declared designs produce one canonical algebraic model", {
    metadata <- data.frame(
        condition = factor(
            rep(c("B", "A"), each = 4L),
            levels = c("B", "A")
        ),
        batch = factor(rep(c("y", "x"), times = 4L)),
        run_order = c(8L, 6L, 4L, 2L, 7L, 5L, 3L, 1L),
        subject = rep(paste0("p", 1:4), 2L),
        acquisition = rep("DIA", 8L),
        row.names = paste0("s", c(8L, 6L, 4L, 2L, 7L, 5L, 3L, 1L))
    )
    design <- missingness_design(
        metadata,
        condition = "condition",
        nuisance = c("run_order", "batch"),
        block = "subject",
        acquisition = "acquisition",
        interactions = list(
            c("condition", "run_order"),
            c("batch", "condition")
        )
    )

    core <- imputefinder:::.new_design_estimability(design)

    expect_identical(
        names(core),
        c("schema", "methods", "model", "rank", "aliasing", "units")
    )
    expect_identical(core$schema, "design_estimability_v1")
    expect_identical(
        core$methods,
        c(
            encoding = "canonical_treatment_contrasts_v1",
            rank = "svd_relative_v1",
            aliasing = "ordered_null_projector_mgs_v1",
            contrast = "row_space_residual_v1",
            units = "declared_block_or_sample_v1"
        )
    )
    expect_identical(
        names(core$model),
        c("matrix", "variables", "terms", "columns")
    )
    expect_identical(
        rownames(core$model$matrix),
        paste0("s", seq_len(8L))
    )
    expect_identical(
        core$model$variables$column,
        c("condition", "batch", "run_order", "subject", "acquisition")
    )
    expect_identical(
        core$model$variables$role,
        c("condition", "nuisance", "nuisance", "block", "acquisition")
    )
    expect_identical(
        core$model$terms$term,
        c(
            "(Intercept)", "condition", "batch", "run_order", "subject",
            "batch:condition", "condition:run_order"
        )
    )
    expect_identical(
        core$model$terms$coefficient_count,
        as.integer(table(factor(
            core$model$columns$term_id,
            levels = core$model$terms$term_id
        )))
    )
    expect_identical(
        colnames(core$model$matrix),
        core$model$columns$coefficient
    )
    expect_true(all(core$model$matrix[, "coef_0001"] == 1))
    expect_identical(
        core$rank$threshold,
        max(dim(core$model$matrix)) * .Machine$double.eps *
            core$rank$singular_values[[1L]]
    )
    expect_identical(names(core$rank$leverage), rownames(core$model$matrix))
    expect_invisible(
        imputefinder:::.validate_design_estimability(core, design)
    )
})

test_that("SVD aliases and row-space contrasts separate crossed from aliased", {
    crossed <- imputefinder:::.new_design_estimability(
        design_core_crossed_fixture()
    )
    aliased <- imputefinder:::.new_design_estimability(
        design_core_crossed_fixture(perfect = TRUE)
    )

    expect_identical(crossed$rank$rank, 3L)
    expect_identical(crossed$rank$nullity, 0L)
    expect_true(crossed$rank$full_column_rank)
    expect_identical(nrow(crossed$aliasing$affected_coefficients), 0L)
    expect_identical(dim(crossed$aliasing$null_basis), c(3L, 0L))

    expect_identical(aliased$rank$rank, 2L)
    expect_identical(aliased$rank$nullity, 1L)
    expect_false(aliased$rank$full_column_rank)
    expect_identical(
        aliased$aliasing$affected_terms$term,
        c("condition", "batch")
    )
    expect_equal(
        unname(aliased$aliasing$null_basis[, 1L]),
        c(0, 1 / sqrt(2), -1 / sqrt(2)),
        tolerance = 1e-12
    )

    estimable <- imputefinder:::.design_contrast_estimability(
        crossed,
        list(condition = c(B = 1, A = -1))
    )
    nonestimable <- imputefinder:::.design_contrast_estimability(
        aliased,
        list(condition = c(A = -1, B = 1))
    )

    expect_identical(
        names(estimable),
        c(
            "schema", "contrast", "coefficient", "residual",
            "contrast_norm", "residual_norm", "tolerance", "estimable",
            "affected_terms"
        )
    )
    expect_identical(estimable$schema, "design_contrast_estimability_v1")
    expect_true(estimable$estimable)
    expect_lte(estimable$residual_norm, estimable$tolerance)
    expect_false(nonestimable$estimable)
    expect_gt(nonestimable$residual_norm, nonestimable$tolerance)
    expect_identical(
        nonestimable$affected_terms$term,
        c("condition", "batch")
    )
    expect_identical(
        nonestimable$contrast$weights,
        c(A = -1, B = 1)
    )

    raw <- stats::setNames(
        c(0, 1, 0),
        colnames(crossed$model$matrix)
    )
    expect_true(
        imputefinder:::.design_contrast_estimability(crossed, raw)$estimable
    )
})

test_that("block units keep paired technical siblings inseparable", {
    core <- imputefinder:::.new_design_estimability(
        design_core_paired_fixture()
    )

    expect_identical(core$units$resampling_role, "block")
    expect_identical(core$units$resampling_column, "subject")
    expect_identical(core$units$independent_unit_count, 3L)
    expect_true(core$units$grouped_technical_siblings)
    expect_true(all(core$units$sample$weight == 1L))
    expect_true(all(core$units$sample$unit_sample_count == 4L))
    expect_true(all(core$units$sample$unit_condition_sample_count == 2L))
    expect_identical(
        core$units$unit$sample_count,
        rep(4L, 3L)
    )
    expect_identical(
        core$units$unit$condition_count,
        rep(2L, 3L)
    )
    expect_identical(
        core$units$condition$sample_count,
        c(6L, 6L)
    )
    expect_identical(
        core$units$condition$unit_count,
        c(3L, 3L)
    )
})

test_that("independent unequal replication retains every named design row", {
    metadata <- data.frame(
        condition = c("C", "B", "A", "B", "A", "B"),
        row.names = paste0("s", c(6L, 5L, 1L, 3L, 2L, 4L))
    )
    design <- missingness_design(metadata, condition = "condition")

    core <- imputefinder:::.new_design_estimability(design)

    expect_identical(core$units$resampling_role, "sample")
    expect_identical(core$units$resampling_column, NA_character_)
    expect_identical(core$units$independent_unit_count, 6L)
    expect_false(core$units$grouped_technical_siblings)
    expect_identical(core$units$sample$sample, paste0("s", 1:6))
    expect_identical(core$units$sample$unit, paste0("s", 1:6))
    expect_identical(core$units$sample$weight, rep(1L, 6L))
    expect_identical(core$units$condition$condition, c("A", "B", "C"))
    expect_identical(core$units$condition$sample_count, c(2L, 3L, 1L))
    expect_identical(core$units$condition$unit_count, c(2L, 3L, 1L))
    expect_identical(
        rownames(core$model$matrix),
        core$units$sample$sample
    )
})

test_that("canonical core ignores sample order and label-preserving encodings", {
    metadata <- data.frame(
        condition = factor(
            c("B", "A", "B", "A", "B", "A"),
            levels = c("B", "A")
        ),
        batch = factor(
            c("y", "x", "x", "y", "y", "x"),
            levels = c("y", "x")
        ),
        run_order = c(6L, 1L, 4L, 3L, 2L, 5L),
        subject = rep(paste0("p", 1:3), 2L),
        row.names = paste0("s", c(6L, 1L, 4L, 3L, 2L, 5L))
    )
    baseline <- missingness_design(
        metadata,
        condition = "condition",
        nuisance = c("run_order", "batch"),
        block = "subject",
        interactions = list(
            c("condition", "run_order"),
            c("batch", "condition")
        )
    )

    reencoded <- metadata[rev(rownames(metadata)), rev(names(metadata)),
        drop = FALSE
    ]
    reencoded$condition <- as.character(reencoded$condition)
    reencoded$batch <- as.character(reencoded$batch)
    reencoded$run_order <- as.double(reencoded$run_order)
    reordered <- missingness_design(
        reencoded,
        condition = "condition",
        nuisance = c("batch", "run_order"),
        block = "subject",
        interactions = rev(list(
            c("run_order", "condition"),
            c("condition", "batch")
        ))
    )

    expect_identical(
        imputefinder:::.new_design_estimability(baseline),
        imputefinder:::.new_design_estimability(reordered)
    )
})

test_that("design-core schemas reject drift and malformed contrasts", {
    design <- design_core_crossed_fixture()
    core <- imputefinder:::.new_design_estimability(design)

    wrong_schema <- core
    wrong_schema$schema <- "design_estimability_v2"
    error <- tryCatch(
        imputefinder:::.validate_design_estimability(wrong_schema, design),
        imputefinder_design_estimability_lifecycle_error = identity
    )
    expect_s3_class(
        error,
        "imputefinder_design_estimability_lifecycle_error"
    )
    expect_identical(error$expected_schema, "design_estimability_v1")
    expect_identical(error$actual_schema, "design_estimability_v2")

    for (mutator in list(
        function(x) { x$rank$rank <- x$rank$rank - 1L; x },
        function(x) { x$model$matrix[[1L]] <- 0; x },
        function(x) { x$units$sample$unit[[1L]] <- "other"; x }
    )) {
        expect_error(
            imputefinder:::.validate_design_estimability(
                mutator(core),
                design
            ),
            class = "imputefinder_design_estimability_error"
        )
    }

    invalid <- list(
        list(condition = c(A = 1, B = 1)),
        list(condition = c(A = -1, absent = 1)),
        list(batch = c(x = -1, y = 1), condition = c(A = -1, B = 1)),
        c(unknown = 1)
    )
    for (contrast in invalid) {
        expect_error(
            imputefinder:::.design_contrast_estimability(core, contrast),
            class = "imputefinder_design_contrast_error"
        )
    }
})

test_that("input-first analysis stores the mandatory core without side effects", {
    fixture <- normative_fixture()
    metadata <- data.frame(
        condition = fixture$group,
        batch = rep(c("x", "y"), 4L),
        row.names = colnames(fixture$x),
        stringsAsFactors = FALSE
    )
    design <- missingness_design(
        metadata[rev(rownames(metadata)), , drop = FALSE],
        condition = "condition",
        nuisance = "batch"
    )
    original_x <- fixture$x
    original_design <- design
    classic_before <- classify_normative_fixture()
    caller_kind <- RNGkind()
    caller_seed <- get0(
        ".Random.seed",
        envir = globalenv(),
        inherits = FALSE,
        ifnotfound = NULL
    )
    caller_options <- options()
    caller_directory <- getwd()

    analysis <- analyze_missingness(
        fixture$x,
        design,
        cutoffs = fixture$cutoffs,
        modules = character()
    )
    classic_after <- classify_normative_fixture()

    expect_identical(fixture$x, original_x)
    expect_identical(design, original_design)
    expect_identical(classic_after, classic_before)
    expect_identical(RNGkind(), caller_kind)
    expect_identical(
        get0(
            ".Random.seed",
            envir = globalenv(),
            inherits = FALSE,
            ifnotfound = NULL
        ),
        caller_seed
    )
    expect_identical(options(), caller_options)
    expect_identical(getwd(), caller_directory)
    expect_null(analysis$sentinel)
    expect_s3_class(analysis$classic, "imputefinder_result")
    expect_identical(
        analysis$design$estimability,
        imputefinder:::.new_design_estimability(
            analysis$design$declared
        )
    )
    expect_invisible(
        imputefinder:::.validate_imputefinder_analysis(analysis)
    )
})
