test_that("unified input routing preserves matrix behavior", {
    fixture <- matrix_input_fixture()

    prepared <- prepare_input(fixture$x, fixture$group)

    expect_identical(
        prepared,
        prepare_matrix_input(fixture$x, fixture$group)
    )
    expect_error(
        prepare_input(as.data.frame(fixture$x), fixture$group),
        "`x` must be an ordinary numeric matrix or a SummarizedExperiment object",
        fixed = TRUE
    )
})

test_that("a single SummarizedExperiment assay is selected implicitly", {
    fixture <- summarized_experiment_input_fixture()
    original <- fixture$se

    prepared <- prepare_input(
        fixture$se,
        group_col = "condition"
    )

    expect_identical(prepared$data, fixture$x)
    expect_identical(
        prepared$groups_by_sample,
        setNames(as.character(fixture$group), colnames(fixture$x))
    )
    expect_identical(prepared$representation, "SummarizedExperiment")
    expect_identical(prepared$assay_index, 1L)
    expect_identical(prepared$assay_name, "intensity")
    expect_identical(fixture$se, original)
})

test_that("matrix and SummarizedExperiment adapters feed one core contract", {
    fixture <- summarized_experiment_input_fixture()

    matrix_prepared <- prepare_input(fixture$x, fixture$group)
    se_prepared <- prepare_input(fixture$se, group_col = "condition")

    expect_identical(se_prepared$data, matrix_prepared$data)
    expect_identical(
        se_prepared$groups_by_sample,
        matrix_prepared$groups_by_sample
    )
})

test_that("SummarizedExperiment assay selection is explicit and unambiguous", {
    fixture <- summarized_experiment_input_fixture(multiple_assays = TRUE)

    expect_error(
        prepare_input(fixture$se, group_col = "condition"),
        "`assay` must name one assay when `x` contains multiple assays",
        fixed = TRUE
    )

    prepared <- prepare_input(
        fixture$se,
        group_col = "condition",
        assay = "auxiliary"
    )
    expect_identical(prepared$data, fixture$x + 100)
    expect_identical(prepared$assay_index, 2L)
    expect_identical(prepared$assay_name, "auxiliary")

    expect_error(
        prepare_input(
            fixture$se,
            group_col = "condition",
            assay = "missing"
        ),
        "`assay` must match exactly one assay name",
        fixed = TRUE
    )
    expect_error(
        prepare_input(fixture$se, group_col = "condition", assay = 1L),
        "`assay` must be NULL or one non-empty assay name",
        fixed = TRUE
    )

    duplicated <- fixture$se
    SummarizedExperiment::assayNames(duplicated) <- c("same", "same")
    expect_error(
        prepare_input(duplicated, group_col = "condition", assay = "same"),
        "`assay` must match exactly one assay name",
        fixed = TRUE
    )
})

test_that("SummarizedExperiment routing validates metadata arguments", {
    fixture <- summarized_experiment_input_fixture()

    expect_error(
        prepare_input(
            fixture$se,
            group = fixture$group,
            group_col = "condition"
        ),
        "`group` must be NULL for SummarizedExperiment input",
        fixed = TRUE
    )
    expect_error(
        prepare_input(fixture$se),
        "`group_col` must name exactly one colData column",
        fixed = TRUE
    )
    expect_error(
        prepare_input(fixture$se, group_col = "missing"),
        "`group_col` must name exactly one colData column",
        fixed = TRUE
    )

    duplicated <- fixture$se
    names(SummarizedExperiment::colData(duplicated)) <- c(
        "condition",
        "condition"
    )
    expect_error(
        prepare_input(duplicated, group_col = "condition"),
        "`group_col` must name exactly one colData column",
        fixed = TRUE
    )

    missing_group <- fixture$se
    missing_group$condition[2] <- NA
    expect_error(
        prepare_input(missing_group, group_col = "condition"),
        "Condition labels must not be missing or empty",
        fixed = TRUE
    )
})

test_that("SummarizedExperiment assays satisfy the matrix value contract", {
    fixture <- summarized_experiment_input_fixture()

    nonnumeric <- fixture$se
    SummarizedExperiment::assay(nonnumeric, "intensity") <- matrix(
        "x",
        nrow = nrow(fixture$x),
        ncol = ncol(fixture$x),
        dimnames = dimnames(fixture$x)
    )
    expect_error(
        prepare_input(nonnumeric, group_col = "condition"),
        "The selected assay must be an ordinary numeric matrix",
        fixed = TRUE
    )

    for (case in list(nan = NaN, infinite = Inf, negative_infinite = -Inf)) {
        invalid <- fixture$se
        assay_data <- fixture$x
        assay_data[1, 1] <- case[[1L]]
        SummarizedExperiment::assay(invalid, "intensity") <- assay_data

        expected <- if (is.nan(case[[1L]])) {
            "The selected assay must not contain NaN values"
        } else {
            "The selected assay must not contain Inf or -Inf values"
        }
        expect_error(
            prepare_input(invalid, group_col = "condition"),
            expected,
            fixed = TRUE
        )
    }

    empty <- SummarizedExperiment::SummarizedExperiment(
        rowData = data.frame(row.names = c("p1", "p2")),
        colData = data.frame(
            condition = c("A", "B"),
            row.names = c("s1", "s2")
        )
    )
    expect_error(
        prepare_input(empty, group_col = "condition"),
        "`x` must contain at least one assay",
        fixed = TRUE
    )
})

test_that("output restoration preserves the selected representation", {
    matrix_fixture <- matrix_input_fixture()
    matrix_prepared <- prepare_input(
        matrix_fixture$x,
        matrix_fixture$group
    )
    expect_identical(
        restore_output_data(
            matrix_fixture$x,
            matrix_prepared,
            matrix_fixture$x
        ),
        matrix_fixture$x
    )

    fixture <- summarized_experiment_input_fixture(multiple_assays = TRUE)
    original <- fixture$se
    prepared <- prepare_input(
        fixture$se,
        group_col = "condition",
        assay = "intensity"
    )
    modified <- fixture$x["protein_b", , drop = FALSE]
    modified[1, 1] <- -10

    restored <- restore_output_data(modified, prepared, fixture$se)

    expect_s4_class(restored, "SummarizedExperiment")
    expect_identical(
        SummarizedExperiment::assay(restored, "intensity"),
        modified
    )
    expect_identical(
        SummarizedExperiment::assay(restored, "auxiliary"),
        SummarizedExperiment::assay(
            fixture$se,
            "auxiliary"
        )[rownames(modified), , drop = FALSE]
    )
    expect_identical(
        as.character(SummarizedExperiment::rowData(restored)$annotation),
        "second"
    )
    expect_identical(
        as.data.frame(SummarizedExperiment::colData(restored)),
        as.data.frame(SummarizedExperiment::colData(fixture$se))
    )
    expect_identical(
        methods::slot(restored, "metadata"),
        methods::slot(fixture$se, "metadata")
    )
    expect_identical(fixture$se, original)

    empty <- restore_output_data(
        fixture$x[integer(), , drop = FALSE],
        prepared,
        fixture$se
    )
    expect_identical(dim(empty), c(0L, ncol(fixture$x)))
    expect_identical(
        dim(SummarizedExperiment::assay(empty, "intensity")),
        c(0L, ncol(fixture$x))
    )
})

test_that("output restoration enforces original feature and sample order", {
    fixture <- summarized_experiment_input_fixture()
    prepared <- prepare_input(fixture$se, group_col = "condition")

    expect_error(
        restore_output_data(
            fixture$x[rev(rownames(fixture$x)), , drop = FALSE],
            prepared,
            fixture$se
        ),
        "Core output feature names must preserve original order",
        fixed = TRUE
    )
    expect_error(
        restore_output_data(
            fixture$x[, rev(colnames(fixture$x)), drop = FALSE],
            prepared,
            fixture$se
        ),
        "Core output sample names must match the input order",
        fixed = TRUE
    )
})
