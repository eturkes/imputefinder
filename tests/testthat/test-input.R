test_that("unnamed matrix groups align positionally", {
    fixture <- matrix_input_fixture()
    original <- fixture$x

    prepared <- prepare_matrix_input(fixture$x, fixture$group)

    expect_identical(prepared$data, fixture$x)
    expect_identical(
        prepared$groups_by_sample,
        setNames(fixture$group, colnames(fixture$x))
    )
    expect_identical(prepared$representation, "matrix")
    expect_identical(fixture$x, original)
})

test_that("named matrix groups align by sample name", {
    fixture <- matrix_input_fixture()
    group <- c(s3 = "B", s2 = "B", s1 = "A")

    prepared <- prepare_matrix_input(fixture$x, group)

    expect_identical(
        prepared$groups_by_sample,
        c(s2 = "B", s1 = "A", s3 = "B")
    )
})

test_that("design rows align by sample name", {
    fixture <- matrix_input_fixture()
    design <- data.frame(
        batch = c(2L, 1L, 1L),
        condition = factor(c("B", "B", "A")),
        row.names = c("s3", "s2", "s1")
    )

    prepared <- prepare_matrix_input(
        fixture$x,
        design,
        group_col = "condition"
    )

    expect_identical(
        prepared$groups_by_sample,
        c(s2 = "B", s1 = "A", s3 = "B")
    )
})

test_that("matrix data require finite numeric values and stable identifiers", {
    fixture <- matrix_input_fixture()

    expect_error(
        prepare_matrix_input(as.data.frame(fixture$x), fixture$group),
        "`x` must be an ordinary numeric matrix",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(matrix("x", 1, 1), "A"),
        "`x` must be an ordinary numeric matrix",
        fixed = TRUE
    )

    invalid <- fixture$x
    invalid[1, 1] <- NaN
    expect_error(
        prepare_matrix_input(invalid, fixture$group),
        "`x` must not contain NaN values",
        fixed = TRUE
    )

    for (value in c(Inf, -Inf)) {
        invalid[1, 1] <- value
        expect_error(
            prepare_matrix_input(invalid, fixture$group),
            "`x` must not contain Inf or -Inf values",
            fixed = TRUE
        )
    }

    invalid <- fixture$x
    rownames(invalid) <- NULL
    expect_error(
        prepare_matrix_input(invalid, fixture$group),
        "`x` must have non-empty, unique feature names",
        fixed = TRUE
    )

    invalid <- fixture$x
    rownames(invalid)[2] <- rownames(invalid)[1]
    expect_error(
        prepare_matrix_input(invalid, fixture$group),
        "`x` must have non-empty, unique feature names",
        fixed = TRUE
    )

    invalid <- fixture$x
    rownames(invalid)[1] <- ""
    expect_error(
        prepare_matrix_input(invalid, fixture$group),
        "`x` must have non-empty, unique feature names",
        fixed = TRUE
    )

    invalid <- fixture$x
    colnames(invalid) <- NULL
    expect_error(
        prepare_matrix_input(invalid, fixture$group),
        "`x` must have non-empty, unique sample names",
        fixed = TRUE
    )

    invalid <- fixture$x
    colnames(invalid)[2] <- colnames(invalid)[1]
    expect_error(
        prepare_matrix_input(invalid, fixture$group),
        "`x` must have non-empty, unique sample names",
        fixed = TRUE
    )

    invalid <- fixture$x
    colnames(invalid)[1] <- ""
    expect_error(
        prepare_matrix_input(invalid, fixture$group),
        "`x` must have non-empty, unique sample names",
        fixed = TRUE
    )

    empty <- matrix(numeric(), nrow = 0L, ncol = 0L)
    expect_error(
        prepare_matrix_input(empty, character()),
        "`x` must have non-empty, unique feature names",
        fixed = TRUE
    )
})

test_that("the public path rejects unsupported special-value inputs", {
    fixture <- matrix_input_fixture()

    for (value in list(NaN, Inf, -Inf)) {
        invalid <- fixture$x
        invalid[1L, 1L] <- value
        expect_error(
            classify_missingness(
                invalid,
                fixture$group,
                cutoffs = c(A = 1, B = 1)
            ),
            if (is.nan(value)) "must not contain NaN" else
                "must not contain Inf or -Inf",
            fixed = TRUE
        )
    }

    no_condition_minimum <- fixture$x
    no_condition_minimum[, fixture$group == "A"] <- NA_real_
    expect_error(
        classify_missingness(
            no_condition_minimum,
            fixture$group,
            cutoffs = c(A = 1, B = 1)
        ),
        "Condition `A` has no finite intensity",
        fixed = TRUE
    )
})

test_that("atomic group validation is explicit", {
    fixture <- matrix_input_fixture()

    expect_error(
        prepare_matrix_input(fixture$x, NULL),
        "`group` is required for matrix input",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(fixture$x, c("A", "B")),
        "Unnamed `group` must have one value per sample",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(
            fixture$x,
            c(s1 = "A", s2 = "B", other = "B")
        ),
        "Named `group` must have unique names matching all sample names exactly",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(
            fixture$x,
            setNames(fixture$group, c("s2", "", "s3"))
        ),
        "Named `group` must have unique names matching all sample names exactly",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(
            fixture$x,
            c(s1 = "A", s1 = "A", s3 = "B")
        ),
        "Named `group` must have unique names matching all sample names exactly",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(fixture$x, c("A", NA, "B")),
        "Condition labels must not be missing or empty",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(fixture$x, c("A", "", "B")),
        "Condition labels must not be missing or empty",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(fixture$x, fixture$group, group_col = "condition"),
        "`group_col` must be NULL when `group` is an atomic vector",
        fixed = TRUE
    )
})

test_that("design group validation is explicit", {
    fixture <- matrix_input_fixture()
    design <- data.frame(
        condition = c("B", "A", "B"),
        row.names = c("s2", "s1", "s3")
    )

    expect_error(
        prepare_matrix_input(fixture$x, design),
        "`group_col` must name exactly one design column",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(fixture$x, design, group_col = "missing"),
        "`group_col` must name exactly one design column",
        fixed = TRUE
    )

    duplicate_columns <- cbind(design, design, deparse.level = 0)
    names(duplicate_columns) <- c("condition", "condition")
    expect_error(
        prepare_matrix_input(
            fixture$x,
            duplicate_columns,
            group_col = "condition"
        ),
        "`group_col` must name exactly one design column",
        fixed = TRUE
    )

    mismatched <- design
    rownames(mismatched)[3] <- "other"
    expect_error(
        prepare_matrix_input(fixture$x, mismatched, group_col = "condition"),
        "Design row names must be unique and match all sample names exactly",
        fixed = TRUE
    )

    duplicated <- design
    attr(duplicated, "row.names") <- c("s2", "s2", "s3")
    expect_error(
        prepare_matrix_input(fixture$x, duplicated, group_col = "condition"),
        "Design row names must be unique and match all sample names exactly",
        fixed = TRUE
    )

    design$condition[2] <- NA
    expect_error(
        prepare_matrix_input(fixture$x, design, group_col = "condition"),
        "Condition labels must not be missing or empty",
        fixed = TRUE
    )

    design$condition <- I(list("B", "A", "B"))
    expect_error(
        prepare_matrix_input(fixture$x, design, group_col = "condition"),
        "The selected design condition column must be an atomic vector",
        fixed = TRUE
    )
})

test_that("matrix routing rejects assay selection and unsupported group objects", {
    fixture <- matrix_input_fixture()

    expect_error(
        prepare_matrix_input(fixture$x, fixture$group, assay = "intensity"),
        "`assay` must be NULL for matrix input",
        fixed = TRUE
    )
    expect_error(
        prepare_matrix_input(fixture$x, list("A", "B", "B")),
        "`group` must be an atomic vector or data frame for matrix input",
        fixed = TRUE
    )
})
