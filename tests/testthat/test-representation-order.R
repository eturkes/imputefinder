expect_core_result_equivalent <- function(matrix_result, experiment_result) {
    components <- c(
        "classifications",
        "groups",
        "feature_status",
        "cutoffs",
        "cutoff_diagnostics",
        "profiles",
        "seed_log",
        "groups_by_sample"
    )
    for (component in components) {
        expect_identical(
            experiment_result[[component]],
            matrix_result[[component]],
            info = component
        )
    }
}

expect_experiment_representation <- function(
    experiment_result,
    matrix_result,
    original
) {
    output <- experiment_result$data
    retained <- matrix_result$feature_status$feature[
        matrix_result$feature_status$retained
    ]

    expect_identical(class(output), class(original))
    expect_identical(
        SummarizedExperiment::assayNames(output),
        SummarizedExperiment::assayNames(original)
    )
    expect_identical(
        SummarizedExperiment::assay(output, "intensity"),
        matrix_result$data
    )
    expect_identical(
        SummarizedExperiment::assay(output, "reference"),
        SummarizedExperiment::assay(original, "reference")[
            retained,
            ,
            drop = FALSE
        ]
    )
    expect_identical(
        SummarizedExperiment::rowData(output),
        SummarizedExperiment::rowData(original)[retained, , drop = FALSE]
    )
    expect_identical(
        SummarizedExperiment::colData(output),
        SummarizedExperiment::colData(original)
    )
    expect_identical(
        methods::slot(output, "metadata"),
        methods::slot(original, "metadata")
    )
    expect_identical(rownames(output), retained)
    expect_identical(colnames(output), colnames(original))
}

test_that("manual matrix and experiment paths preserve ordered representation", {
    fixture <- normative_fixture()
    feature_order <- c(7L, 1L, 5L, 2L, 6L, 4L, 3L)
    sample_order <- c(8L, 1L, 5L, 4L, 7L, 2L, 6L, 3L)
    x <- fixture$x[feature_order, sample_order, drop = FALSE]
    group <- stats::setNames(
        fixture$group[sample_order],
        colnames(x)
    )
    experiment <- rich_summarized_experiment(x, unname(group))
    original_matrix <- x
    original_experiment <- experiment

    matrix_result <- classify_missingness(
        x,
        group,
        cutoffs = c(B = 12, A = 12),
        seed = 17L
    )
    experiment_result <- classify_missingness(
        experiment,
        group_col = "condition",
        assay = "intensity",
        cutoffs = c(B = 12, A = 12),
        seed = 17L
    )

    expect_identical(x, original_matrix)
    expect_identical(experiment, original_experiment)
    expect_core_result_equivalent(matrix_result, experiment_result)
    expect_experiment_representation(
        experiment_result,
        matrix_result,
        original_experiment
    )
    expect_identical(names(matrix_result$groups), c("A", "B"))
    expect_identical(names(matrix_result$cutoffs), c("A", "B"))
    expect_identical(
        rownames(matrix_result$data),
        rownames(x)[matrix_result$feature_status$retained]
    )
    expect_identical(colnames(matrix_result$data), colnames(x))
})

test_that("automatic matrix and experiment paths preserve ordered representation", {
    fixture <- automatic_cutoff_matrix_fixture()
    feature_order <- c(
        seq(2L, nrow(fixture$x), by = 2L),
        seq(1L, nrow(fixture$x), by = 2L)
    )
    sample_order <- c(8L, 1L, 5L, 4L, 7L, 2L, 6L, 3L)
    x <- fixture$x[feature_order, sample_order, drop = FALSE]
    group <- stats::setNames(
        fixture$group[sample_order],
        colnames(x)
    )
    experiment <- rich_summarized_experiment(x, unname(group))
    original_experiment <- experiment

    matrix_result <- classify_missingness(x, group, seed = 29L)
    experiment_result <- classify_missingness(
        experiment,
        group_col = "condition",
        assay = "intensity",
        seed = 29L
    )

    expect_identical(experiment, original_experiment)
    expect_core_result_equivalent(matrix_result, experiment_result)
    expect_experiment_representation(
        experiment_result,
        matrix_result,
        original_experiment
    )
    expect_equal(matrix_result$cutoffs, fixture$expected, tolerance = 1e-12)
    expect_identical(
        vapply(
            matrix_result$cutoff_diagnostics,
            `[[`,
            character(1L),
            "source"
        ),
        c(A = "automatic", B = "automatic")
    )
    expect_identical(
        rownames(matrix_result$data),
        rownames(x)[matrix_result$feature_status$retained]
    )
    expect_identical(colnames(matrix_result$data), colnames(x))
})

test_that("an empty experiment result retains samples and metadata", {
    x <- matrix(
        c(8, NA, 9, NA),
        nrow = 1L,
        dimnames = list("all_mnar", paste0("s", seq_len(4L)))
    )
    group <- c("A", "A", "B", "B")
    experiment <- rich_summarized_experiment(x, group)
    original <- experiment
    matrix_result <- classify_missingness(
        x,
        group,
        cutoffs = c(A = 10, B = 10)
    )
    experiment_result <- classify_missingness(
        experiment,
        group_col = "condition",
        assay = "intensity",
        cutoffs = c(A = 10, B = 10)
    )

    expect_identical(experiment, original)
    expect_core_result_equivalent(matrix_result, experiment_result)
    expect_experiment_representation(
        experiment_result,
        matrix_result,
        original
    )
    expect_identical(dim(experiment_result$data), c(0L, 4L))
    expect_identical(
        vapply(
            SummarizedExperiment::assays(experiment_result$data),
            nrow,
            integer(1L)
        ),
        c(reference = 0L, intensity = 0L)
    )
})
