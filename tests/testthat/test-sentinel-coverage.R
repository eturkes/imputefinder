sentinel_coverage_fixture <- function() {
    x <- rbind(
        f1 = c(10, 11, 12, 13, 14, 15),
        f2 = c(20, 21, NA, 22, NA, NA),
        f3 = c(NA, NA, NA, NA, 30, 31),
        f4 = rep(NA_real_, 6L),
        f5 = c(40, NA, NA, NA, 41, NA)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    metadata <- data.frame(
        condition = factor(
            c("A", "A", "A", "B", "B", "B"),
            levels = c("B", "A")
        ),
        batch = factor(
            c("x", "x", "z", "x", "y", "y"),
            levels = c("z", "y", "x")
        ),
        run_order = c(3L, 1L, 5L, 2L, 6L, 4L),
        subject = c("p1", "p2", "p3", "p1", "p2", "p4"),
        acquisition = c("DDA", "DDA", "DIA", "DDA", "DIA", "DIA"),
        row.names = colnames(x)
    )
    design <- missingness_design(
        metadata,
        condition = "condition",
        nuisance = c("run_order", "batch"),
        block = "subject",
        acquisition = "acquisition"
    )

    list(x = x, metadata = metadata, design = design)
}

sentinel_coverage_analysis <- function(fixture = sentinel_coverage_fixture()) {
    analyze_missingness(
        fixture$x,
        fixture$design,
        cutoffs = c(A = 0, B = 0),
        modules = "sentinel"
    )
}

test_that("static sentinel coverage has one exact schema", {
    analysis <- sentinel_coverage_analysis()
    coverage <- analysis$sentinel$coverage

    expect_identical(
        names(analysis$sentinel),
        c("pre_rescue", "coverage")
    )
    expect_identical(
        names(coverage),
        c(
            "schema", "sample", "condition", "role_level",
            "condition_role", "feature_overlap"
        )
    )
    expect_identical(coverage$schema, "sentinel_static_coverage_v1")
    expect_identical(
        names(coverage$sample),
        c(
            "sample", "condition", "input_feature_count",
            "globally_observable_feature_count", "detected_feature_count",
            "detection_fraction", "mean_intensity", "minimum_intensity",
            "median_intensity", "maximum_intensity"
        )
    )
    expect_identical(
        names(coverage$condition),
        c(
            "condition", "sample_count", "input_feature_count",
            "globally_observable_feature_count", "detected_feature_count",
            "complete_feature_count", "observed_cell_count",
            "eligible_cell_count", "detection_fraction"
        )
    )
    expect_identical(
        names(coverage$role_level),
        c(
            "role", "column", "encoding", "level", "numeric_value",
            "sample_count", "condition_count", "singleton"
        )
    )
    expect_identical(
        names(coverage$condition_role),
        c(
            "role", "column", "encoding", "level", "numeric_value",
            "condition", "sample_count",
            "globally_observable_feature_count", "detected_feature_count",
            "observed_cell_count", "eligible_cell_count",
            "detection_fraction", "empty"
        )
    )
    expect_identical(
        names(coverage$feature_overlap),
        c(
            "condition_left", "condition_right",
            "globally_observable_feature_count", "left_detected_count",
            "right_detected_count", "shared_detected_count",
            "union_detected_count", "left_only_count", "right_only_count",
            "neither_count", "jaccard"
        )
    )
    expect_invisible(imputefinder:::.validate_imputefinder_analysis(analysis))
})

test_that("sample and condition summaries use pre-rescue eligible support", {
    coverage <- sentinel_coverage_analysis()$sentinel$coverage

    expect_identical(coverage$sample$sample, paste0("s", seq_len(6L)))
    expect_identical(coverage$sample$condition, rep(c("A", "B"), each = 3L))
    expect_identical(coverage$sample$input_feature_count, rep(5L, 6L))
    expect_identical(
        coverage$sample$globally_observable_feature_count,
        rep(4L, 6L)
    )
    expect_identical(
        coverage$sample$detected_feature_count,
        c(3L, 2L, 1L, 2L, 3L, 2L)
    )
    expect_identical(
        coverage$sample$detection_fraction,
        c(3, 2, 1, 2, 3, 2) / 4
    )
    expect_equal(
        coverage$sample$mean_intensity,
        c(70 / 3, 16, 12, 17.5, 85 / 3, 23)
    )
    expect_identical(
        coverage$sample$minimum_intensity,
        c(10, 11, 12, 13, 14, 15)
    )
    expect_identical(
        coverage$sample$median_intensity,
        c(20, 16, 12, 17.5, 30, 23)
    )
    expect_identical(
        coverage$sample$maximum_intensity,
        c(40, 21, 12, 22, 41, 31)
    )

    expect_identical(coverage$condition$condition, c("A", "B"))
    expect_identical(coverage$condition$sample_count, c(3L, 3L))
    expect_identical(
        coverage$condition$detected_feature_count,
        c(3L, 4L)
    )
    expect_identical(
        coverage$condition$complete_feature_count,
        c(1L, 1L)
    )
    expect_identical(coverage$condition$observed_cell_count, c(6L, 7L))
    expect_identical(coverage$condition$eligible_cell_count, c(12L, 12L))
    expect_identical(coverage$condition$detection_fraction, c(1 / 2, 7 / 12))
})

test_that("declared role grids retain singleton levels and empty cells", {
    coverage <- sentinel_coverage_analysis()$sentinel$coverage
    levels <- coverage$role_level
    cells <- coverage$condition_role

    expect_identical(
        unique(levels$role),
        c("condition", "nuisance", "block", "acquisition")
    )
    expect_identical(
        unique(levels$column[levels$role == "nuisance"]),
        c("batch", "run_order")
    )
    batch_z <- levels$column == "batch" & levels$level == "z"
    expect_identical(levels$sample_count[batch_z], 1L)
    expect_identical(levels$condition_count[batch_z], 1L)
    expect_true(levels$singleton[batch_z])
    expect_true(all(levels$singleton[levels$column == "run_order"]))
    expect_identical(
        levels$numeric_value[levels$column == "run_order"],
        as.double(seq_len(6L))
    )
    expect_identical(
        levels$level[levels$column == "run_order"],
        as.character(seq_len(6L))
    )

    batch <- cells[cells$column == "batch", , drop = FALSE]
    expect_identical(
        batch[, c("level", "condition")],
        data.frame(
            level = rep(c("x", "y", "z"), each = 2L),
            condition = rep(c("A", "B"), 3L),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(batch$sample_count, c(2L, 1L, 0L, 2L, 1L, 0L))
    expect_identical(batch$empty, batch$sample_count == 0L)
    expect_identical(
        batch$observed_cell_count,
        c(5L, 2L, 0L, 5L, 1L, 0L)
    )
    expect_identical(
        batch$eligible_cell_count,
        c(8L, 4L, 0L, 8L, 4L, 0L)
    )
    expect_identical(
        is.na(batch$detection_fraction),
        batch$empty
    )

    subject <- cells[cells$column == "subject", , drop = FALSE]
    expect_identical(nrow(subject), 8L)
    expect_true(any(subject$condition == "B" & subject$level == "p3" &
        subject$empty))
    expect_true(any(subject$condition == "A" & subject$level == "p4" &
        subject$empty))
})

test_that("feature overlap is pairwise and excludes global absences", {
    overlap <- sentinel_coverage_analysis()$sentinel$coverage$feature_overlap

    expect_identical(nrow(overlap), 1L)
    expect_identical(overlap$condition_left, "A")
    expect_identical(overlap$condition_right, "B")
    expect_identical(overlap$globally_observable_feature_count, 4L)
    expect_identical(overlap$left_detected_count, 3L)
    expect_identical(overlap$right_detected_count, 4L)
    expect_identical(overlap$shared_detected_count, 3L)
    expect_identical(overlap$union_detected_count, 4L)
    expect_identical(overlap$left_only_count, 0L)
    expect_identical(overlap$right_only_count, 1L)
    expect_identical(overlap$neither_count, 0L)
    expect_identical(overlap$jaccard, 0.75)
})

test_that("coverage is invariant to named order and label representation", {
    fixture <- sentinel_coverage_fixture()
    baseline <- sentinel_coverage_analysis(fixture)$sentinel$coverage
    sample_order <- c("s6", "s2", "s4", "s1", "s5", "s3")
    metadata <- fixture$metadata[sample_order, rev(names(fixture$metadata)),
        drop = FALSE
    ]
    metadata$condition <- as.character(metadata$condition)
    metadata$batch <- as.character(metadata$batch)
    metadata$run_order <- as.double(metadata$run_order)
    reordered <- list(
        x = fixture$x[rev(rownames(fixture$x)), sample_order, drop = FALSE],
        design = missingness_design(
            metadata,
            condition = "condition",
            nuisance = c("batch", "run_order"),
            block = "subject",
            acquisition = "acquisition"
        )
    )

    expect_identical(
        sentinel_coverage_analysis(reordered)$sentinel$coverage,
        baseline
    )

    cancellation <- matrix(
        c(1e16, 1, -1e16, -1e16, 2, 1e16),
        nrow = 3L,
        dimnames = list(c("z", "a", "m"), c("left", "right"))
    )
    cancellation_design <- missingness_design(
        data.frame(
            condition = c("A", "B"),
            row.names = colnames(cancellation)
        ),
        condition = "condition"
    )
    cancellation_call <- function(data) {
        analyze_missingness(
            data,
            cancellation_design,
            cutoffs = c(A = 0, B = 0),
            modules = "sentinel"
        )
    }
    original <- cancellation_call(cancellation)
    reversed <- cancellation_call(cancellation[3:1, , drop = FALSE])
    expect_identical(
        original$sentinel$pre_rescue$sample,
        reversed$sentinel$pre_rescue$sample
    )
    expect_identical(
        original$sentinel$coverage,
        reversed$sentinel$coverage
    )
})

test_that("zero-support samples and absent optional roles remain explicit", {
    x <- matrix(
        c(NA, NA, 1, NA, 2, 3),
        nrow = 2L,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )
    design <- missingness_design(
        data.frame(
            condition = c("A", "A", "B"),
            row.names = colnames(x)
        ),
        condition = "condition"
    )
    analysis <- analyze_missingness(
        x,
        design,
        cutoffs = c(A = 0, B = 0),
        modules = "sentinel"
    )
    coverage <- analysis$sentinel$coverage

    expect_identical(coverage$sample$detected_feature_count, c(0L, 1L, 2L))
    expect_identical(coverage$sample$detection_fraction, c(0, 1 / 2, 1))
    expect_true(all(is.na(coverage$sample[1L, c(
        "mean_intensity", "minimum_intensity", "median_intensity",
        "maximum_intensity"
    )])))
    expect_identical(nrow(coverage$condition_role), 0L)
    expect_identical(unique(coverage$role_level$role), "condition")
    expect_invisible(imputefinder:::.validate_imputefinder_analysis(analysis))
})

test_that("globally absent input yields identified zero denominators", {
    x <- matrix(
        NA_real_,
        nrow = 2L,
        ncol = 2L,
        dimnames = list(c("f2", "f1"), c("s2", "s1"))
    )
    design <- missingness_design(
        data.frame(
            condition = c("B", "A"),
            row.names = colnames(x)
        ),
        condition = "condition"
    )
    analysis <- analyze_missingness(
        x,
        design,
        cutoffs = c(A = 0, B = 0),
        modules = "sentinel"
    )
    coverage <- analysis$sentinel$coverage

    expect_s3_class(analysis$classic, "imputefinder_classic_failure")
    expect_identical(
        coverage$sample$globally_observable_feature_count,
        c(0L, 0L)
    )
    expect_identical(coverage$sample$detected_feature_count, c(0L, 0L))
    expect_true(all(is.na(coverage$sample$detection_fraction)))
    expect_true(all(is.na(coverage$condition$detection_fraction)))
    expect_identical(
        coverage$feature_overlap$union_detected_count,
        0L
    )
    expect_identical(coverage$feature_overlap$jaccard, NA_real_)
    expect_invisible(imputefinder:::.validate_imputefinder_analysis(analysis))
})

test_that("coverage lifecycle and internal cross-counts reject corruption", {
    analysis <- sentinel_coverage_analysis()

    wrong_schema <- analysis
    wrong_schema$sentinel$coverage$schema <- "sentinel_static_coverage_v2"
    error <- tryCatch(
        imputefinder:::.validate_imputefinder_analysis(wrong_schema),
        imputefinder_sentinel_coverage_lifecycle_error = identity
    )
    expect_s3_class(
        error,
        "imputefinder_sentinel_coverage_lifecycle_error"
    )
    expect_identical(error$expected_schema, "sentinel_static_coverage_v1")
    expect_identical(error$actual_schema, "sentinel_static_coverage_v2")

    corruptions <- list(
        function(x) {
            x$sentinel$coverage$sample$detected_feature_count[[1L]] <- 0L
            x
        },
        function(x) {
            x$sentinel$coverage$condition_role <-
                x$sentinel$coverage$condition_role[-1L, , drop = FALSE]
            x
        },
        function(x) {
            x$sentinel$coverage$feature_overlap$shared_detected_count[[1L]] <-
                2L
            x
        },
        function(x) {
            x$sentinel$coverage$sample$minimum_intensity[[1L]] <- 50
            x
        }
    )
    for (mutator in corruptions) {
        expect_error(
            imputefinder:::.validate_imputefinder_analysis(mutator(analysis)),
            class = "imputefinder_sentinel_coverage_schema_error"
        )
    }

    legacy <- analysis
    legacy$sentinel$coverage <- NULL
    expect_identical(names(legacy$sentinel), "pre_rescue")
    expect_invisible(imputefinder:::.validate_imputefinder_analysis(legacy))
})
