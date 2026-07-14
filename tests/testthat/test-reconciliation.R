reconcile_normative_states <- function() {
    fixture <- normative_fixture()
    post_seed <- post_seed_statistics(fixture$x, fixture$group)
    resolved <- imputefinder:::.resolve_manual_cutoffs(
        post_seed$statistics,
        fixture$cutoffs
    )
    classified <- imputefinder:::.assign_condition_states(
        post_seed$statistics,
        resolved$cutoffs
    )

    imputefinder:::.reconcile_condition_states(
        classified,
        post_seed$rescued$feature_status
    )
}

three_condition_fixture <- function() {
    x <- rbind(
        complete_all = rep(8, 12L),
        mixed = c(
            8, NA, NA, NA,
            12, 12, NA, 12,
            12, 12, 12, 12
        ),
        all_mnar = rep(c(8, NA, NA, NA), 3L),
        insufficient_multi = c(
            12, NA, NA, NA,
            8, NA, NA, NA,
            12, 12, NA, NA
        ),
        mar_all = rep(c(12, 12, NA, 12), 3L)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    group <- rep(c("A", "B", "C"), each = 4L)
    permutation <- c(9:12, 1:4, 5:8)

    list(
        x = x[, permutation, drop = FALSE],
        group = group[permutation],
        cutoffs = c(C = 10, A = 10, B = 10)
    )
}

test_that("reconciliation applies every drop rule in deterministic precedence", {
    reconciled <- reconcile_normative_states()

    expect_identical(
        reconciled$feature_status,
        data.frame(
            feature = c(
                "on_off",
                "mar_both",
                "sparse_mar",
                "sparse_mnar",
                "all_mnar",
                "globally_absent",
                "complete_low"
            ),
            retained = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE),
            drop_reason = c(
                NA,
                NA,
                "insufficient:A",
                NA,
                "MNAR_all_conditions",
                "all_missing",
                NA
            ),
            stringsAsFactors = FALSE
        )
    )

    sparse_rows <- reconciled$classifications$feature == "sparse_mar"
    mnar_rows <- reconciled$classifications$feature == "all_mnar"
    expect_true(all(!reconciled$classifications$retained[sparse_rows]))
    expect_true(all(
        reconciled$classifications$drop_reason[sparse_rows] ==
            "insufficient:A"
    ))
    expect_true(all(!reconciled$classifications$retained[mnar_rows]))
    expect_true(all(
        reconciled$classifications$drop_reason[mnar_rows] ==
            "MNAR_all_conditions"
    ))
    expect_false(any(
        reconciled$classifications$feature == "globally_absent"
    ))
})

test_that("condition groups contain retained features in original order", {
    reconciled <- reconcile_normative_states()
    groups <- imputefinder:::.retained_condition_groups(
        reconciled$classifications,
        reconciled$feature_status$feature
    )

    expect_identical(
        groups,
        list(
            A = list(
                MNAR = c("on_off", "sparse_mnar"),
                MAR = "mar_both",
                complete = "complete_low",
                MAR_or_complete = c("mar_both", "complete_low")
            ),
            B = list(
                MNAR = character(),
                MAR = "mar_both",
                complete = c("on_off", "sparse_mnar", "complete_low"),
                MAR_or_complete = c(
                    "on_off",
                    "mar_both",
                    "sparse_mnar",
                    "complete_low"
                )
            )
        )
    )

    retained <- reconciled$feature_status$feature[
        reconciled$feature_status$retained
    ]
    for (condition in names(groups)) {
        expect_length(
            intersect(
                groups[[condition]]$MNAR,
                groups[[condition]]$MAR_or_complete
            ),
            0L
        )
        expect_setequal(
            c(
                groups[[condition]]$MNAR,
                groups[[condition]]$MAR_or_complete
            ),
            retained
        )
    }
})

test_that("public reconciliation generalizes to three conditions", {
    fixture <- three_condition_fixture()

    result <- classify_missingness(
        fixture$x,
        fixture$group,
        cutoffs = fixture$cutoffs
    )

    expect_identical(names(result$groups), c("A", "B", "C"))
    expect_identical(
        result$feature_status$drop_reason[
            result$feature_status$feature == "insufficient_multi"
        ],
        "insufficient:A,C"
    )
    expect_identical(
        rownames(result$data),
        c("complete_all", "mixed", "mar_all")
    )
    expect_identical(result$groups$A$MNAR, "mixed")
    expect_identical(result$groups$B$MAR, c("mixed", "mar_all"))
    expect_identical(result$groups$C$complete, c("complete_all", "mixed"))
    expect_false("all_mnar" %in% rownames(result$data))
})

test_that("the public matrix result satisfies the manual M3 contract", {
    fixture <- normative_fixture()
    original <- fixture$x

    result <- classify_missingness(
        fixture$x,
        fixture$group,
        cutoffs = fixture$cutoffs
    )

    expect_s3_class(result, "imputefinder_result")
    expect_identical(fixture$x, original)
    expect_identical(
        rownames(result$data),
        c("on_off", "mar_both", "sparse_mnar", "complete_low")
    )
    expect_identical(colnames(result$data), colnames(fixture$x))
    retained_original <- fixture$x[rownames(result$data), , drop = FALSE]
    originally_observed <- !is.na(retained_original)
    expect_identical(
        result$data[originally_observed],
        retained_original[originally_observed]
    )
    expect_identical(result$cutoffs, fixture$cutoffs)
    expect_identical(
        vapply(result$cutoff_diagnostics, `[[`, character(1L), "source"),
        c(A = "manual", B = "manual")
    )
    expect_identical(
        result$groups_by_sample,
        setNames(fixture$group, colnames(fixture$x))
    )
    expect_true(is.call(result$call))
})

test_that("all-complete conditions run without automatic cutoffs", {
    x <- rbind(
        low = c(1, 1, 2, 2),
        high = c(10, 10, 20, 20)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))

    result <- classify_missingness(x, c("B", "B", "A", "A"))

    expect_identical(result$cutoffs, c(A = NA_real_, B = NA_real_))
    expect_true(all(result$classifications$state == "complete"))
    expect_true(all(result$feature_status$retained))
    expect_identical(rownames(result$data), rownames(x))
    expect_identical(
        vapply(result$cutoff_diagnostics, `[[`, character(1L), "source"),
        c(A = "not_needed", B = "not_needed")
    )
})

test_that("an all-ineligible result preserves samples and audit rows", {
    x <- matrix(
        c(8, NA, 9, NA),
        nrow = 1L,
        dimnames = list("all_mnar", paste0("s", seq_len(4L)))
    )

    result <- classify_missingness(
        x,
        c("A", "A", "B", "B"),
        cutoffs = c(A = 10, B = 10)
    )

    expect_identical(dim(result$data), c(0L, 4L))
    expect_identical(colnames(result$data), colnames(x))
    expect_true(all(lengths(unlist(result$groups, recursive = FALSE)) == 0L))
    expect_identical(result$feature_status$retained, FALSE)
    expect_identical(
        result$feature_status$drop_reason,
        "MNAR_all_conditions"
    )
    expect_true(all(!result$classifications$retained))
})

test_that("public SummarizedExperiment and matrix classifications agree", {
    skip_if_not_installed("SummarizedExperiment")
    fixture <- normative_fixture()
    auxiliary <- fixture$x + 100
    experiment <- SummarizedExperiment::SummarizedExperiment(
        assays = list(intensity = fixture$x, auxiliary = auxiliary),
        rowData = data.frame(
            marker = seq_len(nrow(fixture$x)),
            row.names = rownames(fixture$x)
        ),
        colData = data.frame(
            condition = fixture$group,
            row.names = colnames(fixture$x)
        ),
        metadata = list(source = "M3 integration fixture")
    )
    original <- experiment

    matrix_result <- classify_normative_fixture()
    experiment_result <- classify_missingness(
        experiment,
        group_col = "condition",
        assay = "intensity",
        cutoffs = fixture$cutoffs
    )

    expect_s4_class(experiment_result$data, "SummarizedExperiment")
    expect_identical(experiment, original)
    expect_identical(
        SummarizedExperiment::assay(experiment_result$data, "intensity"),
        matrix_result$data
    )
    expect_identical(
        SummarizedExperiment::assay(experiment_result$data, "auxiliary"),
        auxiliary[rownames(matrix_result$data), , drop = FALSE]
    )
    expect_identical(
        experiment_result$classifications,
        matrix_result$classifications
    )
    expect_identical(experiment_result$groups, matrix_result$groups)
    expect_identical(experiment_result$profiles, matrix_result$profiles)
    expect_identical(
        experiment_result$cutoff_diagnostics,
        matrix_result$cutoff_diagnostics
    )
    expect_identical(
        experiment_result$feature_status,
        matrix_result$feature_status
    )
})
