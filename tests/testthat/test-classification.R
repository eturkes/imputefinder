state_boundary_fixture <- function() {
    x <- rbind(
        complete_low = c(1, 1, 1, 1, 20, 20, 20),
        sparse_mnar = c(9, NA, NA, NA, 9, NA, NA),
        equal_majority = c(10, 10, NA, 10, 10, NA, 10),
        above_majority = c(11, 11, NA, 11, 11, NA, 11),
        exactly_half = c(11, 11, NA, NA, 11, NA, NA),
        arithmetic_mean = c(1, 1, 20, NA, 20, 20, 20)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))

    list(
        x = x,
        group = c(rep("A", 4L), rep("B", 3L)),
        cutoffs = c(A = 10, B = 10)
    )
}

classified_boundaries <- function() {
    fixture <- state_boundary_fixture()
    prepared <- prepare_matrix_input(fixture$x, fixture$group)
    statistics <- imputefinder:::.feature_condition_statistics(
        prepared$data,
        prepared$groups_by_sample,
        imputefinder:::.empty_seed_log()
    )
    resolved <- imputefinder:::.resolve_manual_cutoffs(
        statistics,
        fixture$cutoffs
    )

    imputefinder:::.assign_condition_states(statistics, resolved$cutoffs)
}

classified_state <- function(classifications, feature, condition) {
    classifications$state[
        classifications$feature == feature &
            classifications$condition == condition
    ]
}

test_that("four-state boundaries follow the normative decision order", {
    classified <- classified_boundaries()

    expect_identical(classified_state(classified, "complete_low", "A"), "complete")
    expect_identical(classified_state(classified, "sparse_mnar", "A"), "MNAR")
    expect_identical(classified_state(classified, "equal_majority", "A"), "MAR")
    expect_identical(classified_state(classified, "above_majority", "A"), "MAR")
    expect_identical(classified_state(classified, "exactly_half", "A"), "insufficient")
})

test_that("strict majority is derived independently for unequal group sizes", {
    classified <- classified_boundaries()

    expect_identical(classified_state(classified, "equal_majority", "B"), "MAR")
    expect_identical(classified_state(classified, "exactly_half", "B"), "insufficient")
    expect_identical(
        classified$sample_count[classified$condition == "A"] |> unique(),
        4L
    )
    expect_identical(
        classified$sample_count[classified$condition == "B"] |> unique(),
        3L
    )
})

test_that("one- and two-sample conditions follow the explicit state rules", {
    x <- rbind(
        single_rescue = c(NA, 15, 15),
        low_half = c(5, 8, NA),
        high_half = c(5, 12, NA),
        support = c(5, 5, 5)
    )
    colnames(x) <- c("a1", "b1", "b2")

    result <- classify_missingness(
        x,
        c("A", "B", "B"),
        cutoffs = c(B = 10)
    )

    a_states <- result$classifications$state[
        result$classifications$condition == "A"
    ]
    b_states <- stats::setNames(
        result$classifications$state[
            result$classifications$condition == "B"
        ],
        result$classifications$feature[
            result$classifications$condition == "B"
        ]
    )
    expect_true(all(a_states == "complete"))
    expect_identical(result$cutoffs[["A"]], NA_real_)
    expect_identical(result$cutoff_diagnostics$A$source, "not_needed")
    expect_identical(b_states[["low_half"]], "MNAR")
    expect_identical(b_states[["high_half"]], "insufficient")
    expect_identical(result$seed_log$feature, "single_rescue")
    expect_identical(result$seed_log$condition, "A")
})

test_that("the arithmetic mean drives intensity classification", {
    classified <- classified_boundaries()
    row <- classified[
        classified$feature == "arithmetic_mean" &
            classified$condition == "A",
        ,
        drop = FALSE
    ]

    expect_equal(row$mean_intensity, mean(c(1, 1, 20)))
    expect_gt(row$mean_intensity, 5)
    expect_lt(stats::median(c(1, 1, 20)), 5)

    recut <- imputefinder:::.assign_condition_states(
        classified[, c(
            "feature",
            "condition",
            "sample_count",
            "observed_count",
            "missing_count",
            "missing_fraction",
            "mean_intensity",
            "seeded"
        )],
        c(A = 5, B = 10)
    )
    expect_identical(classified_state(recut, "arithmetic_mean", "A"), "MAR")
})

test_that("state assignment produces the typed classification skeleton", {
    classified <- classified_boundaries()

    expect_identical(names(classified), names(imputefinder:::.empty_classifications()))
    expect_identical(
        vapply(classified, typeof, character(1L)),
        vapply(imputefinder:::.empty_classifications(), typeof, character(1L))
    )
    expect_true(all(is.na(classified$retained)))
    expect_true(all(is.na(classified$drop_reason)))
    expect_setequal(unique(classified$state), imputefinder:::.missingness_states())
})

test_that("incomplete conditions require a resolved cutoff before states", {
    fixture <- state_boundary_fixture()
    prepared <- prepare_matrix_input(fixture$x, fixture$group)
    statistics <- imputefinder:::.feature_condition_statistics(
        prepared$data,
        prepared$groups_by_sample,
        imputefinder:::.empty_seed_log()
    )

    expect_error(
        imputefinder:::.assign_condition_states(
            statistics,
            c(A = 10, B = NA_real_)
        ),
        "Condition `B` has incomplete features but no resolved cutoff.",
        fixed = TRUE
    )
})
