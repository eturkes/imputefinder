test_that("manual cutoff recovers the frozen causal on/off case", {
    simulation <- scientific_routine_fixture("manual_on_off")
    result <- classify_missingness(
        simulation$data,
        simulation$groups,
        cutoffs = c(A = 12, B = 12),
        seed = 1L
    )
    score <- score_scientific_routine(result, simulation)

    expect_gte(score[["mnar_f1"]], 0.60)
    expect_gte(score[["mar_f1"]], 0.60)
    expect_gte(score[["state_macro_f1"]], 0.66)
    expect_gte(score[["retention_f1"]], 0.90)
    expect_identical(score[["on_off_retention_recall"]], 1)
    expect_identical(score[["on_off_seed_recall"]], 1)
})

test_that("automatic cutoff passes the frozen cliff and permutation gates", {
    simulation <- scientific_routine_fixture("automatic_cliff")
    manual <- classify_missingness(
        simulation$data,
        simulation$groups,
        cutoffs = c(A = 12, B = 12),
        seed = 1L
    )
    automatic <- classify_missingness(
        simulation$data,
        simulation$groups,
        seed = 1L
    )
    manual_score <- score_scientific_routine(manual, simulation)
    automatic_score <- score_scientific_routine(automatic, simulation)

    expect_true(all(abs(automatic$cutoffs - 12) <= 1))
    expect_gte(automatic_score[["mnar_f1"]], 0.72)
    expect_gte(automatic_score[["mar_f1"]], 0.72)
    expect_gte(automatic_score[["state_macro_f1"]], 0.76)
    expect_gte(automatic_score[["retention_f1"]], 0.94)
    expect_gte(
        automatic_score[["mnar_f1"]] - manual_score[["mnar_f1"]],
        -0.05
    )
    expect_gte(
        automatic_score[["retention_f1"]] -
            manual_score[["retention_f1"]],
        -0.03
    )
    expect_identical(automatic_score[["on_off_retention_recall"]], 1)
    expect_identical(automatic_score[["on_off_seed_recall"]], 1)

    permutation <- permute_scientific_routine(simulation)
    permuted <- classify_missingness(
        simulation$data[
            permutation$features,
            permutation$samples,
            drop = FALSE
        ],
        stats::setNames(
            factor(
                unname(simulation$groups[permutation$samples]),
                levels = c("B", "A")
            ),
            colnames(simulation$data)[permutation$samples]
        ),
        seed = 1L
    )

    expect_identical(
        canonical_scientific_routine(permuted),
        canonical_scientific_routine(automatic)
    )
})
