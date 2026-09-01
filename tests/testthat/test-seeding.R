seed_missing_conditions <- function(x, group, seed = 1L) {
    prepared <- prepare_matrix_input(x, group)
    imputefinder:::.seed_missing_conditions(
        data = prepared$data,
        groups_by_sample = prepared$groups_by_sample,
        seed = seed
    )
}

rescue_fixture <- function() {
    x <- rbind(
        rescue_a_1 = c(NA, NA, 10, 11),
        rescue_a_2 = c(NA, NA, 12, 13),
        rescue_b = c(5, 6, NA, NA),
        partial = c(2, NA, 7, NA),
        globally_absent = rep(NA, 4)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))

    list(x = x, group = rep(c("A", "B"), each = 2))
}

test_that("globally absent features are audited, excluded, and never seeded", {
    fixture <- normative_fixture()

    rescued <- seed_missing_conditions(fixture$x, fixture$group)
    absent_status <- rescued$feature_status[
        rescued$feature_status$feature == "globally_absent",
        ,
        drop = FALSE
    ]

    expect_identical(
        rownames(rescued$data),
        setdiff(rownames(fixture$x), "globally_absent")
    )
    expect_false(absent_status$retained)
    expect_identical(absent_status$drop_reason, "all_missing")
    expect_false(any(rescued$seed_log$feature == "globally_absent"))
    expect_true(all(is.na(
        rescued$feature_status$retained[
            rescued$feature_status$feature != "globally_absent"
        ]
    )))
})

test_that("fully missing condition blocks receive one pre-seed minimum", {
    fixture <- rescue_fixture()

    rescued <- seed_missing_conditions(fixture$x, fixture$group, seed = 2L)

    expect_identical(rescued$condition_minima, c(A = 2, B = 7))
    expect_identical(
        rescued$seed_log[, c("feature", "condition")],
        data.frame(
            feature = c("rescue_a_1", "rescue_a_2", "rescue_b"),
            condition = c("A", "A", "B"),
            stringsAsFactors = FALSE
        )
    )

    for (row in seq_len(nrow(rescued$seed_log))) {
        entry <- rescued$seed_log[row, , drop = FALSE]
        condition_samples <- colnames(fixture$x)[
            fixture$group == entry$condition
        ]

        expect_true(is.na(fixture$x[entry$feature, entry$sample]))
        expect_true(is.na(entry$old_value))
        expect_identical(
            unname(rescued$data[entry$feature, entry$sample]),
            entry$inserted_value
        )
        expect_identical(
            entry$inserted_value,
            unname(rescued$condition_minima[entry$condition])
        )
        expect_equal(
            sum(!is.na(rescued$data[entry$feature, condition_samples])),
            1L
        )
    }

    a_log <- rescued$seed_log[rescued$seed_log$condition == "A", ]
    expect_identical(a_log$sample[[1L]], a_log$sample[[2L]])
    expect_true(all(rescued$seed_log$seed == 2L))
})

test_that("rescue preserves every originally observed value", {
    fixture <- rescue_fixture()
    original <- fixture$x

    rescued <- seed_missing_conditions(fixture$x, fixture$group)
    surviving_input <- fixture$x[rownames(rescued$data), , drop = FALSE]
    originally_observed <- !is.na(surviving_input)

    expect_identical(
        rescued$data[originally_observed],
        surviving_input[originally_observed]
    )
    expect_identical(colnames(rescued$data), colnames(fixture$x))
    expect_identical(fixture$x, original)
})

test_that("rescue is reproducible for a fixed seed", {
    fixture <- rescue_fixture()

    first <- seed_missing_conditions(fixture$x, fixture$group, seed = 8L)
    second <- seed_missing_conditions(fixture$x, fixture$group, seed = 8L)

    expect_identical(first, second)
})

test_that("rescue preserves an existing caller RNG state", {
    fixture <- rescue_fixture()

    with_preserved_random_state({
        set.seed(90210L)
        caller_seed <- get(
            ".Random.seed",
            envir = globalenv(),
            inherits = FALSE
        )

        rescued <- seed_missing_conditions(fixture$x, fixture$group)

        expect_gt(nrow(rescued$seed_log), 0L)
        expect_identical(
            get(
                ".Random.seed",
                envir = globalenv(),
                inherits = FALSE
            ),
            caller_seed
        )
    })
})

test_that("rescue does not create a caller RNG state", {
    fixture <- rescue_fixture()

    with_preserved_random_state({
        RNGkind("Knuth-TAOCP-2002", "Box-Muller", "Rejection")
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
        }
        caller_kind <- RNGkind()
        expect_false(exists(
            ".Random.seed",
            envir = globalenv(),
            inherits = FALSE
        ))

        rescued <- seed_missing_conditions(fixture$x, fixture$group)

        expect_gt(nrow(rescued$seed_log), 0L)
        expect_false(exists(
            ".Random.seed",
            envir = globalenv(),
            inherits = FALSE
        ))
        expect_identical(RNGkind(), caller_kind)
    })
})

test_that("rescue is independent of the caller RNG configuration", {
    x <- matrix(seq_len(600L), nrow = 60L)
    x[, seq_len(5L)] <- NA_real_
    x[[1L, 1L]] <- 1
    rownames(x) <- paste0("feature_", seq_len(nrow(x)))
    colnames(x) <- paste0("sample_", seq_len(ncol(x)))
    group <- rep(c("A", "B"), each = 5L)

    with_preserved_random_state({
        RNGkind("Mersenne-Twister", "Inversion", "Rejection")
        baseline <- seed_missing_conditions(x, group, seed = 23L)

        RNGkind("Knuth-TAOCP-2002", "Box-Muller", "Rejection")
        set.seed(90210L)
        caller_kind <- RNGkind()
        caller_seed <- get(
            ".Random.seed",
            envir = globalenv(),
            inherits = FALSE
        )
        configured <- seed_missing_conditions(x, group, seed = 23L)

        expect_identical(configured, baseline)
        expect_identical(RNGkind(), caller_kind)
        expect_identical(
            get(
                ".Random.seed",
                envir = globalenv(),
                inherits = FALSE
            ),
            caller_seed
        )
    })
})

test_that("named rescue assignments are invariant to input order", {
    fixture <- rescue_fixture()
    baseline <- seed_missing_conditions(fixture$x, fixture$group, seed = 19L)

    row_order <- c(5L, 3L, 1L, 4L, 2L)
    column_order <- c(4L, 3L, 2L, 1L)
    permuted_x <- fixture$x[row_order, column_order, drop = FALSE]
    permuted_group <- factor(
        fixture$group[column_order],
        levels = c("B", "A")
    )
    permuted <- seed_missing_conditions(
        permuted_x,
        permuted_group,
        seed = 19L
    )

    expect_identical(permuted$condition_minima, baseline$condition_minima)
    expect_identical(permuted$seed_log, baseline$seed_log)
    expect_identical(
        permuted$data[rownames(baseline$data), colnames(baseline$data)],
        baseline$data
    )
    permuted_status <- permuted$feature_status[
        match(
            baseline$feature_status$feature,
            permuted$feature_status$feature
        ),
        ,
        drop = FALSE
    ]
    rownames(permuted_status) <- NULL
    expect_identical(
        permuted_status,
        baseline$feature_status
    )

    expected_rows <- rownames(permuted_x)[
        rownames(permuted_x) != "globally_absent"
    ]
    expect_identical(rownames(permuted$data), expected_rows)
    expect_identical(colnames(permuted$data), colnames(permuted_x))
    expect_identical(permuted$feature_status$feature, rownames(permuted_x))
})

test_that("named group order does not alter rescue assignments", {
    fixture <- rescue_fixture()
    groups_by_sample <- stats::setNames(fixture$group, colnames(fixture$x))

    aligned <- seed_missing_conditions(
        fixture$x,
        groups_by_sample,
        seed = 31L
    )
    shuffled <- seed_missing_conditions(
        fixture$x,
        groups_by_sample[c("s3", "s1", "s4", "s2")],
        seed = 31L
    )

    expect_identical(shuffled, aligned)
})

test_that("every condition requires a finite pre-seed minimum", {
    x <- rbind(
        detected_b = c(NA, NA, 5, 6),
        globally_absent = rep(NA, 4)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    group <- rep(c("A", "B"), each = 2)

    expect_error(
        seed_missing_conditions(x, group),
        "Condition `A` has no finite intensity; a rescue minimum is undefined.",
        fixed = TRUE
    )
})

test_that("rescue seed must be one non-missing integer", {
    fixture <- rescue_fixture()

    for (seed in list(NA_integer_, NaN, Inf, 1.5, c(1L, 2L), "1")) {
        expect_error(
            seed_missing_conditions(fixture$x, fixture$group, seed = seed),
            "`seed` must be one non-missing integer.",
            fixed = TRUE
        )
    }
})
