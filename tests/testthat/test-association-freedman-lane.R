association_fl_detection_matrix <- function(counts, feature_count = max(counts)) {
    samples <- sprintf("sample_%02d", seq_along(counts))
    output <- matrix(
        NA_real_,
        nrow = feature_count,
        ncol = length(counts),
        dimnames = list(
            sprintf("feature_%02d", seq_len(feature_count)),
            samples
        )
    )
    for (index in seq_along(counts)) {
        output[seq_len(counts[[index]]), index] <- 1
    }
    output
}

association_fl_seed_fixture <- function(draw_ids = 1:4) {
    instances <- data.frame(
        candidate_id = "a_fraction_freedman_lane",
        acquisition = "all",
        hypothesis_id = paste0("a_", strrep("1", 64L)),
        stringsAsFactors = FALSE
    )
    imputefinder:::.association_seed_manifest(
        strrep("0", 64L),
        instances,
        as.integer(draw_ids)
    )
}

association_fl_block_swaps <- function(block_count) {
    n <- 2L * block_count
    lapply(seq_len(block_count), function(block) {
        map <- seq_len(n)
        rows <- 2L * block + c(-1L, 0L)
        map[rows] <- rev(rows)
        as.integer(map)
    })
}

association_fl_outcome <- function(artifact, label) {
    index <- which(artifact$hypotheses$label == label)
    artifact$outcomes[[artifact$hypotheses$hypothesis[[index]]]]
}

test_that("permutation map hashes and global seeds match frozen bytes", {
    expect_identical(
        imputefinder:::.association_permutation_map_sha256(as.integer(1:5)),
        "1407b0f356b05b4083391bfd1d4876d4e61f05435460f6b5a61b81b96da22b70"
    )
    expect_identical(
        imputefinder:::.association_permutation_map_sha256(as.integer(5:1)),
        "7537141aed7d7017cf3a200af0fe7bfc13e1d53f1fe0bdd8b6612ea916177e3a"
    )
    expect_identical(
        imputefinder:::.association_permutation_map_sha256(
            as.integer(c(2, 5, 1, 3, 4))
        ),
        "3c8d56eb2105171dd4f165b31736396045028c69124ff0dfdcc4f94980ccadff"
    )

    seeds <- association_fl_seed_fixture()
    expect_identical(seeds$protocol_id, rep(
        "m13_a_association_protocol_v3",
        4L
    ))
    expect_identical(seeds$draw_id, 1:4)
    expect_identical(
        seeds$seed,
        c(89457621L, 112932404L, 177825668L, 33585167L)
    )
    expect_identical(seeds$seed_nonce, rep(0L, 4L))

    reversed <- data.frame(
        candidate_id = rep("a_fraction_freedman_lane", 2L),
        acquisition = c("z", "a"),
        hypothesis_id = c(
            paste0("a_", strrep("2", 64L)),
            paste0("a_", strrep("1", 64L))
        ),
        stringsAsFactors = FALSE
    )
    canonical <- imputefinder:::.association_seed_manifest(
        strrep("f", 64L),
        reversed,
        as.integer(c(2, 1))
    )
    expect_identical(canonical$acquisition, c("a", "a", "z", "z"))
    expect_identical(canonical$draw_id, c(1L, 2L, 1L, 2L))
    expect_identical(anyDuplicated(canonical$seed), 0L)
})

test_that("Monte Carlo maps match frozen streams and restore RNG exactly", {
    seeds <- association_fl_seed_fixture()
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) .Random.seed else NULL
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (had_seed) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
        )) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)

    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(8491)
    before_kind <- RNGkind()
    before_seed <- .Random.seed
    generated <- imputefinder:::.association_monte_carlo_maps(
        seeds,
        "independent",
        9L,
        groups = list(as.integer(1:9))
    )
    expect_identical(RNGkind(), before_kind)
    expect_identical(.Random.seed, before_seed)
    expect_identical(generated$maps, list(
        `1` = c(9L, 3L, 4L, 2L, 8L, 1L, 6L, 5L, 7L),
        `2` = c(8L, 2L, 5L, 1L, 6L, 3L, 7L, 4L, 9L),
        `3` = c(2L, 1L, 8L, 5L, 6L, 9L, 3L, 4L, 7L),
        `4` = c(6L, 4L, 1L, 9L, 5L, 2L, 7L, 3L, 8L)
    ))
    expect_identical(generated$manifest$map_retry, rep(0L, 4L))
    expect_identical(
        generated$manifest$map_sha256,
        c(
            "18ee5a32a3c8345274c5c4200ac257e85192dd2849c4325dba38d47e0562e3cc",
            "04facf4394b59c3bbfa194a42c6d99c96c24dbbfad87d0fc66d043e4a3f0feac",
            "1e4e1989a5079de24b81812468270e57e2dd338549771d258119bd93c2ff5f06",
            "8da4c17b54cddac01b43160ec374ca2abe0571fd6e197a286e1acecf8268b19e"
        )
    )
    wrong_protocol <- seeds
    wrong_protocol$protocol_id <- "wrong_protocol"
    negative_nonce <- seeds
    negative_nonce$seed_nonce[[1L]] <- -1L
    double_draw <- seeds
    double_draw$draw_id <- as.double(double_draw$draw_id)
    for (malformed in list(wrong_protocol, negative_nonce, double_draw)) {
        expect_error(
            imputefinder:::.association_monte_carlo_maps(
                malformed,
                "independent",
                9L,
                groups = list(as.integer(1:9))
            ),
            class = "imputefinder_association_permutation_error"
        )
    }

    rm(".Random.seed", envir = .GlobalEnv)
    absent <- !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    blocked <- imputefinder:::.association_monte_carlo_maps(
        seeds,
        "blocked",
        34L,
        swaps = association_fl_block_swaps(17L)
    )
    expect_true(absent)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    expect_identical(blocked$manifest$map_retry, rep(0L, 4L))
    expect_identical(
        blocked$manifest$map_sha256,
        c(
            "0fa6c9f528db3305ddc8cbaf7dba85f1f19fa5934cb1c08309c6348a7259dfc8",
            "c918ba056b916fe8bdba6d2b1e654319eb3372bbc64cc5f6c3e04eb56073c5cb",
            "0129bc561c8a0577199f05f4cc5454727dc54bc3f4f0aacb90b13544c9ec6e8c",
            "5dc9da3391f7dfffd75daf1c6ade1aede89bf635dd135ceaa0da3a760f55f1bb"
        )
    )
})

test_that("constrained-null geometry and independent groups are canonical", {
    z <- cbind(
        intercept = rep(1, 6L),
        linear = c(-2, -1, 0, 0, 1, 2),
        curved = c(4, 1, 0, 0, 1, 4)
    )
    contrast <- c(0, 1, -0.25)
    response <- c(0.1, 0.3, 0.2, 0.7, 0.8, 0.9)
    restricted <- imputefinder:::.association_constrained_null(
        z,
        contrast,
        response
    )
    expect_true(restricted$ok)
    expect_equal(
        as.vector(crossprod(contrast, restricted$basis)),
        c(0, 0),
        tolerance = 1e-14
    )
    expect_equal(
        unname(crossprod(restricted$basis)),
        diag(2L),
        tolerance = 1e-14
    )
    expect_equal(
        restricted$fitted + restricted$residual,
        response,
        tolerance = 1e-14
    )
    expect_equal(
        as.vector(crossprod(restricted$z0, restricted$residual)),
        c(0, 0),
        tolerance = 1e-13
    )

    z0 <- matrix(c(rep(0, 4L), seq_len(8L)), ncol = 1L)
    plan <- imputefinder:::.association_independent_permutations(z0)
    expect_true(plan$ok)
    expect_identical(plan$mode, "exact")
    expect_identical(plan$allowable_transformations, 24.0)
    expect_length(plan$maps, 24L)
    expect_identical(plan$maps[[1L]], as.integer(seq_len(12L)))
    expect_identical(
        plan$maps[[2L]],
        c(1L, 2L, 4L, 3L, as.integer(5:12))
    )
    expect_identical(
        plan$maps[[24L]],
        c(4L, 3L, 2L, 1L, as.integer(5:12))
    )
    hashes <- vapply(
        plan$maps,
        imputefinder:::.association_permutation_map_sha256,
        character(1L)
    )
    expect_identical(anyDuplicated(hashes), 0L)
    two_groups <- imputefinder:::.association_independent_permutations(
        matrix(c(0, 0, 1, 1), ncol = 1L)
    )
    expect_identical(two_groups$maps, list(
        1:4,
        c(1L, 2L, 4L, 3L),
        c(2L, 1L, 3L, 4L),
        c(2L, 1L, 4L, 3L)
    ))

    tolerance <- sqrt(.Machine$double.eps)
    nontransitive <- matrix(c(0, 0.75, 1.5) * tolerance, ncol = 1L)
    rejected <- imputefinder:::.association_independent_permutations(
        nontransitive
    )
    expect_false(rejected$ok)
    expect_identical(rejected$code, "association_incompatible_permutation")

    rounded <- imputefinder:::.association_independent_permutations(
        matrix(numeric(), nrow = 28L, ncol = 0L),
        enumerate = FALSE
    )
    expect_identical(
        sprintf("%a", rounded$allowable_transformations),
        "0x1.ec92dd23d6967p+97"
    )

    overflow <- imputefinder:::.association_independent_permutations(
        matrix(numeric(), nrow = 171L, ncol = 0L),
        enumerate = FALSE
    )
    expect_true(overflow$ok)
    expect_identical(overflow$mode, "monte_carlo")
    expect_identical(
        overflow$allowable_transformations,
        .Machine$double.xmax
    )
    expect_identical(
        imputefinder:::.association_factorial_product_double(100000L),
        .Machine$double.xmax
    )
})

test_that("blocked maps swap complete bundles in ascending binary masks", {
    samples <- sprintf("sample_%02d", seq_len(8L))
    metadata <- data.frame(
        condition = rep(c("A", "A", "B", "B"), 2L),
        batch = rep(c("x", "y", "x", "y"), 2L),
        subject = rep(c("p1", "p2"), each = 4L),
        row.names = samples,
        stringsAsFactors = FALSE
    )
    data <- association_fl_detection_matrix(seq_len(8L) + 4L, 12L)
    preparation <- imputefinder:::.new_association_preparation(
        data,
        missingness_design(
            metadata,
            condition = "condition",
            nuisance = "batch",
            block = "subject"
        )
    )
    hypothesis <- preparation$hypotheses[
        preparation$hypotheses$label == "condition[B]",
        ,
        drop = FALSE
    ]
    stratum <- preparation$strata[[1L]]
    fit <- imputefinder:::.association_robust_algebra(stratum)
    restriction <- imputefinder:::.association_freedman_lane_restriction(
        fit,
        hypothesis,
        stratum
    )
    plan <- imputefinder:::.association_blocked_permutations(
        stratum,
        hypothesis,
        restriction$z0
    )
    expect_true(plan$ok)
    expect_identical(plan$mode, "low_resolution")
    expect_identical(plan$allowable_transformations, 4.0)
    expect_identical(plan$maps[[1L]], as.integer(1:8))
    expect_identical(plan$maps[[2L]], c(3L, 4L, 1L, 2L, 5L, 6L, 7L, 8L))
    expect_identical(plan$maps[[3L]], c(1L, 2L, 3L, 4L, 7L, 8L, 5L, 6L))
    expect_identical(plan$maps[[4L]], c(3L, 4L, 1L, 2L, 7L, 8L, 5L, 6L))
    expect_true(all(vapply(plan$maps, function(map) {
        max(abs(restriction$z0[map, , drop = FALSE] - restriction$z0)) <=
            restriction$tolerance
    }, logical(1L))))

    incomplete <- stratum
    incomplete$design$sample_data$batch[[8L]] <- "z"
    incomplete$core <- imputefinder:::.new_design_estimability(
        incomplete$design
    )
    incomplete_fit <- imputefinder:::.association_robust_algebra(incomplete)
    incomplete_restriction <-
        imputefinder:::.association_freedman_lane_restriction(
            incomplete_fit,
            hypothesis,
            incomplete
        )
    incomplete_plan <- imputefinder:::.association_blocked_permutations(
        incomplete,
        hypothesis,
        incomplete_restriction$z0
    )
    expect_true(incomplete_plan$ok)
    expect_identical(incomplete_plan$allowable_transformations, 2.0)

    reversed_labels <- metadata
    reversed_labels$subject <- rep(c("z_block", "a_block"), each = 4L)
    reversed_preparation <- imputefinder:::.new_association_preparation(
        data,
        missingness_design(
            reversed_labels,
            condition = "condition",
            nuisance = "batch",
            block = "subject"
        )
    )
    reversed_hypothesis <- reversed_preparation$hypotheses[
        reversed_preparation$hypotheses$label == "condition[B]",
        ,
        drop = FALSE
    ]
    reversed_stratum <- reversed_preparation$strata[[1L]]
    reversed_fit <- imputefinder:::.association_robust_algebra(
        reversed_stratum
    )
    reversed_restriction <-
        imputefinder:::.association_freedman_lane_restriction(
            reversed_fit,
            reversed_hypothesis,
            reversed_stratum
        )
    reversed_plan <- imputefinder:::.association_blocked_permutations(
        reversed_stratum,
        reversed_hypothesis,
        reversed_restriction$z0
    )
    expect_identical(
        reversed_plan$maps[[2L]],
        c(1L, 2L, 3L, 4L, 7L, 8L, 5L, 6L)
    )
})

test_that("blocked pivots leave untested condition levels fixed", {
    samples <- sprintf("sample_%02d", seq_len(6L))
    metadata <- data.frame(
        condition = rep(c("A", "B", "C"), 2L),
        subject = rep(c("p1", "p2"), each = 3L),
        row.names = samples,
        stringsAsFactors = FALSE
    )
    preparation <- imputefinder:::.new_association_preparation(
        association_fl_detection_matrix(seq_len(6L) + 4L, 10L),
        missingness_design(
            metadata,
            condition = "condition",
            block = "subject"
        )
    )
    hypothesis <- preparation$hypotheses[
        preparation$hypotheses$label == "condition[B]",
        ,
        drop = FALSE
    ]
    stratum <- preparation$strata[[1L]]
    fit <- imputefinder:::.association_robust_algebra(stratum)
    restriction <- imputefinder:::.association_freedman_lane_restriction(
        fit,
        hypothesis,
        stratum
    )
    plan <- imputefinder:::.association_blocked_permutations(
        stratum,
        hypothesis,
        restriction$z0
    )
    expect_true(plan$ok)
    expect_identical(plan$allowable_transformations, 4.0)
    expect_true(all(vapply(plan$maps, function(map) {
        identical(map[c(3L, 6L)], c(3L, 6L))
    }, logical(1L))))
})

test_that("cached HC3 and CR2 response refits equal full refits", {
    independent_x <- association_fl_detection_matrix(
        c(8L, 10L, 9L, 12L, 11L, 15L, 13L, 16L, 18L),
        18L
    )
    independent_metadata <- data.frame(
        condition = c(rep("A", 4L), rep("B", 5L)),
        row.names = colnames(independent_x),
        stringsAsFactors = FALSE
    )
    blocked_x <- association_fl_detection_matrix(
        c(8L, 12L, 10L, 13L, 9L, 15L, 12L, 16L,
          11L, 15L, 13L, 18L, 14L, 17L, 15L, 20L),
        20L
    )
    blocked_metadata <- data.frame(
        condition = rep(c("A", "B"), 8L),
        subject = rep(sprintf("p%02d", seq_len(8L)), each = 2L),
        row.names = colnames(blocked_x),
        stringsAsFactors = FALSE
    )
    preparations <- list(
        imputefinder:::.new_association_preparation(
            independent_x,
            missingness_design(
                independent_metadata,
                condition = "condition"
            )
        ),
        imputefinder:::.new_association_preparation(
            blocked_x,
            missingness_design(
                blocked_metadata,
                condition = "condition",
                block = "subject"
            )
        )
    )
    for (preparation in preparations) {
        stratum <- preparation$strata[[1L]]
        fit <- imputefinder:::.association_robust_algebra(stratum)
        response <- as.double(0.1 + 0.75 * rev(stratum$response))
        cached <- imputefinder:::.association_robust_refit(fit, response)
        rebuilt_stratum <- stratum
        rebuilt_stratum$response <- response
        rebuilt <- imputefinder:::.association_robust_algebra(
            rebuilt_stratum
        )
        expect_true(cached$ok)
        expect_true(rebuilt$ok)
        expect_equal(
            cached$coefficients,
            rebuilt$coefficients,
            tolerance = 1e-14
        )
        expect_equal(
            cached$covariance,
            rebuilt$covariance,
            tolerance = 1e-13
        )
        coefficient <- preparation$hypotheses$coefficient[[1L]]
        columns <- colnames(stratum$core$model$matrix)
        cached_result <- imputefinder:::.association_robust_contrast(
            cached,
            coefficient,
            columns
        )
        rebuilt_result <- imputefinder:::.association_robust_contrast(
            rebuilt,
            coefficient,
            columns
        )
        expect_true(cached_result$ok)
        expect_true(rebuilt_result$ok)
        expect_equal(
            unlist(cached_result[c(
                "effect", "standard_error", "statistic", "reference_df"
            )]),
            unlist(rebuilt_result[c(
                "effect", "standard_error", "statistic", "reference_df"
            )]),
            tolerance = 1e-13
        )
    }
})

test_that("orthogonal-complement CR2 reaches blocked Monte Carlo scope", {
    block_count <- 17L
    counts <- as.integer(unlist(lapply(seq_len(block_count), function(block) {
        c(5L + block, 7L + block + block %% 3L)
    })))
    x <- association_fl_detection_matrix(counts, max(counts))
    metadata <- data.frame(
        condition = rep(c("A", "B"), block_count),
        subject = rep(
            sprintf("p%02d", seq_len(block_count)),
            each = 2L
        ),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    preparation <- imputefinder:::.new_association_preparation(
        x,
        missingness_design(
            metadata,
            condition = "condition",
            block = "subject"
        )
    )
    stratum <- preparation$strata[[1L]]
    fit <- imputefinder:::.association_robust_algebra(stratum)
    expect_true(fit$ok)
    hypothesis <- preparation$hypotheses[1L, , drop = FALSE]
    result <- imputefinder:::.association_robust_contrast(
        fit,
        hypothesis$coefficient,
        colnames(stratum$core$model$matrix)
    )
    expect_true(result$ok)
    restriction <- imputefinder:::.association_freedman_lane_restriction(
        fit,
        hypothesis,
        stratum
    )
    plan <- imputefinder:::.association_blocked_permutations(
        stratum,
        hypothesis,
        restriction$z0,
        enumerate = FALSE
    )
    expect_true(plan$ok)
    expect_identical(plan$mode, "monte_carlo")
    expect_identical(plan$allowable_transformations, 131072.0)

    artifact <- imputefinder:::.run_association_freedman_lane(preparation)
    outcome <- artifact$outcomes[[1L]]
    expect_s3_class(outcome, "imputefinder_association")
    expect_identical(outcome$diagnostics$variance_method, "CR2")
    expect_identical(outcome$diagnostics$permutation_mode, "monte_carlo")
    expect_identical(outcome$diagnostics$exceedance_count, 0L)
    expect_identical(outcome$raw_p, 0.0001)
    expect_identical(nrow(artifact$diagnostics$seed_manifest), 9999L)
})

test_that("the exact identity map preserves its literal Wald tie", {
    expect_false(imputefinder:::.association_byte_identical_numeric(0, -0))

    counts <- c(21L, 9L, 31L, 3L, 23L, 43L, 49L, 8L, 50L, 6L)
    x <- association_fl_detection_matrix(counts, 50L)
    metadata <- data.frame(
        condition = c("A", "A", "B", "B", "A", "B", "B", "B", "A", "A"),
        batch = c("c", "a", "b", "c", "c", "a", "c", "b", "c", "b"),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    preparation <- imputefinder:::.new_association_preparation(
        x,
        missingness_design(
            metadata,
            condition = "condition",
            nuisance = "batch"
        )
    )
    index <- which(preparation$hypotheses$label == "condition[B]")
    record <- imputefinder:::.association_freedman_lane_records(
        preparation
    )[[index]]
    expect_identical(record$permutations$allowable_transformations, 1440.0)
    recomposed <- record$restriction$fitted + record$restriction$residual
    expect_equal(
        max(abs(recomposed - record$stratum$response)),
        5.551115123125783e-17,
        tolerance = 0
    )
    drifted_fit <- imputefinder:::.association_robust_refit(
        record$fit,
        as.double(recomposed)
    )
    drifted <- imputefinder:::.association_robust_contrast(
        drifted_fit,
        record$hypothesis$coefficient,
        colnames(record$stratum$core$model$matrix)
    )
    observed_wald <- (record$observed$effect /
        record$observed$standard_error)^2
    drifted_wald <- (drifted$effect / drifted$standard_error)^2
    expect_lt(drifted_wald, observed_wald)

    record$permutations$maps <- list(
        as.integer(seq_along(record$stratum$response))
    )
    record$permutations$allowable_transformations <- 1.0
    evaluated <- imputefinder:::.association_freedman_lane_evaluate(record)
    expect_true(evaluated$ok)
    expect_identical(evaluated$statistic, observed_wald)
    expect_identical(evaluated$exceedance_count, 1L)
})

test_that("nonidentity residual stabilizers preserve literal Wald ties", {
    counts <- c(7L, 26L, 3L, 44L, 45L, 7L, 14L, 11L, 46L, 50L)
    x <- association_fl_detection_matrix(counts, 50L)
    metadata <- data.frame(
        condition = c(rep("A", 5L), rep("B", 5L)),
        batch = c("a", "b", "c", "a", "b", "a", "c", "b", "c", "b"),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    preparation <- imputefinder:::.new_association_preparation(
        x,
        missingness_design(
            metadata,
            condition = "condition",
            nuisance = "batch"
        )
    )
    index <- which(preparation$hypotheses$label == "condition[B]")
    record <- imputefinder:::.association_freedman_lane_records(
        preparation
    )[[index]]
    expect_identical(record$permutations$allowable_transformations, 864.0)
    stabilizer <- c(6L, 2L, 3L, 4L, 5L, 1L, 7L, 8L, 9L, 10L)
    expect_identical(
        record$restriction$residual[stabilizer],
        record$restriction$residual
    )
    recomposed <- record$restriction$fitted +
        record$restriction$residual[stabilizer]
    expect_equal(
        max(abs(recomposed - record$stratum$response)),
        2.7755575615628914e-17,
        tolerance = 0
    )
    drifted_fit <- imputefinder:::.association_robust_refit(
        record$fit,
        as.double(recomposed)
    )
    drifted <- imputefinder:::.association_robust_contrast(
        drifted_fit,
        record$hypothesis$coefficient,
        colnames(record$stratum$core$model$matrix)
    )
    observed_wald <- (record$observed$effect /
        record$observed$standard_error)^2
    drifted_wald <- (drifted$effect / drifted$standard_error)^2
    expect_lt(drifted_wald, observed_wald)

    record$permutations$maps <- list(stabilizer)
    record$permutations$allowable_transformations <- 1.0
    evaluated <- imputefinder:::.association_freedman_lane_evaluate(record)
    expect_true(evaluated$ok)
    expect_identical(evaluated$statistic, observed_wald)
    expect_identical(evaluated$exceedance_count, 1L)
})

test_that("blocked exact Freedman-Lane recomputes CR2 Wald statistics", {
    counts <- c(
        8L, 12L, 10L, 13L, 9L, 15L, 12L, 16L,
        11L, 15L, 13L, 18L, 14L, 17L, 15L, 20L
    )
    x <- association_fl_detection_matrix(counts, 20L)
    metadata <- data.frame(
        condition = rep(c("A", "B"), 8L),
        subject = rep(sprintf("p%02d", seq_len(8L)), each = 2L),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    preparation <- imputefinder:::.new_association_preparation(
        x,
        missingness_design(
            metadata,
            condition = "condition",
            block = "subject"
        )
    )
    robust <- imputefinder:::.run_association_ols_hc3_cr2(preparation)
    artifact <- imputefinder:::.run_association_freedman_lane(preparation)
    outcome <- association_fl_outcome(artifact, "condition[B]")
    robust_outcome <- association_fl_outcome(robust, "condition[B]")

    expect_identical(artifact$candidate, "a_fraction_freedman_lane")
    expect_equal(outcome$effect, robust_outcome$effect, tolerance = 1e-14)
    expect_equal(
        outcome$standard_error,
        robust_outcome$standard_error,
        tolerance = 1e-14
    )
    expect_equal(
        outcome$statistic,
        robust_outcome$statistic^2,
        tolerance = 1e-13
    )
    expect_identical(outcome$diagnostics$variance_method, "CR2")
    expect_identical(outcome$diagnostics$permutation_mode, "exact")
    expect_identical(outcome$diagnostics$allowable_transformations, 256.0)
    expect_identical(outcome$diagnostics$evaluated_transformations, 256L)
    expect_identical(outcome$permutation_count, 256L)
    expect_identical(outcome$diagnostics$exceedance_count, 1L)
    expect_identical(outcome$raw_p, 0.00390625)
    expect_identical(
        artifact$diagnostics$seed_manifest,
        imputefinder:::.empty_association_seed_manifest()
    )
    expect_invisible(
        imputefinder:::.validate_association_candidate_artifact(
            artifact,
            preparation
        )
    )

    coherent_wrong_count <- artifact
    coherent_outcome <- coherent_wrong_count$outcomes[[1L]]
    coherent_outcome$raw_p <- 2 / 128
    coherent_outcome$adjusted_p <- coherent_outcome$raw_p
    coherent_outcome$flag <- TRUE
    coherent_outcome$permutation_count <- 128L
    coherent_outcome$diagnostics$allowable_transformations <- 128.0
    coherent_outcome$diagnostics$evaluated_transformations <- 128L
    coherent_wrong_count$outcomes[[1L]] <- coherent_outcome
    coherent_wrong_exceedance <- artifact
    wrong_exceedance <- coherent_wrong_exceedance$outcomes[[1L]]
    wrong_exceedance$raw_p <- 3 / 256
    wrong_exceedance$adjusted_p <- wrong_exceedance$raw_p
    wrong_exceedance$diagnostics$exceedance_count <- 3L
    coherent_wrong_exceedance$outcomes[[1L]] <- wrong_exceedance
    forged_effect <- artifact
    forged <- forged_effect$outcomes[[1L]]
    shift <- 0.01
    forged$effect <- forged$effect + shift
    forged$conf_low <- forged$conf_low + shift
    forged$conf_high <- forged$conf_high + shift
    forged$statistic <- (forged$effect / forged$standard_error)^2
    forged_effect$outcomes[[1L]] <- forged
    for (corrupted in list(
        coherent_wrong_count,
        coherent_wrong_exceedance,
        forged_effect
    )) {
        expect_error(
            imputefinder:::.validate_association_candidate_artifact(
                corrupted,
                preparation
            ),
            class = "imputefinder_association_artifact_error"
        )
    }
})

test_that("independent Monte Carlo Freedman-Lane binds all 9999 maps", {
    counts <- c(8L, 10L, 9L, 12L, 11L, 15L, 13L, 16L, 18L)
    x <- association_fl_detection_matrix(counts, 18L)
    metadata <- data.frame(
        condition = c(rep("A", 4L), rep("B", 5L)),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    preparation <- imputefinder:::.new_association_preparation(
        x,
        missingness_design(metadata, condition = "condition")
    )
    before_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    before_seed <- if (had_seed) .Random.seed else NULL
    artifact <- imputefinder:::.run_association_freedman_lane(preparation)
    expect_identical(RNGkind(), before_kind)
    expect_identical(
        exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
        had_seed
    )
    if (had_seed) {
        expect_identical(.Random.seed, before_seed)
    }

    outcome <- association_fl_outcome(artifact, "condition[B]")
    seed <- artifact$diagnostics$seed_manifest
    expect_s3_class(outcome, "imputefinder_association")
    expect_identical(outcome$diagnostics$variance_method, "HC3")
    expect_identical(outcome$diagnostics$permutation_mode, "monte_carlo")
    expect_identical(outcome$diagnostics$allowable_transformations, 362880.0)
    expect_identical(outcome$diagnostics$evaluated_transformations, 9999L)
    expect_identical(outcome$permutation_count, 10000L)
    expect_equal(outcome$effect, 0.269444444444445, tolerance = 1e-14)
    expect_equal(outcome$statistic, 8.40923535253228, tolerance = 1e-12)
    expect_identical(outcome$diagnostics$exceedance_count, 317L)
    expect_identical(outcome$raw_p, 0.0318)
    expect_identical(nrow(seed), 9999L)
    expect_identical(seed$draw_id, seq_len(9999L))
    expect_identical(anyDuplicated(seed$seed), 0L)
    expect_identical(anyDuplicated(seed$map_sha256), 0L)
    expect_true(all(grepl("^[0-9a-f]{64}$", seed$map_sha256)))
    expect_identical(seed$seed[c(1L, 9999L)], c(4099601L, 252728501L))
    expect_identical(
        seed$map_sha256[c(1L, 9999L)],
        c(
            "b242994955770d61f6f9f95c2634d848e3870c314941a39f29b88066e1139a2b",
            "043056a2628561e6268cef08c622f19eca5b3b9c7ce77bf23067f4c7a5cb968d"
        )
    )
    expect_invisible(
        imputefinder:::.validate_association_candidate_artifact(
            artifact,
            preparation
        )
    )
    wrong_seed <- artifact
    wrong_seed$diagnostics$seed_manifest$seed[[1L]] <-
        wrong_seed$diagnostics$seed_manifest$seed[[1L]] + 1L
    duplicate_map <- artifact
    duplicate_map$diagnostics$seed_manifest$map_sha256[[2L]] <-
        duplicate_map$diagnostics$seed_manifest$map_sha256[[1L]]
    over_retry <- artifact
    over_retry$diagnostics$seed_manifest$map_retry[[1L]] <- 1000001L
    for (corrupted in list(wrong_seed, duplicate_map, over_retry)) {
        expect_error(
            imputefinder:::.validate_association_candidate_artifact(
                corrupted,
                preparation
            ),
            class = "imputefinder_association_artifact_error"
        )
    }
})

test_that("retry exhaustion omits only the failed Monte Carlo seed rows", {
    counts <- c(4L, 7L, 9L, 11L, 14L, 6L, 8L, 10L, 13L, 16L,
        5L, 12L, 15L, 18L, 20L)
    x <- association_fl_detection_matrix(counts, 20L)
    metadata <- data.frame(
        condition = rep(c("A", "B", "C"), each = 5L),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    preparation <- imputefinder:::.new_association_preparation(
        x,
        missingness_design(metadata, condition = "condition")
    )
    records <- imputefinder:::.association_freedman_lane_records(preparation)
    expect_identical(
        unname(vapply(records, function(record) {
            record$permutations$mode
        }, character(1L))),
        rep("monte_carlo", 2L)
    )
    planned_seeds <- imputefinder:::.association_freedman_lane_seeds(
        records,
        preparation$input_sha256
    )
    failed_index <- match("condition[B]", preparation$hypotheses$label)
    failed_id <- preparation$hypotheses$hypothesis[[failed_index]]
    original_maps <- imputefinder:::.association_monte_carlo_maps
    artifact <- testthat::with_mocked_bindings(
        imputefinder:::.run_association_freedman_lane(preparation),
        .association_monte_carlo_maps = function(
            seed_manifest,
            design,
            n,
            groups = NULL,
            swaps = NULL
        ) {
            if (identical(unique(seed_manifest$hypothesis_id), failed_id)) {
                imputefinder:::.abort_association(
                    "Mocked Monte Carlo retry exhaustion.",
                    "imputefinder_association_permutation_retry_error"
                )
            }
            original_maps(
                seed_manifest,
                design,
                n,
                groups = groups,
                swaps = swaps
            )
        },
        .package = "imputefinder"
    )
    failed <- artifact$outcomes[[failed_id]]
    expect_s3_class(failed, "imputefinder_unavailable")
    expect_identical(failed$code, "association_numerical_failure")

    available_id <- setdiff(preparation$hypotheses$hypothesis, failed_id)
    expect_s3_class(
        artifact$outcomes[[available_id]],
        "imputefinder_association"
    )
    manifest <- artifact$diagnostics$seed_manifest
    expect_identical(nrow(manifest), 9999L)
    expect_identical(unique(manifest$hypothesis_id), available_id)
    expected <- planned_seeds[
        planned_seeds$hypothesis_id == available_id,
        ,
        drop = FALSE
    ]
    row.names(expected) <- NULL
    expect_identical(
        manifest[names(expected)],
        expected
    )
})

test_that("Freedman-Lane stage failures abstain before result allocation", {
    counts <- c(8L, 10L, 9L, 12L, 11L, 15L, 13L, 16L)
    x <- association_fl_detection_matrix(counts, 16L)
    metadata <- data.frame(
        condition = rep(c("A", "B"), 4L),
        run = seq_len(8L),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    preparation <- imputefinder:::.new_association_preparation(
        x,
        missingness_design(
            metadata,
            condition = "condition",
            nuisance = "run"
        )
    )
    artifact <- imputefinder:::.run_association_freedman_lane(preparation)
    condition <- association_fl_outcome(artifact, "condition[B]")
    expect_s3_class(condition, "imputefinder_unavailable")
    expect_identical(
        condition$code,
        "association_low_permutation_resolution"
    )

    blocked_metadata <- data.frame(
        condition = rep(c("A", "B"), 8L),
        run = c(rep(c(0, 1), 4L), rep(c(1, 0), 4L)),
        subject = rep(sprintf("p%02d", seq_len(8L)), each = 2L),
        row.names = sprintf("sample_%02d", seq_len(16L)),
        stringsAsFactors = FALSE
    )
    blocked_x <- association_fl_detection_matrix(
        c(8L, 12L, 10L, 13L, 9L, 15L, 12L, 16L,
          11L, 15L, 13L, 18L, 14L, 17L, 15L, 20L),
        20L
    )
    blocked <- imputefinder:::.new_association_preparation(
        blocked_x,
        missingness_design(
            blocked_metadata,
            condition = "condition",
            nuisance = "run",
            block = "subject"
        )
    )
    blocked_artifact <-
        imputefinder:::.run_association_freedman_lane(blocked)
    blocked_run <- association_fl_outcome(blocked_artifact, "run")
    expect_s3_class(blocked_run, "imputefinder_unavailable")
    expect_identical(
        blocked_run$code,
        "association_incompatible_permutation"
    )

    acquisition_metadata <- data.frame(
        condition = c(rep("A", 4L), rep("B", 5L)),
        acquisition = rep("D|A", 9L),
        row.names = sprintf("sample_%02d", seq_len(9L)),
        stringsAsFactors = FALSE
    )
    acquisition <- imputefinder:::.new_association_preparation(
        association_fl_detection_matrix(
            c(8L, 10L, 9L, 12L, 11L, 15L, 13L, 16L, 18L),
            18L
        ),
        missingness_design(
            acquisition_metadata,
            condition = "condition",
            acquisition = "acquisition"
        )
    )
    acquisition_artifact <-
        imputefinder:::.run_association_freedman_lane(acquisition)
    expect_s3_class(
        acquisition_artifact$outcomes[[1L]],
        "imputefinder_association"
    )
    expect_identical(
        acquisition_artifact$outcomes[[1L]]$diagnostics$permutation_mode,
        "monte_carlo"
    )
    expect_identical(
        unique(acquisition_artifact$diagnostics$seed_manifest$acquisition),
        "D|A"
    )
})
