association_detection_matrix <- function(counts, feature_count = max(counts)) {
    sample <- sprintf("sample_%02d", seq_along(counts))
    x <- matrix(
        NA_real_,
        nrow = feature_count,
        ncol = length(counts),
        dimnames = list(sprintf("feature_%02d", seq_len(feature_count)), sample)
    )
    for (index in seq_along(counts)) {
        x[seq_len(counts[[index]]), index] <- 1
    }
    x
}

association_hc3_fixture <- function() {
    counts <- c(8L, 11L, 10L, 14L, 12L, 15L, 13L, 16L, 15L, 18L, 17L, 20L)
    x <- association_detection_matrix(counts, 20L)
    metadata <- data.frame(
        condition = rep(c("A", "B"), 6L),
        run = seq_len(12L),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    list(
        x = x,
        design = missingness_design(
            metadata,
            condition = "condition",
            nuisance = "run"
        )
    )
}

association_cr2_fixture <- function() {
    counts <- c(8L, 12L, 10L, 13L, 9L, 15L, 12L, 16L, 11L, 15L, 13L, 18L, 14L, 17L, 15L, 20L)
    x <- association_detection_matrix(counts, 20L)
    metadata <- data.frame(
        condition = rep(c("A", "B"), 8L),
        subject = rep(sprintf("p%02d", seq_len(8L)), each = 2L),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    list(
        x = x,
        design = missingness_design(
            metadata,
            condition = "condition",
            block = "subject"
        )
    )
}

association_outcome_by_label <- function(artifact, label) {
    index <- which(artifact$hypotheses$label == label)
    artifact$outcomes[[artifact$hypotheses$hypothesis[[index]]]]
}

test_that("independent OLS uses the frozen HC3 rank-coordinate estimator", {
    fixture <- association_hc3_fixture()
    preparation <- imputefinder:::.new_association_preparation(
        fixture$x,
        fixture$design
    )
    artifact <- imputefinder:::.run_association_ols_hc3_cr2(preparation)
    outcome <- association_outcome_by_label(artifact, "condition[B]")

    expect_s3_class(artifact, "imputefinder_association_candidate_artifact")
    expect_identical(
        names(artifact),
        c(
            "schema", "protocol", "candidate", "input_sha256", "response",
            "hypotheses", "support", "outcomes", "multiplicity",
            "diagnostics"
        )
    )
    expect_identical(artifact$schema, "association_candidate_artifact_v1")
    expect_identical(artifact$candidate, "a_fraction_ols_hc3_cr2")
    expect_identical(artifact$protocol, preparation$protocol)
    expect_identical(artifact$input_sha256, preparation$input_sha256)
    expect_true(all(artifact$hypotheses$family_member))
    expect_identical(
        names(outcome),
        c(
            "status", "quantity", "candidate", "stratum", "hypothesis",
            "effect", "standard_error", "conf_low", "conf_high",
            "statistic", "reference_df", "raw_p", "adjusted_p", "flag",
            "log_odds", "log_odds_standard_error", "log_odds_conf_low",
            "log_odds_conf_high", "permutation_count", "diagnostics"
        )
    )
    expect_identical(class(outcome), "imputefinder_association")
    expect_equal(outcome$effect, 0.11583333333333329, tolerance = 1e-13)
    expect_equal(outcome$standard_error, 0.01543962314284522, tolerance = 1e-13)
    expect_equal(outcome$statistic, 7.502342010660468, tolerance = 1e-12)
    expect_identical(outcome$reference_df, 9.0)
    expect_equal(outcome$raw_p, 3.68367944543356e-05, tolerance = 1e-13)
    expect_equal(outcome$conf_low, 0.08090647924984104, tolerance = 1e-13)
    expect_equal(outcome$conf_high, 0.15076018741682554, tolerance = 1e-13)
    expect_identical(outcome$diagnostics$variance_method, "HC3")
    expect_identical(outcome$diagnostics$residual_df, 9L)
    expect_equal(
        outcome$diagnostics$leverage_max,
        0.34523809523809523,
        tolerance = 1e-13
    )
    raw <- vapply(artifact$outcomes, `[[`, numeric(1L), "raw_p")
    adjusted <- vapply(artifact$outcomes, `[[`, numeric(1L), "adjusted_p")
    expect_identical(adjusted, stats::p.adjust(raw, method = "holm"))
    expect_identical(artifact$multiplicity$family_size, 2L)
    expect_invisible(
        imputefinder:::.validate_association_candidate_artifact(
            artifact,
            preparation
        )
    )

    scaled_artifact <- function(scale) {
        scaled_design <- fixture$design$sample_data
        scaled_design$run <- scaled_design$run * scale
        imputefinder:::.run_association_ols_hc3_cr2(
            imputefinder:::.new_association_preparation(
                fixture$x,
                missingness_design(
                    scaled_design,
                    condition = "condition",
                    nuisance = "run"
                )
            )
        )
    }
    scaled <- scaled_artifact(3e-8)
    scaled_condition <- association_outcome_by_label(scaled, "condition[B]")
    expect_s3_class(scaled_condition, "imputefinder_association")
    expect_equal(scaled_condition$effect, outcome$effect, tolerance = 1e-10)
    expect_equal(
        scaled_condition$standard_error,
        outcome$standard_error,
        tolerance = 1e-10
    )
    expect_equal(
        scaled_condition$statistic,
        outcome$statistic,
        tolerance = 1e-10
    )
    expect_true(all(scaled$hypotheses$family_member))

    below_global_covariance_scale <- scaled_artifact(1e-9)
    below_condition <- association_outcome_by_label(
        below_global_covariance_scale,
        "condition[B]"
    )
    base_run <- association_outcome_by_label(artifact, "run")
    scaled_run <- association_outcome_by_label(
        below_global_covariance_scale,
        "run"
    )
    expect_identical(
        below_condition$code,
        "association_singular_covariance"
    )
    expect_s3_class(scaled_run, "imputefinder_association")
    expect_equal(scaled_run$effect * 1e-9, base_run$effect, tolerance = 1e-10)
    expect_equal(
        scaled_run$standard_error * 1e-9,
        base_run$standard_error,
        tolerance = 1e-10
    )
    expect_equal(scaled_run$statistic, base_run$statistic, tolerance = 1e-10)
    expect_equal(scaled_run$raw_p, base_run$raw_p, tolerance = 1e-10)
})

test_that("blocked OLS uses full-design CR2 and Satterthwaite df", {
    fixture <- association_cr2_fixture()
    preparation <- imputefinder:::.new_association_preparation(
        fixture$x,
        fixture$design
    )
    artifact <- imputefinder:::.run_association_ols_hc3_cr2(preparation)
    outcome <- association_outcome_by_label(artifact, "condition[B]")

    expect_identical(outcome$diagnostics$variance_method, "CR2")
    expect_equal(outcome$effect, 0.2125, tolerance = 1e-13)
    expect_equal(outcome$standard_error, 0.018298126367785, tolerance = 1e-9)
    expect_equal(outcome$statistic, 11.613210868087545, tolerance = 1e-8)
    expect_equal(outcome$reference_df, 7, tolerance = 1e-10)
    expect_equal(outcome$raw_p, 7.919630641550018e-06, tolerance = 1e-9)
    expect_equal(outcome$conf_low, 0.1692318066320262, tolerance = 1e-9)
    expect_equal(outcome$conf_high, 0.2557681933679738, tolerance = 1e-9)
    expect_equal(outcome$diagnostics$satterthwaite_df, 7, tolerance = 1e-10)
    expect_identical(outcome$diagnostics$residual_df, 7L)
    expect_identical(artifact$multiplicity$family_size, 1L)
    expect_invisible(
        imputefinder:::.validate_association_candidate_artifact(
            artifact,
            preparation
        )
    )

    singleton_metadata <- fixture$design$sample_data
    singleton_metadata$batch <- c(
        rep("common", nrow(singleton_metadata) - 1L),
        "singleton"
    )
    singleton <- imputefinder:::.run_association_ols_hc3_cr2(
        imputefinder:::.new_association_preparation(
            fixture$x,
            missingness_design(
                singleton_metadata,
                condition = "condition",
                nuisance = "batch",
                block = "subject"
            )
        )
    )
    singleton_condition <- association_outcome_by_label(
        singleton,
        "condition[B]"
    )
    expect_s3_class(singleton_condition, "imputefinder_association")
    expect_identical(singleton_condition$diagnostics$variance_method, "CR2")
    expect_equal(singleton_condition$diagnostics$leverage_max, 1)
    expect_gte(singleton_condition$reference_df, 5)

    low_reference_x <- association_detection_matrix(
        c(8L, 10L, 9L, 14L, 10L, 11L, 11L, 17L, 12L, 15L, 13L, 20L),
        20L
    )
    low_reference_metadata <- data.frame(
        condition = rep(c("A", "B"), 6L),
        run = c(0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6),
        subject = rep(sprintf("p%02d", seq_len(6L)), each = 2L),
        row.names = colnames(low_reference_x),
        stringsAsFactors = FALSE
    )
    low_reference <- imputefinder:::.run_association_ols_hc3_cr2(
        imputefinder:::.new_association_preparation(
            low_reference_x,
            missingness_design(
                low_reference_metadata,
                condition = "condition",
                nuisance = "run",
                block = "subject"
            )
        )
    )
    low_reference_condition <- which(
        low_reference$hypotheses$label == "condition[B]"
    )
    expect_true(low_reference$hypotheses$estimable[low_reference_condition])
    expect_true(low_reference$support$eligible[low_reference_condition])
    expect_identical(
        low_reference$outcomes[[low_reference_condition]]$code,
        "association_low_reference_df"
    )
})

test_that("declared acquisition strata own separate Holm families", {
    counts <- c(
        8L, 10L, 9L, 11L, 12L, 15L, 13L, 16L,
        9L, 12L, 10L, 14L, 13L, 17L, 15L, 20L
    )
    x <- association_detection_matrix(counts, 20L)
    metadata <- data.frame(
        condition = rep(rep(c("A", "B"), each = 4L), 2L),
        run = rep(seq_len(8L), 2L),
        acquisition = rep(c("DDA", "DIA"), each = 8L),
        row.names = colnames(x),
        stringsAsFactors = FALSE
    )
    preparation <- imputefinder:::.new_association_preparation(
        x,
        missingness_design(
            metadata,
            condition = "condition",
            nuisance = "run",
            acquisition = "acquisition"
        )
    )
    artifact <- imputefinder:::.run_association_ols_hc3_cr2(preparation)

    expect_length(unique(artifact$response$stratum), 2L)
    expect_identical(artifact$multiplicity$family_size, c(2L, 2L))
    expect_true(all(artifact$hypotheses$family_member))
    for (stratum in unique(artifact$hypotheses$stratum)) {
        selected <- artifact$hypotheses$stratum == stratum
        raw <- vapply(
            artifact$outcomes[selected],
            `[[`,
            numeric(1L),
            "raw_p"
        )
        adjusted <- vapply(
            artifact$outcomes[selected],
            `[[`,
            numeric(1L),
            "adjusted_p"
        )
        expect_identical(adjusted, stats::p.adjust(raw, method = "holm"))
    }
    expect_invisible(
        imputefinder:::.validate_association_candidate_artifact(
            artifact,
            preparation
        )
    )
})

test_that("robust candidate enforces unavailable precedence", {
    confounded_x <- association_detection_matrix(
        c(8L, 10L, 12L, 14L, 11L, 13L, 15L, 20L),
        20L
    )
    confounded_metadata <- data.frame(
        condition = rep(c("A", "B"), each = 4L),
        batch = rep(c("x", "y"), each = 4L),
        row.names = colnames(confounded_x),
        stringsAsFactors = FALSE
    )
    confounded <- imputefinder:::.run_association_ols_hc3_cr2(
        imputefinder:::.new_association_preparation(
            confounded_x,
            missingness_design(
                confounded_metadata,
                condition = "condition",
                nuisance = "batch"
            )
        )
    )
    expect_identical(
        unname(vapply(confounded$outcomes, `[[`, character(1L), "code")),
        rep("association_nonestimable", 2L)
    )

    low_support_x <- association_detection_matrix(
        c(8L, 10L, 12L, 14L, 16L, 20L),
        20L
    )
    low_support_metadata <- data.frame(
        condition = rep(c("A", "B"), each = 3L),
        row.names = colnames(low_support_x),
        stringsAsFactors = FALSE
    )
    low_support <- imputefinder:::.run_association_ols_hc3_cr2(
        imputefinder:::.new_association_preparation(
            low_support_x,
            missingness_design(low_support_metadata, condition = "condition")
        )
    )
    expect_identical(
        low_support$outcomes[[1L]]$code,
        "association_low_independent_support"
    )

    low_df_metadata <- expand.grid(
        n3 = c(0, 1),
        n2 = c(0, 1),
        condition = c("A", "B"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    rownames(low_df_metadata) <- sprintf("sample_%02d", seq_len(8L))
    low_df_x <- association_detection_matrix(
        c(8L, 10L, 12L, 14L, 11L, 15L, 17L, 20L),
        20L
    )
    colnames(low_df_x) <- rownames(low_df_metadata)
    low_df <- imputefinder:::.run_association_ols_hc3_cr2(
        imputefinder:::.new_association_preparation(
            low_df_x,
            missingness_design(
                low_df_metadata,
                condition = "condition",
                nuisance = c("n2", "n3"),
                interactions = list(
                    c("condition", "n2"),
                    c("condition", "n3")
                )
            )
        )
    )
    condition <- which(low_df$hypotheses$label == "condition[B]")
    expect_true(low_df$support$eligible[condition])
    expect_true(low_df$hypotheses$estimable[condition])
    expect_identical(
        low_df$outcomes[[condition]]$code,
        "association_low_residual_df"
    )

    degenerate_metadata <- data.frame(
        condition = rep(c("A", "B"), each = 4L),
        row.names = sprintf("sample_%02d", seq_len(8L)),
        stringsAsFactors = FALSE
    )
    degenerate_x <- matrix(
        1,
        nrow = 1L,
        ncol = 8L,
        dimnames = list("feature", rownames(degenerate_metadata))
    )
    degenerate <- imputefinder:::.run_association_ols_hc3_cr2(
        imputefinder:::.new_association_preparation(
            degenerate_x,
            missingness_design(degenerate_metadata, condition = "condition")
        )
    )
    expect_identical(
        degenerate$outcomes[[1L]]$code,
        "association_degenerate_response"
    )
})

test_that("rank deficiency and numerical abstention remain explicit", {
    rank_deficient_x <- association_detection_matrix(
        c(8L, 11L, 10L, 14L, 12L, 15L, 13L, 16L),
        16L
    )
    rank_deficient_metadata <- data.frame(
        condition = rep(c("A", "B"), 4L),
        constant = rep(1, 8L),
        row.names = colnames(rank_deficient_x),
        stringsAsFactors = FALSE
    )
    rank_deficient <- imputefinder:::.run_association_ols_hc3_cr2(
        imputefinder:::.new_association_preparation(
            rank_deficient_x,
            missingness_design(
                rank_deficient_metadata,
                condition = "condition",
                nuisance = "constant"
            )
        )
    )
    condition <- which(rank_deficient$hypotheses$label == "condition[B]")
    constant <- which(rank_deficient$hypotheses$label == "constant")
    expect_s3_class(
        rank_deficient$outcomes[[condition]],
        "imputefinder_association"
    )
    expect_equal(
        rank_deficient$outcomes[[condition]]$effect,
        0.203125,
        tolerance = 1e-13
    )
    expect_identical(
        rank_deficient$outcomes[[constant]]$code,
        "association_nonestimable"
    )

    high_leverage_x <- association_detection_matrix(
        c(8L, 11L, 10L, 14L, 12L, 15L, 13L, 16L, 15L, 18L),
        18L
    )
    high_leverage_metadata <- data.frame(
        condition = rep(c("A", "B"), 5L),
        batch = c(rep("common", 9L), "singleton"),
        row.names = colnames(high_leverage_x),
        stringsAsFactors = FALSE
    )
    high_leverage <- imputefinder:::.run_association_ols_hc3_cr2(
        imputefinder:::.new_association_preparation(
            high_leverage_x,
            missingness_design(
                high_leverage_metadata,
                condition = "condition",
                nuisance = "batch"
            )
        )
    )
    condition <- which(high_leverage$hypotheses$label == "condition[B]")
    singleton <- which(
        high_leverage$hypotheses$label == "batch[singleton]"
    )
    expect_identical(
        high_leverage$outcomes[[condition]]$code,
        "association_numerical_failure"
    )
    expect_identical(
        high_leverage$outcomes[[singleton]]$code,
        "association_low_independent_support"
    )

    singular_x <- association_detection_matrix(
        c(rep(8L, 4L), rep(12L, 4L)),
        12L
    )
    singular_metadata <- data.frame(
        condition = rep(c("A", "B"), each = 4L),
        row.names = colnames(singular_x),
        stringsAsFactors = FALSE
    )
    singular <- imputefinder:::.run_association_ols_hc3_cr2(
        imputefinder:::.new_association_preparation(
            singular_x,
            missingness_design(singular_metadata, condition = "condition")
        )
    )
    expect_identical(
        singular$outcomes[[1L]]$code,
        "association_singular_covariance"
    )
})

test_that("common PSD cleanup clamps only tolerance-scale eigenvalues", {
    tolerance <- 2 * .Machine$double.eps
    clamped <- imputefinder:::.association_clean_psd(diag(c(1, -tolerance)))

    expect_type(clamped, "list")
    expect_equal(clamped$matrix, diag(c(1, 0)), tolerance = 1e-15)
    expect_equal(clamped$tolerance, tolerance, tolerance = 1e-15)
    expect_null(imputefinder:::.association_clean_psd(diag(c(1, -1e-8))))

    angle <- pi / 4
    rotation <- matrix(
        c(cos(angle), sin(angle), -sin(angle), cos(angle)),
        nrow = 2L
    )
    ill_conditioned <- rotation %*% diag(c(4e10, 3e-5)) %*% t(rotation)
    symmetrized <- (ill_conditioned + t(ill_conditioned)) / 2
    unchanged <- imputefinder:::.association_clean_psd(ill_conditioned)
    expect_identical(unchanged$matrix, symmetrized)

    projection <- matrix(c(1, 2, -1, 3, 2, -2, 1, 4), nrow = 4L)
    reference_df <- imputefinder:::.association_cr2_reference_df(projection)
    expect_true(is.finite(reference_df))
    expect_equal(
        imputefinder:::.association_cr2_reference_df(projection * 1e150),
        reference_df,
        tolerance = 1e-14
    )
    expect_equal(
        imputefinder:::.association_cr2_reference_df(projection * 1e-150),
        reference_df,
        tolerance = 1e-14
    )
})

test_that("candidate artifact validation rejects Holm and join drift", {
    fixture <- association_hc3_fixture()
    preparation <- imputefinder:::.new_association_preparation(
        fixture$x,
        fixture$design
    )
    artifact <- imputefinder:::.run_association_ols_hc3_cr2(preparation)
    wrong_adjustment <- artifact
    wrong_adjustment$outcomes[[1L]]$adjusted_p <-
        min(1, wrong_adjustment$outcomes[[1L]]$raw_p + 0.1)
    wrong_quantity <- artifact
    wrong_quantity$outcomes[[1L]]$quantity <- "other"
    wrong_family <- artifact
    wrong_family$hypotheses$family_member[[1L]] <- FALSE
    wrong_hash <- artifact
    wrong_hash$diagnostics$input_sha256 <- paste(rep("0", 64L), collapse = "")
    wrong_rank_join <- artifact
    wrong_rank_join$outcomes[[1L]]$diagnostics$rank <-
        wrong_rank_join$outcomes[[1L]]$diagnostics$rank + 1L
    malformed_diagnostics <- artifact
    malformed_diagnostics$diagnostics <- 1
    detached_hash <- artifact
    detached_hash$input_sha256 <- paste(rep("0", 64L), collapse = "")
    detached_hash$diagnostics$input_sha256 <- detached_hash$input_sha256
    wrong_seed_types <- artifact
    wrong_seed_types$diagnostics$seed_manifest <- as.data.frame(
        stats::setNames(
            rep(list(logical()), 9L),
            names(wrong_seed_types$diagnostics$seed_manifest)
        )
    )

    singular_x <- association_detection_matrix(
        c(rep(8L, 4L), rep(12L, 4L)),
        12L
    )
    singular_metadata <- data.frame(
        condition = rep(c("A", "B"), each = 4L),
        row.names = colnames(singular_x),
        stringsAsFactors = FALSE
    )
    singular_preparation <- imputefinder:::.new_association_preparation(
        singular_x,
        missingness_design(singular_metadata, condition = "condition")
    )
    singular_artifact <- imputefinder:::.run_association_ols_hc3_cr2(
        singular_preparation
    )
    impossible_code <- singular_artifact
    impossible_code$outcomes[[1L]]$code <-
        "association_low_permutation_resolution"
    impossible_code$diagnostics$warnings <-
        "association_low_permutation_resolution"
    impossible_stage <- singular_artifact
    impossible_stage$outcomes[[1L]]$code <-
        "association_low_independent_support"
    impossible_stage$diagnostics$warnings <-
        "association_low_independent_support"

    for (corrupted in list(
        wrong_adjustment,
        wrong_quantity,
        wrong_family,
        wrong_hash,
        wrong_rank_join,
        malformed_diagnostics,
        detached_hash,
        wrong_seed_types
    )) {
        expect_error(
            imputefinder:::.validate_association_candidate_artifact(
                corrupted,
                preparation
            ),
            class = "imputefinder_association_artifact_error"
        )
    }
    expect_error(
        imputefinder:::.validate_association_candidate_artifact(
            impossible_code,
            singular_preparation
        ),
        class = "imputefinder_association_artifact_error"
    )
    expect_error(
        imputefinder:::.validate_association_candidate_artifact(
            impossible_stage,
            singular_preparation
        ),
        class = "imputefinder_association_artifact_error"
    )
    expect_error(
        imputefinder:::.validate_association_candidate_artifact(
            singular_artifact,
            preparation
        ),
        class = "imputefinder_association_artifact_error"
    )
    expect_error(
        imputefinder:::.validate_association_candidate_artifact(artifact),
        class = "imputefinder_association_artifact_error"
    )
})
