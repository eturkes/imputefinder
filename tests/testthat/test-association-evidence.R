association_metric_fixture <- function() {
    run <- imputefinder:::.association_run_grid()
    run$input_sha256 <- strrep("0", 64L)
    run$result_sha256 <- strrep("1", 64L)
    run$result_kind <- "panel"
    runtime <- c(
        a_fraction_ols_hc3_cr2 = 3,
        a_fraction_freedman_lane = 1,
        a_fraction_quasibinomial = 2
    )
    run$elapsed_seconds <- as.double(unname(runtime[run$candidate_id]))
    run$status <- "measured"
    run$failure_code <- NA_character_
    run <- run[imputefinder:::.ASSOCIATION_RUN_BINDING_FIELDS]

    target_grid <- imputefinder:::.association_target_grid()
    opportunity <- target_grid[c(
        "candidate_id", "scenario_id", "replicate"
    )]
    opportunity$stratum <- "all"
    opportunity$hypothesis <- paste0("a_", strrep("2", 64L))
    opportunity$available <- TRUE
    opportunity <- opportunity[
        imputefinder:::.ASSOCIATION_OPPORTUNITY_AUDIT_FIELDS
    ]
    opportunity <- imputefinder:::.association_canonicalize_audit(
        opportunity,
        c(
            "candidate_id", "scenario_id", "replicate", "stratum",
            "hypothesis"
        )
    )

    null <- imputefinder:::.association_null_grid()
    null$available_p_count <- 1L
    null$false_flag <- FALSE
    null <- null[imputefinder:::.ASSOCIATION_NULL_AUDIT_FIELDS]

    target <- target_grid
    target$outcome_status <- "available"
    target$code <- NA_character_
    target$truth_status <- "available"
    target$truth <- 0.20
    target$effect <- 0.20
    target$conf_low <- 0.10
    target$conf_high <- 0.30
    target$adjusted_p <- 0.01
    target <- target[imputefinder:::.ASSOCIATION_TARGET_AUDIT_FIELDS]

    list(run = run, opportunity = opportunity, null = null, target = target)
}

association_absent_study_input <- function(scenario_id, replicate = 1L) {
    manifest <- imputefinder:::.association_study_manifest_row(scenario_id)
    sizes <- imputefinder:::.association_parse_condition_sizes(
        manifest$condition_sizes[[1L]]
    )
    condition <- rep(names(sizes), sizes)
    sample <- sprintf("sample_%03d", seq_along(condition))
    metadata <- data.frame(
        condition = condition,
        row.names = sample,
        stringsAsFactors = FALSE
    )
    nuisance <- manifest$nuisance[[1L]]
    if (identical(nuisance, "run_order")) {
        metadata$run_order <- seq_along(condition)
    } else if (identical(nuisance, "batch")) {
        metadata$batch <- rep(c("batch_1", "batch_2"), length.out = nrow(metadata))
    } else if (identical(nuisance, "technical_replicate")) {
        metadata$technical_replicate <- rep(
            c("tech_01", "tech_02"),
            length.out = nrow(metadata)
        )
    }
    block <- manifest$block[[1L]]
    if (identical(block, "subject")) {
        per_condition <- unname(sizes[[1L]])
        technical_count <- if (identical(
            nuisance,
            "technical_replicate"
        )) 2L else 1L
        subjects <- sprintf(
            "subject_%02d",
            rep(
                seq_len(per_condition / technical_count),
                each = technical_count
            )
        )
        metadata$subject <- rep(subjects, length(sizes))
    }
    metadata$acquisition_mode <- manifest$acquisition[[1L]]
    design <- missingness_design(
        metadata,
        condition = "condition",
        nuisance = if (nzchar(nuisance)) nuisance else NULL,
        block = if (nzchar(block)) block else NULL,
        acquisition = "acquisition_mode"
    )
    features <- sprintf(
        "protein_%05d",
        seq_len(manifest$n_features[[1L]])
    )
    data <- matrix(
        NA_real_,
        nrow = length(features),
        ncol = length(sample),
        dimnames = list(features, sample)
    )
    probability <- matrix(
        1.0,
        nrow = nrow(data),
        ncol = ncol(data),
        dimnames = dimnames(data)
    )
    imputefinder:::.new_association_study_input(
        scenario_id,
        as.integer(replicate),
        data,
        design,
        probability
    )
}

association_absent_resolver <- function(
    implementation_hash = strrep("a", 64L)
) {
    inputs <- new.env(hash = TRUE, parent = emptyenv())
    payload_hashes <- new.env(hash = TRUE, parent = emptyenv())
    abstentions <- new.env(hash = TRUE, parent = emptyenv())
    for (scenario_id in imputefinder:::.ASSOCIATION_STUDY_SCENARIOS) {
        input <- association_absent_study_input(scenario_id)
        assign(scenario_id, input, envir = inputs)
        assign(
            scenario_id,
            imputefinder:::.association_study_input_payload_sha256(input),
            envir = payload_hashes
        )
        assign(
            scenario_id,
            imputefinder:::.new_association_preparation(
                input$data,
                input$design
            ),
            envir = abstentions
        )
    }
    input_reader <- function(scenario_id, replicate) {
        input <- get(scenario_id, envir = inputs, inherits = FALSE)
        input$replicate <- replicate
        input
    }
    result_reader <- function(candidate_id, scenario_id, replicate) {
        input <- get(scenario_id, envir = inputs, inherits = FALSE)
        input$replicate <- replicate
        imputefinder:::.new_association_study_result(
            candidate_id,
            scenario_id,
            replicate,
            imputefinder:::.association_study_input_sha256(
                input,
                get(scenario_id, envir = payload_hashes, inherits = FALSE)
            ),
            get(scenario_id, envir = abstentions, inherits = FALSE),
            0.0
        )
    }
    imputefinder:::.new_association_study_resolver(
        input_reader,
        result_reader,
        implementation_hash
    )
}

test_that("production study constants are exact frozen-contract projections", {
    root <- normalizePath(
        testthat::test_path("..", ".."),
        winslash = "/",
        mustWork = TRUE
    )
    path <- file.path(root, "dev", "m13-association-contract.R")
    if (!file.exists(path)) {
        testthat::skip("repository-only frozen-contract projection rail")
    }
    contract <- new.env(parent = globalenv())
    previous <- getwd()
    on.exit(setwd(previous), add = TRUE)
    setwd(root)
    sys.source(path, envir = contract)

    expected_manifest <- contract$m13a_effective_manifest()
    expect_identical(
        imputefinder:::.ASSOCIATION_STUDY_MANIFEST$scenario_id,
        expected_manifest$scenario_id
    )
    expect_identical(
        imputefinder:::.ASSOCIATION_STUDY_MANIFEST$acquisition,
        expected_manifest$acquisition
    )
    expect_identical(
        imputefinder:::.ASSOCIATION_STUDY_MANIFEST$n_features,
        expected_manifest$n_features
    )
    expect_identical(
        imputefinder:::.ASSOCIATION_STUDY_MANIFEST$condition_sizes,
        expected_manifest$condition_sizes
    )
    expected_nuisance <- vapply(
        expected_manifest$nuisance_role,
        function(value) switch(
            value,
            none = "",
            run_order = "run_order",
            batch_crossed = "batch",
            batch_partial = "batch",
            batch_aliased = "batch",
            technical_replicate = "technical_replicate"
        ),
        character(1L)
    )
    expected_block <- ifelse(
        expected_manifest$block_role == "subject",
        "subject",
        ""
    )
    expect_identical(
        imputefinder:::.ASSOCIATION_STUDY_MANIFEST$nuisance,
        unname(expected_nuisance)
    )
    expect_identical(
        imputefinder:::.ASSOCIATION_STUDY_MANIFEST$block,
        unname(expected_block)
    )

    expected_gates <- contract$m13a_gate_registry()
    expected_gates <- expected_gates[
        expected_gates$stage == "selected_winner_synthetic",
        ,
        drop = FALSE
    ]
    gates <- imputefinder:::.association_study_gate_registry()
    expect_identical(gates$gate_id, expected_gates$gate_id)
    expect_identical(gates$m12_gate_id, expected_gates$m12_gate_id)
    expect_identical(gates$registry_version, expected_gates$registry_version)
    expect_identical(gates$metric, expected_gates$metric)
    expect_identical(gates$operator, expected_gates$threshold_operator)
    expect_identical(gates$threshold, expected_gates$threshold_value)
    expect_identical(
        gates$failure_treatment,
        expected_gates$failure_treatment
    )
    expect_identical(
        gates$gate_binding_hash,
        expected_gates$gate_binding_hash
    )
    runtime <- imputefinder:::.association_load_protocol_generator(root)
    expect_true(is.function(runtime$simulate))
    expect_identical(
        parent.env(environment(runtime$simulate)),
        asNamespace("stats")
    )
    expect_true(exists(
        "ave",
        envir = environment(runtime$simulate),
        inherits = TRUE
    ))
})

test_that("study identities and selected-winner gates exactly mirror v4", {
    expect_identical(
        imputefinder:::.ASSOCIATION_CANDIDATES,
        c(
            "a_fraction_ols_hc3_cr2",
            "a_fraction_freedman_lane",
            "a_fraction_quasibinomial"
        )
    )
    expect_identical(
        imputefinder:::.ASSOCIATION_STUDY_REPLICATES,
        1:32
    )
    expect_identical(
        imputefinder:::.ASSOCIATION_TARGETS$label,
        c(
            "batch[batch_2]", "condition[D]", "batch[batch_3]",
            "condition[B]"
        )
    )
    gates <- imputefinder:::.association_study_gate_registry()
    expect_identical(
        gates$gate_binding_hash,
        c(
            "4f00acdaa6db9e9035e4b4b7ee5980f7b3beb835474ca1c451bd5fae4cda1c02",
            "de518da9a2aed4bdb1cfdbc50831136723332dfa5497f6f7a1af209eebe00dea",
            "0a1fe11d0760a924af5e91785f1843914b9b7c7a3e95b69282d8b092a045bcc7",
            "9c3ba375943078477f491a4825265e783d8426099437da64c6369070a343418b",
            "63df7d7ef92ee9b6c29f74f1c6f6667d064c5138e33360665b352e11de79de58"
        )
    )
    expect_identical(gates$threshold, c(0.15, 0.125, 0.90, 0.70, 0.03))
})

test_that("input and resolver records reject scenario and implementation drift", {
    input <- association_absent_study_input("dda_null_balanced")
    expect_invisible(imputefinder:::.validate_association_study_input(
        input,
        "dda_null_balanced",
        1L
    ))

    corrupted <- input
    corrupted$missing_probability[[1L]] <- 1.1
    expect_error(
        imputefinder:::.validate_association_study_input(
            corrupted,
            "dda_null_balanced",
            1L
        ),
        class = "imputefinder_association_evidence_error"
    )
    expect_error(
        imputefinder:::.validate_association_study_input(
            input,
            "dda_null_balanced",
            2L
        ),
        class = "imputefinder_association_evidence_error"
    )

    baseline_hash <- imputefinder:::.association_study_input_sha256(input)
    probability_changed <- input
    probability_changed$missing_probability[] <- 0.9
    expect_invisible(imputefinder:::.validate_association_study_input(
        probability_changed,
        "dda_null_balanced",
        1L
    ))
    expect_false(identical(
        baseline_hash,
        imputefinder:::.association_study_input_sha256(probability_changed)
    ))
    replicate_changed <- input
    replicate_changed$replicate <- 2L
    expect_false(identical(
        baseline_hash,
        imputefinder:::.association_study_input_sha256(replicate_changed)
    ))

    resolver <- association_absent_resolver()
    expect_invisible(imputefinder:::.validate_association_study_resolver(
        resolver,
        strrep("a", 64L)
    ))
    expect_error(
        imputefinder:::.validate_association_study_resolver(
            resolver,
            strrep("b", 64L)
        ),
        class = "imputefinder_association_evidence_error"
    )
    expect_identical(
        names(formals(
            imputefinder:::.new_association_protocol_study_resolver
        )),
        c("implementation_manifest", "root")
    )
    expect_identical(
        names(formals(
            imputefinder:::.association_execute_protocol_candidate
        )),
        c("candidate_id", "preparation")
    )
    expect_false(any(c(
        "input_provenance", "result_provenance"
    ) %in% imputefinder:::.ASSOCIATION_STUDY_RESOLVER_FIELDS))
})

test_that("production artifact inventory rejects missing locks and extra allocations", {
    layout <- imputefinder:::.association_protocol_artifact_layout(FALSE)
    expect_identical(sum(layout$kind == "directory"), 57L)
    expect_identical(sum(layout$kind == "file"), 1665L)
    expect_false(any(grepl("replicate-033", layout$relative_path)))
    allocation <- layout
    allocation$sha256 <- ifelse(
        allocation$kind == "file",
        strrep("1", 64L),
        NA_character_
    )
    allocation_hash <-
        imputefinder:::.association_protocol_allocation_inventory_hash(
            allocation
        )
    changed_allocation <- allocation
    file_row <- which(changed_allocation$kind == "file")[[1L]]
    changed_allocation$sha256[[file_row]] <- strrep("2", 64L)
    expect_false(identical(
        allocation_hash,
        imputefinder:::.association_protocol_allocation_inventory_hash(
            changed_allocation
        )
    ))
    with_evidence <- rbind(
        allocation,
        data.frame(
            relative_path =
                imputefinder:::.ASSOCIATION_STUDY_EVIDENCE_ARTIFACT,
            kind = "file",
            sha256 = strrep("3", 64L),
            stringsAsFactors = FALSE
        )
    )
    with_evidence <- with_evidence[
        order(with_evidence$relative_path, method = "radix"),
        ,
        drop = FALSE
    ]
    expect_identical(
        imputefinder:::.association_protocol_allocation_inventory_hash(
            with_evidence
        ),
        allocation_hash
    )

    root <- withr::local_tempdir(pattern = "m13a-inventory-")
    artifact_root <- file.path(
        root,
        imputefinder:::.ASSOCIATION_STUDY_ARTIFACT_DIRECTORY
    )
    dir.create(artifact_root, recursive = TRUE)
    implementation_hash <- strrep("a", 64L)
    locked <- list(manifest_hash = implementation_hash)
    expect_error(
        imputefinder:::.new_association_protocol_study_resolver(
            locked,
            root
        ),
        class = "imputefinder_association_evidence_error"
    )
    expect_error(
        imputefinder:::.association_protocol_artifact_inventory(
            normalizePath(artifact_root, winslash = "/"),
            implementation_hash,
            complete = FALSE,
            evidence_policy = "allow"
        ),
        class = "imputefinder_association_evidence_error"
    )

    saveRDS(
        locked,
        file.path(
            artifact_root,
            imputefinder:::.ASSOCIATION_STUDY_MANIFEST_ARTIFACT
        )
    )
    inventory <- imputefinder:::.association_protocol_artifact_inventory(
        normalizePath(artifact_root, winslash = "/"),
        implementation_hash,
        complete = FALSE,
        evidence_policy = "allow"
    )
    expect_identical(
        inventory$relative_path,
        imputefinder:::.ASSOCIATION_STUDY_MANIFEST_ARTIFACT
    )
    expect_error(
        imputefinder:::.association_protocol_artifact_inventory(
            normalizePath(artifact_root, winslash = "/"),
            implementation_hash,
            complete = TRUE,
            evidence_policy = "absent"
        ),
        class = "imputefinder_association_evidence_error"
    )

    extra_directory <- file.path(
        artifact_root,
        "inputs",
        "dda_null_balanced"
    )
    dir.create(extra_directory, recursive = TRUE)
    saveRDS(
        TRUE,
        file.path(extra_directory, "replicate-033.rds")
    )
    expect_error(
        imputefinder:::.association_protocol_artifact_inventory(
            normalizePath(artifact_root, winslash = "/"),
            implementation_hash,
            complete = FALSE,
            evidence_policy = "allow"
        ),
        class = "imputefinder_association_evidence_error"
    )

    detached_root <- withr::local_tempdir(pattern = "m13a-detached-")
    detached_artifact_root <- file.path(
        detached_root,
        imputefinder:::.ASSOCIATION_STUDY_ARTIFACT_DIRECTORY
    )
    dir.create(detached_artifact_root, recursive = TRUE)
    saveRDS(
        list(manifest_hash = strrep("b", 64L)),
        file.path(
            detached_artifact_root,
            imputefinder:::.ASSOCIATION_STUDY_MANIFEST_ARTIFACT
        )
    )
    expect_error(
        imputefinder:::.new_association_protocol_study_resolver(
            locked,
            detached_root
        ),
        class = "imputefinder_association_evidence_error"
    )
})

test_that("oracle response and both projection scales replay from probability", {
    samples <- sprintf("sample_%02d", 1:12)
    metadata <- data.frame(
        condition = rep(c("A", "B"), each = 6L),
        row.names = samples,
        stringsAsFactors = FALSE
    )
    design <- missingness_design(metadata, condition = "condition")
    data <- matrix(
        1.0,
        nrow = 20L,
        ncol = length(samples),
        dimnames = list(sprintf("feature_%02d", 1:20), samples)
    )
    probability <- matrix(
        rep(c(0.4, 0.2), each = 6L),
        nrow = nrow(data),
        ncol = ncol(data),
        byrow = TRUE,
        dimnames = dimnames(data)
    )
    preparation <- imputefinder:::.new_association_preparation(data, design)
    oracle <- imputefinder:::.association_oracle_response(
        list(data = data, missing_probability = probability),
        preparation
    )

    expect_equal(unname(oracle$response[1:6]), rep(0.6, 6L))
    expect_equal(unname(oracle$response[7:12]), rep(0.8, 6L))
    expect_equal(
        imputefinder:::.association_target_truth(
            "a_fraction_ols_hc3_cr2",
            preparation,
            oracle,
            "condition[B]"
        ),
        0.2,
        tolerance = 1e-12
    )
    expect_equal(
        imputefinder:::.association_target_truth(
            "a_fraction_quasibinomial",
            preparation,
            oracle,
            "condition[B]"
        ),
        0.2,
        tolerance = 1e-8
    )

    nuisance_metadata <- data.frame(
        condition = rep(c("A", "B"), each = 6L),
        batch = rep(c("batch_1", "batch_2"), 6L),
        row.names = samples,
        stringsAsFactors = FALSE
    )
    nuisance_design <- missingness_design(
        nuisance_metadata,
        condition = "condition",
        nuisance = "batch"
    )
    nuisance_preparation <- imputefinder:::.new_association_preparation(
        data,
        nuisance_design
    )
    expected_detection <- 0.50 +
        0.20 * (nuisance_metadata$condition == "B") +
        0.10 * (nuisance_metadata$batch == "batch_2")
    nuisance_probability <- matrix(
        rep(1 - expected_detection, each = nrow(data)),
        nrow = nrow(data),
        dimnames = dimnames(data)
    )
    nuisance_oracle <- imputefinder:::.association_oracle_response(
        list(data = data, missing_probability = nuisance_probability),
        nuisance_preparation
    )
    expect_equal(
        imputefinder:::.association_target_truth(
            "a_fraction_ols_hc3_cr2",
            nuisance_preparation,
            nuisance_oracle,
            "condition[B]"
        ),
        0.20,
        tolerance = 1e-12
    )
})

test_that("complete audits derive all screening gates and anchored ranking", {
    fixture <- association_metric_fixture()
    set.seed(7301)
    state <- .Random.seed
    kind <- RNGkind()
    metrics <- imputefinder:::.association_candidate_metrics(
        fixture$run,
        fixture$opportunity,
        fixture$null,
        fixture$target
    )

    expect_identical(.Random.seed, state)
    expect_identical(RNGkind(), kind)
    expect_identical(nrow(metrics$screening_metrics), 18L)
    expect_true(all(metrics$screening_metrics$status == "measured"))
    expect_true(all(metrics$screening_metrics$passed))
    expect_identical(nrow(metrics$full_gate_results), 15L)
    expect_true(all(metrics$full_gate_results$status == "measured"))
    expect_true(all(metrics$full_gate_results$passed))
    expect_identical(
        metrics$candidate_order,
        imputefinder:::.ASSOCIATION_CANDIDATES
    )
    expect_identical(
        metrics$selected_candidate,
        "a_fraction_ols_hc3_cr2"
    )
    expect_identical(metrics$state, "winner_locked")
    runtime_extremes <- metrics$ranking_metrics
    runtime_extremes$median_runtime <- c(1e6, 0, 0)
    expect_identical(
        imputefinder:::.association_rank_candidates(runtime_extremes),
        imputefinder:::.ASSOCIATION_CANDIDATES
    )
    boundary <- metrics$ranking_metrics
    boundary$scope_coverage <- c(0.88, 0.90, 0.50)
    boundary$power <- c(0.78, 0.80, 0.50)
    boundary$absolute_bias <- c(0.03, 0.01, 0.20)
    expect_identical(
        imputefinder:::.association_rank_candidates(boundary)[1:2],
        imputefinder:::.ASSOCIATION_CANDIDATES[1:2]
    )
    outside <- boundary
    outside$scope_coverage[[1L]] <- 0.879999999
    expect_identical(
        imputefinder:::.association_rank_candidates(outside)[1:2],
        imputefinder:::.ASSOCIATION_CANDIDATES[2:1]
    )

    null_upper <- metrics$screening_metrics$estimate[
        metrics$screening_metrics$candidate_id ==
            "a_fraction_ols_hc3_cr2" &
            metrics$screening_metrics$metric ==
                "family_false_flag_clopper_pearson_upper95"
    ]
    expect_identical(
        null_upper,
        imputefinder:::.association_clopper_pearson_upper(0L, 64L)
    )
    expect_true(length(unique(
        metrics$full_gate_results$evidence_hash[1:5]
    )) == 5L)
})

test_that("full gates retain heterogeneous numerators intervals and median bias", {
    fixture <- association_metric_fixture()
    candidate_id <- imputefinder:::.ASSOCIATION_CANDIDATES[[1L]]
    null <- fixture$null[
        fixture$null$candidate_id == candidate_id,
        ,
        drop = FALSE
    ]
    null$false_flag <- (null$acquisition == "DDA" & null$replicate <= 4L) |
        (null$acquisition == "DIA" & null$replicate <= 3L)
    target <- fixture$target[
        fixture$target$candidate_id == candidate_id,
        ,
        drop = FALSE
    ]
    null_gates <- imputefinder:::.association_candidate_full_gates(
        candidate_id,
        null,
        target
    )
    pooled <- null_gates$gate_id == "a_assoc_null_upper_v4"
    expect_identical(null_gates$numerator[pooled], 7)
    expect_identical(null_gates$denominator[pooled], 64)
    expect_identical(
        null_gates$estimate[pooled],
        imputefinder:::.association_clopper_pearson_upper(7L, 64L)
    )
    expect_identical(null_gates$lower[pooled], 0)
    expect_identical(null_gates$upper[pooled], null_gates$estimate[pooled])

    stratum <- null_gates$gate_id == "a_assoc_null_stratum_v4"
    stratum_interval <- imputefinder:::.association_exact_binomial_interval(
        4L,
        32L
    )
    expect_identical(null_gates$numerator[stratum], 4)
    expect_identical(null_gates$denominator[stratum], 32)
    expect_identical(null_gates$estimate[stratum], 4 / 32)
    expect_identical(null_gates$lower[stratum], stratum_interval[[1L]])
    expect_identical(null_gates$upper[stratum], stratum_interval[[2L]])
    expect_match(null_gates$note[stratum], "DDA", fixed = TRUE)

    endpoint_target <- target
    endpoint_target$conf_low[[1L]] <- endpoint_target$truth[[1L]]
    endpoint_target$conf_high[[2L]] <- endpoint_target$truth[[2L]]
    endpoint_gates <- imputefinder:::.association_candidate_full_gates(
        candidate_id,
        fixture$null[fixture$null$candidate_id == candidate_id, , drop = FALSE],
        endpoint_target
    )
    coverage <- endpoint_gates$gate_id == "a_assoc_interval_coverage_v4"
    expect_identical(endpoint_gates$numerator[coverage], 128)

    heterogeneous <- target
    missed <- heterogeneous$replicate <= 2L
    heterogeneous$conf_low[missed] <- 0.21
    heterogeneous$conf_high[missed] <- 0.30
    heterogeneous$adjusted_p[heterogeneous$replicate > 24L] <- 0.50
    heterogeneous$effect <- ifelse(
        heterogeneous$replicate <= 16L,
        0.22,
        0.24
    )
    gates <- imputefinder:::.association_candidate_full_gates(
        candidate_id,
        fixture$null[fixture$null$candidate_id == candidate_id, , drop = FALSE],
        heterogeneous
    )
    coverage <- gates$gate_id == "a_assoc_interval_coverage_v4"
    power <- gates$gate_id == "a_assoc_alternative_power_v4"
    bias <- gates$gate_id == "a_assoc_effect_bias_v4"
    expect_identical(gates$numerator[coverage], 120)
    expect_identical(gates$denominator[coverage], 128)
    expect_identical(gates$estimate[coverage], 120 / 128)
    expect_identical(gates$numerator[power], 96)
    expect_identical(gates$denominator[power], 128)
    expect_identical(gates$estimate[power], 96 / 128)
    expect_equal(gates$numerator[bias], 3.84, tolerance = 1e-14)
    expect_identical(gates$denominator[bias], 128)
    expect_equal(gates$estimate[bias], 0.03, tolerance = 1e-14)
    expect_identical(
        c(gates$lower[coverage], gates$upper[coverage]),
        c(0.890625, 0.9765625)
    )
    expect_identical(
        c(gates$lower[power], gates$upper[power]),
        c(0.671875, 0.8203125)
    )
    expect_equal(
        c(gates$lower[bias], gates$upper[bias]),
        c(0.02, 0.04),
        tolerance = 1e-14
    )
    expect_true(all(nzchar(gates$note)))
})

test_that("candidate-agnostic opportunities reject denominator drift", {
    fixture <- association_metric_fixture()
    drift <- fixture$opportunity[-1L, , drop = FALSE]
    expect_error(
        imputefinder:::.association_candidate_metrics(
            fixture$run,
            drift,
            fixture$null,
            fixture$target
        ),
        class = "imputefinder_association_evidence_error"
    )
})

test_that("target scales reject impossible intervals and probabilities", {
    fixture <- association_metric_fixture()
    corrupted <- fixture$target
    corrupted$conf_low[[1L]] <- 0.4
    corrupted$conf_high[[1L]] <- 0.3
    expect_error(
        imputefinder:::.association_candidate_metrics(
            fixture$run,
            fixture$opportunity,
            fixture$null,
            corrupted
        ),
        class = "imputefinder_association_evidence_error"
    )

    corrupted <- fixture$target
    corrupted$adjusted_p[[1L]] <- -0.1
    expect_error(
        imputefinder:::.association_candidate_metrics(
            fixture$run,
            fixture$opportunity,
            fixture$null,
            corrupted
        ),
        class = "imputefinder_association_evidence_error"
    )

    corrupted <- fixture$target
    corrupted$truth[[1L]] <- 0
    expect_error(
        imputefinder:::.association_candidate_metrics(
            fixture$run,
            fixture$opportunity,
            fixture$null,
            corrupted
        ),
        class = "imputefinder_association_evidence_error"
    )

    corrupted <- fixture$target
    quasi <- which(
        corrupted$candidate_id == "a_fraction_quasibinomial"
    )[[1L]]
    corrupted$effect[[quasi]] <- 1.01
    expect_error(
        imputefinder:::.association_candidate_metrics(
            fixture$run,
            fixture$opportunity,
            fixture$null,
            corrupted
        ),
        class = "imputefinder_association_evidence_error"
    )
})

test_that("zero-p null families cannot improve candidate selection", {
    fixture <- association_metric_fixture()
    selected <- fixture$null$candidate_id == "a_fraction_ols_hc3_cr2"
    fixture$null$available_p_count[which(selected)[[1L]]] <- 0L
    metrics <- imputefinder:::.association_candidate_metrics(
        fixture$run,
        fixture$opportunity,
        fixture$null,
        fixture$target
    )

    robust_null <- metrics$screening_metrics$candidate_id ==
        "a_fraction_ols_hc3_cr2" & metrics$screening_metrics$metric %in% c(
            "family_false_flag_clopper_pearson_upper95",
            "maximum_acquisition_false_flag_fraction"
        )
    expect_true(all(
        metrics$screening_metrics$status[robust_null] == "failed"
    ))
    expect_false("a_fraction_ols_hc3_cr2" %in%
        metrics$full_gate_results$candidate_id)
    expect_identical(
        metrics$selected_candidate,
        "a_fraction_freedman_lane"
    )
})

test_that("execution failures retain dispositions and reject globally", {
    samples <- sprintf("sample_%02d", 1:12)
    metadata <- data.frame(
        condition = rep(c("A", "B"), each = 6L),
        row.names = samples,
        stringsAsFactors = FALSE
    )
    data <- matrix(
        1.0,
        nrow = 20L,
        ncol = length(samples),
        dimnames = list(sprintf("feature_%02d", 1:20), samples)
    )
    preparation <- imputefinder:::.new_association_preparation(
        data,
        missingness_design(metadata, condition = "condition")
    )
    failure <- imputefinder:::.new_association_execution_failure(
        "association_numerical_failure",
        "Candidate engine raised an error under sealed execution."
    )
    record <- imputefinder:::.new_association_study_result(
        "a_fraction_ols_hc3_cr2",
        "dda_mixed_outlier",
        1L,
        strrep("a", 64L),
        failure,
        0.0
    )
    resolved <- imputefinder:::.association_resolve_run(
        record,
        preparation,
        NULL,
        "a_fraction_ols_hc3_cr2",
        "dda_mixed_outlier",
        1L,
        "DDA",
        NULL
    )
    expect_identical(
        nrow(resolved$outcomes),
        nrow(preparation$hypotheses)
    )
    expect_true(all(resolved$outcomes$status == "unavailable"))
    expect_true(all(
        resolved$outcomes$code == "association_numerical_failure"
    ))

    fixture <- association_metric_fixture()
    failed_run <- fixture$run$candidate_id ==
        "a_fraction_ols_hc3_cr2" &
        fixture$run$scenario_id == "dda_mixed_outlier" &
        fixture$run$replicate == 1L
    fixture$run$result_sha256[failed_run] <- NA_character_
    fixture$run$result_kind[failed_run] <- "execution_failure"
    fixture$run$status[failed_run] <- "failed"
    fixture$run$failure_code[failed_run] <-
        "association_numerical_failure"
    metrics <- imputefinder:::.association_candidate_metrics(
        fixture$run,
        fixture$opportunity,
        fixture$null,
        fixture$target
    )
    robust_thresholds <- metrics$screening_metrics$candidate_id ==
        "a_fraction_ols_hc3_cr2" &
        metrics$screening_metrics$metric != "common_opportunity_coverage"
    expect_true(all(
        metrics$screening_metrics$status[robust_thresholds] == "failed"
    ))
    expect_false("a_fraction_ols_hc3_cr2" %in%
        metrics$full_gate_results$candidate_id)
})

test_that("screening scope cannot hide a full-gate paired target miss", {
    fixture <- association_metric_fixture()
    paired <- fixture$target$candidate_id == "a_fraction_quasibinomial" &
        fixture$target$scenario_id == "dia_monotone_paired"
    fixture$target$outcome_status[paired] <- "unavailable"
    fixture$target$code[paired] <- "association_quasibinomial_paired_scope"
    fixture$target[paired, c(
        "effect", "conf_low", "conf_high", "adjusted_p"
    )] <- NA_real_
    fixture$opportunity$available[
        fixture$opportunity$candidate_id == "a_fraction_quasibinomial" &
            fixture$opportunity$scenario_id == "dia_monotone_paired"
    ] <- FALSE
    metrics <- imputefinder:::.association_candidate_metrics(
        fixture$run,
        fixture$opportunity,
        fixture$null,
        fixture$target
    )

    quasi_screen <- metrics$screening_metrics$candidate_id ==
        "a_fraction_quasibinomial" &
        metrics$screening_metrics$metric != "common_opportunity_coverage"
    expect_true(all(metrics$screening_metrics$passed[quasi_screen]))
    quasi_alternative <- metrics$full_gate_results$candidate_id ==
        "a_fraction_quasibinomial" &
        metrics$full_gate_results$gate_id %in% c(
            "a_assoc_interval_coverage_v4",
            "a_assoc_alternative_power_v4",
            "a_assoc_effect_bias_v4"
        )
    expect_true(all(
        metrics$full_gate_results$status[quasi_alternative] == "failed"
    ))
    expect_false("a_fraction_quasibinomial" %in%
        metrics$ranking_metrics$candidate_id)
})

test_that("entirely unavailable evidence resolves no winner", {
    fixture <- association_metric_fixture()
    fixture$run$result_sha256[] <- NA_character_
    fixture$run$result_kind[] <- "unavailable"
    fixture$run$status[] <- "failed"
    fixture$run$failure_code[] <- "association_no_observable_features"
    fixture$opportunity$available[] <- FALSE
    fixture$null$available_p_count[] <- 0L
    fixture$target$outcome_status[] <- "unavailable"
    fixture$target$code[] <- "association_no_observable_features"
    fixture$target[c(
        "effect", "conf_low", "conf_high", "adjusted_p"
    )] <- NA_real_
    metrics <- imputefinder:::.association_candidate_metrics(
        fixture$run,
        fixture$opportunity,
        fixture$null,
        fixture$target
    )

    expect_identical(metrics$state, "no_winner")
    expect_true(is.na(metrics$selected_candidate))
    expect_length(metrics$candidate_order, 0L)
    expect_identical(nrow(metrics$full_gate_results), 0L)
    expect_identical(nrow(metrics$ranking_metrics), 0L)
})

test_that("complete resolver replay owns evidence and self-hash", {
    resolver <- association_absent_resolver()
    evidence <- imputefinder:::.resolve_association_candidate_evidence(
        resolver
    )

    expect_identical(
        class(evidence),
        "imputefinder_association_candidate_evidence_fixture"
    )
    expect_identical(
        evidence$schema,
        imputefinder:::.ASSOCIATION_CANDIDATE_EVIDENCE_FIXTURE_SCHEMA
    )
    expect_true(is.na(evidence$artifact_inventory_hash))
    expect_identical(
        names(evidence),
        imputefinder:::.ASSOCIATION_CANDIDATE_EVIDENCE_FIELDS
    )
    expect_identical(nrow(evidence$run_bindings), 1248L)
    expect_true(all(evidence$run_bindings$result_kind == "unavailable"))
    expect_identical(nrow(evidence$outcome_bindings), 1248L)
    expect_true(all(
        evidence$outcome_bindings$hypothesis ==
            "association_panel_abstention"
    ))
    expect_identical(evidence$state, "no_winner")
    expect_true(is.na(evidence$selected_candidate))
    expect_identical(
        evidence$evidence_hash,
        imputefinder:::.association_candidate_evidence_hash(evidence)
    )
    components <- imputefinder:::.association_candidate_evidence_components(
        evidence
    )
    expect_identical(
        components$scalar_fields$position,
        c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 18L, 19L)
    )
    reordered <- evidence
    reordered$run_bindings <- reordered$run_bindings[
        c(2L, 1L, seq.int(3L, nrow(reordered$run_bindings))),
        ,
        drop = FALSE
    ]
    expect_false(identical(
        evidence$evidence_hash,
        imputefinder:::.association_candidate_evidence_hash(reordered)
    ))
    expect_invisible(
        imputefinder:::.validate_association_candidate_evidence_replay(
            evidence,
            resolver,
            strrep("a", 64L)
        )
    )

    detached <- evidence
    detached$state <- "winner_locked"
    detached$selected_candidate <- "a_fraction_ols_hc3_cr2"
    detached$evidence_hash <-
        imputefinder:::.association_candidate_evidence_hash(detached)
    expect_error(
        imputefinder:::.validate_association_candidate_evidence_replay(
            detached,
            resolver,
            strrep("a", 64L)
        ),
        class = "imputefinder_association_evidence_error"
    )
})

test_that("Clopper-Pearson boundary arithmetic is exact and guarded", {
    expect_identical(
        imputefinder:::.association_clopper_pearson_upper(64L, 64L),
        1.0
    )
    expect_identical(
        imputefinder:::.association_clopper_pearson_upper(3L, 64L),
        0.11671658790335937
    )
    expect_error(
        imputefinder:::.association_clopper_pearson_upper(-1L, 64L),
        class = "imputefinder_association_evidence_error"
    )
    expect_error(
        imputefinder:::.association_clopper_pearson_upper(65L, 64L),
        class = "imputefinder_association_evidence_error"
    )
})

test_that("whole-replicate bootstrap is frozen and restores seed absence", {
    bootstrap_gates <- imputefinder:::.association_study_gate_registry()$gate_id[
        3:5
    ]
    seeds <- unname(unlist(lapply(
        imputefinder:::.ASSOCIATION_CANDIDATES,
        function(candidate_id) vapply(
            bootstrap_gates,
            function(gate_id) imputefinder:::.association_bootstrap_seed(
                candidate_id,
                gate_id
            ),
            integer(1L)
        )
    )))
    expect_identical(
        seeds,
        c(
            192773681L, 114382253L, 190856277L,
            81073899L, 220375734L, 217308341L,
            138631724L, 73949617L, 233448835L
        )
    )
    expect_identical(anyDuplicated(seeds), 0L)
    kind_before <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    seed_before <- if (had_seed) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit({
        do.call(RNGkind, as.list(kind_before))
        if (had_seed) {
            assign(".Random.seed", seed_before, envir = .GlobalEnv)
        } else if (exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
        )) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    if (had_seed) {
        rm(".Random.seed", envir = .GlobalEnv)
    }
    absent_kind <- RNGkind()
    values <- as.double(rep(c(rep(1, 29L), rep(0, 3L)), 4L))
    scenarios <- rep(
        imputefinder:::.ASSOCIATION_TARGETS$scenario_id,
        each = 32L
    )
    interval <- imputefinder:::.association_bootstrap_interval(
        values,
        scenarios,
        rep(1:32, times = 4L),
        "a_fraction_ols_hc3_cr2",
        "a_assoc_interval_coverage_v4",
        "mean"
    )
    expect_identical(interval, c(0.8515625, 0.953125))
    expect_identical(RNGkind(), absent_kind)
    expect_false(exists(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
    ))

    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(8841)
    present_kind <- RNGkind()
    present_seed <- .Random.seed
    stratified_values <- as.double(c(
        rep(0, 32L), rep(c(1, 0), c(8L, 24L)),
        rep(c(1, 0), c(24L, 8L)), rep(1, 32L)
    ))
    stratified_interval <- imputefinder:::.association_bootstrap_interval(
        stratified_values,
        scenarios,
        rep(1:32, times = 4L),
        "a_fraction_ols_hc3_cr2",
        "a_assoc_alternative_power_v4",
        "mean"
    )
    median_values <- as.double(c(
        rep(c(0.01, 0.05), c(20L, 12L)),
        rep(c(0.02, 0.06), c(16L, 16L)),
        rep(c(0.03, 0.07), c(12L, 20L)),
        rep(c(0.04, 0.08), c(8L, 24L))
    ))
    median_interval <- imputefinder:::.association_bootstrap_interval(
        median_values,
        scenarios,
        rep(1:32, times = 4L),
        "a_fraction_ols_hc3_cr2",
        "a_assoc_effect_bias_v4",
        "median"
    )
    expect_identical(stratified_interval, c(0.4453125, 0.5546875))
    expect_false(identical(
        stratified_interval,
        c(0.4140625, 0.5859375)
    ))
    expect_identical(median_interval, c(0.04, 0.06))
    expect_identical(RNGkind(), present_kind)
    expect_identical(.Random.seed, present_seed)
    expect_error(
        imputefinder:::.association_bootstrap_interval(
            stratified_values,
            rev(scenarios),
            rep(1:32, times = 4L),
            "a_fraction_ols_hc3_cr2",
            "a_assoc_alternative_power_v4",
            "mean"
        ),
        class = "imputefinder_association_evidence_error"
    )
    expect_identical(
        imputefinder:::.association_expand_bootstrap_interval(
            c(0.2, 0.8),
            0.1
        ),
        c(0.1, 0.8)
    )
    expect_identical(
        imputefinder:::.association_expand_bootstrap_interval(
            c(0.2, 0.8),
            0.9
        ),
        c(0.2, 0.9)
    )
})
