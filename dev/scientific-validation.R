#!/usr/bin/env Rscript

# M7 scientific-regression protocol + benchmark. M7a freezes simulation truth,
# metrics, gates, and the routine/long-run split before results are inspected.
# Run from the repository root:
# Rscript --vanilla dev/scientific-validation.R --verify
# Rscript --vanilla dev/scientific-validation.R --benchmark

.scientific_conditions <- c("A", "B")
.scientific_states <- c("complete", "MNAR", "MAR", "insufficient")

.scientific_scenario <- function(
    family,
    n_features,
    samples_per_condition,
    cutoff,
    mar_rate,
    automatic_gate = "required",
    mnar_width = 1.2,
    mnar_max_probability = 0.8,
    on_off_fraction = 0.04,
    mean_center_offset = 1,
    mean_sd = 2.2,
    condition_sd = 0.35,
    measurement_sd = 0.18
) {
    list(
        family = family,
        n_features = as.integer(n_features),
        samples_per_condition = as.integer(samples_per_condition),
        cutoff = as.numeric(cutoff),
        mar_rate = as.numeric(mar_rate),
        automatic_gate = automatic_gate,
        mnar_width = as.numeric(mnar_width),
        mnar_max_probability = as.numeric(mnar_max_probability),
        on_off_fraction = as.numeric(on_off_fraction),
        mean_center_offset = as.numeric(mean_center_offset),
        mean_sd = as.numeric(mean_sd),
        condition_sd = as.numeric(condition_sd),
        measurement_sd = as.numeric(measurement_sd)
    )
}

scientific_validation_catalog <- function() {
    group_grid <- expand.grid(
        samples_per_condition = c(4L, 8L, 20L),
        mar_rate = c(0.05, 0.25),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    group_scenarios <- lapply(
        seq_len(nrow(group_grid)),
        function(index) {
            row <- group_grid[index, , drop = FALSE]
            .scientific_scenario(
                family = "group_size_mar",
                n_features = 2000L,
                samples_per_condition = row$samples_per_condition,
                cutoff = 12,
                mar_rate = row$mar_rate,
                automatic_gate = if (row$mar_rate != 0.25) {
                    "required"
                } else if (row$samples_per_condition == 20L) {
                    "evidence_sensitivity"
                } else {
                    "stress"
                }
            )
        }
    )
    names(group_scenarios) <- sprintf(
        "group_n%02d_mar%02d",
        group_grid$samples_per_condition,
        as.integer(group_grid$mar_rate * 100)
    )

    sweep_cutoffs <- 8:14
    sweep_scenarios <- lapply(
        sweep_cutoffs,
        function(cutoff) {
            .scientific_scenario(
                family = "cutoff_sweep",
                n_features = 1600L,
                samples_per_condition = 8L,
                cutoff = cutoff,
                mar_rate = 0.05
            )
        }
    )
    names(sweep_scenarios) <- sprintf("sweep_c%02d", sweep_cutoffs)

    c(group_scenarios, sweep_scenarios)
}

scientific_validation_seeds <- function() {
    c(
        104677L, 130387L, 155939L, 181123L,
        206407L, 232019L, 257707L, 283397L
    )
}

scientific_routine_catalog <- function() {
    list(
        manual_on_off = .scientific_scenario(
            family = "routine_manual",
            n_features = 320L,
            samples_per_condition = 4L,
            cutoff = 12,
            mar_rate = 0.25,
            automatic_gate = "not_run",
            on_off_fraction = 0.05
        ),
        automatic_cliff = .scientific_scenario(
            family = "routine_automatic",
            n_features = 800L,
            samples_per_condition = 8L,
            cutoff = 12,
            mar_rate = 0.05,
            automatic_gate = "required",
            on_off_fraction = 0.05
        )
    )
}

scientific_validation_protocol <- function() {
    list(
        protocol_version = "1",
        conditions = .scientific_conditions,
        simulation_seeds = scientific_validation_seeds(),
        rescue_seeds = c(1L, 7L, 29L),
        catalog = scientific_validation_catalog(),
        routine = scientific_routine_catalog(),
        permutation_scenarios = c(
            "group_n04_mar25",
            "group_n08_mar05",
            "sweep_c08",
            "sweep_c14"
        ),
        metrics = c(
            "state_accuracy",
            "state_macro_f1",
            "mnar_precision",
            "mnar_recall",
            "mnar_f1",
            "mar_precision",
            "mar_recall",
            "mar_f1",
            "insufficient_f1",
            "retention_precision",
            "retention_recall",
            "retention_f1",
            "on_off_retention_recall",
            "on_off_seed_recall"
        ),
        gates = list(
            manual = data.frame(
                mar_rate = c(0.05, 0.25),
                q10_mnar_f1 = c(0.72, 0.60),
                q10_mar_f1 = c(0.72, 0.60),
                q10_state_macro_f1 = c(0.76, 0.66),
                q10_retention_f1 = c(0.94, 0.90),
                stringsAsFactors = FALSE
            ),
            on_off_retention_recall = 1,
            on_off_seed_recall = 1,
            automatic = list(
                required_success_rate = 0.875,
                stress_success_rate = 0.75,
                median_absolute_cutoff_error = 0.75,
                q90_absolute_cutoff_error = 1.25,
                minimum_median_signed_error = -0.35,
                maximum_median_signed_error = 0.75,
                q10_mnar_f1_delta_from_manual = -0.05,
                q10_retention_f1_delta_from_manual = -0.03,
                finite_cutoff = "strictly_inside_observed_mean_range"
            ),
            evidence_sensitivity = list(
                eligible_decision = "score_against_automatic_gates",
                ineligible_decision = "structured_failure",
                maximum_finite_cutoff_error = 1.25
            ),
            cutoff_sweep = list(
                state_and_retention = "exact",
                missing_masks = "exact",
                maximum_shift_error = 1e-10
            ),
            routine = list(
                manual = c(
                    mnar_f1 = 0.60,
                    mar_f1 = 0.60,
                    state_macro_f1 = 0.66,
                    retention_f1 = 0.90,
                    on_off_retention_recall = 1,
                    on_off_seed_recall = 1
                ),
                automatic = c(
                    mnar_f1 = 0.72,
                    mar_f1 = 0.72,
                    state_macro_f1 = 0.76,
                    retention_f1 = 0.94,
                    cutoff_error = 1,
                    mnar_f1_delta_from_manual = -0.05,
                    retention_f1_delta_from_manual = -0.03,
                    on_off_retention_recall = 1,
                    on_off_seed_recall = 1
                )
            ),
            permutation = "exact_named_result",
            rng_state = "exact",
            observed_values = "exact_except_logged_seeds"
        )
    )
}

.with_scientific_seed <- function(seed, code) {
    old_kind <- RNGkind()
    caller_has_seed <- exists(
        ".Random.seed",
        envir = globalenv(),
        inherits = FALSE
    )
    if (caller_has_seed) {
        caller_seed <- get(
            ".Random.seed",
            envir = globalenv(),
            inherits = FALSE
        )
    }

    on.exit(
        {
            do.call(RNGkind, as.list(old_kind))
            if (caller_has_seed) {
                assign(".Random.seed", caller_seed, envir = globalenv())
            } else if (exists(
                ".Random.seed",
                envir = globalenv(),
                inherits = FALSE
            )) {
                rm(".Random.seed", envir = globalenv())
            }
        },
        add = TRUE
    )

    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    force(code)
}

.normalise_scientific_scenario <- function(scenario, catalog) {
    if (is.character(scenario) && length(scenario) == 1L && !is.na(scenario)) {
        if (!scenario %in% names(catalog)) {
            stop("Unknown scientific-validation scenario.", call. = FALSE)
        }
        return(catalog[[scenario]])
    }
    if (!is.list(scenario)) {
        stop("`scenario` must be a catalog name or scenario list.", call. = FALSE)
    }
    scenario
}

simulate_scientific_scenario <- function(
    scenario,
    seed,
    catalog = scientific_validation_catalog()
) {
    scenario <- .normalise_scientific_scenario(scenario, catalog)
    .validate_scientific_scenario(scenario)
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        !is.finite(seed) || seed != trunc(seed) ||
        seed < -.Machine$integer.max || seed > .Machine$integer.max) {
        stop("`seed` must be one finite integer.", call. = FALSE)
    }

    .with_scientific_seed(as.integer(seed), {
        .generate_scientific_scenario(scenario, as.integer(seed))
    })
}

.validate_scientific_scenario <- function(scenario) {
    required <- c(
        "family", "n_features", "samples_per_condition", "cutoff",
        "mar_rate", "automatic_gate", "mnar_width",
        "mnar_max_probability", "on_off_fraction", "mean_center_offset",
        "mean_sd", "condition_sd", "measurement_sd"
    )
    if (!is.list(scenario) || !all(required %in% names(scenario))) {
        stop("Scientific-validation scenario is malformed.", call. = FALSE)
    }
    valid <- scenario$n_features >= 80L &&
        scenario$samples_per_condition >= 1L &&
        is.finite(scenario$cutoff) &&
        scenario$mar_rate >= 0 && scenario$mar_rate < 1 &&
        scenario$mnar_width > 0 &&
        scenario$mnar_max_probability > 0 &&
        scenario$mnar_max_probability <= 1 &&
        scenario$on_off_fraction > 0 && scenario$on_off_fraction < 0.5 &&
        scenario$mean_sd > 0 && scenario$condition_sd >= 0 &&
        scenario$measurement_sd >= 0 &&
        scenario$automatic_gate %in% c(
            "required", "stress", "evidence_sensitivity", "not_run"
        )
    if (!valid) {
        stop("Scientific-validation scenario has invalid values.", call. = FALSE)
    }
    invisible(scenario)
}

.generate_scientific_scenario <- function(scenario, seed) {
    n <- scenario$n_features
    samples <- scenario$samples_per_condition
    cutoff <- scenario$cutoff
    features <- sprintf("protein_%05d", seq_len(n))
    sample_names <- unlist(
        lapply(
            .scientific_conditions,
            function(condition) {
                sprintf("%s_%02d", condition, seq_len(samples))
            }
        ),
        use.names = FALSE
    )
    groups <- stats::setNames(
        rep(.scientific_conditions, each = samples),
        sample_names
    )

    probability <- (seq_len(n) - 0.5) / n
    relative_mean <- stats::qnorm(
        probability,
        mean = scenario$mean_center_offset,
        sd = scenario$mean_sd
    )
    relative_mean <- sample(relative_mean, n, replace = FALSE)
    condition_effect <- stats::rnorm(n, sd = scenario$condition_sd)
    condition_means <- cbind(
        A = cutoff + relative_mean + condition_effect / 2,
        B = cutoff + relative_mean - condition_effect / 2
    )

    on_off_count <- max(1L, as.integer(floor(n * scenario$on_off_fraction)))
    on_off_features <- sample(seq_len(n), 2L * on_off_count, replace = FALSE)
    off_a <- on_off_features[seq_len(on_off_count)]
    off_b <- on_off_features[on_off_count + seq_len(on_off_count)]
    condition_means[off_a, "A"] <- cutoff - 2.5
    condition_means[off_a, "B"] <- cutoff + 2.5
    condition_means[off_b, "A"] <- cutoff + 2.5
    condition_means[off_b, "B"] <- cutoff - 2.5

    complete_data <- do.call(
        cbind,
        lapply(
            .scientific_conditions,
            function(condition) {
                matrix(
                    rep(condition_means[, condition], samples),
                    nrow = n,
                    ncol = samples
                )
            }
        )
    )
    complete_data <- complete_data + matrix(
        stats::rnorm(length(complete_data), sd = scenario$measurement_sd),
        nrow = n
    )
    dimnames(complete_data) <- list(features, sample_names)

    mnar_probability <- scenario$mnar_max_probability * pmin(
        pmax((cutoff - complete_data) / scenario$mnar_width, 0),
        1
    )
    mnar_probability[off_a, groups == "A"] <- 1
    mnar_probability[off_b, groups == "B"] <- 1
    mar_mask <- matrix(
        stats::runif(length(complete_data)) < scenario$mar_rate,
        nrow = n,
        dimnames = dimnames(complete_data)
    )
    mnar_mask <- matrix(
        stats::runif(length(complete_data)) < mnar_probability,
        nrow = n,
        dimnames = dimnames(complete_data)
    )
    missing_mask <- mar_mask | mnar_mask
    data <- complete_data
    data[missing_mask] <- NA_real_

    feature_truth <- data.frame(
        feature = features,
        off_condition = NA_character_,
        on_condition = NA_character_,
        mean_A = condition_means[, "A"],
        mean_B = condition_means[, "B"],
        stringsAsFactors = FALSE
    )
    feature_truth$off_condition[off_a] <- "A"
    feature_truth$on_condition[off_a] <- "B"
    feature_truth$off_condition[off_b] <- "B"
    feature_truth$on_condition[off_b] <- "A"

    oracle <- .scientific_oracle(
        data = data,
        groups = groups,
        mnar_mask = mnar_mask,
        feature_truth = feature_truth
    )

    list(
        scenario = scenario,
        simulation_seed = seed,
        data = data,
        groups = groups,
        complete_data = complete_data,
        masks = list(
            MAR = mar_mask,
            MNAR = mnar_mask,
            missing = missing_mask
        ),
        mnar_probability = mnar_probability,
        feature_truth = feature_truth,
        oracle = oracle
    )
}

.scientific_oracle <- function(data, groups, mnar_mask, feature_truth) {
    features <- rownames(data)
    conditions <- sort(unique(unname(groups)), method = "radix")
    globally_absent <- rowSums(!is.na(data)) == 0L
    surviving_features <- features[!globally_absent]
    condition_count <- length(conditions)
    row_count <- length(surviving_features) * condition_count

    classifications <- data.frame(
        feature = rep(surviving_features, each = condition_count),
        condition = rep(conditions, times = length(surviving_features)),
        state = rep(NA_character_, row_count),
        stringsAsFactors = FALSE
    )
    for (condition_index in seq_along(conditions)) {
        condition <- conditions[[condition_index]]
        columns <- groups == condition
        block <- data[surviving_features, columns, drop = FALSE]
        causal_mnar <- rowSums(
            mnar_mask[surviving_features, columns, drop = FALSE]
        ) > 0L
        observed <- rowSums(!is.na(block))
        missing <- ncol(block) - observed
        state <- rep("insufficient", length(surviving_features))
        state[missing == 0L] <- "complete"
        state[missing > 0L & causal_mnar] <- "MNAR"
        mar_candidate <- missing > 0L & !causal_mnar
        state[mar_candidate & observed > ncol(block) / 2] <- "MAR"
        rows <- seq.int(condition_index, row_count, by = condition_count)
        classifications$state[rows] <- state
    }

    feature_status <- data.frame(
        feature = features,
        retained = FALSE,
        drop_reason = NA_character_,
        stringsAsFactors = FALSE
    )
    feature_status$drop_reason[globally_absent] <- "all_missing"
    if (length(surviving_features) > 0L) {
        states <- matrix(
            classifications$state,
            nrow = length(surviving_features),
            ncol = condition_count,
            byrow = TRUE,
            dimnames = list(surviving_features, conditions)
        )
        surviving_status <- match(surviving_features, features)
        for (feature_index in seq_along(surviving_features)) {
            feature_states <- states[feature_index, ]
            insufficient <- conditions[feature_states == "insufficient"]
            if (length(insufficient) > 0L) {
                feature_status$drop_reason[surviving_status[[feature_index]]] <-
                    paste0(
                        "insufficient:",
                        paste(sort(insufficient, method = "radix"), collapse = ",")
                    )
            } else if (all(feature_states == "MNAR")) {
                feature_status$drop_reason[surviving_status[[feature_index]]] <-
                    "MNAR_all_conditions"
            } else {
                feature_status$retained[surviving_status[[feature_index]]] <- TRUE
            }
        }
    }

    feature_truth$oracle_retained <- feature_status$retained[
        match(feature_truth$feature, feature_status$feature)
    ]
    list(
        classifications = classifications,
        feature_status = feature_status,
        feature_truth = feature_truth
    )
}

.scientific_binary_metrics <- function(truth, predicted) {
    true_positive <- sum(truth & predicted)
    false_positive <- sum(!truth & predicted)
    false_negative <- sum(truth & !predicted)
    precision <- if (true_positive + false_positive == 0L) {
        NA_real_
    } else {
        true_positive / (true_positive + false_positive)
    }
    recall <- if (true_positive + false_negative == 0L) {
        NA_real_
    } else {
        true_positive / (true_positive + false_negative)
    }
    f1 <- if (!is.finite(precision) || !is.finite(recall) ||
        precision + recall == 0) {
        NA_real_
    } else {
        2 * precision * recall / (precision + recall)
    }

    c(precision = precision, recall = recall, f1 = f1)
}

.scientific_state_f1 <- function(truth, predicted, state) {
    unname(.scientific_binary_metrics(truth == state, predicted == state)[["f1"]])
}

score_scientific_result <- function(result, simulation) {
    required <- c("classifications", "feature_status", "seed_log")
    if (!is.list(result) || !all(required %in% names(result))) {
        stop("Scientific score requires a complete result schema.", call. = FALSE)
    }
    oracle <- simulation$oracle
    truth_keys <- paste(
        oracle$classifications$feature,
        oracle$classifications$condition,
        sep = "\r"
    )
    predicted_keys <- paste(
        result$classifications$feature,
        result$classifications$condition,
        sep = "\r"
    )
    if (anyDuplicated(predicted_keys) || !setequal(predicted_keys, truth_keys)) {
        stop("Predicted classifications do not align with oracle blocks.", call. = FALSE)
    }
    predicted_state <- result$classifications$state[
        match(truth_keys, predicted_keys)
    ]
    truth_state <- oracle$classifications$state
    if (!is.character(predicted_state) || anyNA(predicted_state) ||
        any(!predicted_state %in% .scientific_states)) {
        stop("Predicted states are invalid.", call. = FALSE)
    }

    mnar <- .scientific_binary_metrics(
        truth_state == "MNAR",
        predicted_state == "MNAR"
    )
    mar <- .scientific_binary_metrics(
        truth_state == "MAR",
        predicted_state == "MAR"
    )
    state_f1 <- vapply(
        .scientific_states,
        function(state) .scientific_state_f1(truth_state, predicted_state, state),
        numeric(1L)
    )

    truth_status <- oracle$feature_status
    predicted_status <- result$feature_status
    if (anyDuplicated(predicted_status$feature) ||
        !setequal(predicted_status$feature, truth_status$feature)) {
        stop("Predicted feature status does not align with oracle features.", call. = FALSE)
    }
    predicted_retained <- predicted_status$retained[
        match(truth_status$feature, predicted_status$feature)
    ]
    if (!is.logical(predicted_retained) || anyNA(predicted_retained)) {
        stop("Predicted retention values are invalid.", call. = FALSE)
    }
    retention <- .scientific_binary_metrics(
        truth_status$retained,
        predicted_retained
    )

    feature_truth <- oracle$feature_truth
    on_off <- !is.na(feature_truth$off_condition)
    eligible_on_off <- on_off & feature_truth$oracle_retained
    on_off_retention <- if (!any(eligible_on_off)) {
        NA_real_
    } else {
        mean(predicted_retained[eligible_on_off])
    }
    not_globally_absent <- is.na(truth_status$drop_reason) |
        truth_status$drop_reason != "all_missing"
    expected_on_off_seed <- on_off & not_globally_absent
    expected_seed <- paste(
        feature_truth$feature[expected_on_off_seed],
        feature_truth$off_condition[expected_on_off_seed],
        sep = "\r"
    )
    actual_seed <- paste(result$seed_log$feature, result$seed_log$condition, sep = "\r")
    on_off_seed_recall <- if (length(expected_seed) == 0L) {
        NA_real_
    } else {
        mean(expected_seed %in% actual_seed)
    }

    c(
        state_accuracy = mean(predicted_state == truth_state),
        state_macro_f1 = mean(state_f1, na.rm = TRUE),
        mnar_precision = unname(mnar[["precision"]]),
        mnar_recall = unname(mnar[["recall"]]),
        mnar_f1 = unname(mnar[["f1"]]),
        mar_precision = unname(mar[["precision"]]),
        mar_recall = unname(mar[["recall"]]),
        mar_f1 = unname(mar[["f1"]]),
        insufficient_f1 = unname(state_f1[["insufficient"]]),
        retention_precision = unname(retention[["precision"]]),
        retention_recall = unname(retention[["recall"]]),
        retention_f1 = unname(retention[["f1"]]),
        on_off_retention_recall = on_off_retention,
        on_off_seed_recall = on_off_seed_recall
    )
}

scientific_validation_permutation <- function(simulation) {
    feature_count <- nrow(simulation$data)
    feature_order <- c(
        seq.int(2L, feature_count, by = 2L),
        seq.int(1L, feature_count, by = 2L)
    )
    conditions <- rev(sort(unique(unname(simulation$groups)), method = "radix"))
    sample_order <- unlist(
        lapply(
            conditions,
            function(condition) {
                rev(which(simulation$groups == condition))
            }
        ),
        use.names = FALSE
    )
    list(features = feature_order, samples = sample_order)
}

scientific_validation_manifest <- function(
    catalog = scientific_validation_catalog()
) {
    rows <- lapply(
        names(catalog),
        function(name) {
            scenario <- catalog[[name]]
            data.frame(
                scenario = name,
                family = scenario$family,
                features = scenario$n_features,
                samples_per_condition = scenario$samples_per_condition,
                mar_rate = scenario$mar_rate,
                cutoff = scenario$cutoff,
                on_off_features = 2L * max(
                    1L,
                    as.integer(floor(
                        scenario$n_features * scenario$on_off_fraction
                    ))
                ),
                automatic_gate = scenario$automatic_gate,
                stringsAsFactors = FALSE
            )
        }
    )
    manifest <- do.call(rbind, rows)
    rownames(manifest) <- NULL
    manifest
}

scientific_validation_protocol_hash <- function(
    protocol = scientific_validation_protocol()
) {
    path <- tempfile("imputefinder-scientific-protocol-", fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(protocol, path, version = 2, compress = FALSE)
    unname(tools::md5sum(path))
}

.validate_scientific_protocol <- function(protocol) {
    manifest <- scientific_validation_manifest(protocol$catalog)
    expected_group <- as.vector(outer(
        c(4L, 8L, 20L),
        c(5L, 25L),
        function(n, mar) sprintf("group_n%02d_mar%02d", n, mar)
    ))
    expected_sweep <- sprintf("sweep_c%02d", 8:14)
    score_simulation <- simulate_scientific_scenario(
        protocol$routine$manual_on_off,
        protocol$simulation_seeds[[1L]],
        catalog = protocol$routine
    )
    score_names <- names(score_scientific_result(
        .perfect_scientific_result(score_simulation),
        score_simulation
    ))
    stopifnot(
        identical(protocol$protocol_version, "1"),
        identical(protocol$conditions, c("A", "B")),
        length(protocol$simulation_seeds) == 8L,
        !anyDuplicated(protocol$simulation_seeds),
        identical(protocol$rescue_seeds, c(1L, 7L, 29L)),
        setequal(names(protocol$catalog), c(expected_group, expected_sweep)),
        identical(
            sort(unique(manifest$samples_per_condition)),
            c(4L, 8L, 20L)
        ),
        identical(sort(unique(manifest$mar_rate)), c(0.05, 0.25)),
        identical(
            sort(manifest$cutoff[manifest$family == "cutoff_sweep"]),
            as.numeric(8:14)
        ),
        all(manifest$on_off_features > 0L),
        identical(protocol$metrics, score_names)
    )
    invisible(manifest)
}

.perfect_scientific_result <- function(simulation) {
    feature_truth <- simulation$oracle$feature_truth
    status <- simulation$oracle$feature_status
    drop_reason <- status$drop_reason[
        match(feature_truth$feature, status$feature)
    ]
    eligible <- !is.na(feature_truth$off_condition) &
        (is.na(drop_reason) | drop_reason != "all_missing")
    seed_log <- data.frame(
        feature = feature_truth$feature[eligible],
        condition = feature_truth$off_condition[eligible],
        stringsAsFactors = FALSE
    )
    list(
        classifications = simulation$oracle$classifications,
        feature_status = simulation$oracle$feature_status,
        seed_log = seed_log
    )
}

verify_scientific_validation <- function() {
    caller_has_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
    if (caller_has_seed) {
        caller_seed <- get(".Random.seed", globalenv(), inherits = FALSE)
    }
    caller_kind <- RNGkind()
    on.exit(
        {
            do.call(RNGkind, as.list(caller_kind))
            if (caller_has_seed) {
                assign(".Random.seed", caller_seed, envir = globalenv())
            } else if (exists(".Random.seed", globalenv(), inherits = FALSE)) {
                rm(".Random.seed", envir = globalenv())
            }
        },
        add = TRUE
    )

    protocol <- scientific_validation_protocol()
    manifest <- .validate_scientific_protocol(protocol)
    seed <- protocol$simulation_seeds[[1L]]

    set.seed(991L)
    present_seed <- get(".Random.seed", globalenv(), inherits = FALSE)
    simulation <- simulate_scientific_scenario(
        protocol$routine$automatic_cliff,
        seed,
        catalog = protocol$routine
    )
    stopifnot(identical(
        get(".Random.seed", globalenv(), inherits = FALSE),
        present_seed
    ))
    repeated <- simulate_scientific_scenario(
        protocol$routine$automatic_cliff,
        seed,
        catalog = protocol$routine
    )
    stopifnot(identical(simulation, repeated))

    rm(".Random.seed", envir = globalenv())
    absent_seed_simulation <- simulate_scientific_scenario(
        protocol$routine$manual_on_off,
        seed,
        catalog = protocol$routine
    )
    stopifnot(!exists(".Random.seed", globalenv(), inherits = FALSE))

    stopifnot(
        identical(simulation$masks$missing, simulation$masks$MAR |
            simulation$masks$MNAR),
        identical(is.na(simulation$data), simulation$masks$missing),
        all(is.finite(simulation$complete_data)),
        identical(names(simulation$groups), colnames(simulation$data)),
        identical(rownames(simulation$data), simulation$feature_truth$feature),
        all(simulation$mnar_probability >= 0 &
            simulation$mnar_probability <= 1),
        all(simulation$mnar_probability[
            simulation$complete_data >= simulation$scenario$cutoff &
                simulation$mnar_probability < 1
        ] == 0)
    )

    on_off <- !is.na(simulation$feature_truth$off_condition)
    for (index in which(on_off)) {
        feature <- simulation$feature_truth$feature[[index]]
        condition <- simulation$feature_truth$off_condition[[index]]
        stopifnot(all(simulation$masks$MNAR[
            feature,
            simulation$groups == condition
        ]))
    }

    perfect <- score_scientific_result(
        .perfect_scientific_result(absent_seed_simulation),
        absent_seed_simulation
    )
    stopifnot(
        identical(names(perfect), protocol$metrics),
        all(perfect[is.finite(perfect)] == 1)
    )

    sweep_8 <- simulate_scientific_scenario("sweep_c08", seed)
    sweep_14 <- simulate_scientific_scenario("sweep_c14", seed)
    stopifnot(
        identical(sweep_8$masks, sweep_14$masks),
        isTRUE(all.equal(
            sweep_8$complete_data + 6,
            sweep_14$complete_data,
            tolerance = 1e-14,
            check.attributes = TRUE
        )),
        identical(
            sweep_8$oracle$classifications$state,
            sweep_14$oracle$classifications$state
        ),
        identical(
            sweep_8$oracle$feature_status,
            sweep_14$oracle$feature_status
        )
    )

    permutation <- scientific_validation_permutation(simulation)
    stopifnot(
        identical(sort(permutation$features), seq_len(nrow(simulation$data))),
        identical(sort(permutation$samples), seq_len(ncol(simulation$data))),
        identical(
            unique(unname(simulation$groups[permutation$samples])),
            c("B", "A")
        )
    )

    list(
        manifest = manifest,
        routine = scientific_validation_manifest(protocol$routine),
        hash = scientific_validation_protocol_hash(protocol)
    )
}

.scientific_public_classifier <- function() {
    if (!requireNamespace("imputefinder", quietly = TRUE)) {
        stop(
            paste0(
                "M7b requires an installed imputefinder package. Install the ",
                "current source tree into the active library first."
            ),
            call. = FALSE
        )
    }

    getExportedValue("imputefinder", "classify_missingness")
}

.scientific_rng_snapshot <- function() {
    list(
        kind = RNGkind(),
        seed = if (exists(".Random.seed", globalenv(), inherits = FALSE)) {
            get(".Random.seed", globalenv(), inherits = FALSE)
        } else {
            NULL
        }
    )
}

.restore_scientific_rng <- function(snapshot) {
    do.call(RNGkind, as.list(snapshot$kind))
    if (is.null(snapshot$seed)) {
        if (exists(".Random.seed", globalenv(), inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
        }
    } else {
        assign(".Random.seed", snapshot$seed, envir = globalenv())
    }
    invisible(snapshot)
}

.scientific_output_audit <- function(result, simulation, rescue_seed) {
    fail <- function(...) {
        list(ok = FALSE, detail = paste0(...))
    }
    required <- c("data", "feature_status", "seed_log", "groups_by_sample")
    if (!is.list(result) || !all(required %in% names(result)) ||
        !is.matrix(result$data) || !is.data.frame(result$feature_status) ||
        !is.data.frame(result$seed_log)) {
        return(fail("result schema is incomplete"))
    }

    input <- simulation$data
    output <- result$data
    status <- result$feature_status
    retained <- status$feature[status$retained]
    if (!identical(rownames(output), retained) ||
        !identical(colnames(output), colnames(input)) ||
        !identical(
            result$groups_by_sample,
            stats::setNames(
                as.character(simulation$groups),
                names(simulation$groups)
            )
        )) {
        return(fail("output axes, retained order, or aligned groups changed"))
    }

    seed_log <- result$seed_log
    seed_columns <- c(
        "feature", "condition", "sample", "old_value", "inserted_value", "seed"
    )
    if (!identical(names(seed_log), seed_columns)) {
        return(fail("seed log schema changed"))
    }
    seed_keys <- paste(seed_log$feature, seed_log$sample, sep = "\r")
    if (anyDuplicated(seed_keys)) {
        return(fail("seed log contains duplicate feature-sample cells"))
    }

    conditions <- sort(unique(unname(simulation$groups)), method = "radix")
    condition_minima <- stats::setNames(
        vapply(
            conditions,
            function(condition) {
                min(input[, simulation$groups == condition, drop = FALSE], na.rm = TRUE)
            },
            numeric(1L)
        ),
        conditions
    )
    if (nrow(seed_log) > 0L) {
        feature_index <- match(seed_log$feature, rownames(input))
        sample_index <- match(seed_log$sample, colnames(input))
        valid_reference <- !anyNA(feature_index) && !anyNA(sample_index)
        original_values <- if (valid_reference) {
            input[cbind(feature_index, sample_index)]
        } else {
            rep(NA_real_, nrow(seed_log))
        }
        expected_condition <- if (valid_reference) {
            as.character(unname(simulation$groups[seed_log$sample]))
        } else {
            rep(NA_character_, nrow(seed_log))
        }
        expected_minimum <- unname(condition_minima[seed_log$condition])
        valid_log <- valid_reference &&
            all(is.na(original_values)) &&
            all(is.na(seed_log$old_value)) &&
            identical(seed_log$condition, expected_condition) &&
            all(is.finite(seed_log$inserted_value)) &&
            identical(seed_log$inserted_value, expected_minimum) &&
            identical(seed_log$seed, rep(as.integer(rescue_seed), nrow(seed_log)))
        if (!valid_log) {
            return(fail("seed log does not match original cells or condition minima"))
        }
    }

    if (nrow(output) > 0L) {
        original <- input[rownames(output), colnames(output), drop = FALSE]
        originally_observed <- !is.na(original)
        if (!identical(
            unname(output[originally_observed]),
            unname(original[originally_observed])
        )) {
            return(fail("an original observed value changed"))
        }

        inserted <- is.na(original) & !is.na(output)
        output_keys <- as.vector(outer(
            rownames(output),
            colnames(output),
            paste,
            sep = "\r"
        ))
        inserted_keys <- output_keys[as.vector(inserted)]
        retained_seed <- seed_log$feature %in% rownames(output)
        if (!setequal(inserted_keys, seed_keys[retained_seed])) {
            return(fail("returned inserted cells and seed log disagree"))
        }
        if (any(inserted)) {
            inserted_values <- stats::setNames(
                as.vector(output)[as.vector(inserted)],
                inserted_keys
            )
            logged_values <- stats::setNames(
                seed_log$inserted_value[retained_seed],
                seed_keys[retained_seed]
            )
            if (!identical(
                unname(inserted_values[names(logged_values)]),
                unname(logged_values)
            )) {
                return(fail("returned seed values differ from seed provenance"))
            }
        }
    }

    list(ok = TRUE, detail = "")
}

.invoke_scientific_classifier <- function(
    classifier,
    simulation,
    cutoffs,
    rescue_seed
) {
    original <- simulation$data
    before_rng <- .scientific_rng_snapshot()
    value <- tryCatch(
        classifier(
            x = simulation$data,
            group = simulation$groups,
            cutoffs = cutoffs,
            seed = rescue_seed
        ),
        error = identity
    )
    after_rng <- .scientific_rng_snapshot()
    input_unchanged <- identical(simulation$data, original)
    rng_unchanged <- identical(after_rng, before_rng)
    if (!rng_unchanged) {
        .restore_scientific_rng(before_rng)
    }
    output_audit <- if (inherits(value, "error")) {
        list(ok = inherits(value, "imputefinder_cutoff_error"), detail = "")
    } else {
        .scientific_output_audit(value, simulation, rescue_seed)
    }

    list(
        value = value,
        audit = list(
            input_unchanged = input_unchanged,
            rng_unchanged = rng_unchanged,
            output_valid = output_audit$ok,
            detail = output_audit$detail
        )
    )
}

.scientific_score_condition <- function(result, simulation, condition) {
    selected_result <- result
    selected_result$classifications <- result$classifications[
        result$classifications$condition == condition,
        ,
        drop = FALSE
    ]
    selected_simulation <- simulation
    selected_simulation$oracle$classifications <-
        simulation$oracle$classifications[
            simulation$oracle$classifications$condition == condition,
            ,
            drop = FALSE
        ]

    score_scientific_result(selected_result, selected_simulation)
}

.canonical_scientific_frame <- function(x, columns) {
    if (nrow(x) == 0L) {
        rownames(x) <- NULL
        return(x)
    }
    order_arguments <- c(
        lapply(columns, function(column) x[[column]]),
        list(method = "radix")
    )
    ordered <- x[do.call(order, order_arguments), , drop = FALSE]
    rownames(ordered) <- NULL
    ordered
}

.canonical_scientific_profile <- function(profile) {
    list(
        raw = .canonical_scientific_frame(profile$raw, c("feature", "condition")),
        grid = profile$grid,
        metadata = profile$metadata
    )
}

.canonical_scientific_groups <- function(groups) {
    lapply(
        groups[sort(names(groups), method = "radix")],
        function(condition_groups) {
            lapply(
                condition_groups,
                sort,
                method = "radix"
            )
        }
    )
}

.canonical_scientific_result <- function(
    result,
    include_data = FALSE,
    include_seed_assignment = FALSE
) {
    profiles <- lapply(
        result$profiles[sort(names(result$profiles), method = "radix")],
        .canonical_scientific_profile
    )
    seed_columns <- if (include_seed_assignment) {
        names(result$seed_log)
    } else {
        setdiff(names(result$seed_log), c("sample", "seed"))
    }
    seed_log <- result$seed_log[, seed_columns, drop = FALSE]
    seed_log <- .canonical_scientific_frame(seed_log, c("feature", "condition"))
    canonical <- list(
        classifications = .canonical_scientific_frame(
            result$classifications,
            c("feature", "condition")
        ),
        groups = .canonical_scientific_groups(result$groups),
        feature_status = .canonical_scientific_frame(
            result$feature_status,
            "feature"
        ),
        cutoffs = result$cutoffs[sort(names(result$cutoffs), method = "radix")],
        cutoff_diagnostics = result$cutoff_diagnostics[
            sort(names(result$cutoff_diagnostics), method = "radix")
        ],
        profiles = profiles,
        seed_log = seed_log,
        groups_by_sample = result$groups_by_sample[
            sort(names(result$groups_by_sample), method = "radix")
        ]
    )
    if (include_data) {
        canonical$data <- result$data[
            sort(rownames(result$data), method = "radix"),
            sort(colnames(result$data), method = "radix"),
            drop = FALSE
        ]
    }
    canonical
}

.canonical_scientific_failure <- function(error) {
    list(
        class = class(error),
        condition = error$condition,
        reason = error$reason,
        diagnostic = error$diagnostic,
        profile = .canonical_scientific_profile(error$profile)
    )
}

.canonical_scientific_outcome <- function(
    outcome,
    include_data = FALSE,
    include_seed_assignment = FALSE
) {
    value <- outcome$value
    if (inherits(value, "imputefinder_cutoff_error")) {
        return(list(status = "unidentifiable", value =
            .canonical_scientific_failure(value)))
    }
    if (inherits(value, "error")) {
        return(list(
            status = "error",
            value = list(class = class(value), message = conditionMessage(value))
        ))
    }
    list(
        status = "ok",
        value = .canonical_scientific_result(
            value,
            include_data = include_data,
            include_seed_assignment = include_seed_assignment
        )
    )
}

.scientific_audit_row <- function(
    scenario,
    simulation_seed,
    path,
    rescue_seed,
    invocation
) {
    data.frame(
        scenario = scenario,
        simulation_seed = as.integer(simulation_seed),
        path = path,
        rescue_seed = as.integer(rescue_seed),
        input_unchanged = invocation$audit$input_unchanged,
        rng_unchanged = invocation$audit$rng_unchanged,
        output_valid = invocation$audit$output_valid,
        detail = invocation$audit$detail,
        stringsAsFactors = FALSE
    )
}

.run_scientific_automatic <- function(
    classifier,
    simulation,
    rescue_seed
) {
    conditions <- sort(unique(unname(simulation$groups)), method = "radix")
    joint <- .invoke_scientific_classifier(
        classifier,
        simulation,
        cutoffs = NULL,
        rescue_seed = rescue_seed
    )
    invocations <- list(joint = joint)
    outcomes <- stats::setNames(vector("list", length(conditions)), conditions)
    if (!inherits(joint$value, "error")) {
        for (condition in conditions) {
            outcomes[[condition]] <- joint
        }
        return(list(outcomes = outcomes, invocations = invocations))
    }

    failed_condition <- if (
        inherits(joint$value, "imputefinder_cutoff_error") &&
            is.character(joint$value$condition) &&
            length(joint$value$condition) == 1L &&
            joint$value$condition %in% conditions
    ) {
        joint$value$condition
    } else {
        NA_character_
    }
    if (!is.na(failed_condition)) {
        outcomes[[failed_condition]] <- joint
    }

    unresolved <- conditions[vapply(outcomes, is.null, logical(1L))]
    for (condition in unresolved) {
        manual_conditions <- setdiff(conditions, condition)
        manual <- stats::setNames(
            rep(simulation$scenario$cutoff, length(manual_conditions)),
            manual_conditions
        )
        isolated <- .invoke_scientific_classifier(
            classifier,
            simulation,
            cutoffs = manual,
            rescue_seed = rescue_seed
        )
        invocations[[paste0("isolated_", condition)]] <- isolated
        outcomes[[condition]] <- isolated
    }

    list(outcomes = outcomes, invocations = invocations)
}

.scientific_score_values <- function(score, metric_names) {
    values <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)
    if (is.numeric(score) && all(metric_names %in% names(score))) {
        values[] <- score[metric_names]
    }
    values
}

.scientific_manual_row <- function(
    scenario_name,
    simulation,
    condition,
    result,
    metric_names
) {
    score <- if (inherits(result, "error")) {
        NULL
    } else {
        .scientific_score_condition(result, simulation, condition)
    }
    metrics <- .scientific_score_values(score, metric_names)
    data.frame(
        scenario = scenario_name,
        family = simulation$scenario$family,
        simulation_seed = simulation$simulation_seed,
        condition = condition,
        samples_per_condition = simulation$scenario$samples_per_condition,
        mar_rate = simulation$scenario$mar_rate,
        cutoff = simulation$scenario$cutoff,
        as.list(metrics),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

.scientific_automatic_evidence <- function(value, condition) {
    if (inherits(value, "imputefinder_cutoff_error")) {
        diagnostic <- value$diagnostic
        profile <- value$profile
    } else if (!inherits(value, "error")) {
        diagnostic <- value$cutoff_diagnostics[[condition]]
        profile <- value$profiles[[condition]]
    } else {
        return(list(
            diagnostic = NULL,
            profile = NULL,
            feature_count = NA_integer_,
            missing_count = NA_integer_,
            complete_count = NA_integer_,
            minimum_features = NA_integer_,
            minimum_class_features = NA_integer_,
            observed_range = c(minimum = NA_real_, maximum = NA_real_)
        ))
    }

    evidence <- diagnostic$quality$evidence
    feature_count <- evidence$feature_count
    class_counts <- evidence$class_counts
    observed_range <- diagnostic$quality$support$observed_mean_range
    if (is.null(observed_range) && is.list(profile)) {
        observed_range <- profile$metadata$observed_mean_range
    }
    list(
        diagnostic = diagnostic,
        profile = profile,
        feature_count = as.integer(feature_count),
        missing_count = as.integer(unname(class_counts[["missing"]])),
        complete_count = as.integer(unname(class_counts[["complete"]])),
        minimum_features = as.integer(evidence$minimum_feature_count),
        minimum_class_features =
            as.integer(evidence$minimum_class_feature_count),
        observed_range = stats::setNames(
            as.numeric(observed_range),
            c("minimum", "maximum")
        )
    )
}

.scientific_automatic_row <- function(
    scenario_name,
    simulation,
    condition,
    outcome,
    manual_result,
    metric_names
) {
    value <- outcome$value
    status <- if (inherits(value, "imputefinder_cutoff_error")) {
        "unidentifiable"
    } else if (inherits(value, "error")) {
        "error"
    } else {
        "ok"
    }
    score <- if (identical(status, "ok")) {
        .scientific_score_condition(value, simulation, condition)
    } else {
        NULL
    }
    metrics <- .scientific_score_values(score, metric_names)
    manual_score <- if (inherits(manual_result, "error")) {
        NULL
    } else {
        .scientific_score_condition(manual_result, simulation, condition)
    }
    manual_metrics <- .scientific_score_values(manual_score, metric_names)
    evidence <- .scientific_automatic_evidence(value, condition)
    cutoff <- if (identical(status, "ok")) {
        unname(value$cutoffs[[condition]])
    } else {
        NA_real_
    }
    signed_error <- cutoff - simulation$scenario$cutoff
    eligible <- is.finite(evidence$feature_count) &&
        is.finite(evidence$missing_count) &&
        is.finite(evidence$complete_count) &&
        is.finite(evidence$minimum_features) &&
        is.finite(evidence$minimum_class_features) &&
        evidence$feature_count >= evidence$minimum_features &&
        evidence$missing_count >= evidence$minimum_class_features &&
        evidence$complete_count >= evidence$minimum_class_features
    inside_support <- is.finite(cutoff) &&
        is.finite(evidence$observed_range[["minimum"]]) &&
        is.finite(evidence$observed_range[["maximum"]]) &&
        cutoff > evidence$observed_range[["minimum"]] &&
        cutoff < evidence$observed_range[["maximum"]]
    reason <- if (inherits(value, "imputefinder_cutoff_error")) {
        value$reason
    } else if (inherits(value, "error")) {
        conditionMessage(value)
    } else {
        ""
    }
    structured_condition_exact <- if (
        inherits(value, "imputefinder_cutoff_error")
    ) {
        identical(value$condition, condition)
    } else {
        TRUE
    }
    structured_profile_available <- if (
        inherits(value, "imputefinder_cutoff_error")
    ) {
        is.list(value$diagnostic) && is.list(value$profile) &&
            identical(value$diagnostic$profile, value$profile$metadata)
    } else {
        TRUE
    }

    data.frame(
        scenario = scenario_name,
        family = simulation$scenario$family,
        automatic_gate = simulation$scenario$automatic_gate,
        simulation_seed = simulation$simulation_seed,
        condition = condition,
        samples_per_condition = simulation$scenario$samples_per_condition,
        mar_rate = simulation$scenario$mar_rate,
        true_cutoff = simulation$scenario$cutoff,
        status = status,
        structured_condition_exact = structured_condition_exact,
        structured_profile_available = structured_profile_available,
        cutoff = cutoff,
        signed_cutoff_error = signed_error,
        absolute_cutoff_error = abs(signed_error),
        inside_observed_support = inside_support,
        evidence_eligible = eligible,
        feature_count = evidence$feature_count,
        missing_feature_count = evidence$missing_count,
        complete_feature_count = evidence$complete_count,
        reason = reason,
        as.list(metrics),
        mnar_f1_delta_from_manual = metrics[["mnar_f1"]] -
            manual_metrics[["mnar_f1"]],
        retention_f1_delta_from_manual = metrics[["retention_f1"]] -
            manual_metrics[["retention_f1"]],
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

.scientific_outcome_condition_signature <- function(outcome, condition) {
    value <- outcome$value
    if (inherits(value, "imputefinder_cutoff_error")) {
        return(.canonical_scientific_outcome(outcome))
    }
    if (inherits(value, "error")) {
        return(.canonical_scientific_outcome(outcome))
    }
    list(
        status = "ok",
        cutoff = unname(value$cutoffs[[condition]]),
        diagnostic = value$cutoff_diagnostics[[condition]],
        profile = .canonical_scientific_profile(value$profiles[[condition]]),
        classification = .canonical_scientific_frame(
            value$classifications[
                value$classifications$condition == condition,
                ,
                drop = FALSE
            ],
            c("feature", "condition")
        ),
        feature_status = .canonical_scientific_frame(
            value$feature_status,
            "feature"
        )
    )
}

.scientific_comparison_row <- function(
    audit,
    scenario,
    simulation_seed,
    condition = NA_character_,
    rescue_seed = NA_integer_,
    passed,
    detail = ""
) {
    data.frame(
        audit = audit,
        scenario = scenario,
        simulation_seed = as.integer(simulation_seed),
        condition = condition,
        rescue_seed = as.integer(rescue_seed),
        passed = isTRUE(passed),
        detail = detail,
        stringsAsFactors = FALSE
    )
}

.scientific_sweep_snapshot <- function(
    simulation,
    manual_result,
    automatic
) {
    manual <- if (inherits(manual_result, "error")) {
        list(status = "error", message = conditionMessage(manual_result))
    } else {
        list(
            status = "ok",
            classifications = .canonical_scientific_frame(
                manual_result$classifications[
                    , c("feature", "condition", "state", "retained", "drop_reason"),
                    drop = FALSE
                ],
                c("feature", "condition")
            ),
            feature_status = .canonical_scientific_frame(
                manual_result$feature_status,
                "feature"
            )
        )
    }
    conditions <- names(automatic$outcomes)
    automatic_snapshot <- lapply(
        conditions,
        function(condition) {
            value <- automatic$outcomes[[condition]]$value
            if (inherits(value, "imputefinder_cutoff_error")) {
                return(list(status = "unidentifiable"))
            }
            if (inherits(value, "error")) {
                return(list(status = "error", message = conditionMessage(value)))
            }
            list(
                status = "ok",
                cutoff_error = unname(value$cutoffs[[condition]]) -
                    simulation$scenario$cutoff,
                classifications = .canonical_scientific_frame(
                    value$classifications[
                        value$classifications$condition == condition,
                        c("feature", "condition", "state", "retained", "drop_reason"),
                        drop = FALSE
                    ],
                    c("feature", "condition")
                ),
                feature_status = .canonical_scientific_frame(
                    value$feature_status,
                    "feature"
                )
            )
        }
    )
    names(automatic_snapshot) <- conditions

    list(
        missing_mask = simulation$masks$missing,
        manual = manual,
        automatic = automatic_snapshot
    )
}

.compare_scientific_sweep <- function(
    reference,
    candidate,
    scenario,
    simulation_seed,
    tolerance
) {
    rows <- list(
        .scientific_comparison_row(
            "sweep_missing_mask",
            scenario,
            simulation_seed,
            passed = identical(reference$missing_mask, candidate$missing_mask)
        ),
        .scientific_comparison_row(
            "sweep_manual_state_retention",
            scenario,
            simulation_seed,
            passed = identical(reference$manual, candidate$manual)
        )
    )
    next_row <- length(rows) + 1L
    for (condition in names(reference$automatic)) {
        expected <- reference$automatic[[condition]]
        observed <- candidate$automatic[[condition]]
        status_exact <- identical(expected$status, observed$status)
        state_exact <- status_exact
        shift_exact <- status_exact
        if (status_exact && identical(expected$status, "ok")) {
            state_exact <- identical(
                expected[c("classifications", "feature_status")],
                observed[c("classifications", "feature_status")]
            )
            shift_exact <- isTRUE(
                abs(expected$cutoff_error - observed$cutoff_error) <= tolerance
            )
        }
        rows[[next_row]] <- .scientific_comparison_row(
            "sweep_automatic_status",
            scenario,
            simulation_seed,
            condition,
            passed = status_exact
        )
        rows[[next_row + 1L]] <- .scientific_comparison_row(
            "sweep_automatic_state_retention",
            scenario,
            simulation_seed,
            condition,
            passed = state_exact
        )
        rows[[next_row + 2L]] <- .scientific_comparison_row(
            "sweep_automatic_shift_error",
            scenario,
            simulation_seed,
            condition,
            passed = shift_exact,
            detail = if (
                status_exact && identical(expected$status, "ok")
            ) {
                format(
                    abs(expected$cutoff_error - observed$cutoff_error),
                    digits = 17L
                )
            } else {
                ""
            }
        )
        next_row <- next_row + 3L
    }
    do.call(rbind, rows)
}

.scientific_permutation_audit <- function(
    classifier,
    protocol,
    scenario_name
) {
    seed <- protocol$simulation_seeds[[1L]]
    rescue_seed <- protocol$rescue_seeds[[1L]]
    simulation <- simulate_scientific_scenario(
        scenario_name,
        seed,
        catalog = protocol$catalog
    )
    permutation <- scientific_validation_permutation(simulation)
    permuted <- simulation
    permuted$data <- simulation$data[
        permutation$features,
        permutation$samples,
        drop = FALSE
    ]
    permuted$groups <- stats::setNames(
        factor(
            unname(simulation$groups[permutation$samples]),
            levels = c("B", "A")
        ),
        colnames(permuted$data)
    )

    manual_cutoffs <- stats::setNames(
        rep(simulation$scenario$cutoff, length(protocol$conditions)),
        rev(protocol$conditions)
    )
    baseline_manual <- .invoke_scientific_classifier(
        classifier,
        simulation,
        manual_cutoffs,
        rescue_seed
    )
    permuted_manual <- .invoke_scientific_classifier(
        classifier,
        permuted,
        manual_cutoffs,
        rescue_seed
    )
    manual_exact <- !inherits(baseline_manual$value, "error") &&
        !inherits(permuted_manual$value, "error") &&
        identical(
            .canonical_scientific_result(
                baseline_manual$value,
                include_data = TRUE,
                include_seed_assignment = TRUE
            ),
            .canonical_scientific_result(
                permuted_manual$value,
                include_data = TRUE,
                include_seed_assignment = TRUE
            )
        )

    baseline_automatic <- .run_scientific_automatic(
        classifier,
        simulation,
        rescue_seed
    )
    permuted_automatic <- .run_scientific_automatic(
        classifier,
        permuted,
        rescue_seed
    )
    automatic_rows <- lapply(
        protocol$conditions,
        function(condition) {
            baseline <- .canonical_scientific_outcome(
                baseline_automatic$outcomes[[condition]],
                include_data = TRUE,
                include_seed_assignment = TRUE
            )
            candidate <- .canonical_scientific_outcome(
                permuted_automatic$outcomes[[condition]],
                include_data = TRUE,
                include_seed_assignment = TRUE
            )
            .scientific_comparison_row(
                "permutation_automatic_exact",
                scenario_name,
                seed,
                condition,
                rescue_seed,
                identical(baseline, candidate)
            )
        }
    )

    comparison <- do.call(
        rbind,
        c(
            list(.scientific_comparison_row(
                "permutation_manual_exact",
                scenario_name,
                seed,
                rescue_seed = rescue_seed,
                passed = manual_exact
            )),
            automatic_rows
        )
    )
    invocations <- c(
        list(
            baseline_manual = baseline_manual,
            permuted_manual = permuted_manual
        ),
        stats::setNames(
            baseline_automatic$invocations,
            paste0("baseline_automatic_", names(baseline_automatic$invocations))
        ),
        stats::setNames(
            permuted_automatic$invocations,
            paste0("permuted_automatic_", names(permuted_automatic$invocations))
        )
    )
    audit <- do.call(
        rbind,
        lapply(
            names(invocations),
            function(path) {
                .scientific_audit_row(
                    scenario_name,
                    seed,
                    paste0("permutation_", path),
                    rescue_seed,
                    invocations[[path]]
                )
            }
        )
    )

    list(comparison = comparison, invocations = audit)
}

.scientific_object_hash <- function(x) {
    path <- tempfile("imputefinder-scientific-result-", fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(x, path, version = 2, compress = FALSE)
    unname(tools::md5sum(path))
}

.scientific_record_automatic_audits <- function(
    automatic,
    scenario_name,
    simulation_seed,
    rescue_seed,
    prefix = "automatic"
) {
    do.call(
        rbind,
        lapply(
            names(automatic$invocations),
            function(path) {
                .scientific_audit_row(
                    scenario_name,
                    simulation_seed,
                    paste(prefix, path, sep = "_"),
                    rescue_seed,
                    automatic$invocations[[path]]
                )
            }
        )
    )
}

run_scientific_validation_benchmark <- function(progress = TRUE) {
    protocol <- scientific_validation_protocol()
    manifest <- .validate_scientific_protocol(protocol)
    protocol_hash <- scientific_validation_protocol_hash(protocol)
    frozen_hash <- "4011e381bba2d0d747e91d277a45de5e"
    if (!identical(protocol_hash, frozen_hash)) {
        stop(
            "M7b refuses to run because the frozen protocol-v1 hash changed.",
            call. = FALSE
        )
    }
    classifier <- .scientific_public_classifier()
    caller_rng <- .scientific_rng_snapshot()
    on.exit(.restore_scientific_rng(caller_rng), add = TRUE)
    set.seed(7102026L)
    started <- proc.time()[["elapsed"]]

    manual_rows <- automatic_rows <- invocation_rows <- list()
    comparison_rows <- list()
    manual_index <- automatic_index <- invocation_index <- comparison_index <- 1L
    sweep_reference <- list()

    for (scenario_name in names(protocol$catalog)) {
        if (progress) {
            message(sprintf("M7b benchmark: %s", scenario_name))
        }
        for (simulation_seed in protocol$simulation_seeds) {
            simulation <- simulate_scientific_scenario(
                scenario_name,
                simulation_seed,
                catalog = protocol$catalog
            )
            rescue_seed <- protocol$rescue_seeds[[1L]]
            manual_cutoffs <- stats::setNames(
                rep(simulation$scenario$cutoff, length(protocol$conditions)),
                protocol$conditions
            )
            manual <- .invoke_scientific_classifier(
                classifier,
                simulation,
                manual_cutoffs,
                rescue_seed
            )
            invocation_rows[[invocation_index]] <- .scientific_audit_row(
                scenario_name,
                simulation_seed,
                "manual",
                rescue_seed,
                manual
            )
            invocation_index <- invocation_index + 1L
            for (condition in protocol$conditions) {
                manual_rows[[manual_index]] <- .scientific_manual_row(
                    scenario_name,
                    simulation,
                    condition,
                    manual$value,
                    protocol$metrics
                )
                manual_index <- manual_index + 1L
            }

            automatic <- .run_scientific_automatic(
                classifier,
                simulation,
                rescue_seed
            )
            invocation_rows[[invocation_index]] <-
                .scientific_record_automatic_audits(
                    automatic,
                    scenario_name,
                    simulation_seed,
                    rescue_seed
                )
            invocation_index <- invocation_index + 1L
            for (condition in protocol$conditions) {
                automatic_rows[[automatic_index]] <- .scientific_automatic_row(
                    scenario_name,
                    simulation,
                    condition,
                    automatic$outcomes[[condition]],
                    manual$value,
                    protocol$metrics
                )
                automatic_index <- automatic_index + 1L
            }

            for (alternative_seed in protocol$rescue_seeds[-1L]) {
                alternative_manual <- .invoke_scientific_classifier(
                    classifier,
                    simulation,
                    manual_cutoffs,
                    alternative_seed
                )
                invocation_rows[[invocation_index]] <- .scientific_audit_row(
                    scenario_name,
                    simulation_seed,
                    "manual_alternative_seed",
                    alternative_seed,
                    alternative_manual
                )
                invocation_index <- invocation_index + 1L
                comparison_rows[[comparison_index]] <-
                    .scientific_comparison_row(
                        "rescue_seed_manual_invariant",
                        scenario_name,
                        simulation_seed,
                        rescue_seed = alternative_seed,
                        passed = identical(
                            .canonical_scientific_outcome(manual),
                            .canonical_scientific_outcome(alternative_manual)
                        )
                    )
                comparison_index <- comparison_index + 1L

                alternative_automatic <- .run_scientific_automatic(
                    classifier,
                    simulation,
                    alternative_seed
                )
                invocation_rows[[invocation_index]] <-
                    .scientific_record_automatic_audits(
                        alternative_automatic,
                        scenario_name,
                        simulation_seed,
                        alternative_seed,
                        prefix = "automatic_alternative_seed"
                    )
                invocation_index <- invocation_index + 1L
                for (condition in protocol$conditions) {
                    comparison_rows[[comparison_index]] <-
                        .scientific_comparison_row(
                            "rescue_seed_automatic_invariant",
                            scenario_name,
                            simulation_seed,
                            condition,
                            alternative_seed,
                            identical(
                                .scientific_outcome_condition_signature(
                                    automatic$outcomes[[condition]],
                                    condition
                                ),
                                .scientific_outcome_condition_signature(
                                    alternative_automatic$outcomes[[condition]],
                                    condition
                                )
                            )
                        )
                    comparison_index <- comparison_index + 1L
                }
            }

            if (identical(simulation$scenario$family, "cutoff_sweep")) {
                snapshot <- .scientific_sweep_snapshot(
                    simulation,
                    manual$value,
                    automatic
                )
                seed_key <- as.character(simulation_seed)
                if (identical(scenario_name, "sweep_c08")) {
                    sweep_reference[[seed_key]] <- snapshot
                } else {
                    comparison_rows[[comparison_index]] <-
                        .compare_scientific_sweep(
                            sweep_reference[[seed_key]],
                            snapshot,
                            scenario_name,
                            simulation_seed,
                            protocol$gates$cutoff_sweep$maximum_shift_error
                        )
                    comparison_index <- comparison_index + 1L
                }
            }
        }
    }

    for (scenario_name in protocol$permutation_scenarios) {
        if (progress) {
            message(sprintf("M7b permutation audit: %s", scenario_name))
        }
        permutation <- .scientific_permutation_audit(
            classifier,
            protocol,
            scenario_name
        )
        comparison_rows[[comparison_index]] <- permutation$comparison
        comparison_index <- comparison_index + 1L
        invocation_rows[[invocation_index]] <- permutation$invocations
        invocation_index <- invocation_index + 1L
    }

    raw <- list(
        manual = do.call(rbind, manual_rows),
        automatic = do.call(rbind, automatic_rows),
        invocations = do.call(rbind, invocation_rows),
        comparisons = do.call(rbind, comparison_rows)
    )
    rownames(raw$manual) <- rownames(raw$automatic) <- NULL
    rownames(raw$invocations) <- rownames(raw$comparisons) <- NULL
    summaries <- .summarise_scientific_benchmark(raw, protocol)
    gates <- .evaluate_scientific_benchmark(raw, summaries, protocol)
    deterministic <- list(
        protocol_hash = protocol_hash,
        manifest = manifest,
        raw = raw,
        summaries = summaries,
        gates = gates
    )
    context <- list(
        r = R.version.string,
        platform = R.version$platform,
        imputefinder_version = as.character(
            utils::packageVersion("imputefinder")
        ),
        package_library = dirname(find.package("imputefinder")),
        source_md5 = tools::md5sum(sort(list.files("R", full.names = TRUE))),
        harness_md5 = unname(tools::md5sum("dev/scientific-validation.R")),
        elapsed_seconds = unname(proc.time()[["elapsed"]] - started)
    )

    list(
        protocol_hash = protocol_hash,
        result_hash = .scientific_object_hash(deterministic),
        raw = raw,
        summaries = summaries,
        gates = gates,
        context = context
    )
}

.scientific_quantile <- function(x, probability) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    unname(stats::quantile(x, probability, names = FALSE, type = 7L))
}

.summarise_scientific_manual <- function(rows) {
    keys <- unique(rows[, c("scenario", "condition"), drop = FALSE])
    summaries <- lapply(
        seq_len(nrow(keys)),
        function(index) {
            key <- keys[index, , drop = FALSE]
            selected <- rows$scenario == key$scenario &
                rows$condition == key$condition
            x <- rows[selected, , drop = FALSE]
            data.frame(
                scenario = key$scenario,
                condition = key$condition,
                runs = nrow(x),
                mar_rate = x$mar_rate[[1L]],
                q10_mnar_f1 = .scientific_quantile(x$mnar_f1, 0.1),
                q10_mar_f1 = .scientific_quantile(x$mar_f1, 0.1),
                q10_state_macro_f1 =
                    .scientific_quantile(x$state_macro_f1, 0.1),
                q10_retention_f1 =
                    .scientific_quantile(x$retention_f1, 0.1),
                stringsAsFactors = FALSE
            )
        }
    )
    do.call(rbind, summaries)
}

.summarise_scientific_automatic <- function(rows) {
    keys <- unique(rows[, c("scenario", "condition"), drop = FALSE])
    summaries <- lapply(
        seq_len(nrow(keys)),
        function(index) {
            key <- keys[index, , drop = FALSE]
            selected <- rows$scenario == key$scenario &
                rows$condition == key$condition
            x <- rows[selected, , drop = FALSE]
            success <- x$status == "ok"
            data.frame(
                scenario = key$scenario,
                condition = key$condition,
                automatic_gate = x$automatic_gate[[1L]],
                runs = nrow(x),
                eligible_runs = sum(x$evidence_eligible),
                successes = sum(success),
                success_rate = mean(success),
                median_absolute_cutoff_error = .scientific_quantile(
                    x$absolute_cutoff_error,
                    0.5
                ),
                q90_absolute_cutoff_error = .scientific_quantile(
                    x$absolute_cutoff_error,
                    0.9
                ),
                median_signed_cutoff_error = .scientific_quantile(
                    x$signed_cutoff_error,
                    0.5
                ),
                q10_mnar_f1_delta = .scientific_quantile(
                    x$mnar_f1_delta_from_manual,
                    0.1
                ),
                q10_retention_f1_delta = .scientific_quantile(
                    x$retention_f1_delta_from_manual,
                    0.1
                ),
                minimum_missing_features = min(x$missing_feature_count),
                minimum_complete_features = min(x$complete_feature_count),
                stringsAsFactors = FALSE
            )
        }
    )
    do.call(rbind, summaries)
}

.summarise_scientific_benchmark <- function(raw, protocol) {
    comparison <- aggregate(
        raw$comparisons$passed,
        list(audit = raw$comparisons$audit),
        function(x) c(passed = sum(x), total = length(x))
    )
    comparison <- data.frame(
        audit = comparison$audit,
        passed = comparison$x[, "passed"],
        total = comparison$x[, "total"],
        stringsAsFactors = FALSE
    )
    invocation <- data.frame(
        audit = c("input_unchanged", "rng_unchanged", "output_valid"),
        passed = c(
            sum(raw$invocations$input_unchanged),
            sum(raw$invocations$rng_unchanged),
            sum(raw$invocations$output_valid)
        ),
        total = rep(nrow(raw$invocations), 3L),
        stringsAsFactors = FALSE
    )

    list(
        manual = .summarise_scientific_manual(raw$manual),
        automatic = .summarise_scientific_automatic(raw$automatic),
        audits = rbind(invocation, comparison),
        automatic_failures = {
            failures <- as.data.frame(table(
                scenario = raw$automatic$scenario[
                    raw$automatic$status != "ok"
                ],
                condition = raw$automatic$condition[
                    raw$automatic$status != "ok"
                ],
                reason = raw$automatic$reason[
                    raw$automatic$status != "ok"
                ]
            ), stringsAsFactors = FALSE)
            failures[failures$Freq > 0L, , drop = FALSE]
        }
    )
}

.scientific_gate_row <- function(
    gate,
    scope,
    passed,
    observed,
    required,
    detail = ""
) {
    data.frame(
        gate = gate,
        scope = scope,
        passed = isTRUE(passed),
        observed = as.character(observed),
        required = as.character(required),
        detail = detail,
        stringsAsFactors = FALSE
    )
}

.scientific_format_number <- function(x) {
    if (length(x) != 1L || !is.finite(x)) {
        return("NA")
    }
    format(x, digits = 6L, trim = TRUE, scientific = FALSE)
}

.scientific_automatic_metric_gates <- function(
    summary,
    specification,
    success_rate
) {
    scope <- paste(summary$scenario, summary$condition, sep = ":")
    values <- list(
        .scientific_gate_row(
            "automatic_success_rate",
            scope,
            summary$success_rate >= success_rate,
            .scientific_format_number(summary$success_rate),
            paste0(">=", .scientific_format_number(success_rate))
        ),
        .scientific_gate_row(
            "automatic_median_absolute_cutoff_error",
            scope,
            summary$median_absolute_cutoff_error <=
                specification$median_absolute_cutoff_error,
            .scientific_format_number(summary$median_absolute_cutoff_error),
            paste0(
                "<=",
                .scientific_format_number(
                    specification$median_absolute_cutoff_error
                )
            )
        ),
        .scientific_gate_row(
            "automatic_q90_absolute_cutoff_error",
            scope,
            summary$q90_absolute_cutoff_error <=
                specification$q90_absolute_cutoff_error,
            .scientific_format_number(summary$q90_absolute_cutoff_error),
            paste0(
                "<=",
                .scientific_format_number(
                    specification$q90_absolute_cutoff_error
                )
            )
        ),
        .scientific_gate_row(
            "automatic_median_signed_cutoff_error",
            scope,
            summary$median_signed_cutoff_error >=
                specification$minimum_median_signed_error &&
                summary$median_signed_cutoff_error <=
                    specification$maximum_median_signed_error,
            .scientific_format_number(summary$median_signed_cutoff_error),
            sprintf(
                "[%s,%s]",
                .scientific_format_number(
                    specification$minimum_median_signed_error
                ),
                .scientific_format_number(
                    specification$maximum_median_signed_error
                )
            )
        ),
        .scientific_gate_row(
            "automatic_q10_mnar_f1_delta",
            scope,
            summary$q10_mnar_f1_delta >=
                specification$q10_mnar_f1_delta_from_manual,
            .scientific_format_number(summary$q10_mnar_f1_delta),
            paste0(
                ">=",
                .scientific_format_number(
                    specification$q10_mnar_f1_delta_from_manual
                )
            )
        ),
        .scientific_gate_row(
            "automatic_q10_retention_f1_delta",
            scope,
            summary$q10_retention_f1_delta >=
                specification$q10_retention_f1_delta_from_manual,
            .scientific_format_number(summary$q10_retention_f1_delta),
            paste0(
                ">=",
                .scientific_format_number(
                    specification$q10_retention_f1_delta_from_manual
                )
            )
        )
    )
    do.call(rbind, values)
}

.evaluate_scientific_benchmark <- function(raw, summaries, protocol) {
    gates <- list()
    next_gate <- 1L
    manual_specification <- protocol$gates$manual
    for (index in seq_len(nrow(summaries$manual))) {
        summary <- summaries$manual[index, , drop = FALSE]
        specification <- manual_specification[
            manual_specification$mar_rate == summary$mar_rate,
            ,
            drop = FALSE
        ]
        scope <- paste(summary$scenario, summary$condition, sep = ":")
        for (metric in c(
            "q10_mnar_f1", "q10_mar_f1", "q10_state_macro_f1",
            "q10_retention_f1"
        )) {
            observed <- summary[[metric]]
            required <- specification[[metric]]
            gates[[next_gate]] <- .scientific_gate_row(
                paste0("manual_", metric),
                scope,
                observed >= required,
                .scientific_format_number(observed),
                paste0(">=", .scientific_format_number(required))
            )
            next_gate <- next_gate + 1L
        }
    }

    for (metric in c("on_off_retention_recall", "on_off_seed_recall")) {
        values <- raw$manual[[metric]]
        failed <- paste(
            paste(
                raw$manual$scenario[values != 1],
                raw$manual$simulation_seed[values != 1],
                raw$manual$condition[values != 1],
                sep = ":"
            ),
            collapse = ","
        )
        gates[[next_gate]] <- .scientific_gate_row(
            paste0("manual_individual_", metric),
            "all",
            length(values) > 0L && all(values == 1),
            .scientific_format_number(min(values)),
            "=1",
            failed
        )
        next_gate <- next_gate + 1L
    }

    automatic_specification <- protocol$gates$automatic
    for (index in seq_len(nrow(summaries$automatic))) {
        summary <- summaries$automatic[index, , drop = FALSE]
        gate_type <- summary$automatic_gate
        if (identical(gate_type, "evidence_sensitivity")) {
            rows <- raw$automatic$scenario == summary$scenario &
                raw$automatic$condition == summary$condition
            selected <- raw$automatic[rows, , drop = FALSE]
            ineligible <- !selected$evidence_eligible
            ineligible_pass <- all(
                selected$status[ineligible] == "unidentifiable"
            )
            gates[[next_gate]] <- .scientific_gate_row(
                "evidence_ineligible_structured_failure",
                paste(summary$scenario, summary$condition, sep = ":"),
                ineligible_pass,
                sprintf(
                    "%d/%d",
                    sum(selected$status[ineligible] == "unidentifiable"),
                    sum(ineligible)
                ),
                "all"
            )
            next_gate <- next_gate + 1L
            finite <- is.finite(selected$absolute_cutoff_error)
            maximum_error <- if (any(finite)) {
                max(selected$absolute_cutoff_error[finite])
            } else {
                NA_real_
            }
            gates[[next_gate]] <- .scientific_gate_row(
                "evidence_finite_cutoff_error",
                paste(summary$scenario, summary$condition, sep = ":"),
                !any(finite) || maximum_error <=
                    protocol$gates$evidence_sensitivity$
                        maximum_finite_cutoff_error,
                .scientific_format_number(maximum_error),
                paste0(
                    "<=",
                    .scientific_format_number(
                        protocol$gates$evidence_sensitivity$
                            maximum_finite_cutoff_error
                    )
                )
            )
            next_gate <- next_gate + 1L
            eligible <- selected[selected$evidence_eligible, , drop = FALSE]
            if (nrow(eligible) == 0L) {
                next
            }
            eligible_summary <- .summarise_scientific_automatic(eligible)
            gates[[next_gate]] <- .scientific_automatic_metric_gates(
                eligible_summary,
                automatic_specification,
                automatic_specification$required_success_rate
            )
            next_gate <- next_gate + 1L
            next
        }

        required_rate <- if (identical(gate_type, "stress")) {
            automatic_specification$stress_success_rate
        } else {
            automatic_specification$required_success_rate
        }
        gates[[next_gate]] <- .scientific_automatic_metric_gates(
            summary,
            automatic_specification,
            required_rate
        )
        next_gate <- next_gate + 1L
    }

    successes <- raw$automatic$status == "ok"
    gates[[next_gate]] <- .scientific_gate_row(
        "automatic_structured_status",
        "all",
        all(raw$automatic$status %in% c("ok", "unidentifiable")),
        paste(sort(unique(raw$automatic$status)), collapse = ","),
        "ok|unidentifiable"
    )
    next_gate <- next_gate + 1L
    gates[[next_gate]] <- .scientific_gate_row(
        "automatic_cutoff_inside_observed_support",
        "all successful runs",
        all(raw$automatic$inside_observed_support[successes]),
        sprintf(
            "%d/%d",
            sum(raw$automatic$inside_observed_support[successes]),
            sum(successes)
        ),
        "all"
    )
    next_gate <- next_gate + 1L
    structured <- raw$automatic$status == "unidentifiable"
    gates[[next_gate]] <- .scientific_gate_row(
        "automatic_structured_condition",
        "all structured failures",
        all(raw$automatic$structured_condition_exact[structured]),
        sprintf(
            "%d/%d",
            sum(raw$automatic$structured_condition_exact[structured]),
            sum(structured)
        ),
        "all"
    )
    next_gate <- next_gate + 1L
    gates[[next_gate]] <- .scientific_gate_row(
        "automatic_structured_profile",
        "all structured failures",
        all(raw$automatic$structured_profile_available[structured]),
        sprintf(
            "%d/%d",
            sum(raw$automatic$structured_profile_available[structured]),
            sum(structured)
        ),
        "all"
    )
    next_gate <- next_gate + 1L

    invocation_metrics <- c(
        input_unchanged = "input exact",
        rng_unchanged = "caller RNG exact",
        output_valid = "observed values + seed provenance exact"
    )
    for (metric in names(invocation_metrics)) {
        values <- raw$invocations[[metric]]
        gates[[next_gate]] <- .scientific_gate_row(
            paste0("invocation_", metric),
            "all public calls",
            all(values),
            sprintf("%d/%d", sum(values), length(values)),
            "all",
            paste(raw$invocations$detail[!values], collapse = ";")
        )
        next_gate <- next_gate + 1L
    }
    for (audit in unique(raw$comparisons$audit)) {
        selected <- raw$comparisons$audit == audit
        values <- raw$comparisons$passed[selected]
        gates[[next_gate]] <- .scientific_gate_row(
            audit,
            "all",
            all(values),
            sprintf("%d/%d", sum(values), length(values)),
            "all",
            paste(raw$comparisons$detail[selected & !raw$comparisons$passed],
                collapse = ";")
        )
        next_gate <- next_gate + 1L
    }

    result <- do.call(rbind, gates)
    rownames(result) <- NULL
    result
}

print_scientific_validation_benchmark <- function(benchmark) {
    cat("\nManual true-cutoff q10 summaries:\n")
    print(benchmark$summaries$manual, row.names = FALSE, digits = 4)
    cat("\nAutomatic per-condition summaries:\n")
    print(benchmark$summaries$automatic, row.names = FALSE, digits = 4)
    cat("\nExact audit summaries:\n")
    print(benchmark$summaries$audits, row.names = FALSE)
    failures <- benchmark$gates[!benchmark$gates$passed, , drop = FALSE]
    cat(sprintf(
        "\nGates: %d/%d passed; %d failed.\n",
        sum(benchmark$gates$passed),
        nrow(benchmark$gates),
        nrow(failures)
    ))
    if (nrow(failures) > 0L) {
        print(failures, row.names = FALSE)
    }
    cat(sprintf(
        "\nProtocol MD5: %s\nResult MD5: %s\nElapsed: %.3f seconds\n",
        benchmark$protocol_hash,
        benchmark$result_hash,
        benchmark$context$elapsed_seconds
    ))
    invisible(benchmark)
}

.scientific_validation_main <- function(
    arguments = commandArgs(trailingOnly = TRUE)
) {
    if (identical(arguments, "--verify")) {
        verified <- verify_scientific_validation()
        print(verified$manifest, row.names = FALSE, digits = 4)
        cat("\nRoutine subset:\n")
        print(verified$routine, row.names = FALSE, digits = 4)
        message(
            sprintf(
                "M7 protocol verification passed: %d long scenarios; MD5 %s.",
                nrow(verified$manifest),
                verified$hash
            )
        )
        return(invisible(verified))
    }
    if (identical(arguments, "--benchmark")) {
        benchmark <- run_scientific_validation_benchmark()
        print_scientific_validation_benchmark(benchmark)
        failed <- sum(!benchmark$gates$passed)
        if (failed > 0L) {
            stop(
                sprintf("M7b scientific benchmark failed %d frozen gates.", failed),
                call. = FALSE
            )
        }
        message("M7b scientific benchmark passed every frozen gate.")
        return(invisible(benchmark))
    }

    stop(
        paste0(
            "Usage: Rscript --vanilla dev/scientific-validation.R ",
            "--verify|--benchmark"
        ),
        call. = FALSE
    )
}

if (sys.nframe() == 0L) {
    .scientific_validation_main()
}
