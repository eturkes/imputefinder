#!/usr/bin/env Rscript

# M7 scientific-regression protocol. M7a freezes simulation truth, metrics,
# gates, and the routine/long-run split before classifier results are inspected.
# Run from the repository root:
# Rscript --vanilla dev/scientific-validation.R --verify

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

.scientific_validation_main <- function(
    arguments = commandArgs(trailingOnly = TRUE)
) {
    if (!identical(arguments, "--verify")) {
        stop(
            paste0(
                "Usage: Rscript --vanilla dev/scientific-validation.R ",
                "--verify (benchmark execution opens after the M7a freeze)"
            ),
            call. = FALSE
        )
    }
    verified <- verify_scientific_validation()
    print(verified$manifest, row.names = FALSE, digits = 4)
    cat("\nRoutine subset:\n")
    print(verified$routine, row.names = FALSE, digits = 4)
    message(
        sprintf(
            "M7a verification passed: %d long scenarios; protocol MD5 %s.",
            nrow(verified$manifest),
            verified$hash
        )
    )
    invisible(verified)
}

if (sys.nframe() == 0L) {
    .scientific_validation_main()
}
