#!/usr/bin/env Rscript

# M5 automatic-cutoff benchmark harness. Run from the repository root:
# Rscript --vanilla dev/cutoff-validation.R --verify

.cutoff_validation_core <- local({
    required <- c(
        "R/result.R",
        "R/seed_missing_conditions.R",
        "R/profile_missingness.R",
        "R/reconcile.R"
    )
    missing <- required[!file.exists(required)]
    if (length(missing) > 0L) {
        stop(
            "Run the cutoff validation harness from the repository root.",
            call. = FALSE
        )
    }

    core <- new.env(parent = baseenv())
    for (file in required) {
        sys.source(file, envir = core)
    }
    core
})

.validation_transition <- function(left, right, height) {
    stopifnot(
        is.finite(left),
        is.finite(right),
        left < right,
        is.finite(height),
        height > 0,
        height <= 1
    )
    data.frame(left = left, right = right, height = height)
}

.validation_condition <- function(
    transitions = data.frame(left = numeric(), right = numeric(), height = numeric()),
    mar_rate = 0,
    complete = FALSE
) {
    stopifnot(
        is.data.frame(transitions),
        identical(names(transitions), c("left", "right", "height")),
        is.numeric(mar_rate),
        length(mar_rate) == 1L,
        is.finite(mar_rate),
        mar_rate >= 0,
        mar_rate < 1,
        is.logical(complete),
        length(complete) == 1L
    )
    list(
        transitions = transitions,
        mar_rate = mar_rate,
        complete = complete
    )
}

.validation_targets <- function(
    condition,
    identifiable,
    right_boundary,
    tolerance,
    required_success,
    max_left_bias,
    bootstrap_iqr_limit
) {
    data.frame(
        condition = condition,
        identifiable = identifiable,
        right_boundary = right_boundary,
        tolerance = tolerance,
        required_success = required_success,
        max_left_bias = max_left_bias,
        bootstrap_iqr_limit = bootstrap_iqr_limit,
        stringsAsFactors = FALSE
    )
}

.single_target_scenario <- function(
    description,
    n_features,
    samples,
    transitions,
    mar_rate,
    right_boundary,
    tolerance,
    required_success = 0.875,
    max_left_bias = 0.25,
    bootstrap_iqr_limit = 0.5,
    center = 13,
    scale = 2.2,
    measurement_sd = 0.2,
    right_tail_fraction = 0,
    duplicate_step = NA_real_
) {
    identifiable <- is.finite(right_boundary)
    list(
        description = description,
        n_features = as.integer(n_features),
        samples = as.integer(samples),
        intensity = list(
            center = center,
            scale = scale,
            right_tail_fraction = right_tail_fraction,
            duplicate_step = duplicate_step
        ),
        measurement_sd = measurement_sd,
        conditions = list(
            case = .validation_condition(transitions, mar_rate),
            anchor = .validation_condition(complete = TRUE)
        ),
        targets = .validation_targets(
            condition = "case",
            identifiable = identifiable,
            right_boundary = right_boundary,
            tolerance = tolerance,
            required_success = if (identifiable) required_success else 0,
            max_left_bias = if (identifiable) max_left_bias else NA_real_,
            bootstrap_iqr_limit = if (identifiable) {
                bootstrap_iqr_limit
            } else {
                NA_real_
            }
        )
    )
}

cutoff_validation_catalog <- function() {
    sharp <- .validation_transition(11.4, 12, 0.72)
    broad <- .validation_transition(9.5, 12, 0.72)
    standard <- .validation_transition(11, 12, 0.68)

    scenarios <- list(
        sharp_cliff = .single_target_scenario(
            "High-amplitude narrow transition.",
            1600L, 8L, sharp, 0.05, 12, 0.45
        ),
        broad_cliff = .single_target_scenario(
            "High-amplitude broad transition; target remains its right edge.",
            2000L, 8L, broad, 0.05, 12, 0.75,
            max_left_bias = 0.35,
            bootstrap_iqr_limit = 0.65
        ),
        low_amplitude_cliff = .single_target_scenario(
            "Weak intensity-dependent component above sparse MAR.",
            2800L, 8L, .validation_transition(11, 12, 0.18), 0.05,
            12, 0.9, required_success = 0.75,
            bootstrap_iqr_limit = 0.8
        ),
        noisy_right_tail = .single_target_scenario(
            "Sparse high-intensity tail makes the MAR-only profile tail noisy.",
            1800L, 8L, standard, 0.08, 12, 0.9,
            right_tail_fraction = 0.08,
            bootstrap_iqr_limit = 0.75
        ),
        heavy_mar_25 = .single_target_scenario(
            "Twenty-five-percent cell-level intensity-independent MAR background.",
            4000L, 8L, standard, 0.25, 12, 1,
            required_success = 0.75,
            bootstrap_iqr_limit = 0.85
        ),
        sparse_mar_05 = .single_target_scenario(
            "Five-percent cell-level intensity-independent MAR background.",
            1600L, 8L, standard, 0.05, 12, 0.6
        ),
        flat_no_cliff = .single_target_scenario(
            "Intensity-independent missingness only; automatic cutoff is invalid.",
            2000L, 8L,
            data.frame(left = numeric(), right = numeric(), height = numeric()),
            0.12, NA_real_, NA_real_, required_success = 0
        ),
        two_transitions = .single_target_scenario(
            "Dominant low-intensity transition plus a smaller right transition.",
            2600L, 8L,
            rbind(
                .validation_transition(10.8, 11.6, 0.52),
                .validation_transition(13.3, 14.2, 0.16)
            ),
            0.05, 11.6, 0.75,
            bootstrap_iqr_limit = 0.7
        ),
        unequal_class_counts = .single_target_scenario(
            "High-centered intensities leave missing feature-blocks in the minority.",
            2200L, 8L, .validation_transition(11.2, 12, 0.35), 0.02,
            12, 0.75, center = 14.4, scale = 1.8
        ),
        duplicate_means = .single_target_scenario(
            "Rounded latent means and zero measurement noise create exact ties.",
            1800L, 8L, standard, 0.05, 12, 0.5,
            measurement_sd = 0, duplicate_step = 0.25
        ),
        small_feature_count = .single_target_scenario(
            "Ninety-six features, near the predeclared automatic-evidence floor.",
            96L, 8L, sharp, 0.05, 12, 1,
            required_success = 0.75,
            bootstrap_iqr_limit = 1
        )
    )

    scenarios$condition_specific <- list(
        description = "Two independently identifiable cliffs plus a complete anchor.",
        n_features = 2400L,
        samples = 8L,
        intensity = list(
            center = 12.5,
            scale = 2.4,
            right_tail_fraction = 0,
            duplicate_step = NA_real_
        ),
        measurement_sd = 0.2,
        conditions = list(
            A = .validation_condition(
                .validation_transition(10.2, 11, 0.68),
                0.05
            ),
            B = .validation_condition(
                .validation_transition(13.2, 14, 0.68),
                0.05
            ),
            anchor = .validation_condition(complete = TRUE)
        ),
        targets = .validation_targets(
            condition = c("A", "B"),
            identifiable = c(TRUE, TRUE),
            right_boundary = c(11, 14),
            tolerance = c(0.6, 0.6),
            required_success = c(0.875, 0.875),
            max_left_bias = c(0.25, 0.25),
            bootstrap_iqr_limit = c(0.55, 0.55)
        )
    )

    expected <- c(
        "sharp_cliff",
        "broad_cliff",
        "low_amplitude_cliff",
        "noisy_right_tail",
        "heavy_mar_25",
        "sparse_mar_05",
        "flat_no_cliff",
        "two_transitions",
        "unequal_class_counts",
        "duplicate_means",
        "condition_specific",
        "small_feature_count"
    )
    scenarios <- scenarios[expected]
    for (index in seq_along(scenarios)) {
        scenarios[[index]]$seed_offset <- as.integer(index * 10000L)
    }
    scenarios
}

cutoff_validation_seeds <- function() {
    c(104729L, 130363L, 155921L, 181081L, 206369L, 232003L, 257681L, 283369L)
}

.with_validation_seed <- function(seed, code) {
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

.validation_latent_means <- function(scenario) {
    count <- scenario$n_features
    probability <- (seq_len(count) - 0.5) / count
    intensity <- scenario$intensity$center +
        scenario$intensity$scale * stats::qnorm(probability)

    tail_count <- as.integer(round(
        count * scenario$intensity$right_tail_fraction
    ))
    if (tail_count > 0L) {
        tail_rows <- seq.int(count - tail_count + 1L, count)
        intensity[tail_rows] <- seq(16, 21, length.out = tail_count)
    }
    if (is.finite(scenario$intensity$duplicate_step)) {
        step <- scenario$intensity$duplicate_step
        intensity <- round(intensity / step) * step
    } else {
        intensity <- intensity + stats::rnorm(count, sd = 0.03)
    }

    sample(intensity, length(intensity), replace = FALSE)
}

.validation_mnar_probability <- function(intensity, transitions) {
    probability <- numeric(length(intensity))
    if (nrow(transitions) == 0L) {
        return(probability)
    }

    for (index in seq_len(nrow(transitions))) {
        left <- transitions$left[[index]]
        right <- transitions$right[[index]]
        height <- transitions$height[[index]]
        component <- ifelse(
            intensity <= left,
            height,
            ifelse(
                intensity >= right,
                0,
                height * (right - intensity) / (right - left)
            )
        )
        probability <- probability + component
    }

    pmin(probability, 0.98)
}

.simulate_validation_condition <- function(
    latent_mean,
    condition,
    config,
    samples,
    measurement_sd,
    feature_names
) {
    value_count <- length(latent_mean) * samples
    values <- matrix(
        rep(latent_mean, times = samples) +
            stats::rnorm(value_count, sd = measurement_sd),
        nrow = length(latent_mean),
        dimnames = list(
            feature_names,
            sprintf("%s_%02d", condition, seq_len(samples))
        )
    )
    if (config$complete) {
        mar_mask <- mnar_mask <- matrix(FALSE, nrow(values), ncol(values))
    } else {
        mar_mask <- matrix(
            stats::runif(value_count) < config$mar_rate,
            nrow = nrow(values)
        )
        mnar_probability <- .validation_mnar_probability(
            as.numeric(values),
            config$transitions
        )
        mnar_mask <- matrix(
            stats::runif(value_count) < mnar_probability,
            nrow = nrow(values)
        )
    }
    missing_mask <- mar_mask | mnar_mask
    values[missing_mask] <- NA_real_

    has_mnar <- rowSums(mnar_mask) > 0L
    has_mar <- rowSums(mar_mask) > 0L
    has_missing <- rowSums(missing_mask) > 0L
    mechanism <- ifelse(
        !has_missing,
        "complete",
        ifelse(has_mnar, "MNAR", "MAR")
    )
    truth <- data.frame(
        feature = feature_names,
        condition = condition,
        latent_mean = latent_mean,
        generated_missing_count = as.integer(rowSums(missing_mask)),
        generated_mnar_count = as.integer(rowSums(mnar_mask)),
        generated_mar_count = as.integer(rowSums(mar_mask)),
        mechanism = mechanism,
        stringsAsFactors = FALSE
    )

    list(data = values, truth = truth)
}

.validation_pipeline <- function(data, groups_by_sample, rescue_seed = 1L) {
    core <- .cutoff_validation_core
    rescued <- core$.seed_missing_conditions(
        data,
        groups_by_sample,
        rescue_seed
    )
    statistics <- core$.feature_condition_statistics(
        rescued$data,
        groups_by_sample,
        rescued$seed_log
    )
    list(
        rescued = rescued,
        statistics = statistics,
        profiles = core$.build_missingness_profiles(statistics)
    )
}

simulate_cutoff_scenario <- function(
    name,
    seed = cutoff_validation_seeds()[[1L]],
    catalog = cutoff_validation_catalog()
) {
    if (!is.character(name) || length(name) != 1L || !name %in% names(catalog)) {
        stop("`name` must identify one cutoff-validation scenario.", call. = FALSE)
    }
    scenario <- catalog[[name]]
    actual_seed <- seed + scenario$seed_offset
    valid_seed <- is.numeric(actual_seed) && length(actual_seed) == 1L &&
        is.finite(actual_seed) && actual_seed == trunc(actual_seed) &&
        actual_seed >= 0 && actual_seed <= .Machine$integer.max
    if (!valid_seed) {
        stop("Scenario seed + offset must be a non-negative integer.", call. = FALSE)
    }
    actual_seed <- as.integer(actual_seed)

    .with_validation_seed(actual_seed, {
        feature_names <- sprintf("feature_%05d", seq_len(scenario$n_features))
        latent_mean <- .validation_latent_means(scenario)
        generated <- lapply(
            names(scenario$conditions),
            function(condition) {
                .simulate_validation_condition(
                    latent_mean,
                    condition,
                    scenario$conditions[[condition]],
                    scenario$samples,
                    scenario$measurement_sd,
                    feature_names
                )
            }
        )
        names(generated) <- names(scenario$conditions)
        data <- do.call(cbind, lapply(generated, `[[`, "data"))
        groups_by_sample <- stats::setNames(
            rep(names(generated), each = scenario$samples),
            colnames(data)
        )
        generated_truth <- do.call(
            rbind,
            lapply(generated, `[[`, "truth")
        )
        rownames(generated_truth) <- NULL
        pipeline <- .validation_pipeline(data, groups_by_sample)

        key <- paste(
            pipeline$statistics$feature,
            pipeline$statistics$condition,
            sep = "\r"
        )
        truth_key <- paste(
            generated_truth$feature,
            generated_truth$condition,
            sep = "\r"
        )
        truth <- generated_truth[match(key, truth_key), , drop = FALSE]
        statistics <- pipeline$statistics
        truth$sample_count <- statistics$sample_count
        truth$observed_count <- statistics$observed_count
        truth$missing_count <- statistics$missing_count
        truth$seeded <- statistics$seeded
        truth$state <- ifelse(
            statistics$missing_count == 0L,
            "complete",
            ifelse(
                truth$mechanism == "MNAR",
                "MNAR",
                ifelse(
                    statistics$observed_count > statistics$sample_count / 2,
                    "MAR",
                    "insufficient"
                )
            )
        )
        rownames(truth) <- NULL

        list(
            name = name,
            seed = actual_seed,
            scenario = scenario,
            data = data,
            groups_by_sample = groups_by_sample,
            rescued = pipeline$rescued,
            statistics = statistics,
            profiles = pipeline$profiles,
            truth = truth
        )
    })
}

.binary_metrics <- function(predicted, truth) {
    stopifnot(
        is.logical(predicted),
        is.logical(truth),
        length(predicted) == length(truth),
        !anyNA(predicted),
        !anyNA(truth)
    )
    true_positive <- sum(predicted & truth)
    false_positive <- sum(predicted & !truth)
    false_negative <- sum(!predicted & truth)
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

.condition_cutoff_metrics <- function(simulation, condition, cutoff) {
    rows <- simulation$statistics$condition == condition &
        simulation$statistics$missing_count > 0L
    predicted_mnar <- simulation$statistics$mean_intensity[rows] < cutoff
    true_mnar <- simulation$truth$mechanism[rows] == "MNAR"
    .binary_metrics(predicted_mnar, true_mnar)
}

.truth_retained <- function(simulation) {
    features <- rownames(simulation$data)
    conditions <- sort(unique(simulation$truth$condition), method = "radix")
    state <- matrix(
        simulation$truth$state,
        nrow = length(features),
        ncol = length(conditions),
        byrow = TRUE,
        dimnames = list(features, conditions)
    )
    if (!identical(simulation$truth$feature, rep(features, each = length(conditions))) ||
        !identical(simulation$truth$condition, rep(conditions, times = length(features)))) {
        state[cbind(
            match(simulation$truth$feature, features),
            match(simulation$truth$condition, conditions)
        )] <- simulation$truth$state
    }

    !apply(state == "insufficient", 1L, any) &
        !apply(state == "MNAR", 1L, all)
}

.predicted_retention_metrics <- function(simulation, target_cutoffs) {
    conditions <- sort(unique(simulation$statistics$condition), method = "radix")
    cutoffs <- stats::setNames(rep(NA_real_, length(conditions)), conditions)
    cutoffs[names(target_cutoffs)] <- target_cutoffs
    classified <- .cutoff_validation_core$.assign_condition_states(
        simulation$statistics,
        cutoffs
    )
    reconciled <- .cutoff_validation_core$.reconcile_condition_states(
        classified,
        simulation$rescued$feature_status
    )
    predicted <- stats::setNames(
        reconciled$feature_status$retained,
        reconciled$feature_status$feature
    )
    truth <- .truth_retained(simulation)
    names(truth) <- rownames(simulation$data)
    .binary_metrics(predicted[names(truth)], truth)
}

.normalise_detector_result <- function(result) {
    if (!is.list(result) || !is.character(result$status) ||
        length(result$status) != 1L || is.na(result$status) ||
        !result$status %in% c("ok", "unidentifiable")) {
        return(list(
            status = "invalid",
            cutoff = NA_real_,
            reason = "Detector result must use status `ok` or `unidentifiable`."
        ))
    }
    if (identical(result$status, "ok")) {
        valid <- is.numeric(result$cutoff) && length(result$cutoff) == 1L &&
            is.finite(result$cutoff)
        if (!valid) {
            return(list(
                status = "invalid",
                cutoff = NA_real_,
                reason = "An `ok` detector result requires one finite cutoff."
            ))
        }
        return(list(status = "ok", cutoff = as.numeric(result$cutoff), reason = ""))
    }

    reason <- result$reason
    valid_reason <- is.character(reason) && length(reason) == 1L &&
        !is.na(reason) && nzchar(reason)
    valid_cutoff <- is.null(result$cutoff) ||
        (is.numeric(result$cutoff) && length(result$cutoff) == 1L &&
            is.na(result$cutoff))
    if (!valid_reason || !valid_cutoff) {
        return(list(
            status = "invalid",
            cutoff = NA_real_,
            reason = paste0(
                "An `unidentifiable` result requires an NA cutoff and a ",
                "non-empty reason."
            )
        ))
    }
    list(status = "unidentifiable", cutoff = NA_real_, reason = reason)
}

.invoke_detector <- function(detector, profile) {
    started <- proc.time()[["elapsed"]]
    invocation <- tryCatch(
        list(value = detector(profile), error = NULL),
        error = function(error) {
            list(value = NULL, error = conditionMessage(error))
        }
    )
    elapsed_ms <- 1000 * (proc.time()[["elapsed"]] - started)
    normalised <- if (is.null(invocation$error)) {
        .normalise_detector_result(invocation$value)
    } else {
        list(
            status = "invalid",
            cutoff = NA_real_,
            reason = invocation$error
        )
    }
    normalised$elapsed_ms <- elapsed_ms
    normalised
}

.same_detector_decision <- function(left, right) {
    identical(left$status, right$status) &&
        identical(left$cutoff, right$cutoff)
}

benchmark_cutoff_detectors <- function(
    detectors,
    seeds = cutoff_validation_seeds(),
    catalog = cutoff_validation_catalog()
) {
    valid <- is.list(detectors) && length(detectors) > 0L &&
        !is.null(names(detectors)) && all(nzchar(names(detectors))) &&
        !anyDuplicated(names(detectors)) &&
        all(vapply(detectors, is.function, logical(1L)))
    if (!valid) {
        stop("`detectors` must be a uniquely named list of functions.", call. = FALSE)
    }

    rows <- list()
    next_row <- 1L
    for (scenario_name in names(catalog)) {
        targets <- catalog[[scenario_name]]$targets
        for (seed in seeds) {
            simulation <- simulate_cutoff_scenario(
                scenario_name,
                seed,
                catalog
            )
            oracle_retention <- if (all(targets$identifiable)) {
                .predicted_retention_metrics(
                    simulation,
                    stats::setNames(
                        targets$right_boundary,
                        targets$condition
                    )
                )
            } else {
                c(precision = NA_real_, recall = NA_real_, f1 = NA_real_)
            }

            for (detector_name in names(detectors)) {
                detector <- detectors[[detector_name]]
                decisions <- vector("list", nrow(targets))
                for (target_index in seq_len(nrow(targets))) {
                    target <- targets[target_index, , drop = FALSE]
                    condition <- target$condition[[1L]]
                    profile <- simulation$profiles[[condition]]
                    decision <- .invoke_detector(detector, profile)
                    repeated <- .invoke_detector(detector, profile)
                    permuted_profile <- .cutoff_validation_core$.condition_missingness_profile(
                        profile$raw[rev(seq_len(nrow(profile$raw))), , drop = FALSE]
                    )
                    permuted <- .invoke_detector(detector, permuted_profile)
                    decision$repeat_stable <- .same_detector_decision(
                        decision,
                        repeated
                    )
                    decision$order_stable <- .same_detector_decision(
                        decision,
                        permuted
                    )
                    decisions[[target_index]] <- decision
                }

                all_ok <- all(vapply(
                    decisions,
                    function(decision) identical(decision$status, "ok"),
                    logical(1L)
                ))
                candidate_retention <- if (all_ok) {
                    .predicted_retention_metrics(
                        simulation,
                        stats::setNames(
                            vapply(decisions, `[[`, numeric(1L), "cutoff"),
                            targets$condition
                        )
                    )
                } else {
                    c(precision = NA_real_, recall = NA_real_, f1 = NA_real_)
                }

                for (target_index in seq_len(nrow(targets))) {
                    target <- targets[target_index, , drop = FALSE]
                    decision <- decisions[[target_index]]
                    condition <- target$condition[[1L]]
                    identifiable <- target$identifiable[[1L]]
                    oracle_mnar <- if (identifiable) {
                        .condition_cutoff_metrics(
                            simulation,
                            condition,
                            target$right_boundary[[1L]]
                        )
                    } else {
                        c(precision = NA_real_, recall = NA_real_, f1 = NA_real_)
                    }
                    candidate_mnar <- if (identical(decision$status, "ok")) {
                        .condition_cutoff_metrics(
                            simulation,
                            condition,
                            decision$cutoff
                        )
                    } else {
                        c(precision = NA_real_, recall = NA_real_, f1 = NA_real_)
                    }
                    signed_error <- if (identifiable &&
                        identical(decision$status, "ok")) {
                        decision$cutoff - target$right_boundary[[1L]]
                    } else {
                        NA_real_
                    }

                    rows[[next_row]] <- data.frame(
                        detector = detector_name,
                        scenario = scenario_name,
                        seed = simulation$seed,
                        condition = condition,
                        identifiable = identifiable,
                        target_cutoff = target$right_boundary[[1L]],
                        tolerance = target$tolerance[[1L]],
                        status = decision$status,
                        cutoff = decision$cutoff,
                        reason = decision$reason,
                        signed_error = signed_error,
                        absolute_error = abs(signed_error),
                        within_tolerance = identifiable &&
                            is.finite(signed_error) &&
                            abs(signed_error) <= target$tolerance[[1L]],
                        false_confidence = !identifiable &&
                            identical(decision$status, "ok"),
                        structured_failure = !identifiable &&
                            identical(decision$status, "unidentifiable"),
                        mnar_precision = candidate_mnar[["precision"]],
                        mnar_recall = candidate_mnar[["recall"]],
                        mnar_f1 = candidate_mnar[["f1"]],
                        oracle_mnar_precision = oracle_mnar[["precision"]],
                        oracle_mnar_recall = oracle_mnar[["recall"]],
                        oracle_mnar_f1 = oracle_mnar[["f1"]],
                        retention_precision = candidate_retention[["precision"]],
                        retention_recall = candidate_retention[["recall"]],
                        retention_f1 = candidate_retention[["f1"]],
                        oracle_retention_f1 = oracle_retention[["f1"]],
                        repeat_stable = decision$repeat_stable,
                        order_stable = decision$order_stable,
                        elapsed_ms = decision$elapsed_ms,
                        stringsAsFactors = FALSE
                    )
                    next_row <- next_row + 1L
                }
            }
        }
    }

    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

.safe_quantile <- function(x, probability) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    unname(stats::quantile(x, probability, names = FALSE, type = 8L))
}

summarise_cutoff_benchmark <- function(results) {
    required <- c(
        "detector", "scenario", "condition", "identifiable", "status",
        "signed_error", "absolute_error", "false_confidence",
        "structured_failure", "mnar_f1", "oracle_mnar_f1", "retention_f1",
        "oracle_retention_f1", "repeat_stable", "order_stable", "elapsed_ms"
    )
    if (!is.data.frame(results) || !all(required %in% names(results))) {
        stop("`results` is not a cutoff benchmark result.", call. = FALSE)
    }

    groups <- split(
        seq_len(nrow(results)),
        interaction(
            results$detector,
            results$scenario,
            results$condition,
            drop = TRUE,
            lex.order = TRUE
        )
    )
    summaries <- lapply(
        groups,
        function(rows) {
            block <- results[rows, , drop = FALSE]
            data.frame(
                detector = block$detector[[1L]],
                scenario = block$scenario[[1L]],
                condition = block$condition[[1L]],
                identifiable = block$identifiable[[1L]],
                success_rate = mean(block$status == "ok"),
                structured_failure_rate = mean(block$structured_failure),
                false_confidence_rate = mean(block$false_confidence),
                median_signed_error = stats::median(
                    block$signed_error,
                    na.rm = TRUE
                ),
                median_absolute_error = stats::median(
                    block$absolute_error,
                    na.rm = TRUE
                ),
                q90_absolute_error = .safe_quantile(
                    block$absolute_error,
                    0.9
                ),
                q10_mnar_f1_delta = .safe_quantile(
                    block$mnar_f1 - block$oracle_mnar_f1,
                    0.1
                ),
                q10_retention_f1_delta = .safe_quantile(
                    block$retention_f1 - block$oracle_retention_f1,
                    0.1
                ),
                repeat_stable = all(block$repeat_stable),
                order_stable = all(block$order_stable),
                median_elapsed_ms = stats::median(block$elapsed_ms),
                p95_elapsed_ms = .safe_quantile(block$elapsed_ms, 0.95),
                stringsAsFactors = FALSE
            )
        }
    )
    summary <- do.call(rbind, summaries)
    rownames(summary) <- NULL
    summary[order(summary$detector, summary$scenario, summary$condition), ]
}

assess_cutoff_benchmark <- function(
    summary,
    catalog = cutoff_validation_catalog()
) {
    required <- c(
        "detector", "scenario", "condition", "identifiable", "success_rate",
        "structured_failure_rate", "false_confidence_rate",
        "median_signed_error", "median_absolute_error", "q90_absolute_error",
        "q10_mnar_f1_delta", "q10_retention_f1_delta", "repeat_stable",
        "order_stable", "median_elapsed_ms", "p95_elapsed_ms"
    )
    if (!is.data.frame(summary) || !all(required %in% names(summary))) {
        stop("`summary` is not a cutoff benchmark summary.", call. = FALSE)
    }

    assessed <- summary
    assessed$failed_gates <- NA_character_
    assessed$passed <- FALSE
    for (row in seq_len(nrow(assessed))) {
        scenario <- assessed$scenario[[row]]
        condition <- assessed$condition[[row]]
        if (!scenario %in% names(catalog)) {
            stop("Benchmark summary contains an unknown scenario.", call. = FALSE)
        }
        targets <- catalog[[scenario]]$targets
        target_row <- match(condition, targets$condition)
        if (is.na(target_row)) {
            stop("Benchmark summary contains an unknown target condition.", call. = FALSE)
        }
        target <- targets[target_row, , drop = FALSE]
        failed <- character()
        if (target$identifiable[[1L]]) {
            checks <- c(
                success_rate = assessed$success_rate[[row]] >=
                    target$required_success[[1L]],
                median_error = is.finite(assessed$median_absolute_error[[row]]) &&
                    assessed$median_absolute_error[[row]] <=
                        target$tolerance[[1L]],
                q90_error = is.finite(assessed$q90_absolute_error[[row]]) &&
                    assessed$q90_absolute_error[[row]] <=
                        1.5 * target$tolerance[[1L]],
                right_boundary = is.finite(assessed$median_signed_error[[row]]) &&
                    assessed$median_signed_error[[row]] >=
                        -target$max_left_bias[[1L]] &&
                    assessed$median_signed_error[[row]] <=
                        target$tolerance[[1L]],
                mnar_f1 = is.finite(assessed$q10_mnar_f1_delta[[row]]) &&
                    assessed$q10_mnar_f1_delta[[row]] >= -0.05,
                retention_f1 = is.finite(
                    assessed$q10_retention_f1_delta[[row]]
                ) && assessed$q10_retention_f1_delta[[row]] >= -0.03
            )
        } else {
            checks <- c(
                no_finite_cutoff = assessed$success_rate[[row]] == 0,
                structured_failure = assessed$structured_failure_rate[[row]] == 1,
                false_confidence = assessed$false_confidence_rate[[row]] == 0
            )
        }
        checks <- c(
            checks,
            repeat_stable = isTRUE(assessed$repeat_stable[[row]]),
            order_stable = isTRUE(assessed$order_stable[[row]]),
            median_runtime = is.finite(assessed$median_elapsed_ms[[row]]) &&
                assessed$median_elapsed_ms[[row]] <= 250,
            p95_runtime = is.finite(assessed$p95_elapsed_ms[[row]]) &&
                assessed$p95_elapsed_ms[[row]] <= 1000
        )
        failed <- names(checks)[!checks]
        assessed$failed_gates[[row]] <- paste(failed, collapse = ",")
        assessed$passed[[row]] <- length(failed) == 0L
    }

    assessed
}

bootstrap_cutoff_detector <- function(
    detector,
    simulation,
    condition,
    replicates = 200L,
    seed = 424242L
) {
    if (!is.function(detector) || !condition %in% names(simulation$profiles) ||
        length(replicates) != 1L || replicates < 1L) {
        stop("Invalid bootstrap benchmark input.", call. = FALSE)
    }
    raw <- simulation$profiles[[condition]]$raw
    classes <- split(seq_len(nrow(raw)), raw$has_missing)

    .with_validation_seed(as.integer(seed), {
        rows <- vector("list", as.integer(replicates))
        for (replicate_index in seq_len(replicates)) {
            sampled <- unlist(
                lapply(
                    classes,
                    function(rows) sample(rows, length(rows), replace = TRUE)
                ),
                use.names = FALSE
            )
            boot <- raw[sampled, , drop = FALSE]
            boot$feature <- sprintf(
                "bootstrap_%04d_%05d",
                replicate_index,
                seq_len(nrow(boot))
            )
            profile <- .cutoff_validation_core$.condition_missingness_profile(
                boot
            )
            decision <- .invoke_detector(detector, profile)
            rows[[replicate_index]] <- data.frame(
                replicate = replicate_index,
                status = decision$status,
                cutoff = decision$cutoff,
                reason = decision$reason,
                stringsAsFactors = FALSE
            )
        }
        do.call(rbind, rows)
    })
}

assess_cutoff_bootstrap <- function(results, target) {
    valid_target <- is.data.frame(target) && nrow(target) == 1L &&
        all(c(
            "identifiable", "right_boundary", "tolerance",
            "required_success", "bootstrap_iqr_limit"
        ) %in% names(target))
    if (!is.data.frame(results) ||
        !all(c("status", "cutoff") %in% names(results)) || !valid_target) {
        stop("Invalid bootstrap assessment input.", call. = FALSE)
    }
    successful <- results$status == "ok" & is.finite(results$cutoff)
    success_rate <- mean(successful)
    if (target$identifiable[[1L]]) {
        error <- abs(results$cutoff[successful] - target$right_boundary[[1L]])
        cutoff_iqr <- if (any(successful)) {
            unname(stats::IQR(results$cutoff[successful], type = 8L))
        } else {
            NA_real_
        }
        q90_error <- .safe_quantile(error, 0.9)
        checks <- c(
            success_rate = success_rate >= target$required_success[[1L]],
            cutoff_iqr = is.finite(cutoff_iqr) &&
                cutoff_iqr <= target$bootstrap_iqr_limit[[1L]],
            q90_error = is.finite(q90_error) &&
                q90_error <= 1.5 * target$tolerance[[1L]]
        )
    } else {
        cutoff_iqr <- q90_error <- NA_real_
        checks <- c(
            no_finite_cutoff = success_rate == 0,
            structured_failure = all(results$status == "unidentifiable")
        )
    }

    data.frame(
        success_rate = success_rate,
        cutoff_iqr = cutoff_iqr,
        q90_absolute_error = q90_error,
        passed = all(checks),
        failed_gates = paste(names(checks)[!checks], collapse = ","),
        stringsAsFactors = FALSE
    )
}

cutoff_validation_manifest <- function(simulations) {
    rows <- list()
    next_row <- 1L
    for (simulation in simulations) {
        targets <- simulation$scenario$targets
        for (target_index in seq_len(nrow(targets))) {
            target <- targets[target_index, , drop = FALSE]
            condition <- target$condition[[1L]]
            profile <- simulation$profiles[[condition]]
            rows[[next_row]] <- data.frame(
                scenario = simulation$name,
                condition = condition,
                identifiable = target$identifiable[[1L]],
                target = target$right_boundary[[1L]],
                tolerance = target$tolerance[[1L]],
                features = profile$metadata$feature_count,
                samples = simulation$scenario$samples,
                missing = unname(profile$metadata$class_counts[["missing"]]),
                complete = unname(profile$metadata$class_counts[["complete"]]),
                seeded = profile$metadata$seeded_feature_count,
                mean_min = unname(profile$metadata$observed_mean_range[["minimum"]]),
                mean_max = unname(profile$metadata$observed_mean_range[["maximum"]]),
                stringsAsFactors = FALSE
            )
            next_row <- next_row + 1L
        }
    }
    manifest <- do.call(rbind, rows)
    rownames(manifest) <- NULL
    manifest
}

.canonical_seed_log <- function(seed_log) {
    seed_log <- seed_log[order(
        seed_log$condition,
        seed_log$feature,
        seed_log$sample
    ), , drop = FALSE]
    rownames(seed_log) <- NULL
    seed_log
}

verify_cutoff_validation <- function() {
    catalog <- cutoff_validation_catalog()
    expected <- c(
        "sharp_cliff", "broad_cliff", "low_amplitude_cliff",
        "noisy_right_tail", "heavy_mar_25", "sparse_mar_05",
        "flat_no_cliff", "two_transitions", "unequal_class_counts",
        "duplicate_means", "condition_specific", "small_feature_count"
    )
    stopifnot(identical(names(catalog), expected))

    caller_has_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
    if (caller_has_seed) {
        caller_seed <- get(".Random.seed", globalenv(), inherits = FALSE)
    }
    simulations <- lapply(
        names(catalog),
        simulate_cutoff_scenario,
        seed = cutoff_validation_seeds()[[1L]],
        catalog = catalog
    )
    names(simulations) <- names(catalog)
    if (caller_has_seed) {
        stopifnot(identical(
            get(".Random.seed", globalenv(), inherits = FALSE),
            caller_seed
        ))
    } else {
        stopifnot(!exists(".Random.seed", globalenv(), inherits = FALSE))
    }

    repeated <- simulate_cutoff_scenario(
        "sharp_cliff",
        cutoff_validation_seeds()[[1L]],
        catalog
    )
    stopifnot(identical(simulations$sharp_cliff, repeated))
    manifest <- cutoff_validation_manifest(simulations)
    stopifnot(
        nrow(manifest) == 13L,
        all(manifest$missing >= 12L),
        all(manifest$complete >= 12L),
        all(manifest$identifiable == is.finite(manifest$target)),
        all(
            !manifest$identifiable |
                (manifest$target > manifest$mean_min &
                    manifest$target < manifest$mean_max)
        )
    )

    flat_truth <- simulations$flat_no_cliff$truth
    stopifnot(!any(flat_truth$mechanism == "MNAR"))
    duplicate_raw <- simulations$duplicate_means$profiles$case$raw
    stopifnot(anyDuplicated(duplicate_raw$mean_intensity) > 0L)
    stopifnot(
        nrow(catalog$two_transitions$conditions$case$transitions) == 2L,
        identical(
            catalog$condition_specific$targets$right_boundary,
            c(11, 14)
        ),
        catalog$small_feature_count$n_features == 96L,
        manifest$mean_max[manifest$scenario == "noisy_right_tail"] > 20
    )
    unequal <- manifest[manifest$scenario == "unequal_class_counts", ]
    stopifnot(unequal$complete / unequal$missing > 2)

    baseline <- simulations$sharp_cliff
    reversed_profiles <- .cutoff_validation_core$.build_missingness_profiles(
        baseline$statistics[rev(seq_len(nrow(baseline$statistics))), ]
    )
    for (condition in names(baseline$profiles)) {
        stopifnot(
            identical(
                baseline$profiles[[condition]]$grid,
                reversed_profiles[[condition]]$grid
            ),
            identical(
                baseline$profiles[[condition]]$metadata,
                reversed_profiles[[condition]]$metadata
            )
        )
    }
    permutation <- rev(seq_len(ncol(baseline$data)))
    feature_permutation <- rev(seq_len(nrow(baseline$data)))
    permuted_pipeline <- .validation_pipeline(
        baseline$data[feature_permutation, permutation, drop = FALSE],
        baseline$groups_by_sample[permutation]
    )
    stopifnot(identical(
        .canonical_seed_log(baseline$rescued$seed_log),
        .canonical_seed_log(permuted_pipeline$rescued$seed_log)
    ))
    for (condition in names(baseline$profiles)) {
        stopifnot(identical(
            baseline$profiles[[condition]]$grid,
            permuted_pipeline$profiles[[condition]]$grid
        ))
    }

    constant <- function(profile) {
        list(status = "ok", cutoff = 12)
    }
    structured_failure <- function(profile) {
        list(
            status = "unidentifiable",
            cutoff = NA_real_,
            reason = "No transition in harness self-check."
        )
    }
    invalid_check <- .invoke_detector(
        function(profile) 12,
        simulations$sharp_cliff$profiles$case
    )
    sharp_check <- benchmark_cutoff_detectors(
        list(constant = constant),
        seeds = cutoff_validation_seeds()[[1L]],
        catalog = catalog["sharp_cliff"]
    )
    flat_check <- benchmark_cutoff_detectors(
        list(structured_failure = structured_failure),
        seeds = cutoff_validation_seeds()[[1L]],
        catalog = catalog["flat_no_cliff"]
    )
    stopifnot(
        sharp_check$status == "ok",
        sharp_check$absolute_error == 0,
        sharp_check$repeat_stable,
        sharp_check$order_stable,
        flat_check$status == "unidentifiable",
        flat_check$structured_failure,
        !flat_check$false_confidence,
        invalid_check$status == "invalid",
        nrow(summarise_cutoff_benchmark(rbind(sharp_check, flat_check))) == 2L,
        assess_cutoff_benchmark(
            summarise_cutoff_benchmark(rbind(sharp_check, flat_check)),
            catalog
        )$passed |> all()
    )

    manifest
}

.cutoff_validation_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
    if (!identical(arguments, "--verify")) {
        stop(
            "Usage: Rscript --vanilla dev/cutoff-validation.R --verify",
            call. = FALSE
        )
    }
    manifest <- verify_cutoff_validation()
    print(manifest, row.names = FALSE, digits = 4)
    message(
        sprintf(
            "M5a verification passed: %d scenarios, %d target profiles.",
            length(unique(manifest$scenario)),
            nrow(manifest)
        )
    )
    invisible(manifest)
}

if (sys.nframe() == 0L) {
    .cutoff_validation_main()
}
