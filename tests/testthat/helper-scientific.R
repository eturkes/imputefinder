# Compact protocol-v1 mirror for routine M7 regressions. The independent long
# benchmark and normative protocol remain in dev/scientific-validation.R.

.routine_scientific_conditions <- c("A", "B")
.routine_scientific_states <- c("complete", "MNAR", "MAR", "insufficient")

.scientific_routine_scenario <- function(name) {
    scenario <- switch(
        name,
        manual_on_off = list(
            n_features = 320L,
            samples_per_condition = 4L,
            mar_rate = 0.25
        ),
        automatic_cliff = list(
            n_features = 800L,
            samples_per_condition = 8L,
            mar_rate = 0.05
        ),
        NULL
    )
    if (is.null(scenario)) {
        stop("Unknown routine scientific scenario.", call. = FALSE)
    }

    c(
        scenario,
        list(
            cutoff = 12,
            mnar_width = 1.2,
            mnar_max_probability = 0.8,
            on_off_fraction = 0.05,
            mean_center_offset = 1,
            mean_sd = 2.2,
            condition_sd = 0.35,
            measurement_sd = 0.18
        )
    )
}

.with_routine_scientific_seed <- function(seed, code) {
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

.scientific_routine_oracle <- function(
    data,
    groups,
    mnar_mask,
    feature_truth
) {
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
            reason <- if (length(insufficient) > 0L) {
                paste0(
                    "insufficient:",
                    paste(sort(insufficient, method = "radix"), collapse = ",")
                )
            } else if (all(feature_states == "MNAR")) {
                "MNAR_all_conditions"
            } else {
                NA_character_
            }
            status_index <- surviving_status[[feature_index]]
            feature_status$drop_reason[[status_index]] <- reason
            feature_status$retained[[status_index]] <- is.na(reason)
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

scientific_routine_fixture <- function(name, seed = 104677L) {
    scenario <- .scientific_routine_scenario(name)
    .with_routine_scientific_seed(seed, {
        n <- scenario$n_features
        samples <- scenario$samples_per_condition
        cutoff <- scenario$cutoff
        features <- sprintf("protein_%05d", seq_len(n))
        sample_names <- unlist(
            lapply(
                .routine_scientific_conditions,
                function(condition) {
                    sprintf("%s_%02d", condition, seq_len(samples))
                }
            ),
            use.names = FALSE
        )
        groups <- stats::setNames(
            rep(.routine_scientific_conditions, each = samples),
            sample_names
        )
        relative_mean <- stats::qnorm(
            (seq_len(n) - 0.5) / n,
            mean = scenario$mean_center_offset,
            sd = scenario$mean_sd
        )
        relative_mean <- sample(relative_mean, n, replace = FALSE)
        condition_effect <- stats::rnorm(n, sd = scenario$condition_sd)
        condition_means <- cbind(
            A = cutoff + relative_mean + condition_effect / 2,
            B = cutoff + relative_mean - condition_effect / 2
        )

        on_off_count <- max(
            1L,
            as.integer(floor(n * scenario$on_off_fraction))
        )
        on_off_features <- sample(
            seq_len(n),
            2L * on_off_count,
            replace = FALSE
        )
        off_a <- on_off_features[seq_len(on_off_count)]
        off_b <- on_off_features[on_off_count + seq_len(on_off_count)]
        condition_means[off_a, "A"] <- cutoff - 2.5
        condition_means[off_a, "B"] <- cutoff + 2.5
        condition_means[off_b, "A"] <- cutoff + 2.5
        condition_means[off_b, "B"] <- cutoff - 2.5

        complete_data <- do.call(
            cbind,
            lapply(
                .routine_scientific_conditions,
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
            stats::rnorm(
                length(complete_data),
                sd = scenario$measurement_sd
            ),
            nrow = n
        )
        dimnames(complete_data) <- list(features, sample_names)

        mnar_probability <- scenario$mnar_max_probability * pmin(
            pmax(
                (cutoff - complete_data) / scenario$mnar_width,
                0
            ),
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
        data <- complete_data
        data[mar_mask | mnar_mask] <- NA_real_

        feature_truth <- data.frame(
            feature = features,
            off_condition = NA_character_,
            stringsAsFactors = FALSE
        )
        feature_truth$off_condition[off_a] <- "A"
        feature_truth$off_condition[off_b] <- "B"

        list(
            scenario = scenario,
            simulation_seed = as.integer(seed),
            data = data,
            groups = groups,
            oracle = .scientific_routine_oracle(
                data,
                groups,
                mnar_mask,
                feature_truth
            )
        )
    })
}

.scientific_routine_binary_metrics <- function(truth, predicted) {
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

score_scientific_routine <- function(result, simulation) {
    truth <- simulation$oracle$classifications
    truth_keys <- paste(truth$feature, truth$condition, sep = "\r")
    predicted_keys <- paste(
        result$classifications$feature,
        result$classifications$condition,
        sep = "\r"
    )
    stopifnot(
        !anyDuplicated(predicted_keys),
        setequal(predicted_keys, truth_keys)
    )
    predicted_state <- result$classifications$state[
        match(truth_keys, predicted_keys)
    ]
    stopifnot(
        is.character(predicted_state),
        !anyNA(predicted_state),
        all(predicted_state %in% .routine_scientific_states)
    )

    state_f1 <- vapply(
        .routine_scientific_states,
        function(state) {
            unname(.scientific_routine_binary_metrics(
                truth$state == state,
                predicted_state == state
            )[["f1"]])
        },
        numeric(1L)
    )
    mnar <- .scientific_routine_binary_metrics(
        truth$state == "MNAR",
        predicted_state == "MNAR"
    )
    mar <- .scientific_routine_binary_metrics(
        truth$state == "MAR",
        predicted_state == "MAR"
    )

    truth_status <- simulation$oracle$feature_status
    status_match <- match(truth_status$feature, result$feature_status$feature)
    stopifnot(
        !anyDuplicated(result$feature_status$feature),
        setequal(result$feature_status$feature, truth_status$feature),
        !anyNA(status_match)
    )
    predicted_retained <- result$feature_status$retained[status_match]
    stopifnot(is.logical(predicted_retained), !anyNA(predicted_retained))
    retention <- .scientific_routine_binary_metrics(
        truth_status$retained,
        predicted_retained
    )
    feature_truth <- simulation$oracle$feature_truth
    on_off <- !is.na(feature_truth$off_condition)
    eligible_on_off <- on_off & feature_truth$oracle_retained
    not_globally_absent <- is.na(truth_status$drop_reason) |
        truth_status$drop_reason != "all_missing"
    expected_on_off_seed <- on_off & not_globally_absent
    expected_seed <- paste(
        feature_truth$feature[expected_on_off_seed],
        feature_truth$off_condition[expected_on_off_seed],
        sep = "\r"
    )
    actual_seed <- paste(
        result$seed_log$feature,
        result$seed_log$condition,
        sep = "\r"
    )

    c(
        state_macro_f1 = mean(state_f1, na.rm = TRUE),
        mnar_f1 = unname(mnar[["f1"]]),
        mar_f1 = unname(mar[["f1"]]),
        retention_f1 = unname(retention[["f1"]]),
        on_off_retention_recall = mean(predicted_retained[eligible_on_off]),
        on_off_seed_recall = mean(expected_seed %in% actual_seed)
    )
}

permute_scientific_routine <- function(simulation) {
    feature_count <- nrow(simulation$data)
    conditions <- rev(sort(unique(unname(simulation$groups)), method = "radix"))
    list(
        features = c(
            seq.int(2L, feature_count, by = 2L),
            seq.int(1L, feature_count, by = 2L)
        ),
        samples = unlist(
            lapply(
                conditions,
                function(condition) {
                    rev(which(simulation$groups == condition))
                }
            ),
            use.names = FALSE
        )
    )
}

.canonical_routine_frame <- function(x, columns) {
    if (nrow(x) == 0L) {
        rownames(x) <- NULL
        return(x)
    }
    arguments <- c(
        lapply(columns, function(column) x[[column]]),
        list(method = "radix")
    )
    ordered <- x[do.call(order, arguments), , drop = FALSE]
    rownames(ordered) <- NULL
    ordered
}

.canonical_routine_profile <- function(profile) {
    list(
        raw = .canonical_routine_frame(
            profile$raw,
            c("feature", "condition")
        ),
        grid = profile$grid,
        metadata = profile$metadata
    )
}

canonical_scientific_routine <- function(result) {
    condition_names <- sort(names(result$groups), method = "radix")
    profiles <- lapply(
        result$profiles[condition_names],
        .canonical_routine_profile
    )
    groups <- lapply(
        result$groups[condition_names],
        function(condition_groups) {
            lapply(condition_groups, sort, method = "radix")
        }
    )
    canonical <- list(
        classifications = .canonical_routine_frame(
            result$classifications,
            c("feature", "condition")
        ),
        groups = groups,
        feature_status = .canonical_routine_frame(
            result$feature_status,
            "feature"
        ),
        cutoffs = result$cutoffs[condition_names],
        cutoff_diagnostics = result$cutoff_diagnostics[condition_names],
        profiles = profiles,
        seed_log = .canonical_routine_frame(
            result$seed_log,
            c("feature", "condition")
        ),
        groups_by_sample = result$groups_by_sample[
            sort(names(result$groups_by_sample), method = "radix")
        ]
    )
    canonical$data <- result$data[
        sort(rownames(result$data), method = "radix"),
        sort(colnames(result$data), method = "radix"),
        drop = FALSE
    ]
    canonical
}
