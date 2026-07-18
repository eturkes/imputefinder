#!/usr/bin/env Rscript

# M12b deterministic Tier 1-2 generator + stable-v1 expansion harness.
# Candidate A-C estimators never run here. External result artifacts stay closed.
# Run from the repository root:
# Rscript --vanilla dev/m12-generator-validation.R --verify
# Rscript --vanilla dev/m12-generator-validation.R --v1

.m12b_contract <- local({
    path <- "dev/m12-validation-contract.R"
    if (!file.exists(path)) {
        stop("Run the M12b harness from the repository root.", call. = FALSE)
    }
    environment <- new.env(parent = baseenv())
    sys.source(path, envir = environment)
    environment
})

.M12B_PROTOCOL_ID <- "m12_generator_protocol_v1"
.M12B_DEVELOPMENT_REPLICATES <- 64L
.M12B_AUDIT_REPLICATES <- c(1L, 64L)
.M12B_RESCUE_SEEDS <- c(1L, 7L, 29L)

m12b_acquisition_profiles <- function() {
    data.frame(
        acquisition = c("DDA", "DIA"),
        profile_id = c("synthetic_dda_v1", "synthetic_dia_v1"),
        intensity_center = c(13.0, 13.2),
        feature_sd = c(2.2, 2.0),
        detection_midpoint = c(12.0, 11.8),
        detection_slope = c(1.25, 1.50),
        measurement_sd = c(0.22, 0.16),
        stringsAsFactors = FALSE
    )
}

m12b_scenario_parameters <- function() {
    data.frame(
        scenario_id = c(
            "dda_null_balanced", "dia_null_unequal",
            "dda_monotone_unequal", "dia_monotone_paired",
            "dda_mixed_outlier", "dia_blockwise_paired",
            "dda_batch_crossed", "dia_batch_partial",
            "dda_batch_perfect", "dia_structural_off",
            "dda_no_cliff", "dia_no_cliff_low_support",
            "dda_grouped_leakage_trap"
        ),
        alternative_fraction = c(
            0, 0, 0.24, 0.24, 0.20, 0.20, 0, 0.20, 0.20, 0.12,
            0, 0, 0.20
        ),
        condition_effect = c(
            0, 0, 1.20, 1.00, 1.00, 0.80, 0, 0.80, 0.90, 0.80,
            0, 0, 0.80
        ),
        module_size = rep(25L, 13L),
        module_sd = c(
            0.25, 0.22, 0.18, 0.18, 0.22, 0.20, 0.20, 0.18,
            0.20, 0.18, 0.20, 0.18, 0.28
        ),
        subject_sd = c(
            0, 0, 0, 0.45, 0, 0.50, 0, 0, 0, 0, 0, 0, 0.70
        ),
        run_order_effect = c(
            0, 0, 0.20, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0
        ),
        batch_intensity_effect = c(
            0, 0, 0, 0, 0, 0, 0.35, 0.30, 0.35, 0, 0, 0, 0
        ),
        independent_missing = c(
            0.12, 0.08, 0.03, 0.02, 0.08, 0, 0, 0.04, 0, 0,
            0.18, 0.35, 0.06
        ),
        monotone_max = c(
            0, 0, 0.82, 0.72, 0.82, 0, 0, 0.70, 0, 0, 0, 0, 0
        ),
        batch_dropout_max = c(
            0, 0, 0, 0, 0, 0, 0.72, 0.60, 0.75, 0, 0, 0, 0
        ),
        block_dropout_max = c(
            0, 0, 0, 0, 0, 0.55, 0, 0, 0, 0, 0, 0, 0.65
        ),
        structural_fraction = c(
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0.12, 0, 0, 0
        ),
        heteroscedasticity = c(
            0, 0.35, 0.35, 0, 0.35, 0, 0, 0, 0, 0, 0, 0.35, 0
        ),
        outlier_fraction = c(
            0, 0, 0, 0, 0.004, 0, 0, 0, 0, 0, 0, 0, 0
        ),
        outlier_shift = c(
            0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0
        ),
        stringsAsFactors = FALSE
    )
}

.m12b_seed_from_id <- function(seed_id) {
    digest <- unname(tools::sha256sum(
        bytes = charToRaw(enc2utf8(paste0("m12_sim_seed_stream_v1\r", seed_id)))
    ))
    as.integer(strtoi(substr(digest, 1L, 7L), base = 16L) + 1L)
}

m12b_seed_manifest <- function() {
    scenario_ids <- .m12b_contract$.m12_generator_manifest()$scenario_id
    rows <- do.call(
        rbind,
        lapply(scenario_ids, function(scenario_id) {
            replicate <- seq_len(.M12B_DEVELOPMENT_REPLICATES)
            seed_id <- sprintf("%s__r%03d", scenario_id, replicate)
            data.frame(
                seed_id = seed_id,
                scenario_id = scenario_id,
                replicate = as.integer(replicate),
                simulation_seed = vapply(
                    seed_id,
                    .m12b_seed_from_id,
                    integer(1L)
                ),
                rescue_seed = rep(
                    .M12B_RESCUE_SEEDS,
                    length.out = length(replicate)
                ),
                role = rep("development", length(replicate)),
                stringsAsFactors = FALSE,
                row.names = NULL
            )
        })
    )
    row.names(rows) <- NULL
    rows
}

m12b_generator_protocol <- function() {
    list(
        descriptor = data.frame(
            protocol_id = .M12B_PROTOCOL_ID,
            generator_version = "m12_generator_v1",
            seed_stream = "m12_sim_seed_stream_v1",
            rng_kind = "Mersenne-Twister",
            normal_kind = "Inversion",
            sample_kind = "Rejection",
            development_replicates = .M12B_DEVELOPMENT_REPLICATES,
            audit_replicates = paste(.M12B_AUDIT_REPLICATES, collapse = ";"),
            causal_scope = "generator_specific_latent_and_realized_masks",
            acquisition_scope = "separate_synthetic_parameter_strata_not_physics",
            stringsAsFactors = FALSE
        ),
        acquisition_profiles = m12b_acquisition_profiles(),
        scenario_parameters = m12b_scenario_parameters(),
        seed_manifest = m12b_seed_manifest()
    )
}

.m12b_fail <- function(...) {
    stop(..., call. = FALSE)
}

.m12b_tokens <- function(value) {
    strsplit(value, ";", fixed = TRUE)[[1L]]
}

.m12b_validate_protocol <- function(protocol = m12b_generator_protocol()) {
    expected <- c(
        "descriptor", "acquisition_profiles", "scenario_parameters",
        "seed_manifest"
    )
    if (!identical(names(protocol), expected)) {
        .m12b_fail("M12b protocol components/order differ from v1.")
    }
    descriptor <- protocol$descriptor
    if (nrow(descriptor) != 1L ||
        !identical(descriptor$protocol_id, .M12B_PROTOCOL_ID) ||
        !is.integer(descriptor$development_replicates) ||
        descriptor$development_replicates != .M12B_DEVELOPMENT_REPLICATES) {
        .m12b_fail("M12b protocol descriptor is malformed.")
    }

    profiles <- protocol$acquisition_profiles
    expected_profile_fields <- c(
        "acquisition", "profile_id", "intensity_center", "feature_sd",
        "detection_midpoint", "detection_slope", "measurement_sd"
    )
    numeric_profile_fields <- expected_profile_fields[3:7]
    if (!identical(names(profiles), expected_profile_fields) ||
        !identical(profiles$acquisition, c("DDA", "DIA")) ||
        any(!vapply(profiles[numeric_profile_fields], is.numeric, logical(1L))) ||
        any(!is.finite(as.matrix(profiles[numeric_profile_fields]))) ||
        any(profiles$feature_sd <= 0) || any(profiles$detection_slope <= 0) ||
        any(profiles$measurement_sd <= 0)) {
        .m12b_fail("M12b acquisition profiles are malformed.")
    }

    parameters <- protocol$scenario_parameters
    manifest <- .m12b_contract$.m12_generator_manifest()
    expected_parameter_fields <- c(
        "scenario_id", "alternative_fraction", "condition_effect",
        "module_size", "module_sd", "subject_sd", "run_order_effect",
        "batch_intensity_effect", "independent_missing", "monotone_max",
        "batch_dropout_max", "block_dropout_max", "structural_fraction",
        "heteroscedasticity", "outlier_fraction", "outlier_shift"
    )
    if (!identical(names(parameters), expected_parameter_fields) ||
        !identical(parameters$scenario_id, manifest$scenario_id) ||
        !is.integer(parameters$module_size) || any(parameters$module_size <= 0L)) {
        .m12b_fail("M12b scenario parameter schema/order differs from manifest.")
    }
    numeric_fields <- setdiff(names(parameters), c("scenario_id", "module_size"))
    if (any(!vapply(parameters[numeric_fields], is.numeric, logical(1L))) ||
        any(!is.finite(as.matrix(parameters[numeric_fields]))) ||
        any(parameters$alternative_fraction < 0 |
            parameters$alternative_fraction >= 0.5) ||
        any(parameters$condition_effect < 0) ||
        any(parameters$module_sd < 0) || any(parameters$subject_sd < 0) ||
        any(parameters$run_order_effect < 0) ||
        any(parameters$batch_intensity_effect < 0) ||
        any(parameters$independent_missing < 0 |
            parameters$independent_missing >= 1) ||
        any(parameters$monotone_max < 0 | parameters$monotone_max >= 1) ||
        any(parameters$batch_dropout_max < 0 |
            parameters$batch_dropout_max >= 1) ||
        any(parameters$block_dropout_max < 0 |
            parameters$block_dropout_max >= 1) ||
        any(parameters$structural_fraction < 0 |
            parameters$structural_fraction >= 0.5) ||
        any(parameters$heteroscedasticity < 0) ||
        any(parameters$outlier_fraction < 0 | parameters$outlier_fraction >= 0.1) ||
        any(parameters$outlier_shift < 0)) {
        .m12b_fail("M12b scenario parameters leave their numeric rails.")
    }

    has_token <- function(column, token) {
        vapply(
            column,
            function(value) token %in% .m12b_tokens(value),
            logical(1L),
            USE.NAMES = FALSE
        )
    }
    patterns <- manifest$missingness_patterns
    stresses <- manifest$intensity_stress
    checks <- c(
        monotone = identical(
            parameters$monotone_max > 0,
            has_token(patterns, "monotone_abundance") |
                has_token(patterns, "mixed")
        ),
        batch = identical(
            parameters$batch_dropout_max > 0,
            has_token(patterns, "batch_associated")
        ),
        block = identical(
            parameters$block_dropout_max > 0,
            has_token(patterns, "blockwise")
        ),
        structural = identical(
            parameters$structural_fraction > 0,
            has_token(patterns, "structural_off_compatible")
        ),
        heteroscedastic = identical(
            parameters$heteroscedasticity > 0,
            has_token(stresses, "heteroscedastic")
        ),
        outlier = identical(
            parameters$outlier_fraction > 0,
            has_token(stresses, "outliers")
        ),
        null = all(parameters$alternative_fraction[
            manifest$negative_controls == "association_null"
        ] == 0)
    )
    if (!all(checks)) {
        .m12b_fail(
            "M12b parameter/manifest mismatches: ",
            paste(names(checks)[!checks], collapse = ", ")
        )
    }

    seeds <- protocol$seed_manifest
    expected_seed_fields <- c(
        "seed_id", "scenario_id", "replicate", "simulation_seed",
        "rescue_seed", "role"
    )
    expected_rows <- nrow(manifest) * .M12B_DEVELOPMENT_REPLICATES
    if (!identical(names(seeds), expected_seed_fields) ||
        nrow(seeds) != expected_rows || anyDuplicated(seeds$seed_id) ||
        anyDuplicated(seeds$simulation_seed) ||
        !is.integer(seeds$replicate) || !is.integer(seeds$simulation_seed) ||
        !is.integer(seeds$rescue_seed) || any(seeds$simulation_seed <= 0L) ||
        !all(seeds$rescue_seed %in% .M12B_RESCUE_SEEDS) ||
        !all(seeds$role == "development") ||
        !all(table(seeds$scenario_id) == .M12B_DEVELOPMENT_REPLICATES) ||
        !identical(unique(seeds$scenario_id), manifest$scenario_id)) {
        .m12b_fail("M12b seed manifest is malformed or non-unique.")
    }
    invisible(TRUE)
}

.m12b_with_seed <- function(seed, code) {
    old_kind <- RNGkind()
    caller_has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (caller_has_seed) {
        caller_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    }
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (caller_has_seed) {
            assign(".Random.seed", caller_seed, envir = globalenv())
        } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
        }
    }, add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    force(code)
}

.m12b_parse_sizes <- function(value) {
    pieces <- strsplit(value, ";", fixed = TRUE)[[1L]]
    split <- strsplit(pieces, "=", fixed = TRUE)
    stats::setNames(
        as.integer(vapply(split, `[[`, character(1L), 2L)),
        vapply(split, `[[`, character(1L), 1L)
    )
}

.m12b_make_design <- function(manifest) {
    sizes <- .m12b_parse_sizes(manifest$condition_sizes)
    conditions <- names(sizes)
    technical <- identical(manifest$nuisance_role, "technical_replicate")
    if (identical(manifest$sampling_design, "paired")) {
        technical_replicates <- if (technical) 2L else 1L
        if (any(sizes %% technical_replicates != 0L) ||
            length(unique(sizes / technical_replicates)) != 1L) {
            .m12b_fail("Paired M12b sizes must define complete subject blocks.")
        }
        subject_count <- sizes[[1L]] / technical_replicates
        rows <- do.call(rbind, lapply(conditions, function(condition) {
            expand <- expand.grid(
                subject_index = seq_len(subject_count),
                technical_index = seq_len(technical_replicates),
                KEEP.OUT.ATTRS = FALSE,
                stringsAsFactors = FALSE
            )
            expand$condition <- condition
            expand
        }))
        rows$subject <- sprintf("subject_%02d", rows$subject_index)
        rows$sample <- if (technical) {
            sprintf(
                "%s_%s_tech_%02d",
                rows$condition,
                rows$subject,
                rows$technical_index
            )
        } else {
            sprintf("%s_%s", rows$condition, rows$subject)
        }
    } else {
        rows <- do.call(rbind, lapply(conditions, function(condition) {
            count <- sizes[[condition]]
            data.frame(
                subject_index = seq_len(count),
                technical_index = rep(1L, count),
                condition = rep(condition, count),
                subject = rep(NA_character_, count),
                sample = sprintf("%s_%02d", condition, seq_len(count)),
                stringsAsFactors = FALSE
            )
        }))
    }
    row.names(rows) <- NULL

    condition_index <- match(rows$condition, conditions)
    within_condition <- ave(
        seq_len(nrow(rows)),
        rows$condition,
        FUN = seq_along
    )
    rows$run_order <- rank(
        (within_condition - 1L) * length(conditions) + condition_index,
        ties.method = "first"
    )
    rows$batch <- "batch_1"
    if (identical(manifest$nuisance_role, "batch_crossed")) {
        rows$batch <- paste0("batch_", 1L + (within_condition - 1L) %% 2L)
    } else if (identical(manifest$nuisance_role, "batch_partial")) {
        rows$batch <- unlist(lapply(seq_along(conditions), function(index) {
            count <- sizes[[index]]
            if (index == 1L) {
                c(rep("batch_1", ceiling(0.75 * count)),
                  rep("batch_2", count - ceiling(0.75 * count)))
            } else if (index == length(conditions)) {
                c(rep("batch_2", max(1L, floor(0.25 * count))),
                  rep("batch_3", count - max(1L, floor(0.25 * count))))
            } else {
                rep(c("batch_1", "batch_2", "batch_3"), length.out = count)
            }
        }), use.names = FALSE)
    } else if (identical(manifest$nuisance_role, "batch_aliased")) {
        rows$batch <- paste0("batch_", condition_index)
    }
    rows$technical_replicate <- if (technical) {
        sprintf("tech_%02d", rows$technical_index)
    } else {
        "none"
    }
    rows$resampling_unit <- if (identical(manifest$sampling_design, "paired")) {
        rows$subject
    } else {
        rows$sample
    }
    rows$acquisition_mode <- manifest$acquisition
    out <- rows[c(
        "sample", "condition", "batch", "run_order", "subject",
        "technical_replicate", "resampling_unit", "acquisition_mode"
    )]
    row.names(out) <- out$sample
    out
}

.m12b_component_draw <- function(probability) {
    draw <- matrix(
        stats::runif(length(probability)),
        nrow = nrow(probability),
        dimnames = dimnames(probability)
    )
    draw < probability
}

.m12b_group_dropout <- function(
    complete_data,
    groups,
    maximum,
    profile,
    group_weights = NULL
) {
    probability <- matrix(0, nrow(complete_data), ncol(complete_data),
                          dimnames = dimnames(complete_data))
    mask <- probability > 0
    if (maximum <= 0) {
        return(list(probability = probability, mask = mask))
    }
    levels <- sort(unique(groups), method = "radix")
    if (is.null(group_weights)) {
        group_weights <- rep(1, length(levels))
    }
    for (index in seq_along(levels)) {
        columns <- groups == levels[[index]]
        group_mean <- rowMeans(complete_data[, columns, drop = FALSE])
        group_probability <- maximum * group_weights[[index]] * stats::plogis(
            (profile$detection_midpoint - group_mean) *
                profile$detection_slope
        )
        group_probability <- pmin(group_probability, 0.95)
        group_draw <- stats::runif(nrow(complete_data)) < group_probability
        probability[, columns] <- group_probability
        mask[, columns] <- group_draw
    }
    list(probability = probability, mask = mask)
}

.m12b_rule_oracle <- function(data, design, cutoffs) {
    features <- rownames(data)
    conditions <- sort(unique(design$condition), method = "radix")
    rows <- lapply(conditions, function(condition) {
        columns <- design$condition == condition
        block <- data[, columns, drop = FALSE]
        observed <- rowSums(!is.na(block))
        missing <- ncol(block) - observed
        latent_mean <- rowMeans(block, na.rm = TRUE)
        latent_mean[observed == 0L] <- -Inf
        state <- rep("insufficient", length(features))
        state[missing == 0L] <- "complete"
        state[missing > 0L & latent_mean < cutoffs[[condition]]] <- "MNAR"
        mar <- missing > 0L & latent_mean >= cutoffs[[condition]] &
            observed > ncol(block) / 2
        state[mar] <- "MAR"
        data.frame(
            feature = features,
            condition = condition,
            observed_count = as.integer(observed),
            missing_count = as.integer(missing),
            state = state,
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out <- out[order(match(out$feature, features), out$condition), , drop = FALSE]
    row.names(out) <- NULL
    out
}

simulate_m12b_scenario <- function(
    scenario_id,
    replicate = 1L,
    protocol = m12b_generator_protocol()
) {
    .m12b_validate_protocol(protocol)
    if (!is.character(scenario_id) || length(scenario_id) != 1L ||
        is.na(scenario_id) || !nzchar(scenario_id) ||
        !is.integer(replicate) || length(replicate) != 1L ||
        is.na(replicate) || replicate < 1L ||
        replicate > .M12B_DEVELOPMENT_REPLICATES) {
        .m12b_fail("M12b scenario ID/replicate is invalid.")
    }
    manifest <- .m12b_contract$.m12_generator_manifest()
    scenario_index <- match(scenario_id, manifest$scenario_id)
    if (is.na(scenario_index)) {
        .m12b_fail("Unknown M12b scenario: ", scenario_id)
    }
    manifest <- manifest[scenario_index, , drop = FALSE]
    parameters <- protocol$scenario_parameters[scenario_index, , drop = FALSE]
    profile <- protocol$acquisition_profiles[
        protocol$acquisition_profiles$acquisition == manifest$acquisition,
        ,
        drop = FALSE
    ]
    seed_row <- protocol$seed_manifest[
        protocol$seed_manifest$scenario_id == scenario_id &
            protocol$seed_manifest$replicate == replicate,
        ,
        drop = FALSE
    ]
    if (nrow(seed_row) != 1L) {
        .m12b_fail("M12b seed lookup is not one-to-one.")
    }

    .m12b_with_seed(seed_row$simulation_seed, {
        design <- .m12b_make_design(manifest)
        n_features <- manifest$n_features
        n_samples <- nrow(design)
        features <- sprintf("protein_%05d", seq_len(n_features))
        samples <- design$sample
        conditions <- sort(unique(design$condition), method = "radix")

        quantiles <- (seq_len(n_features) - 0.5) / n_features
        baseline <- stats::qnorm(
            quantiles,
            mean = profile$intensity_center,
            sd = profile$feature_sd
        )
        baseline <- sample(baseline, n_features, replace = FALSE)
        modules <- rep(
            seq_len(ceiling(n_features / parameters$module_size)),
            each = parameters$module_size,
            length.out = n_features
        )
        modules <- sample(modules, n_features, replace = FALSE)

        available <- seq_len(n_features)
        structural_count <- as.integer(floor(
            n_features * parameters$structural_fraction
        ))
        structural_features <- if (structural_count) {
            sample(available, structural_count, replace = FALSE)
        } else {
            integer()
        }
        available <- setdiff(available, structural_features)
        alternative_count <- as.integer(floor(
            n_features * parameters$alternative_fraction
        ))
        alternative_features <- if (alternative_count) {
            sample(available, alternative_count, replace = FALSE)
        } else {
            integer()
        }
        effect_sign <- rep(0, n_features)
        if (alternative_count) {
            signs <- rep(c(-1, 1), length.out = alternative_count)
            effect_sign[alternative_features] <- sample(
                signs,
                alternative_count,
                replace = FALSE
            )
        }
        condition_scores <- seq(-1, 1, length.out = length(conditions))
        condition_effect <- matrix(
            0,
            nrow = n_features,
            ncol = length(conditions),
            dimnames = list(features, conditions)
        )
        if (alternative_count) {
            condition_effect[alternative_features, ] <-
                effect_sign[alternative_features] *
                parameters$condition_effect *
                rep(condition_scores, each = alternative_count)
        }

        structural_condition <- rep(NA_character_, n_features)
        if (structural_count) {
            structural_condition[structural_features] <- rep(
                conditions,
                length.out = structural_count
            )
            for (condition in conditions) {
                indices <- which(structural_condition == condition)
                condition_effect[indices, condition] <-
                    profile$detection_midpoint - 4 - baseline[indices]
            }
        }
        condition_mean <- baseline + condition_effect

        sample_condition <- match(design$condition, conditions)
        complete_data <- condition_mean[, sample_condition, drop = FALSE]
        colnames(complete_data) <- samples

        sample_effect <- numeric(n_samples)
        subjects <- sort(unique(stats::na.omit(design$subject)), method = "radix")
        if (length(subjects) && parameters$subject_sd > 0) {
            subject_effect <- stats::setNames(
                stats::rnorm(length(subjects), sd = parameters$subject_sd),
                subjects
            )
            sample_effect <- sample_effect + unname(subject_effect[design$subject])
        }
        batches <- sort(unique(design$batch), method = "radix")
        if (length(batches) > 1L && parameters$batch_intensity_effect > 0) {
            batch_effect <- stats::setNames(
                seq(-1, 1, length.out = length(batches)) *
                    parameters$batch_intensity_effect,
                batches
            )
            sample_effect <- sample_effect + unname(batch_effect[design$batch])
        }
        if (parameters$run_order_effect > 0) {
            centered_run <- if (n_samples > 1L) {
                2 * (design$run_order - 1) / (n_samples - 1) - 1
            } else {
                0
            }
            sample_effect <- sample_effect +
                centered_run * parameters$run_order_effect
        }
        complete_data <- sweep(complete_data, 2L, sample_effect, "+")

        module_noise <- matrix(
            stats::rnorm(
                max(modules) * n_samples,
                sd = parameters$module_sd
            ),
            nrow = max(modules),
            ncol = n_samples
        )
        complete_data <- complete_data + module_noise[modules, , drop = FALSE]
        residual_scale <- profile$measurement_sd * exp(
            parameters$heteroscedasticity *
                (profile$intensity_center - baseline) / 2
        )
        residual_scale <- pmin(pmax(residual_scale, profile$measurement_sd / 2),
                               profile$measurement_sd * 3)
        residual <- matrix(
            stats::rnorm(n_features * n_samples),
            nrow = n_features,
            ncol = n_samples
        ) * residual_scale
        complete_data <- complete_data + residual
        dimnames(complete_data) <- list(features, samples)

        outlier_mask <- matrix(
            FALSE,
            nrow = n_features,
            ncol = n_samples,
            dimnames = dimnames(complete_data)
        )
        outlier_count <- as.integer(floor(
            length(complete_data) * parameters$outlier_fraction
        ))
        if (outlier_count) {
            cells <- sample(seq_along(complete_data), outlier_count, replace = FALSE)
            signs <- sample(rep(c(-1, 1), length.out = outlier_count))
            complete_data[cells] <- complete_data[cells] +
                signs * parameters$outlier_shift
            outlier_mask[cells] <- TRUE
        }

        zero_probability <- matrix(
            0,
            nrow = n_features,
            ncol = n_samples,
            dimnames = dimnames(complete_data)
        )
        monotone_probability <- zero_probability
        if (parameters$monotone_max > 0) {
            monotone_probability <- parameters$monotone_max * stats::plogis(
                (profile$detection_midpoint - complete_data) *
                    profile$detection_slope
            )
        }
        monotone_mask <- .m12b_component_draw(monotone_probability)

        independent_probability <- zero_probability +
            parameters$independent_missing
        independent_mask <- .m12b_component_draw(independent_probability)

        batch_weights <- if (length(batches) > 1L) {
            seq(0.35, 1, length.out = length(batches))
        } else {
            1
        }
        batch <- .m12b_group_dropout(
            complete_data,
            design$batch,
            parameters$batch_dropout_max,
            profile,
            batch_weights
        )
        block <- .m12b_group_dropout(
            complete_data,
            ifelse(is.na(design$subject), design$sample, design$subject),
            parameters$block_dropout_max,
            profile
        )

        structural_mask <- zero_probability > 0
        if (structural_count) {
            for (condition in conditions) {
                structural_mask[
                    structural_condition == condition,
                    design$condition == condition
                ] <- TRUE
            }
        }
        missing_mask <- monotone_mask | independent_mask | batch$mask |
            block$mask | structural_mask
        missing_probability <- 1 -
            (1 - monotone_probability) *
            (1 - independent_probability) *
            (1 - batch$probability) *
            (1 - block$probability)
        missing_probability[structural_mask] <- 1
        data <- complete_data
        data[missing_mask] <- NA_real_

        marginal_detection <- do.call(rbind, lapply(conditions, function(condition) {
            columns <- design$condition == condition
            data.frame(
                feature = features,
                condition = condition,
                direct_condition_effect = condition_effect[, condition],
                latent_condition_mean = condition_mean[, condition],
                marginal_detection_probability = rowMeans(
                    1 - missing_probability[, columns, drop = FALSE]
                ),
                stringsAsFactors = FALSE
            )
        }))
        marginal_detection <- marginal_detection[
            order(match(marginal_detection$feature, features),
                  marginal_detection$condition),
            ,
            drop = FALSE
        ]
        row.names(marginal_detection) <- NULL

        feature_truth <- data.frame(
            feature = features,
            baseline_intensity = baseline,
            module = as.integer(modules),
            alternative = seq_len(n_features) %in% alternative_features,
            effect_sign = effect_sign,
            structural_off_compatible_condition = structural_condition,
            stringsAsFactors = FALSE
        )
        model_data <- data.frame(
            condition = factor(design$condition),
            batch = factor(design$batch)
        )
        design_terms <- stats::model.matrix(~ 0 + condition, data = model_data)
        if (nlevels(model_data$batch) > 1L) {
            batch_terms <- stats::model.matrix(~ batch, data = model_data)[
                ,
                -1L,
                drop = FALSE
            ]
            design_terms <- cbind(design_terms, batch_terms)
        }
        cutoff <- stats::setNames(
            rep(profile$detection_midpoint, length(conditions)),
            conditions
        )

        simulation <- list(
            protocol_id = .M12B_PROTOCOL_ID,
            scenario_id = scenario_id,
            replicate = replicate,
            simulation_seed = seed_row$simulation_seed,
            rescue_seed = seed_row$rescue_seed,
            manifest = manifest,
            parameters = parameters,
            acquisition_profile = profile,
            design = design,
            complete_data = complete_data,
            data = data,
            masks = list(
                monotone = monotone_mask,
                intensity_independent = independent_mask,
                batch = batch$mask,
                block = block$mask,
                structural_off_compatible = structural_mask,
                missing = missing_mask,
                outlier = outlier_mask
            ),
            missing_probability = missing_probability,
            feature_truth = feature_truth,
            condition_truth = marginal_detection,
            design_truth = list(
                model_columns = ncol(design_terms),
                model_rank = qr(design_terms)$rank,
                condition_nuisance_separable = manifest$confounding != "perfect",
                resampling_unit = if (manifest$block_role == "subject") {
                    "subject"
                } else {
                    "sample"
                }
            ),
            manual_cutoffs = cutoff,
            rule_oracle = .m12b_rule_oracle(data, design, cutoff)
        )
        .m12b_validate_simulation(simulation)
        simulation
    })
}

.m12b_validate_simulation <- function(simulation) {
    required <- c(
        "protocol_id", "scenario_id", "replicate", "simulation_seed",
        "rescue_seed", "manifest", "parameters", "acquisition_profile",
        "design", "complete_data", "data", "masks", "missing_probability",
        "feature_truth", "condition_truth", "design_truth", "manual_cutoffs",
        "rule_oracle"
    )
    if (!is.list(simulation) || !identical(names(simulation), required)) {
        .m12b_fail("M12b simulation schema/order differs from v1.")
    }
    data <- simulation$data
    complete <- simulation$complete_data
    design <- simulation$design
    masks <- simulation$masks
    manifest <- simulation$manifest
    parameters <- simulation$parameters
    expected_masks <- c(
        "monotone", "intensity_independent", "batch", "block",
        "structural_off_compatible", "missing", "outlier"
    )
    if (!is.matrix(data) || !is.numeric(data) ||
        !identical(dim(data), c(manifest$n_features, nrow(design))) ||
        !identical(dimnames(data), dimnames(complete)) ||
        !identical(colnames(data), row.names(design)) ||
        any(!is.finite(complete)) || any(is.nan(data)) ||
        any(is.infinite(data)) || !identical(names(masks), expected_masks) ||
        any(!vapply(masks, is.logical, logical(1L))) ||
        any(!vapply(masks, function(mask) identical(dim(mask), dim(data)),
                    logical(1L)))) {
        .m12b_fail("M12b generated matrix/mask schema is malformed.")
    }
    union_mask <- masks$monotone | masks$intensity_independent | masks$batch |
        masks$block | masks$structural_off_compatible
    if (!identical(union_mask, masks$missing) ||
        !identical(is.na(data), masks$missing) ||
        !identical(data[!masks$missing], complete[!masks$missing]) ||
        any(!is.finite(simulation$missing_probability)) ||
        any(simulation$missing_probability < 0 |
            simulation$missing_probability > 1)) {
        .m12b_fail("M12b realized masks/probabilities disagree.")
    }
    condition_counts <- table(design$condition)
    if (!identical(
        as.integer(condition_counts[names(.m12b_parse_sizes(
            manifest$condition_sizes
        ))]),
        unname(.m12b_parse_sizes(manifest$condition_sizes))
    ) || any(colSums(!is.na(data)) == 0L)) {
        .m12b_fail("M12b design counts or finite sample support disagree.")
    }
    if (manifest$sampling_design == "paired") {
        coverage <- table(design$subject, design$condition)
        if (any(coverage == 0L) || length(unique(as.vector(coverage))) != 1L ||
            !identical(design$resampling_unit, design$subject)) {
            .m12b_fail("M12b paired/block design is incomplete.")
        }
    }
    if (manifest$confounding == "perfect") {
        condition_batch <- table(design$condition, design$batch) > 0L
        if (any(rowSums(condition_batch) != 1L) ||
            any(colSums(condition_batch) != 1L) ||
            simulation$design_truth$model_rank >=
                simulation$design_truth$model_columns ||
            simulation$design_truth$condition_nuisance_separable) {
            .m12b_fail("M12b perfect-confounding truth is not constructed.")
        }
    }
    if (manifest$confounding == "partial") {
        condition_batch <- table(design$condition, design$batch) > 0L
        if (!any(!condition_batch) || any(rowSums(condition_batch) < 2L) ||
            simulation$design_truth$model_rank !=
                simulation$design_truth$model_columns) {
            .m12b_fail("M12b partial-confounding truth is not constructed.")
        }
    }
    if (manifest$nuisance_role == "batch_crossed" &&
        any(table(design$condition, design$batch) == 0L)) {
        .m12b_fail("M12b crossed batch design has an empty cell.")
    }
    if (manifest$negative_controls == "association_null" &&
        any(simulation$feature_truth$alternative)) {
        .m12b_fail("M12b null control contains condition alternatives.")
    }
    if ("no_cliff" %in% .m12b_tokens(manifest$missingness_patterns) &&
        (parameters$monotone_max > 0 || parameters$batch_dropout_max > 0 ||
         parameters$block_dropout_max > 0 ||
         parameters$structural_fraction > 0)) {
        .m12b_fail("M12b no-cliff control contains a structured cliff.")
    }
    structural <- simulation$feature_truth$structural_off_compatible_condition
    if (any(!is.na(structural))) {
        for (condition in unique(stats::na.omit(structural))) {
            if (!all(masks$structural_off_compatible[
                which(structural == condition),
                design$condition == condition,
                drop = FALSE
            ])) {
                .m12b_fail("M12b structural-off-compatible blocks are incomplete.")
            }
        }
    }
    if (manifest$negative_controls == "leakage_trap" &&
        (length(unique(design$resampling_unit)) >= nrow(design) ||
         !anyDuplicated(design$resampling_unit))) {
        .m12b_fail("M12b leakage trap lacks grouped technical siblings.")
    }
    invisible(TRUE)
}

.m12b_hash_object <- function(value) {
    unname(tools::sha256sum(bytes = serialize(value, NULL, version = 3L)))
}

.m12b_hash_frames <- function(frames) {
    vapply(frames, .m12b_contract$.m12_hash_frame, character(1L))
}

m12b_protocol_hashes <- function(protocol = m12b_generator_protocol()) {
    .m12b_validate_protocol(protocol)
    hashes <- .m12b_hash_frames(protocol)
    aggregate <- paste(names(hashes), hashes, sep = "\t", collapse = "\n")
    c(
        hashes,
        protocol = unname(tools::sha256sum(
            bytes = charToRaw(enc2utf8(aggregate))
        ))
    )
}

m12b_generator_audit <- function() {
    protocol <- m12b_generator_protocol()
    scenario_ids <- protocol$scenario_parameters$scenario_id
    audit_grid <- expand.grid(
        scenario_id = scenario_ids,
        replicate = .M12B_AUDIT_REPLICATES,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    audit_grid <- audit_grid[order(
        match(audit_grid$scenario_id, scenario_ids),
        audit_grid$replicate
    ), , drop = FALSE]
    rows <- lapply(seq_len(nrow(audit_grid)), function(index) {
        scenario_id <- audit_grid$scenario_id[[index]]
        replicate <- audit_grid$replicate[[index]]
        caller_kind <- RNGkind()
        caller_has_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
        if (caller_has_seed) {
            caller_seed <- get(".Random.seed", globalenv(), inherits = FALSE)
        }
        simulation <- simulate_m12b_scenario(
            scenario_id,
            replicate = as.integer(replicate),
            protocol = protocol
        )
        if (!identical(RNGkind(), caller_kind) ||
            !identical(
                exists(".Random.seed", globalenv(), inherits = FALSE),
                caller_has_seed
            ) || (caller_has_seed && !identical(
                get(".Random.seed", globalenv(), inherits = FALSE),
                caller_seed
            ))) {
            .m12b_fail("M12b generation changed caller RNG state.")
        }
        masks <- simulation$masks
        data.frame(
            scenario_id = scenario_id,
            replicate = as.integer(replicate),
            simulation_seed = simulation$simulation_seed,
            n_features = nrow(simulation$data),
            n_samples = ncol(simulation$data),
            missing_cells = sum(masks$missing),
            missing_fraction = mean(masks$missing),
            monotone_cells = sum(masks$monotone),
            independent_cells = sum(masks$intensity_independent),
            batch_cells = sum(masks$batch),
            block_cells = sum(masks$block),
            structural_cells = sum(masks$structural_off_compatible),
            all_missing_features = sum(rowSums(!is.na(simulation$data)) == 0L),
            simulation_sha256 = .m12b_hash_object(list(
                design = simulation$design,
                complete_data = simulation$complete_data,
                data = simulation$data,
                masks = simulation$masks,
                missing_probability = simulation$missing_probability,
                feature_truth = simulation$feature_truth,
                condition_truth = simulation$condition_truth,
                design_truth = simulation$design_truth,
                rule_oracle = simulation$rule_oracle
            )),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    row.names(out) <- NULL
    out
}

.m12b_without_finite_doubles <- function(value) {
    if (is.double(value)) {
        value[is.finite(value)] <- 0
    } else if (is.list(value)) {
        value[] <- lapply(value, .m12b_without_finite_doubles)
    }
    value
}

.m12b_classic_matrix <- function(data) {
    if (is.matrix(data)) {
        return(data)
    }
    SummarizedExperiment::assay(data, "intensity")
}

.m12b_canonical_classic <- function(result) {
    classifications <- result$classifications
    classifications <- classifications[order(
        classifications$feature,
        classifications$condition,
        method = "radix"
    ), , drop = FALSE]
    row.names(classifications) <- NULL
    classifications <- .m12b_without_finite_doubles(classifications)

    status <- result$feature_status[
        order(result$feature_status$feature, method = "radix"),
        ,
        drop = FALSE
    ]
    row.names(status) <- NULL
    seed_log <- result$seed_log
    seed_log <- seed_log[order(
        seed_log$feature,
        seed_log$condition,
        seed_log$sample,
        method = "radix"
    ), c("feature", "condition", "sample", "seed"), drop = FALSE]
    row.names(seed_log) <- NULL
    groups <- result$groups[sort(names(result$groups), method = "radix")]
    groups <- lapply(groups, function(condition) {
        condition <- condition[sort(names(condition), method = "radix")]
        lapply(condition, sort, method = "radix")
    })
    data <- .m12b_classic_matrix(result$data)
    data <- data[
        sort(rownames(data), method = "radix"),
        sort(colnames(data), method = "radix"),
        drop = FALSE
    ]
    list(
        schema = c(class = class(result), fields = names(result)),
        data = data,
        classifications = classifications,
        groups = groups,
        feature_status = status,
        seed_log = seed_log,
        cutoffs = result$cutoffs[sort(names(result$cutoffs), method = "radix")]
    )
}

.m12b_observed_value_gate <- function(result, original) {
    returned <- .m12b_classic_matrix(result$data)
    original <- original[rownames(returned), colnames(returned), drop = FALSE]
    observed <- is.finite(original)
    if (!identical(returned[observed], original[observed])) {
        return(FALSE)
    }
    changed <- which(is.na(original) & is.finite(returned), arr.ind = TRUE)
    changed_keys <- if (nrow(changed)) {
        sort(paste(
            rownames(returned)[changed[, "row"]],
            colnames(returned)[changed[, "col"]],
            sep = "\r"
        ), method = "radix")
    } else {
        character()
    }
    retained_log <- result$seed_log[
        result$seed_log$feature %in% rownames(returned),
        ,
        drop = FALSE
    ]
    log_keys <- if (nrow(retained_log)) {
        sort(paste(
            retained_log$feature,
            retained_log$sample,
            sep = "\r"
        ), method = "radix")
    } else {
        character()
    }
    identical(changed_keys, log_keys)
}

m12b_v1_audit <- function() {
    if (!requireNamespace("imputefinder", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
        !requireNamespace("S4Vectors", quietly = TRUE)) {
        .m12b_fail("Install the current package and imports before --v1.")
    }
    classifier <- getExportedValue("imputefinder", "classify_missingness")
    scenario_ids <- m12b_scenario_parameters()$scenario_id
    rows <- lapply(scenario_ids, function(scenario_id) {
        simulation <- simulate_m12b_scenario(scenario_id, 1L)
        x <- simulation$data
        original <- x
        group <- stats::setNames(simulation$design$condition, colnames(x))
        cutoff <- simulation$manual_cutoffs
        caller_kind <- RNGkind()
        caller_has_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
        if (caller_has_seed) {
            caller_seed <- get(".Random.seed", globalenv(), inherits = FALSE)
        }
        result <- classifier(
            x,
            group = group,
            cutoffs = cutoff,
            seed = simulation$rescue_seed
        )
        repeated <- classifier(
            x,
            group = group,
            cutoffs = cutoff,
            seed = simulation$rescue_seed
        )
        if (!identical(x, original) || !identical(result, repeated) ||
            !identical(RNGkind(), caller_kind) ||
            !identical(
                exists(".Random.seed", globalenv(), inherits = FALSE),
                caller_has_seed
            ) || (caller_has_seed && !identical(
                get(".Random.seed", globalenv(), inherits = FALSE),
                caller_seed
            )) || !.m12b_observed_value_gate(result, original)) {
            .m12b_fail("Stable v1 repeat/input/RNG gate failed for ", scenario_id)
        }

        row_order <- rev(seq_len(nrow(x)))
        column_order <- rev(seq_len(ncol(x)))
        permuted <- classifier(
            x[row_order, column_order, drop = FALSE],
            group = group[colnames(x)[column_order]],
            cutoffs = cutoff[rev(names(cutoff))],
            seed = simulation$rescue_seed
        )
        canonical <- .m12b_canonical_classic(result)
        if (!identical(
            .m12b_hash_object(canonical),
            .m12b_hash_object(.m12b_canonical_classic(permuted))
        )) {
            .m12b_fail("Stable v1 named-order gate failed for ", scenario_id)
        }

        experiment <- SummarizedExperiment::SummarizedExperiment(
            assays = list(intensity = x),
            colData = S4Vectors::DataFrame(
                condition = simulation$design$condition,
                row.names = colnames(x)
            )
        )
        experiment_original <- experiment
        experiment_result <- classifier(
            experiment,
            group_col = "condition",
            assay = "intensity",
            cutoffs = cutoff,
            seed = simulation$rescue_seed
        )
        if (!identical(experiment, experiment_original) || !identical(
            .m12b_hash_object(canonical),
            .m12b_hash_object(.m12b_canonical_classic(experiment_result))
        )) {
            .m12b_fail("Stable v1 matrix/SE gate failed for ", scenario_id)
        }

        data.frame(
            scenario_id = scenario_id,
            retained_features = nrow(.m12b_classic_matrix(result$data)),
            seeded_blocks = nrow(result$seed_log),
            all_missing_features = sum(
                result$feature_status$drop_reason == "all_missing",
                na.rm = TRUE
            ),
            exact_sha256 = .m12b_hash_object(canonical),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    row.names(out) <- NULL
    out
}

.M12B_EXPECTED_PROTOCOL_HASHES <- c(
    descriptor = "8076ad5a990841a1fcd37d5145b3f01eaae9cb9f9dd2c922712c150bc221d132",
    acquisition_profiles = "d868a3bdf9d2daaf2c9ed482833416cbbd07b83956a205056bc2d371af731fbe",
    scenario_parameters = "019a81a0438958e7080243955e1e6673999a29e82fbe2ed70b5d8fc0d100b86e",
    seed_manifest = "c3f550c378b3d679e44b472ac01178776ef7e3f386ee4d5a43dde928d752de79",
    protocol = "cdea1bf874152e63fba08c49e390da1dd18690105e43a2dc7a03fc6866d0d080"
)
.M12B_EXPECTED_GENERATOR_AUDIT_HASH <-
    "4d216966c8cc54b1b23685d8801c7c65b2afcc5879dd15430d57a88f765e8a00"
.M12B_EXPECTED_V1_AUDIT_HASH <-
    "8463f4306fe9b6192336487e4ca365c8c7e3525dc06aeca5b66b1a99b1a10250"

.m12b_self_tests <- function() {
    protocol <- m12b_generator_protocol()
    bad_scope <- protocol
    bad_scope$scenario_parameters$structural_fraction[[1L]] <- 0.1
    duplicate_seed <- protocol
    duplicate_seed$seed_manifest$simulation_seed[[2L]] <-
        duplicate_seed$seed_manifest$simulation_seed[[1L]]
    bad_pair <- protocol
    bad_pair$scenario_parameters$block_dropout_max[[1L]] <- 0.2
    c(
        pattern_parameter_rail = inherits(
            try(.m12b_validate_protocol(bad_scope), silent = TRUE),
            "try-error"
        ),
        unique_seed_rail = inherits(
            try(.m12b_validate_protocol(duplicate_seed), silent = TRUE),
            "try-error"
        ),
        block_parameter_rail = inherits(
            try(.m12b_validate_protocol(bad_pair), silent = TRUE),
            "try-error"
        )
    )
}

.m12b_main <- function(args) {
    if (length(args) != 1L || !args %in% c("--verify", "--v1")) {
        .m12b_fail(
            "usage: Rscript --vanilla dev/m12-generator-validation.R ",
            "--verify|--v1"
        )
    }
    protocol_hashes <- m12b_protocol_hashes()
    self_tests <- .m12b_self_tests()
    if (!all(self_tests)) {
        .m12b_fail(
            "M12b self-test failures: ",
            paste(names(self_tests)[!self_tests], collapse = ", ")
        )
    }
    generator_audit <- m12b_generator_audit()
    generator_hash <- .m12b_contract$.m12_hash_frame(generator_audit)
    if (!anyNA(.M12B_EXPECTED_PROTOCOL_HASHES) &&
        !identical(protocol_hashes, .M12B_EXPECTED_PROTOCOL_HASHES)) {
        mismatch <- names(protocol_hashes)[
            protocol_hashes != .M12B_EXPECTED_PROTOCOL_HASHES
        ]
        .m12b_fail("M12b protocol hash mismatch: ", paste(mismatch, collapse = ", "))
    }
    if (!is.na(.M12B_EXPECTED_GENERATOR_AUDIT_HASH) &&
        !identical(generator_hash, .M12B_EXPECTED_GENERATOR_AUDIT_HASH)) {
        .m12b_fail("M12b generator audit hash mismatch.")
    }

    cat("protocol_id: ", .M12B_PROTOCOL_ID, "\n", sep = "")
    cat("self_tests: ", paste(names(self_tests), self_tests,
                              sep = "=", collapse = "; "), "\n", sep = "")
    cat(paste(names(protocol_hashes), protocol_hashes, sep = ": "), sep = "\n")
    cat("\ngenerator_audit: ", generator_hash, "\n", sep = "")
    print(generator_audit, row.names = FALSE)

    if (identical(args, "--v1")) {
        v1_audit <- m12b_v1_audit()
        v1_hash <- .m12b_contract$.m12_hash_frame(v1_audit)
        if (!is.na(.M12B_EXPECTED_V1_AUDIT_HASH) &&
            !identical(v1_hash, .M12B_EXPECTED_V1_AUDIT_HASH)) {
            .m12b_fail("M12b stable-v1 audit hash mismatch.")
        }
        cat("v1_audit: ", v1_hash, "\n", sep = "")
        print(v1_audit, row.names = FALSE)
    }
    if (anyNA(.M12B_EXPECTED_PROTOCOL_HASHES) ||
        is.na(.M12B_EXPECTED_GENERATOR_AUDIT_HASH) ||
        (identical(args, "--v1") && is.na(.M12B_EXPECTED_V1_AUDIT_HASH))) {
        cat("hash_state: bootstrap_unfrozen\n")
    }
    invisible(TRUE)
}

if (sys.nframe() == 0L) {
    .m12b_main(commandArgs(trailingOnly = TRUE))
}
