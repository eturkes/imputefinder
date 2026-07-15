.resolve_manual_cutoffs <- function(statistics, cutoffs, conditions = NULL) {
    .validate_cutoff_statistics(statistics)
    if (is.null(conditions)) {
        conditions <- sort(unique(statistics$condition), method = "radix")
    } else {
        conditions <- .normalise_cutoff_conditions(conditions)
    }
    unknown_statistics <- setdiff(unique(statistics$condition), conditions)
    if (length(unknown_statistics) > 0L) {
        stop(
            "Cutoff conditions must include every statistics condition.",
            call. = FALSE
        )
    }
    manual <- .normalise_manual_cutoffs(cutoffs, conditions)

    resolved <- stats::setNames(rep(NA_real_, length(conditions)), conditions)
    diagnostics <- stats::setNames(
        vector("list", length(conditions)),
        conditions
    )

    for (condition_index in seq_along(conditions)) {
        condition <- conditions[[condition_index]]
        condition_statistics <- statistics[
            statistics$condition == condition,
            ,
            drop = FALSE
        ]
        has_missing <- any(condition_statistics$missing_count > 0L)

        if (!has_missing) {
            diagnostics[[condition_index]] <- .not_needed_cutoff_diagnostic()
        } else if (condition %in% names(manual)) {
            cutoff <- unname(manual[[condition]])
            resolved[[condition]] <- cutoff
            diagnostics[[condition_index]] <- .manual_cutoff_diagnostic(
                condition,
                cutoff,
                condition_statistics$mean_intensity
            )
        }
    }

    list(cutoffs = resolved, diagnostics = diagnostics)
}

.resolve_cutoffs <- function(
    statistics,
    profiles,
    cutoffs,
    conditions = NULL
) {
    resolved <- .resolve_manual_cutoffs(statistics, cutoffs, conditions)
    conditions <- names(resolved$cutoffs)
    valid_profiles <- is.list(profiles) &&
        !is.null(names(profiles)) &&
        !anyNA(names(profiles)) &&
        all(nzchar(names(profiles))) &&
        !anyDuplicated(names(profiles)) &&
        setequal(names(profiles), conditions)
    if (!valid_profiles) {
        stop(
            "`profiles` must contain exactly one named profile per condition.",
            call. = FALSE
        )
    }

    for (condition in conditions) {
        if (!is.null(resolved$diagnostics[[condition]])) {
            next
        }
        profile <- profiles[[condition]]
        decision <- .detect_automatic_cutoff(profile)
        diagnostic <- .automatic_cutoff_diagnostic(decision)
        if (!identical(decision$status, "ok")) {
            .abort_automatic_cutoff(
                condition,
                decision,
                diagnostic,
                profile
            )
        }
        resolved$cutoffs[[condition]] <- decision$cutoff
        resolved$diagnostics[[condition]] <- diagnostic
    }

    resolved
}

.normalise_cutoff_conditions <- function(conditions) {
    valid <- is.character(conditions) &&
        is.null(dim(conditions)) &&
        !anyNA(conditions) &&
        all(nzchar(conditions)) &&
        !anyDuplicated(conditions)
    if (!valid) {
        stop(
            "Cutoff conditions must be unique, non-empty labels.",
            call. = FALSE
        )
    }

    sort(conditions, method = "radix")
}

.normalise_manual_cutoffs <- function(cutoffs, conditions) {
    if (is.null(cutoffs)) {
        return(stats::setNames(numeric(), character()))
    }

    valid_shape <- is.numeric(cutoffs) &&
        is.null(dim(cutoffs)) &&
        length(cutoffs) > 0L &&
        !is.null(names(cutoffs))
    if (!valid_shape) {
        stop(
            "`cutoffs` must be NULL or a non-empty named numeric vector.",
            call. = FALSE
        )
    }

    cutoff_names <- names(cutoffs)
    valid_names <- length(cutoff_names) == length(cutoffs) &&
        !anyNA(cutoff_names) &&
        all(nzchar(cutoff_names)) &&
        !anyDuplicated(cutoff_names)
    if (!valid_names) {
        stop(
            "Manual cutoff names must be unique, non-empty condition labels.",
            call. = FALSE
        )
    }
    if (any(!is.finite(cutoffs))) {
        stop(
            "Manual cutoffs must contain only finite values.",
            call. = FALSE
        )
    }

    unknown <- sort(setdiff(cutoff_names, conditions), method = "radix")
    if (length(unknown) > 0L) {
        stop(
            sprintf(
                "Manual cutoff names must match condition labels exactly: %s.",
                paste(sprintf("`%s`", unknown), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    stats::setNames(as.numeric(cutoffs), cutoff_names)
}

.manual_cutoff_diagnostic <- function(condition, cutoff, mean_intensity) {
    observed_range <- range(mean_intensity)
    warnings <- character()
    if (cutoff < observed_range[[1L]] || cutoff > observed_range[[2L]]) {
        warnings <- sprintf(
            paste0(
                "Manual cutoff %s is outside the observed feature-mean ",
                "range [%s, %s] for condition `%s`."
            ),
            format(cutoff, trim = TRUE),
            format(observed_range[[1L]], trim = TRUE),
            format(observed_range[[2L]], trim = TRUE),
            condition
        )
    }

    list(
        source = "manual",
        method = "manual",
        method_version = NA_character_,
        quality = list(
            observed_mean_range = stats::setNames(
                observed_range,
                c("minimum", "maximum")
            )
        ),
        warnings = warnings
    )
}

.not_needed_cutoff_diagnostic <- function() {
    list(
        source = "not_needed",
        method = NA_character_,
        method_version = NA_character_,
        quality = list(),
        warnings = character()
    )
}

.automatic_cutoff_specification <- function() {
    list(
        method = "derivative_boundary",
        method_version = "1",
        minimum_features = 80L,
        minimum_class_features = 12L,
        trend_alpha = 0.001,
        minimum_drop = 0.05,
        density_floor = 0.1,
        minimum_supported_grid_points = 25L,
        polynomial_degree = 3L,
        half_window_bandwidths = 0.5,
        recovery_fraction = 0.5,
        boundary_correction_bandwidths = 1,
        maximum_peak_offset_bandwidths = 0.5,
        minimum_peak_to_noise = 1.5
    )
}

.empty_automatic_cutoff_quality <- function(specification) {
    list(
        evidence = list(
            feature_count = NA_integer_,
            class_counts = c(missing = NA_integer_, complete = NA_integer_),
            minimum_feature_count = specification$minimum_features,
            minimum_class_feature_count =
                specification$minimum_class_features
        ),
        support = list(
            observed_mean_range = c(minimum = NA_real_, maximum = NA_real_),
            bandwidth = NA_real_,
            supported_grid_points = NA_integer_,
            density_floor = specification$density_floor,
            density_supported_grid_points = NA_integer_,
            minimum_density_supported_grid_points =
                specification$minimum_supported_grid_points,
            density_supported_range = c(
                minimum = NA_real_,
                maximum = NA_real_
            )
        ),
        trend = list(
            descending = FALSE,
            credible = FALSE,
            slope = NA_real_,
            likelihood_ratio = NA_real_,
            p_value = NA_real_,
            warnings = character(),
            maximum_p_value = specification$trend_alpha
        ),
        derivative = list(
            polynomial_degree = specification$polynomial_degree,
            half_window_points = NA_integer_,
            half_window_bandwidths =
                specification$half_window_bandwidths,
            recovery_fraction = specification$recovery_fraction,
            lobe_count = NA_integer_,
            credible_lobe_count = NA_integer_,
            minimum_drop = specification$minimum_drop,
            minimum_peak_to_noise = specification$minimum_peak_to_noise,
            boundary_correction_bandwidths =
                specification$boundary_correction_bandwidths,
            maximum_peak_offset_bandwidths =
                specification$maximum_peak_offset_bandwidths
        ),
        selected_transition = NULL
    )
}

.automatic_cutoff_result <- function(
    status,
    cutoff,
    reason,
    quality,
    specification,
    warnings = character()
) {
    list(
        status = status,
        cutoff = cutoff,
        reason = reason,
        method = specification$method,
        method_version = specification$method_version,
        quality = quality,
        warnings = warnings
    )
}

.automatic_cutoff_failure <- function(
    reason,
    quality,
    specification,
    warnings = character()
) {
    .automatic_cutoff_result(
        status = "unidentifiable",
        cutoff = NA_real_,
        reason = reason,
        quality = quality,
        specification = specification,
        warnings = warnings
    )
}

.automatic_cutoff_input_failure <- function(
    reason,
    quality,
    warnings = character()
) {
    list(
        ok = FALSE,
        reason = reason,
        quality = quality,
        warnings = warnings
    )
}

.validate_automatic_cutoff_profile <- function(profile, quality) {
    valid_container <- is.list(profile) &&
        is.data.frame(profile$raw) &&
        is.data.frame(profile$grid) &&
        is.list(profile$metadata)
    if (!valid_container) {
        return(.automatic_cutoff_input_failure(
            "The stored profile is malformed.",
            quality
        ))
    }

    required_raw <- c("mean_intensity", "has_missing")
    required_grid <- c(
        "intensity",
        "weighted_missing_density",
        "weighted_complete_density",
        "missing_proportion",
        "supported"
    )
    if (!all(required_raw %in% names(profile$raw)) ||
        !all(required_grid %in% names(profile$grid))) {
        return(.automatic_cutoff_input_failure(
            "The stored profile lacks required evidence fields.",
            quality
        ))
    }

    list(ok = TRUE, quality = quality, warnings = character())
}

.automatic_cutoff_evidence <- function(profile, quality, specification) {
    raw <- profile$raw
    valid <- is.numeric(raw$mean_intensity) &&
        is.logical(raw$has_missing) &&
        !anyNA(raw$has_missing) &&
        all(is.finite(raw$mean_intensity))
    if (!valid) {
        return(.automatic_cutoff_input_failure(
            "Raw profile evidence contains invalid values.",
            quality
        ))
    }

    class_counts <- c(
        missing = as.integer(sum(raw$has_missing)),
        complete = as.integer(sum(!raw$has_missing))
    )
    quality$evidence$feature_count <- as.integer(nrow(raw))
    quality$evidence$class_counts <- class_counts
    if (nrow(raw) > 0L) {
        quality$support$observed_mean_range <- stats::setNames(
            range(raw$mean_intensity),
            c("minimum", "maximum")
        )
    }
    if (nrow(raw) < specification$minimum_features) {
        reason <- sprintf(
            "Automatic cutoff requires at least %d feature blocks.",
            specification$minimum_features
        )
        return(.automatic_cutoff_input_failure(reason, quality))
    }
    if (any(class_counts < specification$minimum_class_features)) {
        reason <- sprintf(
            paste0(
                "Automatic cutoff requires at least %d missing and %d ",
                "complete feature blocks."
            ),
            specification$minimum_class_features,
            specification$minimum_class_features
        )
        return(.automatic_cutoff_input_failure(reason, quality))
    }

    list(ok = TRUE, raw = raw, quality = quality, warnings = character())
}

.valid_automatic_cutoff_grid <- function(grid) {
    is.numeric(grid$intensity) &&
        is.numeric(grid$weighted_missing_density) &&
        is.numeric(grid$weighted_complete_density) &&
        is.numeric(grid$missing_proportion) &&
        is.logical(grid$supported)
}

.automatic_density_keep <- function(total_density, finite, density_floor) {
    relative_density <- total_density / max(total_density[finite])
    keep <- finite & relative_density >= density_floor
    runs <- rle(keep)
    run_end <- cumsum(runs$lengths)
    run_start <- run_end - runs$lengths + 1L
    eligible_runs <- which(runs$values)
    peak <- which.max(ifelse(finite, total_density, -Inf))
    peak_run <- eligible_runs[
        run_start[eligible_runs] <= peak & run_end[eligible_runs] >= peak
    ]
    if (length(peak_run) == 1L) {
        keep[] <- FALSE
        keep[seq.int(run_start[[peak_run]], run_end[[peak_run]])] <- TRUE
    }

    list(relative_density = relative_density, keep = keep)
}

.finite_automatic_cutoff_grid <- function(
    grid,
    total_density,
    observed_range
) {
    is.finite(grid$intensity) &
        is.finite(grid$missing_proportion) &
        is.finite(total_density) &
        total_density > 0 &
        !is.na(grid$supported) &
        grid$supported &
        grid$intensity >= observed_range[["minimum"]] &
        grid$intensity <= observed_range[["maximum"]]
}

.record_density_support <- function(quality, grid, keep) {
    quality$support$density_supported_grid_points <- as.integer(sum(keep))
    if (sum(keep) > 0L) {
        quality$support$density_supported_range <- stats::setNames(
            range(grid$intensity[keep]),
            c("minimum", "maximum")
        )
    }
    quality
}

.automatic_cutoff_grid_support <- function(
    profile,
    quality,
    specification
) {
    grid <- profile$grid
    if (!.valid_automatic_cutoff_grid(grid)) {
        return(.automatic_cutoff_input_failure(
            "The profile grid or smoothing metadata is invalid.",
            quality
        ))
    }
    total_density <- grid$weighted_missing_density +
        grid$weighted_complete_density
    observed_range <- quality$support$observed_mean_range
    finite <- .finite_automatic_cutoff_grid(
        grid,
        total_density,
        observed_range
    )
    quality$support$supported_grid_points <- as.integer(sum(finite))
    if (!any(finite)) {
        return(.automatic_cutoff_input_failure(
            "The profile has no supported observed-intensity grid.",
            quality
        ))
    }
    density <- .automatic_density_keep(
        total_density,
        finite,
        specification$density_floor
    )
    keep <- density$keep
    quality <- .record_density_support(quality, grid, keep)
    if (sum(keep) < specification$minimum_supported_grid_points) {
        return(.automatic_cutoff_input_failure(
            "The profile has too little density-supported grid coverage.",
            quality
        ))
    }

    list(
        ok = TRUE,
        grid = grid,
        keep = keep,
        relative_density = density$relative_density,
        quality = quality,
        warnings = character()
    )
}

.valid_automatic_cutoff_numeric <- function(
    x,
    y,
    weight,
    spacing,
    bandwidth
) {
    is.numeric(bandwidth) &&
        length(bandwidth) == 1L &&
        all(is.finite(c(x, y, weight, spacing, bandwidth))) &&
        all(weight > 0) &&
        bandwidth > 0 &&
        all(spacing > 0) &&
        max(abs(spacing - stats::median(spacing))) <=
            100 * .Machine$double.eps * max(1, max(abs(x)))
}

.automatic_cutoff_trend_quality <- function(
    raw,
    quality,
    specification
) {
    trend <- .automatic_cutoff_trend(
        raw$mean_intensity,
        raw$has_missing
    )
    trend$credible <- trend$descending &&
        is.finite(trend$p_value) &&
        trend$p_value <= specification$trend_alpha
    trend$maximum_p_value <- specification$trend_alpha
    quality$trend <- trend
    if (length(trend$warnings) > 0L) {
        return(.automatic_cutoff_input_failure(
            "The missingness-trend model emitted numerical warnings.",
            quality,
            trend$warnings
        ))
    }
    if (!trend$credible) {
        return(.automatic_cutoff_input_failure(
            paste0(
                "No statistically credible descending missingness trend ",
                "was detected."
            ),
            quality
        ))
    }

    list(ok = TRUE, quality = quality, warnings = character())
}

.automatic_cutoff_numeric_candidate <- function(
    profile,
    raw,
    support,
    specification
) {
    keep <- support$keep
    x <- support$grid$intensity[keep]
    y <- support$grid$missing_proportion[keep]
    weight <- support$relative_density[keep]
    spacing <- diff(x)
    bandwidth <- profile$metadata$bandwidth
    if (!.valid_automatic_cutoff_numeric(
        x, y, weight, spacing, bandwidth
    )) {
        return(.automatic_cutoff_input_failure(
            "The profile grid or smoothing metadata is invalid.",
            support$quality
        ))
    }

    quality <- support$quality
    quality$support$bandwidth <- unname(as.numeric(bandwidth))
    trend <- .automatic_cutoff_trend_quality(raw, quality, specification)
    if (!trend$ok) {
        return(trend)
    }

    list(
        ok = TRUE,
        x = x,
        y = y,
        weight = weight,
        bandwidth = unname(as.numeric(bandwidth)),
        spacing = stats::median(spacing),
        quality = trend$quality,
        warnings = character()
    )
}

.prepare_automatic_cutoff_input <- function(profile, specification) {
    quality <- .empty_automatic_cutoff_quality(specification)
    validated <- .validate_automatic_cutoff_profile(profile, quality)
    if (!validated$ok) {
        return(validated)
    }
    evidence <- .automatic_cutoff_evidence(
        profile,
        validated$quality,
        specification
    )
    if (!evidence$ok) {
        return(evidence)
    }
    support <- .automatic_cutoff_grid_support(
        profile,
        evidence$quality,
        specification
    )
    if (!support$ok) {
        return(support)
    }
    .automatic_cutoff_numeric_candidate(
        profile,
        evidence$raw,
        support,
        specification
    )
}

.fit_automatic_cutoff_trend <- function(scaled, has_missing) {
    captured <- new.env(parent = emptyenv())
    captured$warnings <- character()
    fitted <- tryCatch(
        withCallingHandlers(
            stats::glm.fit(
                x = cbind(intercept = 1, intensity = scaled),
                y = as.integer(has_missing),
                family = stats::binomial()
            ),
            warning = function(condition) {
                captured$warnings <- c(
                    captured$warnings,
                    conditionMessage(condition)
                )
                invokeRestart("muffleWarning")
            }
        ),
        error = function(...) NULL
    )

    list(fitted = fitted, warnings = unique(captured$warnings))
}

.automatic_cutoff_trend <- function(intensity, has_missing) {
    stable_order <- order(intensity, has_missing, method = "radix")
    intensity <- intensity[stable_order]
    has_missing <- has_missing[stable_order]
    spread <- stats::sd(intensity)
    if (!is.finite(spread) || spread <= 0) {
        return(list(
            descending = FALSE,
            credible = FALSE,
            slope = NA_real_,
            likelihood_ratio = NA_real_,
            p_value = NA_real_,
            warnings = character()
        ))
    }

    scaled <- (intensity - mean(intensity)) / spread
    fit <- .fit_automatic_cutoff_trend(scaled, has_missing)
    if (is.null(fit$fitted)) {
        return(list(
            descending = FALSE,
            credible = FALSE,
            slope = NA_real_,
            likelihood_ratio = NA_real_,
            p_value = NA_real_,
            warnings = fit$warnings
        ))
    }

    fitted <- fit$fitted
    likelihood_ratio <- fitted$null.deviance - fitted$deviance
    p_value <- if (is.finite(likelihood_ratio) && likelihood_ratio >= 0) {
        stats::pchisq(likelihood_ratio, df = 1, lower.tail = FALSE)
    } else {
        NA_real_
    }
    slope <- unname(fitted$coefficients[["intensity"]])
    list(
        descending = is.finite(slope) && slope < 0,
        credible = FALSE,
        slope = slope,
        likelihood_ratio = likelihood_ratio,
        p_value = p_value,
        warnings = fit$warnings
    )
}

.automatic_local_polynomial <- function(y, spacing, half_window, degree) {
    offsets <- seq.int(-half_window, half_window)
    design <- outer(offsets, 0:degree, `^`)
    projector <- solve(crossprod(design), t(design))
    smooth_coefficient <- projector[1L, ]
    derivative_coefficient <- projector[2L, ] / spacing
    indices <- seq.int(half_window + 1L, length(y) - half_window)
    smooth <- derivative <- rep(NA_real_, length(y))
    for (index in indices) {
        window <- seq.int(index - half_window, index + half_window)
        smooth[[index]] <- sum(smooth_coefficient * y[window])
        derivative[[index]] <- sum(derivative_coefficient * y[window])
    }

    list(smooth = smooth, derivative = derivative, indices = indices)
}

.automatic_derivative_crossing <- function(
    x,
    derivative,
    left,
    right,
    level
) {
    x_pair <- x[c(left, right)]
    y_pair <- derivative[c(left, right)]
    if (any(!is.finite(c(x_pair, y_pair))) ||
        y_pair[[1L]] == y_pair[[2L]]) {
        return(mean(x_pair))
    }

    x_pair[[1L]] +
        (level - y_pair[[1L]]) * diff(x_pair) / diff(y_pair)
}

.automatic_derivative_minima <- function(filtered) {
    valid <- filtered$indices
    valid[
        valid > min(valid) &
            valid < max(valid) &
            filtered$derivative[valid] <=
                filtered$derivative[valid - 1L] &
            filtered$derivative[valid] <
                filtered$derivative[valid + 1L] &
            filtered$derivative[valid] < 0
    ]
}

.automatic_derivative_lobe <- function(peak, filtered, specification) {
    valid <- filtered$indices
    peak_value <- filtered$derivative[[peak]]
    recovery_level <- specification$recovery_fraction * peak_value
    right_candidates <- valid[
        valid > peak & filtered$derivative[valid] >= recovery_level
    ]
    left_candidates <- valid[
        valid < peak & filtered$derivative[valid] >= recovery_level
    ]
    if (length(right_candidates) == 0L ||
        length(left_candidates) == 0L) {
        return(NULL)
    }

    right <- right_candidates[[1L]]
    left <- utils::tail(left_candidates, 1L)
    outside <- valid[valid < left | valid > right]
    noise <- stats::mad(
        filtered$derivative[outside],
        center = stats::median(filtered$derivative[outside]),
        constant = 1.4826,
        na.rm = TRUE
    )
    list(
        peak = peak,
        peak_value = peak_value,
        recovery_level = recovery_level,
        left = left,
        right = right,
        drop = filtered$smooth[[left]] - filtered$smooth[[right]],
        peak_to_noise = abs(peak_value) /
            max(noise, .Machine$double.eps)
    )
}

.credible_automatic_derivative_lobe <- function(lobe, specification) {
    is.finite(lobe$drop) &&
        lobe$drop >= specification$minimum_drop &&
        is.finite(lobe$peak_to_noise) &&
        lobe$peak_to_noise >= specification$minimum_peak_to_noise
}

.automatic_derivative_half_window <- function(candidate, specification) {
    max(
        specification$polynomial_degree,
        as.integer(round(
            specification$half_window_bandwidths *
                candidate$bandwidth / candidate$spacing
        ))
    )
}

.automatic_cutoff_derivative <- function(candidate, specification) {
    quality <- candidate$quality
    half_window <- .automatic_derivative_half_window(candidate, specification)
    quality$derivative$half_window_points <- half_window
    if (2L * half_window + 1L > length(candidate$y)) {
        return(.automatic_cutoff_input_failure(
            "The supported profile is too short for derivative smoothing.",
            quality
        ))
    }
    filtered <- .automatic_local_polynomial(
        candidate$y,
        candidate$spacing,
        half_window,
        specification$polynomial_degree
    )
    lobes <- lapply(
        .automatic_derivative_minima(filtered),
        .automatic_derivative_lobe,
        filtered = filtered,
        specification = specification
    )
    lobes <- Filter(Negate(is.null), lobes)
    quality$derivative$lobe_count <- as.integer(length(lobes))
    credible <- vapply(
        lobes,
        .credible_automatic_derivative_lobe,
        logical(1L),
        specification = specification
    )
    lobes <- lobes[credible]
    quality$derivative$credible_lobe_count <- as.integer(length(lobes))
    if (length(lobes) == 0L) {
        return(.automatic_cutoff_input_failure(
            paste0(
                "No credible descending derivative lobe has an interior ",
                "boundary."
            ),
            quality
        ))
    }

    list(
        ok = TRUE,
        quality = quality,
        filtered = filtered,
        lobes = lobes,
        warnings = character()
    )
}

.select_automatic_cutoff_transition <- function(
    candidate,
    filtered,
    lobes,
    specification
) {
    peak_intensity <- vapply(
        lobes,
        function(lobe) candidate$x[[lobe$peak]],
        numeric(1L)
    )
    lobe <- lobes[[which.min(peak_intensity)]]
    left_crossing <- .automatic_derivative_crossing(
        candidate$x,
        filtered$derivative,
        lobe$left,
        lobe$left + 1L,
        lobe$recovery_level
    )
    right_crossing <- .automatic_derivative_crossing(
        candidate$x,
        filtered$derivative,
        lobe$right - 1L,
        lobe$right,
        lobe$recovery_level
    )
    corrected_boundary <- right_crossing -
        specification$boundary_correction_bandwidths * candidate$bandwidth
    peak_offset_boundary <- candidate$x[[lobe$peak]] +
        specification$maximum_peak_offset_bandwidths * candidate$bandwidth
    list(
        cutoff = min(corrected_boundary, peak_offset_boundary),
        quality = list(
            peak_intensity = candidate$x[[lobe$peak]],
            left_half_depth_intensity = left_crossing,
            right_half_depth_intensity = right_crossing,
            peak_derivative = lobe$peak_value,
            half_depth_derivative = lobe$recovery_level,
            drop = lobe$drop,
            peak_to_noise = lobe$peak_to_noise,
            kde_corrected_boundary = corrected_boundary,
            peak_offset_boundary = peak_offset_boundary
        )
    )
}

.detect_automatic_cutoff <- function(profile) {
    specification <- .automatic_cutoff_specification()
    candidate <- .prepare_automatic_cutoff_input(profile, specification)
    if (!candidate$ok) {
        return(.automatic_cutoff_failure(
            candidate$reason,
            candidate$quality,
            specification,
            candidate$warnings
        ))
    }
    derivative <- .automatic_cutoff_derivative(candidate, specification)
    if (!derivative$ok) {
        return(.automatic_cutoff_failure(
            derivative$reason,
            derivative$quality,
            specification,
            derivative$warnings
        ))
    }

    quality <- derivative$quality
    filtered <- derivative$filtered
    lobes <- derivative$lobes
    transition <- .select_automatic_cutoff_transition(
        candidate,
        filtered,
        lobes,
        specification
    )
    cutoff <- transition$cutoff
    quality$selected_transition <- transition$quality
    if (!is.finite(cutoff) ||
        cutoff <= min(candidate$x) ||
        cutoff >= max(candidate$x)) {
        return(.automatic_cutoff_failure(
            "The derivative boundary is not an interior finite cutoff.",
            quality,
            specification
        ))
    }

    .automatic_cutoff_result(
        status = "ok",
        cutoff = unname(as.numeric(cutoff)),
        reason = "",
        quality = quality,
        specification = specification
    )
}

.automatic_cutoff_diagnostic <- function(decision) {
    list(
        source = "automatic",
        method = decision$method,
        method_version = decision$method_version,
        quality = decision$quality,
        warnings = decision$warnings
    )
}

.abort_automatic_cutoff <- function(
    condition,
    decision,
    diagnostic,
    profile
) {
    diagnostic$profile <- if (is.list(profile)) profile$metadata else NULL
    message <- sprintf(
        paste0(
            "Automatic cutoff is unidentifiable for condition `%s`: %s ",
            "Inspect `profile` on this error and supply a manual cutoff ",
            "for `%s`."
        ),
        condition,
        decision$reason,
        condition
    )
    error <- structure(
        list(
            message = message,
            call = NULL,
            condition = condition,
            reason = decision$reason,
            diagnostic = diagnostic,
            profile = profile
        ),
        class = c(
            "imputefinder_cutoff_unidentifiable",
            "imputefinder_cutoff_error",
            "error",
            "condition"
        )
    )

    stop(error)
}

.validate_cutoff_statistics <- function(statistics) {
    required <- c("condition", "missing_count", "mean_intensity")
    if (!is.data.frame(statistics) || !all(required %in% names(statistics))) {
        stop(
            "`statistics` must contain condition, missing_count, and ",
            "mean_intensity.",
            call. = FALSE
        )
    }
    if (anyNA(statistics$condition) || any(!nzchar(statistics$condition)) ||
        any(!is.finite(statistics$mean_intensity))) {
        stop("`statistics` contains invalid cutoff inputs.", call. = FALSE)
    }

    invisible(statistics)
}
