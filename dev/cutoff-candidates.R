# Fixed automatic-cutoff candidates for the frozen M5 validation protocol.
# Development-only until M5c promotes one detector into R/detect_cutoff.R.

.cutoff_candidate_constants <- list(
    minimum_features = 80L,
    minimum_class_features = 12L,
    trend_alpha = 0.001,
    minimum_drop = 0.05,
    density_floor = 0.1,
    segmented = list(
        minimum_fit_gain = 0.35,
        width_bandwidths = seq(0.5, 8, by = 0.25),
        center_stride = 4L,
        boundary_correction_bandwidths = 1
    ),
    derivative = list(
        polynomial_degree = 3L,
        half_window_bandwidths = 0.5,
        recovery_fraction = 0.5,
        boundary_correction_bandwidths = 1,
        maximum_peak_offset_bandwidths = 0.5,
        minimum_peak_to_noise = 1.5
    )
)

.cutoff_candidate_failure <- function(reason) {
    list(status = "unidentifiable", cutoff = NA_real_, reason = reason)
}

.cutoff_candidate_input <- function(profile) {
    failure <- function(reason) list(ok = FALSE, reason = reason)
    if (!is.list(profile) || !is.data.frame(profile$raw) ||
        !is.data.frame(profile$grid) || !is.list(profile$metadata)) {
        return(failure("The stored profile is malformed."))
    }
    required_raw <- c("mean_intensity", "has_missing")
    required_grid <- c(
        "intensity", "weighted_missing_density",
        "weighted_complete_density", "missing_proportion", "supported"
    )
    if (!all(required_raw %in% names(profile$raw)) ||
        !all(required_grid %in% names(profile$grid))) {
        return(failure("The stored profile lacks required evidence fields."))
    }

    constants <- .cutoff_candidate_constants
    raw <- profile$raw
    classes <- table(factor(raw$has_missing, levels = c(FALSE, TRUE)))
    if (nrow(raw) < constants$minimum_features) {
        return(failure(sprintf(
            "Automatic cutoff requires at least %d feature blocks.",
            constants$minimum_features
        )))
    }
    if (any(classes < constants$minimum_class_features)) {
        return(failure(sprintf(
            paste0(
                "Automatic cutoff requires at least %d missing and %d ",
                "complete feature blocks."
            ),
            constants$minimum_class_features,
            constants$minimum_class_features
        )))
    }
    if (any(!is.finite(raw$mean_intensity)) || anyNA(raw$has_missing)) {
        return(failure("Raw profile evidence contains invalid values."))
    }

    grid <- profile$grid
    total_density <- grid$weighted_missing_density +
        grid$weighted_complete_density
    observed_range <- range(raw$mean_intensity)
    finite <- is.finite(grid$intensity) &
        is.finite(grid$missing_proportion) &
        is.finite(total_density) & total_density > 0 &
        !is.na(grid$supported) & grid$supported &
        grid$intensity >= observed_range[[1L]] &
        grid$intensity <= observed_range[[2L]]
    if (!any(finite)) {
        return(failure("The profile has no supported observed-intensity grid."))
    }
    relative_density <- total_density / max(total_density[finite])
    keep <- finite & relative_density >= constants$density_floor
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
    if (sum(keep) < 25L) {
        return(failure("The profile has too little density-supported grid coverage."))
    }

    x <- grid$intensity[keep]
    y <- grid$missing_proportion[keep]
    weight <- relative_density[keep]
    spacing <- diff(x)
    bandwidth <- profile$metadata$bandwidth
    valid_numeric <- all(is.finite(c(x, y, weight, spacing, bandwidth))) &&
        all(weight > 0) && bandwidth > 0 && all(spacing > 0) &&
        max(abs(spacing - stats::median(spacing))) <=
            100 * .Machine$double.eps * max(1, max(abs(x)))
    if (!valid_numeric) {
        return(failure("The profile grid or smoothing metadata is invalid."))
    }

    trend <- .cutoff_candidate_trend(raw$mean_intensity, raw$has_missing)
    if (!trend$descending || !is.finite(trend$p_value) ||
        trend$p_value > constants$trend_alpha) {
        return(failure(
            "No statistically credible descending missingness trend was detected."
        ))
    }

    list(
        ok = TRUE,
        x = x,
        y = y,
        weight = weight,
        bandwidth = as.numeric(bandwidth),
        spacing = stats::median(spacing),
        trend = trend
    )
}

.cutoff_candidate_trend <- function(intensity, has_missing) {
    spread <- stats::sd(intensity)
    if (!is.finite(spread) || spread <= 0) {
        return(list(
            descending = FALSE,
            slope = NA_real_,
            likelihood_ratio = NA_real_,
            p_value = NA_real_
        ))
    }
    scaled <- (intensity - mean(intensity)) / spread
    fitted <- tryCatch(
        suppressWarnings(stats::glm.fit(
            x = cbind(intercept = 1, intensity = scaled),
            y = as.integer(has_missing),
            family = stats::binomial()
        )),
        error = function(...) NULL
    )
    if (is.null(fitted)) {
        return(list(
            descending = FALSE,
            slope = NA_real_,
            likelihood_ratio = NA_real_,
            p_value = NA_real_
        ))
    }
    statistic <- fitted$null.deviance - fitted$deviance
    p_value <- if (is.finite(statistic) && statistic >= 0) {
        stats::pchisq(statistic, df = 1, lower.tail = FALSE)
    } else {
        NA_real_
    }
    slope <- unname(fitted$coefficients[["intensity"]])
    list(
        descending = is.finite(slope) && slope < 0,
        slope = slope,
        likelihood_ratio = statistic,
        p_value = p_value
    )
}

.weighted_ramp_fit <- function(x, y, weight, left, right) {
    if (!is.finite(left) || !is.finite(right) || left >= right) {
        return(NULL)
    }
    transition <- pmin(pmax((x - left) / (right - left), 0), 1)
    total_weight <- sum(weight)
    mean_transition <- sum(weight * transition) / total_weight
    mean_y <- sum(weight * y) / total_weight
    centered_transition <- transition - mean_transition
    denominator <- sum(weight * centered_transition^2)
    if (!is.finite(denominator) || denominator <= 0) {
        return(NULL)
    }
    slope <- sum(weight * centered_transition * (y - mean_y)) / denominator
    intercept <- mean_y - slope * mean_transition
    residual <- y - intercept - slope * transition
    list(
        sse = sum(weight * residual^2),
        high = intercept,
        low = intercept + slope,
        drop = -slope,
        left = left,
        right = right
    )
}

.segmented_plateau_fit <- function(candidate) {
    constants <- .cutoff_candidate_constants$segmented
    x <- candidate$x
    y <- candidate$y
    weight <- candidate$weight
    bandwidth <- candidate$bandwidth
    centers <- x[seq.int(1L, length(x), by = constants$center_stride)]
    widths <- bandwidth * constants$width_bandwidths
    best <- NULL
    for (width in widths) {
        valid_centers <- centers[
            centers - width / 2 > min(x) & centers + width / 2 < max(x)
        ]
        for (center in valid_centers) {
            fitted <- .weighted_ramp_fit(
                x, y, weight,
                center - width / 2,
                center + width / 2
            )
            if (!is.null(fitted) && fitted$drop > 0 &&
                (is.null(best) || fitted$sse < best$sse)) {
                best <- fitted
            }
        }
    }
    best
}

# Two-plateau continuous segmented regression; returns the fitted ramp's
# density-deblurred right breakpoint.
segmented_plateau_candidate <- function(profile) {
    candidate <- .cutoff_candidate_input(profile)
    if (!candidate$ok) {
        return(.cutoff_candidate_failure(candidate$reason))
    }
    fitted <- .segmented_plateau_fit(candidate)
    if (is.null(fitted)) {
        return(.cutoff_candidate_failure(
            "No descending two-breakpoint segmented fit was available."
        ))
    }
    flat_mean <- stats::weighted.mean(candidate$y, candidate$weight)
    flat_sse <- sum(candidate$weight * (candidate$y - flat_mean)^2)
    fit_gain <- 1 - fitted$sse / flat_sse
    constants <- .cutoff_candidate_constants
    if (!is.finite(fit_gain) ||
        fit_gain < constants$segmented$minimum_fit_gain ||
        fitted$drop < constants$minimum_drop) {
        return(.cutoff_candidate_failure(
            "The segmented transition is too weak relative to a flat profile."
        ))
    }
    cutoff <- fitted$right -
        constants$segmented$boundary_correction_bandwidths *
            candidate$bandwidth
    if (!is.finite(cutoff) || cutoff <= min(candidate$x) ||
        cutoff >= max(candidate$x)) {
        return(.cutoff_candidate_failure(
            "The segmented boundary is not an interior finite cutoff."
        ))
    }
    list(status = "ok", cutoff = cutoff)
}

.savitzky_golay_profile <- function(y, spacing, half_window, degree) {
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

.interpolate_derivative_crossing <- function(x, derivative, left, right, level) {
    x_pair <- x[c(left, right)]
    y_pair <- derivative[c(left, right)]
    if (any(!is.finite(c(x_pair, y_pair))) || y_pair[[1L]] == y_pair[[2L]]) {
        return(mean(x_pair))
    }
    x_pair[[1L]] + (level - y_pair[[1L]]) * diff(x_pair) / diff(y_pair)
}

# Savitzky-Golay local-polynomial derivative; takes the right half-depth edge of
# the dominant descending derivative lobe and corrects once for KDE broadening.
derivative_boundary_candidate <- function(profile) {
    candidate <- .cutoff_candidate_input(profile)
    if (!candidate$ok) {
        return(.cutoff_candidate_failure(candidate$reason))
    }
    constants <- .cutoff_candidate_constants
    derivative_constants <- constants$derivative
    half_window <- max(
        derivative_constants$polynomial_degree,
        as.integer(round(
            derivative_constants$half_window_bandwidths *
                candidate$bandwidth / candidate$spacing
        ))
    )
    if (2L * half_window + 1L > length(candidate$y)) {
        return(.cutoff_candidate_failure(
            "The supported profile is too short for derivative smoothing."
        ))
    }
    filtered <- .savitzky_golay_profile(
        candidate$y,
        candidate$spacing,
        half_window,
        derivative_constants$polynomial_degree
    )
    valid <- filtered$indices
    local_minimum <- valid[
        valid > min(valid) & valid < max(valid) &
            filtered$derivative[valid] <= filtered$derivative[valid - 1L] &
            filtered$derivative[valid] < filtered$derivative[valid + 1L] &
            filtered$derivative[valid] < 0
    ]
    lobes <- lapply(
        local_minimum,
        function(peak) {
            peak_value <- filtered$derivative[[peak]]
            recovery_level <- derivative_constants$recovery_fraction *
                peak_value
            right_candidates <- valid[
                valid > peak &
                    filtered$derivative[valid] >= recovery_level
            ]
            left_candidates <- valid[
                valid < peak &
                    filtered$derivative[valid] >= recovery_level
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
                transition_drop = filtered$smooth[[left]] -
                    filtered$smooth[[right]],
                peak_to_noise = abs(peak_value) /
                    max(noise, .Machine$double.eps)
            )
        }
    )
    lobes <- Filter(Negate(is.null), lobes)
    eligible <- vapply(
        lobes,
        function(lobe) {
            is.finite(lobe$transition_drop) &&
                lobe$transition_drop >= constants$minimum_drop &&
                is.finite(lobe$peak_to_noise) &&
                lobe$peak_to_noise >=
                    derivative_constants$minimum_peak_to_noise
        },
        logical(1L)
    )
    lobes <- lobes[eligible]
    if (length(lobes) == 0L) {
        return(.cutoff_candidate_failure(
            "No credible descending derivative lobe has an interior boundary."
        ))
    }
    peak_intensity <- vapply(
        lobes,
        function(lobe) candidate$x[[lobe$peak]],
        numeric(1L)
    )
    lobe <- lobes[[which.min(peak_intensity)]]
    right_crossing <- .interpolate_derivative_crossing(
        candidate$x,
        filtered$derivative,
        lobe$right - 1L,
        lobe$right,
        lobe$recovery_level
    )
    cutoff <- min(
        right_crossing -
            derivative_constants$boundary_correction_bandwidths *
                candidate$bandwidth,
        candidate$x[[lobe$peak]] +
            derivative_constants$maximum_peak_offset_bandwidths *
                candidate$bandwidth
    )
    if (!is.finite(cutoff) || cutoff <= min(candidate$x) ||
        cutoff >= max(candidate$x)) {
        return(.cutoff_candidate_failure(
            "The derivative boundary is not an interior finite cutoff."
        ))
    }
    list(status = "ok", cutoff = cutoff)
}

cutoff_candidate_detectors <- function() {
    list(
        segmented_plateau_v1 = segmented_plateau_candidate,
        derivative_boundary_v1 = derivative_boundary_candidate
    )
}
