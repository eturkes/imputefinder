.feature_condition_statistics <- function(
    data,
    groups_by_sample,
    seed_log = .empty_seed_log()
) {
    groups_by_sample <- .validate_result_groups(data, groups_by_sample)
    features <- rownames(data)
    conditions <- sort(unique(unname(groups_by_sample)), method = "radix")
    feature_count <- nrow(data)
    condition_count <- length(conditions)
    row_count <- feature_count * condition_count

    statistics <- data.frame(
        feature = rep(features, each = condition_count),
        condition = rep(conditions, times = feature_count),
        sample_count = integer(row_count),
        observed_count = integer(row_count),
        missing_count = integer(row_count),
        missing_fraction = numeric(row_count),
        mean_intensity = numeric(row_count),
        seeded = logical(row_count),
        stringsAsFactors = FALSE
    )
    if (row_count == 0L) {
        return(statistics)
    }

    seeded <- .seeded_feature_condition_matrix(
        features,
        conditions,
        seed_log
    )
    statistics$seeded <- as.vector(t(seeded))

    for (condition_index in seq_along(conditions)) {
        condition <- conditions[[condition_index]]
        condition_data <- data[
            ,
            groups_by_sample == condition,
            drop = FALSE
        ]
        observed_count <- as.integer(rowSums(!is.na(condition_data)))
        missing_count <- as.integer(ncol(condition_data) - observed_count)
        mean_intensity <- rowMeans(condition_data, na.rm = TRUE)

        if (any(!is.finite(mean_intensity))) {
            stop(
                paste0(
                    "Every surviving feature-condition block must have a ",
                    "finite mean; run condition rescue before statistics."
                ),
                call. = FALSE
            )
        }

        rows <- seq.int(condition_index, row_count, by = condition_count)
        statistics$sample_count[rows] <- as.integer(ncol(condition_data))
        statistics$observed_count[rows] <- observed_count
        statistics$missing_count[rows] <- missing_count
        statistics$missing_fraction[rows] <-
            missing_count / ncol(condition_data)
        statistics$mean_intensity[rows] <- mean_intensity
    }

    statistics
}

.seeded_feature_condition_matrix <- function(features, conditions, seed_log) {
    seeded <- matrix(
        FALSE,
        nrow = length(features),
        ncol = length(conditions),
        dimnames = list(features, conditions)
    )
    if (!is.data.frame(seed_log) ||
        !all(c("feature", "condition") %in% names(seed_log))) {
        stop("`seed_log` must contain feature and condition columns.", call. = FALSE)
    }
    if (nrow(seed_log) == 0L) {
        return(seeded)
    }

    feature_index <- match(seed_log$feature, features)
    condition_index <- match(seed_log$condition, conditions)
    if (anyNA(feature_index) || anyNA(condition_index)) {
        stop(
            "`seed_log` must reference surviving features and conditions.",
            call. = FALSE
        )
    }

    seeded[cbind(feature_index, condition_index)] <- TRUE
    seeded
}

.build_missingness_profiles <- function(statistics) {
    .validate_profile_statistics(statistics)
    conditions <- sort(unique(statistics$condition), method = "radix")
    profiles <- lapply(
        conditions,
        function(condition) {
            condition_statistics <- statistics[
                statistics$condition == condition,
                ,
                drop = FALSE
            ]
            .condition_missingness_profile(condition_statistics)
        }
    )

    stats::setNames(profiles, conditions)
}

.condition_missingness_profile <- function(statistics) {
    raw <- statistics
    raw$has_missing <- raw$missing_count > 0L
    means <- raw$mean_intensity
    class_counts <- c(
        missing = as.integer(sum(raw$has_missing)),
        complete = as.integer(sum(!raw$has_missing))
    )
    bandwidth <- .profile_bandwidth(means)
    grid_points <- 512L
    observed_range <- stats::setNames(
        range(means),
        c("minimum", "maximum")
    )
    grid_range <- .profile_grid_range(observed_range, bandwidth$value)

    missing_density <- .profile_density(
        means[raw$has_missing],
        bandwidth$value,
        grid_points,
        grid_range
    )
    complete_density <- .profile_density(
        means[!raw$has_missing],
        bandwidth$value,
        grid_points,
        grid_range
    )
    intensity <- if (class_counts[["missing"]] > 0L) {
        missing_density$x
    } else {
        complete_density$x
    }
    if (class_counts[["missing"]] > 0L &&
        class_counts[["complete"]] > 0L &&
        !identical(missing_density$x, complete_density$x)) {
        stop("Profile classes must share one intensity grid.", call. = FALSE)
    }

    weighted_missing <- class_counts[["missing"]] * missing_density$y
    weighted_complete <- class_counts[["complete"]] * complete_density$y
    total_density <- weighted_missing + weighted_complete
    supported <- is.finite(total_density) & total_density > 0
    missing_proportion <- rep(NA_real_, grid_points)
    missing_proportion[supported] <-
        weighted_missing[supported] / total_density[supported]

    grid <- data.frame(
        intensity = as.numeric(intensity),
        missing_density = as.numeric(missing_density$y),
        complete_density = as.numeric(complete_density$y),
        weighted_missing_density = as.numeric(weighted_missing),
        weighted_complete_density = as.numeric(weighted_complete),
        missing_proportion = as.numeric(missing_proportion),
        supported = as.logical(supported),
        stringsAsFactors = FALSE
    )
    metadata <- list(
        feature_count = as.integer(nrow(raw)),
        class_counts = class_counts,
        profile_type = if (class_counts[["missing"]] == 0L) {
            "complete_only"
        } else if (class_counts[["complete"]] == 0L) {
            "missing_only"
        } else {
            "mixed"
        },
        seeded_feature_count = as.integer(sum(raw$seeded)),
        grid_points = grid_points,
        density_method = "stats::density",
        kernel = "gaussian",
        bandwidth = bandwidth$value,
        bandwidth_method = bandwidth$method,
        support_extension_bandwidths = 3,
        observed_mean_range = observed_range,
        grid_range = grid_range,
        unsupported_grid_points = as.integer(sum(!supported)),
        warnings = if (all(supported)) {
            character()
        } else {
            sprintf(
                "%d profile grid points have no numerical density support.",
                sum(!supported)
            )
        }
    )

    list(raw = raw, grid = grid, metadata = metadata)
}

.profile_bandwidth <- function(values) {
    values <- sort(values, method = "radix")
    bandwidth_values <- if (length(values) == 1L) {
        rep(values, 2L)
    } else {
        values
    }
    bandwidth <- suppressWarnings(stats::bw.nrd0(bandwidth_values))
    method <- if (length(values) == 1L) {
        "pooled_nrd0_singleton_replication"
    } else {
        "pooled_nrd0"
    }

    if (!is.finite(bandwidth) || bandwidth <= 0) {
        bandwidth <- max(abs(values), 1) * sqrt(.Machine$double.eps)
        method <- "pooled_nrd0_numeric_fallback"
    }
    if (!is.finite(bandwidth) || bandwidth <= 0) {
        stop("Profile bandwidth must be finite and positive.", call. = FALSE)
    }

    list(value = unname(as.numeric(bandwidth)), method = method)
}

.profile_grid_range <- function(observed_range, bandwidth) {
    padding <- 3 * bandwidth
    grid_range <- c(
        minimum = observed_range[["minimum"]] - padding,
        maximum = observed_range[["maximum"]] + padding
    )
    if (any(!is.finite(grid_range)) ||
        grid_range[["minimum"]] >= grid_range[["maximum"]]) {
        stop(
            "Feature means are outside the numerically stable profile range.",
            call. = FALSE
        )
    }

    grid_range
}

.profile_density <- function(values, bandwidth, grid_points, grid_range) {
    if (length(values) == 0L) {
        return(list(
            x = seq(
                grid_range[["minimum"]],
                grid_range[["maximum"]],
                length.out = grid_points
            ),
            y = numeric(grid_points)
        ))
    }

    values <- sort(values, method = "radix")
    estimated <- stats::density(
        values,
        bw = bandwidth,
        kernel = "gaussian",
        n = grid_points,
        from = grid_range[["minimum"]],
        to = grid_range[["maximum"]]
    )
    density <- pmax(as.numeric(estimated$y), 0)
    if (length(density) != grid_points || any(!is.finite(density))) {
        stop("Profile density estimation produced invalid values.", call. = FALSE)
    }

    list(x = as.numeric(estimated$x), y = density)
}

.validate_profile_statistics <- function(statistics) {
    required <- c(
        "feature",
        "condition",
        "sample_count",
        "observed_count",
        "missing_count",
        "missing_fraction",
        "mean_intensity",
        "seeded"
    )
    valid_schema <- is.data.frame(statistics) &&
        all(required %in% names(statistics)) &&
        nrow(statistics) > 0L &&
        is.character(statistics$feature) &&
        is.character(statistics$condition) &&
        is.integer(statistics$sample_count) &&
        is.integer(statistics$observed_count) &&
        is.integer(statistics$missing_count) &&
        is.numeric(statistics$missing_fraction) &&
        is.numeric(statistics$mean_intensity) &&
        is.logical(statistics$seeded)
    if (!valid_schema || anyNA(statistics[, required, drop = FALSE])) {
        stop("`statistics` does not satisfy the profile-input schema.", call. = FALSE)
    }

    valid_names <- all(nzchar(statistics$feature)) &&
        all(nzchar(statistics$condition)) &&
        !anyDuplicated(
            statistics[, c("feature", "condition"), drop = FALSE]
        )
    valid_counts <- statistics$sample_count > 0L &
        statistics$observed_count > 0L &
        statistics$missing_count >= 0L &
        statistics$observed_count + statistics$missing_count ==
            statistics$sample_count
    valid_fraction <- is.finite(statistics$missing_fraction) &
        statistics$missing_fraction >= 0 &
        statistics$missing_fraction <= 1 &
        abs(
            statistics$missing_fraction -
                statistics$missing_count / statistics$sample_count
        ) <= sqrt(.Machine$double.eps)
    valid_means <- is.finite(statistics$mean_intensity)
    if (!all(valid_names) || !all(valid_counts) || !all(valid_fraction) ||
        !all(valid_means)) {
        stop("`statistics` contains invalid profile inputs.", call. = FALSE)
    }

    invisible(statistics)
}
