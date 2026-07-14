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
