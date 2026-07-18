.PRE_RESCUE_EVIDENCE_SCHEMA <- "pre_rescue_evidence_v1"
.PRE_RESCUE_EVIDENCE_FIELDS <- c(
    "schema",
    "feature_condition",
    "sample"
)
.PRE_RESCUE_FEATURE_CONDITION_FIELDS <- c(
    "feature",
    "condition",
    "sample_count",
    "observed_count",
    "missing_count",
    "missing_fraction",
    "mean_intensity"
)
.PRE_RESCUE_SAMPLE_FIELDS <- c(
    "sample",
    "condition",
    "feature_count",
    "observed_count",
    "missing_count",
    "missing_fraction",
    "mean_intensity"
)

.pre_rescue_observed_summary <- function(data, margin) {
    observed <- !is.na(data)
    count <- if (margin == 1L) {
        as.integer(rowSums(observed))
    } else {
        as.integer(colSums(observed))
    }
    mean_intensity <- if (margin == 1L) {
        rowMeans(data, na.rm = TRUE)
    } else {
        colMeans(data, na.rm = TRUE)
    }
    mean_intensity[count == 0L] <- NA_real_

    list(count = count, mean_intensity = unname(mean_intensity))
}

.pre_rescue_feature_condition <- function(data, groups_by_sample) {
    features <- rownames(data)
    conditions <- sort(unique(unname(groups_by_sample)), method = "radix")
    condition_count <- length(conditions)
    row_count <- length(features) * condition_count
    output <- data.frame(
        feature = rep(features, each = condition_count),
        condition = rep(conditions, times = length(features)),
        sample_count = integer(row_count),
        observed_count = integer(row_count),
        missing_count = integer(row_count),
        missing_fraction = numeric(row_count),
        mean_intensity = numeric(row_count),
        stringsAsFactors = FALSE
    )

    for (condition_index in seq_along(conditions)) {
        selected <- groups_by_sample == conditions[[condition_index]]
        block <- data[, selected, drop = FALSE]
        summary <- .pre_rescue_observed_summary(block, 1L)
        rows <- seq.int(condition_index, row_count, by = condition_count)
        sample_count <- as.integer(ncol(block))
        missing_count <- as.integer(sample_count - summary$count)

        output$sample_count[rows] <- sample_count
        output$observed_count[rows] <- summary$count
        output$missing_count[rows] <- missing_count
        output$missing_fraction[rows] <- missing_count / sample_count
        output$mean_intensity[rows] <- summary$mean_intensity
    }

    output
}

.pre_rescue_sample <- function(data, groups_by_sample) {
    summary <- .pre_rescue_observed_summary(data, 2L)
    feature_count <- as.integer(nrow(data))
    missing_count <- as.integer(feature_count - summary$count)

    data.frame(
        sample = colnames(data),
        condition = unname(groups_by_sample),
        feature_count = rep(feature_count, ncol(data)),
        observed_count = summary$count,
        missing_count = missing_count,
        missing_fraction = missing_count / feature_count,
        mean_intensity = summary$mean_intensity,
        stringsAsFactors = FALSE
    )
}

.new_pre_rescue_evidence <- function(data, groups_by_sample) {
    .validate_matrix_data(data)
    groups_by_sample <- .validate_result_groups(data, groups_by_sample)

    list(
        schema = .PRE_RESCUE_EVIDENCE_SCHEMA,
        feature_condition = .pre_rescue_feature_condition(
            data,
            groups_by_sample
        ),
        sample = .pre_rescue_sample(data, groups_by_sample)
    )
}

.pre_rescue_schema_value <- function(evidence) {
    if (is.list(evidence) && "schema" %in% names(evidence)) {
        .sidecar_scalar_character(evidence$schema)
    } else {
        NA_character_
    }
}

.validate_pre_rescue_lifecycle <- function(evidence) {
    actual_schema <- .pre_rescue_schema_value(evidence)
    if (!identical(actual_schema, .PRE_RESCUE_EVIDENCE_SCHEMA)) {
        .abort_sidecar(
            paste0(
                "Unsupported pre-rescue evidence schema; reconstruct the ",
                "analysis with the current `analyze_missingness()`."
            ),
            "imputefinder_pre_rescue_lifecycle_error",
            expected_schema = .PRE_RESCUE_EVIDENCE_SCHEMA,
            actual_schema = actual_schema
        )
    }

    invisible(evidence)
}

.valid_pre_rescue_means <- function(mean_intensity, observed_count) {
    is.numeric(mean_intensity) &&
        length(mean_intensity) == length(observed_count) &&
        !any(is.nan(mean_intensity)) &&
        !any(is.infinite(mean_intensity)) &&
        identical(is.na(mean_intensity), observed_count == 0L)
}

.valid_pre_rescue_counts <- function(
    observed_count,
    missing_count,
    total_count,
    missing_fraction
) {
    is.integer(observed_count) &&
        is.integer(missing_count) &&
        is.integer(total_count) &&
        is.numeric(missing_fraction) &&
        length(observed_count) == length(total_count) &&
        length(missing_count) == length(total_count) &&
        length(missing_fraction) == length(total_count) &&
        !anyNA(observed_count) &&
        !anyNA(missing_count) &&
        !anyNA(total_count) &&
        !anyNA(missing_fraction) &&
        all(total_count > 0L) &&
        all(observed_count >= 0L) &&
        all(missing_count >= 0L) &&
        identical(observed_count + missing_count, total_count) &&
        identical(missing_fraction, missing_count / total_count)
}

.expected_pre_rescue_groups <- function(design, input) {
    condition_role <- design$declared$roles$condition
    stats::setNames(
        design$declared$sample_data[[condition_role]],
        input$sample_names
    )
}

.valid_pre_rescue_feature_condition <- function(
    feature_condition,
    input,
    groups_by_sample
) {
    conditions <- sort(unique(unname(groups_by_sample)), method = "radix")
    expected_feature <- rep(
        input$feature_names,
        each = length(conditions)
    )
    expected_condition <- rep(
        conditions,
        times = length(input$feature_names)
    )
    condition_sizes <- as.integer(table(factor(
        groups_by_sample,
        levels = conditions
    )))
    names(condition_sizes) <- conditions
    expected_sample_count <- unname(
        condition_sizes[expected_condition]
    )

    is.data.frame(feature_condition) &&
        identical(
            names(feature_condition),
            .PRE_RESCUE_FEATURE_CONDITION_FIELDS
        ) &&
        identical(feature_condition$feature, expected_feature) &&
        identical(feature_condition$condition, expected_condition) &&
        identical(
            feature_condition$sample_count,
            expected_sample_count
        ) &&
        .valid_pre_rescue_counts(
            feature_condition$observed_count,
            feature_condition$missing_count,
            feature_condition$sample_count,
            feature_condition$missing_fraction
        ) &&
        .valid_pre_rescue_means(
            feature_condition$mean_intensity,
            feature_condition$observed_count
        )
}

.valid_pre_rescue_sample <- function(sample, input, groups_by_sample) {
    expected_feature_count <- rep(
        input$dimensions[["features"]],
        input$dimensions[["samples"]]
    )

    is.data.frame(sample) &&
        identical(names(sample), .PRE_RESCUE_SAMPLE_FIELDS) &&
        identical(sample$sample, input$sample_names) &&
        identical(sample$condition, unname(groups_by_sample)) &&
        identical(sample$feature_count, expected_feature_count) &&
        .valid_pre_rescue_counts(
            sample$observed_count,
            sample$missing_count,
            sample$feature_count,
            sample$missing_fraction
        ) &&
        .valid_pre_rescue_means(
            sample$mean_intensity,
            sample$observed_count
        )
}

.valid_pre_rescue_cross_counts <- function(
    feature_condition,
    sample,
    groups_by_sample
) {
    conditions <- sort(unique(unname(groups_by_sample)), method = "radix")
    vapply(
        conditions,
        function(condition) {
            feature_rows <- feature_condition$condition == condition
            sample_rows <- sample$condition == condition
            identical(
                sum(feature_condition$observed_count[feature_rows]),
                sum(sample$observed_count[sample_rows])
            ) && identical(
                sum(feature_condition$missing_count[feature_rows]),
                sum(sample$missing_count[sample_rows])
            )
        },
        logical(1L)
    ) |> all()
}

.validate_pre_rescue_evidence <- function(evidence, input, design) {
    .validate_pre_rescue_lifecycle(evidence)
    groups_by_sample <- .expected_pre_rescue_groups(design, input)
    valid <- is.list(evidence) &&
        identical(names(evidence), .PRE_RESCUE_EVIDENCE_FIELDS) &&
        .valid_pre_rescue_feature_condition(
            evidence$feature_condition,
            input,
            groups_by_sample
        ) &&
        .valid_pre_rescue_sample(
            evidence$sample,
            input,
            groups_by_sample
        ) &&
        .valid_pre_rescue_cross_counts(
            evidence$feature_condition,
            evidence$sample,
            groups_by_sample
        )
    if (!valid) {
        .abort_sidecar(
            paste0(
                "Stored pre-rescue evidence is malformed or internally ",
                "inconsistent."
            ),
            "imputefinder_pre_rescue_schema_error",
            field = "sentinel.pre_rescue"
        )
    }

    invisible(evidence)
}

.new_sidecar_sentinel <- function(data, groups_by_sample) {
    list(
        pre_rescue = .new_pre_rescue_evidence(data, groups_by_sample)
    )
}

.validate_sidecar_sentinel <- function(sentinel, input, design) {
    if (is.null(sentinel)) {
        return(invisible(sentinel))
    }
    if (!is.list(sentinel) || !identical(names(sentinel), "pre_rescue")) {
        .abort_sidecar(
            "Stored sentinel output does not satisfy its current schema.",
            "imputefinder_analysis_schema_error",
            field = "sentinel"
        )
    }
    .validate_pre_rescue_evidence(sentinel$pre_rescue, input, design)

    invisible(sentinel)
}
