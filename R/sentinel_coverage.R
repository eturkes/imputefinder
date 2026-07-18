.SENTINEL_COVERAGE_SCHEMA <- "sentinel_static_coverage_v1"
.SENTINEL_COVERAGE_FIELDS <- c(
    "schema",
    "sample",
    "condition",
    "role_level",
    "condition_role",
    "feature_overlap"
)
.SENTINEL_COVERAGE_SAMPLE_FIELDS <- c(
    "sample",
    "condition",
    "input_feature_count",
    "globally_observable_feature_count",
    "detected_feature_count",
    "detection_fraction",
    "mean_intensity",
    "minimum_intensity",
    "median_intensity",
    "maximum_intensity"
)
.SENTINEL_COVERAGE_CONDITION_FIELDS <- c(
    "condition",
    "sample_count",
    "input_feature_count",
    "globally_observable_feature_count",
    "detected_feature_count",
    "complete_feature_count",
    "observed_cell_count",
    "eligible_cell_count",
    "detection_fraction"
)
.SENTINEL_COVERAGE_ROLE_LEVEL_FIELDS <- c(
    "role",
    "column",
    "encoding",
    "level",
    "numeric_value",
    "sample_count",
    "condition_count",
    "singleton"
)
.SENTINEL_COVERAGE_CONDITION_ROLE_FIELDS <- c(
    "role",
    "column",
    "encoding",
    "level",
    "numeric_value",
    "condition",
    "sample_count",
    "globally_observable_feature_count",
    "detected_feature_count",
    "observed_cell_count",
    "eligible_cell_count",
    "detection_fraction",
    "empty"
)
.SENTINEL_COVERAGE_OVERLAP_FIELDS <- c(
    "condition_left",
    "condition_right",
    "globally_observable_feature_count",
    "left_detected_count",
    "right_detected_count",
    "shared_detected_count",
    "union_detected_count",
    "left_only_count",
    "right_only_count",
    "neither_count",
    "jaccard"
)

.sentinel_fraction <- function(numerator, denominator) {
    if (denominator == 0L) NA_real_ else numerator / denominator
}

.sentinel_coverage_context <- function(data, design) {
    .validate_matrix_data(data)
    aligned <- .align_missingness_design(design, colnames(data))
    canonical_features <- .canonical_sidecar_names(rownames(data))
    canonical_names <- .canonical_sidecar_names(colnames(data))
    feature_order <- order(canonical_features, method = "radix")
    sample_order <- order(canonical_names, method = "radix")
    data <- data[feature_order, sample_order, drop = FALSE]
    rownames(data) <- canonical_features[feature_order]
    colnames(data) <- canonical_names[sample_order]
    aligned <- .align_missingness_design(aligned, colnames(data))
    condition_column <- aligned$roles$condition
    condition <- .canonical_design_text(
        aligned$sample_data[[condition_column]]
    )
    observed <- !is.na(data)

    list(
        data = data,
        design = aligned,
        sample = colnames(data),
        condition = condition,
        conditions = sort(unique(condition), method = "radix"),
        observed = observed,
        globally_observable = rowSums(observed) > 0L
    )
}

.sentinel_intensity_range <- function(data) {
    output <- matrix(
        NA_real_,
        nrow = ncol(data),
        ncol = 3L,
        dimnames = list(NULL, c("minimum", "median", "maximum"))
    )
    for (index in seq_len(ncol(data))) {
        values <- data[, index]
        values <- values[!is.na(values)]
        if (length(values)) {
            output[index, ] <- c(
                min(values),
                stats::median(values),
                max(values)
            )
        }
    }
    output
}

.sentinel_sample_coverage <- function(context) {
    input_count <- as.integer(nrow(context$data))
    eligible_count <- as.integer(sum(context$globally_observable))
    detected <- as.integer(colSums(
        context$observed[context$globally_observable, , drop = FALSE]
    ))
    means <- .pre_rescue_observed_summary(context$data, 2L)$mean_intensity
    intensity <- .sentinel_intensity_range(context$data)

    data.frame(
        sample = context$sample,
        condition = context$condition,
        input_feature_count = rep(input_count, length(context$sample)),
        globally_observable_feature_count = rep(
            eligible_count,
            length(context$sample)
        ),
        detected_feature_count = detected,
        detection_fraction = vapply(
            detected,
            .sentinel_fraction,
            numeric(1L),
            denominator = eligible_count
        ),
        mean_intensity = means,
        minimum_intensity = intensity[, "minimum"],
        median_intensity = intensity[, "median"],
        maximum_intensity = intensity[, "maximum"],
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

.sentinel_condition_row <- function(condition, context) {
    selected <- context$condition == condition
    sample_count <- as.integer(sum(selected))
    eligible_count <- as.integer(sum(context$globally_observable))
    block <- context$observed[
        context$globally_observable,
        selected,
        drop = FALSE
    ]
    observed_per_feature <- rowSums(block)
    observed_cells <- as.integer(sum(block))
    eligible_cells <- as.integer(eligible_count * sample_count)

    data.frame(
        condition = condition,
        sample_count = sample_count,
        input_feature_count = as.integer(nrow(context$data)),
        globally_observable_feature_count = eligible_count,
        detected_feature_count = as.integer(sum(observed_per_feature > 0L)),
        complete_feature_count = as.integer(
            sum(observed_per_feature == sample_count)
        ),
        observed_cell_count = observed_cells,
        eligible_cell_count = eligible_cells,
        detection_fraction = .sentinel_fraction(
            observed_cells,
            eligible_cells
        ),
        stringsAsFactors = FALSE
    )
}

.sentinel_condition_coverage <- function(context) {
    output <- lapply(
        context$conditions,
        .sentinel_condition_row,
        context = context
    )
    do.call(rbind, unname(output)) |>
        structure(row.names = seq_along(context$conditions))
}

.sentinel_numeric_levels <- function(values) {
    values <- as.double(values)
    values[values == 0] <- 0
    levels <- sort(unique(values))
    list(
        encoding = "numeric",
        level = sprintf("%.17g", levels),
        numeric_value = levels,
        selectors = lapply(levels, function(level) values == level)
    )
}

.sentinel_categorical_levels <- function(values) {
    values <- .canonical_design_text(values)
    levels <- sort(unique(values), method = "radix")
    list(
        encoding = "treatment",
        level = levels,
        numeric_value = rep(NA_real_, length(levels)),
        selectors = lapply(levels, function(level) values == level)
    )
}

.sentinel_coverage_variable <- function(design, column, role) {
    values <- design$sample_data[[column]]
    levels <- if (!is.factor(values) && is.numeric(values)) {
        .sentinel_numeric_levels(values)
    } else {
        .sentinel_categorical_levels(values)
    }
    levels$role <- role
    levels$column <- .canonical_design_text(column)
    levels
}

.sentinel_coverage_variables <- function(context) {
    roles <- context$design$roles
    columns <- unname(unlist(roles, use.names = FALSE))
    role <- rep(names(roles), lengths(roles))
    variables <- lapply(seq_along(columns), function(index) {
        .sentinel_coverage_variable(
            context$design,
            columns[[index]],
            role[[index]]
        )
    })
    names(variables) <- .canonical_design_text(columns)
    variables
}

.empty_sentinel_role_level <- function() {
    data.frame(
        role = character(),
        column = character(),
        encoding = character(),
        level = character(),
        numeric_value = numeric(),
        sample_count = integer(),
        condition_count = integer(),
        singleton = logical(),
        stringsAsFactors = FALSE
    )
}

.sentinel_role_level_rows <- function(variable, condition) {
    rows <- lapply(seq_along(variable$level), function(index) {
        selected <- variable$selectors[[index]]
        count <- as.integer(sum(selected))
        data.frame(
            role = variable$role,
            column = variable$column,
            encoding = variable$encoding,
            level = variable$level[[index]],
            numeric_value = variable$numeric_value[[index]],
            sample_count = count,
            condition_count = as.integer(length(unique(condition[selected]))),
            singleton = count == 1L,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, unname(rows))
}

.sentinel_role_level_coverage <- function(context, variables) {
    if (!length(variables)) {
        return(.empty_sentinel_role_level())
    }
    rows <- lapply(
        variables,
        .sentinel_role_level_rows,
        condition = context$condition
    )
    output <- do.call(rbind, unname(rows))
    row.names(output) <- NULL
    output
}

.empty_sentinel_condition_role <- function() {
    data.frame(
        role = character(),
        column = character(),
        encoding = character(),
        level = character(),
        numeric_value = numeric(),
        condition = character(),
        sample_count = integer(),
        globally_observable_feature_count = integer(),
        detected_feature_count = integer(),
        observed_cell_count = integer(),
        eligible_cell_count = integer(),
        detection_fraction = numeric(),
        empty = logical(),
        stringsAsFactors = FALSE
    )
}

.sentinel_condition_role_cell <- function(
    context,
    variable,
    level_index,
    condition
) {
    selected <- variable$selectors[[level_index]] &
        context$condition == condition
    sample_count <- as.integer(sum(selected))
    eligible_count <- as.integer(sum(context$globally_observable))
    block <- context$observed[
        context$globally_observable,
        selected,
        drop = FALSE
    ]
    observed_cells <- as.integer(sum(block))
    eligible_cells <- as.integer(eligible_count * sample_count)

    data.frame(
        role = variable$role,
        column = variable$column,
        encoding = variable$encoding,
        level = variable$level[[level_index]],
        numeric_value = variable$numeric_value[[level_index]],
        condition = condition,
        sample_count = sample_count,
        globally_observable_feature_count = eligible_count,
        detected_feature_count = as.integer(sum(rowSums(block) > 0L)),
        observed_cell_count = observed_cells,
        eligible_cell_count = eligible_cells,
        detection_fraction = .sentinel_fraction(
            observed_cells,
            eligible_cells
        ),
        empty = sample_count == 0L,
        stringsAsFactors = FALSE
    )
}

.sentinel_condition_role_rows <- function(variable, context) {
    rows <- list()
    index <- 0L
    for (level_index in seq_along(variable$level)) {
        for (condition in context$conditions) {
            index <- index + 1L
            rows[[index]] <- .sentinel_condition_role_cell(
                context,
                variable,
                level_index,
                condition
            )
        }
    }
    do.call(rbind, rows)
}

.sentinel_condition_role_coverage <- function(context, variables) {
    variables <- variables[vapply(
        variables,
        function(variable) variable$role != "condition",
        logical(1L)
    )]
    if (!length(variables)) {
        return(.empty_sentinel_condition_role())
    }
    rows <- lapply(variables, .sentinel_condition_role_rows, context = context)
    output <- do.call(rbind, unname(rows))
    row.names(output) <- NULL
    output
}

.empty_sentinel_feature_overlap <- function() {
    data.frame(
        condition_left = character(),
        condition_right = character(),
        globally_observable_feature_count = integer(),
        left_detected_count = integer(),
        right_detected_count = integer(),
        shared_detected_count = integer(),
        union_detected_count = integer(),
        left_only_count = integer(),
        right_only_count = integer(),
        neither_count = integer(),
        jaccard = numeric(),
        stringsAsFactors = FALSE
    )
}

.sentinel_detected_by_condition <- function(context, condition) {
    selected <- context$condition == condition
    rowSums(context$observed[
        context$globally_observable,
        selected,
        drop = FALSE
    ]) > 0L
}

.sentinel_feature_overlap_row <- function(pair, context) {
    left <- .sentinel_detected_by_condition(context, pair[[1L]])
    right <- .sentinel_detected_by_condition(context, pair[[2L]])
    shared <- as.integer(sum(left & right))
    union <- as.integer(sum(left | right))
    eligible_count <- as.integer(sum(context$globally_observable))

    data.frame(
        condition_left = pair[[1L]],
        condition_right = pair[[2L]],
        globally_observable_feature_count = eligible_count,
        left_detected_count = as.integer(sum(left)),
        right_detected_count = as.integer(sum(right)),
        shared_detected_count = shared,
        union_detected_count = union,
        left_only_count = as.integer(sum(left & !right)),
        right_only_count = as.integer(sum(!left & right)),
        neither_count = as.integer(sum(!left & !right)),
        jaccard = .sentinel_fraction(shared, union),
        stringsAsFactors = FALSE
    )
}

.sentinel_feature_overlap <- function(context) {
    if (length(context$conditions) < 2L) {
        return(.empty_sentinel_feature_overlap())
    }
    pairs <- utils::combn(context$conditions, 2L, simplify = FALSE)
    rows <- lapply(pairs, .sentinel_feature_overlap_row, context = context)
    output <- do.call(rbind, unname(rows))
    row.names(output) <- NULL
    output
}

.new_sentinel_coverage <- function(data, design) {
    context <- .sentinel_coverage_context(data, design)
    variables <- .sentinel_coverage_variables(context)

    list(
        schema = .SENTINEL_COVERAGE_SCHEMA,
        sample = .sentinel_sample_coverage(context),
        condition = .sentinel_condition_coverage(context),
        role_level = .sentinel_role_level_coverage(context, variables),
        condition_role = .sentinel_condition_role_coverage(
            context,
            variables
        ),
        feature_overlap = .sentinel_feature_overlap(context)
    )
}

.sentinel_coverage_schema_value <- function(coverage) {
    if (is.list(coverage) && "schema" %in% names(coverage)) {
        .sidecar_scalar_character(coverage$schema)
    } else {
        NA_character_
    }
}

.validate_sentinel_coverage_lifecycle <- function(coverage) {
    actual_schema <- .sentinel_coverage_schema_value(coverage)
    if (!identical(actual_schema, .SENTINEL_COVERAGE_SCHEMA)) {
        .abort_sidecar(
            paste0(
                "Unsupported sentinel coverage schema; reconstruct the ",
                "analysis with the current `analyze_missingness()`."
            ),
            "imputefinder_sentinel_coverage_lifecycle_error",
            expected_schema = .SENTINEL_COVERAGE_SCHEMA,
            actual_schema = actual_schema
        )
    }
    invisible(coverage)
}

.valid_sentinel_coverage_shape <- function(coverage) {
    is.list(coverage) &&
        identical(names(coverage), .SENTINEL_COVERAGE_FIELDS) &&
        is.data.frame(coverage$sample) &&
        identical(
            names(coverage$sample),
            .SENTINEL_COVERAGE_SAMPLE_FIELDS
        ) &&
        is.data.frame(coverage$condition) &&
        identical(
            names(coverage$condition),
            .SENTINEL_COVERAGE_CONDITION_FIELDS
        ) &&
        is.data.frame(coverage$role_level) &&
        identical(
            names(coverage$role_level),
            .SENTINEL_COVERAGE_ROLE_LEVEL_FIELDS
        ) &&
        is.data.frame(coverage$condition_role) &&
        identical(
            names(coverage$condition_role),
            .SENTINEL_COVERAGE_CONDITION_ROLE_FIELDS
        ) &&
        is.data.frame(coverage$feature_overlap) &&
        identical(
            names(coverage$feature_overlap),
            .SENTINEL_COVERAGE_OVERLAP_FIELDS
        )
}

.sentinel_coverage_mask_data <- function(input) {
    missing <- .unpack_original_mask(
        input$original_mask,
        input$dimensions,
        input$feature_names,
        input$sample_names
    )
    data <- matrix(
        0,
        nrow = nrow(missing),
        ncol = ncol(missing),
        dimnames = dimnames(missing)
    )
    data[missing] <- NA_real_
    data
}

.valid_sentinel_intensity_support <- function(sample, pre_rescue) {
    intensity_fields <- c(
        "mean_intensity",
        "minimum_intensity",
        "median_intensity",
        "maximum_intensity"
    )
    intensity <- sample[intensity_fields]
    numeric_fields <- vapply(intensity, is.double, logical(1L)) |> all()
    finite_or_missing <- vapply(
        intensity,
        function(value) !any(is.nan(value)) && !any(is.infinite(value)),
        logical(1L)
    ) |> all()
    pre_order <- match(sample$sample, pre_rescue$sample$sample)
    mean_matches <- !anyNA(pre_order) && identical(
        sample$mean_intensity,
        pre_rescue$sample$mean_intensity[pre_order]
    )
    missing <- sample$detected_feature_count == 0L
    all_missing <- rowSums(is.na(intensity)) == length(intensity_fields)
    complete <- rowSums(is.na(intensity)) == 0L
    ordered <- missing | (
        sample$minimum_intensity <= sample$median_intensity &
            sample$median_intensity <= sample$maximum_intensity &
            sample$minimum_intensity <= sample$mean_intensity &
            sample$mean_intensity <= sample$maximum_intensity
    )
    singleton <- sample$detected_feature_count != 1L | (
        sample$minimum_intensity == sample$median_intensity &
            sample$median_intensity == sample$maximum_intensity &
            sample$maximum_intensity == sample$mean_intensity
    )

    numeric_fields && finite_or_missing && mean_matches &&
        identical(all_missing, missing) &&
        identical(complete, !missing) &&
        all(ordered) && all(singleton)
}

.sentinel_coverage_structural_sample <- function(sample) {
    sample[seq_len(6L)]
}

.valid_sentinel_coverage_contents <- function(
    coverage,
    input,
    design,
    pre_rescue
) {
    expected <- .new_sentinel_coverage(
        .sentinel_coverage_mask_data(input),
        design$declared
    )
    identical(
        .sentinel_coverage_structural_sample(coverage$sample),
        .sentinel_coverage_structural_sample(expected$sample)
    ) &&
        .valid_sentinel_intensity_support(coverage$sample, pre_rescue) &&
        identical(coverage$condition, expected$condition) &&
        identical(coverage$role_level, expected$role_level) &&
        identical(coverage$condition_role, expected$condition_role) &&
        identical(coverage$feature_overlap, expected$feature_overlap)
}

.validate_sentinel_coverage <- function(
    coverage,
    input,
    design,
    pre_rescue
) {
    .validate_sentinel_coverage_lifecycle(coverage)
    valid <- .valid_sentinel_coverage_shape(coverage) &&
        .valid_sentinel_coverage_contents(
            coverage,
            input,
            design,
            pre_rescue
        )
    if (!valid) {
        .abort_sidecar(
            paste0(
                "Stored sentinel coverage is malformed or internally ",
                "inconsistent."
            ),
            "imputefinder_sentinel_coverage_schema_error",
            field = "sentinel.coverage"
        )
    }
    invisible(coverage)
}
