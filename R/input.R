.prepare_input <- function(
    x,
    group = NULL,
    group_col = NULL,
    assay = NULL
) {
    if (is.matrix(x)) {
        return(.prepare_matrix_input(x, group, group_col, assay))
    }

    if (methods::is(x, "SummarizedExperiment")) {
        return(
            .prepare_summarized_experiment_input(
                x,
                group,
                group_col,
                assay
            )
        )
    }

    stop(
        "`x` must be an ordinary numeric matrix or a ",
        "SummarizedExperiment object.",
        call. = FALSE
    )
}

.prepare_matrix_input <- function(x, group, group_col = NULL, assay = NULL) {
    .validate_matrix_data(x)

    if (!is.null(assay)) {
        stop("`assay` must be NULL for matrix input.", call. = FALSE)
    }
    if (is.null(group)) {
        stop("`group` is required for matrix input.", call. = FALSE)
    }

    if (is.data.frame(group)) {
        groups_by_sample <- .groups_from_design(
            design = group,
            group_col = group_col,
            sample_names = colnames(x)
        )
    } else if (is.atomic(group) && is.null(dim(group))) {
        groups_by_sample <- .groups_from_vector(
            group = group,
            group_col = group_col,
            sample_names = colnames(x)
        )
    } else {
        stop(
            "`group` must be an atomic vector or data frame for matrix input.",
            call. = FALSE
        )
    }

    list(
        data = x,
        groups_by_sample = groups_by_sample,
        representation = "matrix"
    )
}

.prepare_summarized_experiment_input <- function(
    x,
    group,
    group_col,
    assay
) {
    if (!is.null(group)) {
        stop(
            "`group` must be NULL for SummarizedExperiment input.",
            call. = FALSE
        )
    }

    assay_count <- length(SummarizedExperiment::assays(x))
    if (assay_count == 0L) {
        stop("`x` must contain at least one assay.", call. = FALSE)
    }

    assay_names <- SummarizedExperiment::assayNames(x)
    assay_index <- .resolve_assay_index(assay, assay_names, assay_count)
    data <- SummarizedExperiment::assay(x, assay_index)
    .validate_matrix_data(data, subject = "The selected assay")

    column_index <- .resolve_group_column(
        group_col,
        names(SummarizedExperiment::colData(x)),
        source = "colData"
    )
    condition <- SummarizedExperiment::colData(x)[[column_index]]
    if (!is.atomic(condition) || !is.null(dim(condition))) {
        stop(
            "The selected colData condition column must be an atomic vector.",
            call. = FALSE
        )
    }

    list(
        data = data,
        groups_by_sample = .normalise_condition_labels(
            condition,
            colnames(data)
        ),
        representation = "SummarizedExperiment",
        assay_index = assay_index,
        assay_name = .selected_assay_name(assay_names, assay_index)
    )
}

.selected_assay_name <- function(assay_names, assay_index) {
    named <- length(assay_names) >= assay_index &&
        !is.na(assay_names[[assay_index]]) &&
        nzchar(assay_names[[assay_index]])
    if (named) assay_names[[assay_index]] else NA_character_
}

.resolve_assay_index <- function(assay, assay_names, assay_count) {
    if (is.null(assay)) {
        if (assay_count == 1L) {
            return(1L)
        }

        stop(
            "`assay` must name one assay when `x` contains multiple assays.",
            call. = FALSE
        )
    }

    valid_assay <- is.character(assay) &&
        length(assay) == 1L &&
        !is.na(assay) &&
        nzchar(assay)
    if (!valid_assay) {
        stop(
            "`assay` must be NULL or one non-empty assay name.",
            call. = FALSE
        )
    }

    matches <- which(!is.na(assay_names) & assay_names == assay)
    if (length(matches) != 1L) {
        stop("`assay` must match exactly one assay name.", call. = FALSE)
    }

    as.integer(matches)
}

.resolve_group_column <- function(group_col, column_names, source) {
    valid_name <- is.character(group_col) &&
        length(group_col) == 1L &&
        !is.na(group_col) &&
        nzchar(group_col)
    matches <- if (valid_name) {
        which(column_names == group_col)
    } else {
        integer()
    }

    if (length(matches) != 1L) {
        stop(
            sprintf(
                "`group_col` must name exactly one %s column.",
                source
            ),
            call. = FALSE
        )
    }

    matches
}

.validate_matrix_data <- function(x, subject = "`x`") {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop(
            sprintf("%s must be an ordinary numeric matrix.", subject),
            call. = FALSE
        )
    }

    .validate_axis_names(rownames(x), nrow(x), "feature", subject)
    .validate_axis_names(colnames(x), ncol(x), "sample", subject)

    if (any(is.nan(x))) {
        stop(
            sprintf("%s must not contain NaN values.", subject),
            call. = FALSE
        )
    }
    if (any(is.infinite(x))) {
        stop(
            sprintf("%s must not contain Inf or -Inf values.", subject),
            call. = FALSE
        )
    }

    invisible(x)
}

.validate_axis_names <- function(x, expected_length, axis, subject = "`x`") {
    invalid <- expected_length == 0L ||
        is.null(x) ||
        length(x) != expected_length ||
        anyNA(x) ||
        any(!nzchar(x)) ||
        anyDuplicated(x)

    if (invalid) {
        stop(
            sprintf(
                "%s must have non-empty, unique %s names.",
                subject,
                axis
            ),
            call. = FALSE
        )
    }

    invisible(x)
}

.groups_from_vector <- function(group, group_col, sample_names) {
    if (!is.null(group_col)) {
        stop(
            "`group_col` must be NULL when `group` is an atomic vector.",
            call. = FALSE
        )
    }

    group_names <- names(group)
    if (is.null(group_names)) {
        if (length(group) != length(sample_names)) {
            stop(
                "Unnamed `group` must have one value per sample.",
                call. = FALSE
            )
        }
    } else {
        valid_names <- length(group) == length(sample_names) &&
            !anyNA(group_names) &&
            all(nzchar(group_names)) &&
            !anyDuplicated(group_names) &&
            all(sample_names %in% group_names)

        if (!valid_names) {
            stop(
                "Named `group` must have unique names matching all ",
                "sample names exactly.",
                call. = FALSE
            )
        }
        group <- group[sample_names]
    }

    .normalise_condition_labels(group, sample_names)
}

.groups_from_design <- function(design, group_col, sample_names) {
    column_index <- .resolve_group_column(
        group_col,
        names(design),
        source = "design"
    )

    design_names <- rownames(design)
    valid_names <- length(design_names) == length(sample_names) &&
        !anyNA(design_names) &&
        all(nzchar(design_names)) &&
        !anyDuplicated(design_names) &&
        all(sample_names %in% design_names)

    if (!valid_names) {
        stop(
            "Design row names must be unique and match all sample ",
            "names exactly.",
            call. = FALSE
        )
    }

    condition <- design[[column_index]][match(sample_names, design_names)]
    if (!is.atomic(condition) || !is.null(dim(condition))) {
        stop(
            "The selected design condition column must be an atomic vector.",
            call. = FALSE
        )
    }

    .normalise_condition_labels(condition, sample_names)
}

.normalise_condition_labels <- function(group, sample_names) {
    if (length(group) != length(sample_names) || anyNA(group)) {
        stop("Condition labels must not be missing or empty.", call. = FALSE)
    }

    labels <- as.character(group)
    if (anyNA(labels) || any(!nzchar(labels))) {
        stop("Condition labels must not be missing or empty.", call. = FALSE)
    }

    stats::setNames(labels, sample_names)
}

.restore_output_data <- function(data, prepared, original) {
    .validate_core_output(data, prepared$data)

    if (identical(prepared$representation, "matrix")) {
        return(data)
    }

    if (!identical(prepared$representation, "SummarizedExperiment") ||
        !methods::is(original, "SummarizedExperiment")) {
        stop("Unsupported input representation provenance.", call. = FALSE)
    }

    row_index <- if (nrow(data) == 0L) integer() else rownames(data)
    output <- original[row_index, , drop = FALSE]
    SummarizedExperiment::assay(output, prepared$assay_index) <- data
    output
}

.validate_core_output <- function(data, source_data) {
    if (!is.matrix(data) || !is.numeric(data)) {
        stop("Core output data must be a numeric matrix.", call. = FALSE)
    }
    if (any(is.nan(data)) || any(is.infinite(data))) {
        stop(
            "Core output data must contain only finite values or NA.",
            call. = FALSE
        )
    }

    source_features <- rownames(source_data)
    output_features <- rownames(data)
    empty_output <- nrow(data) == 0L && is.null(output_features)
    valid_feature_names <- empty_output || (
        !is.null(output_features) &&
            !anyNA(output_features) &&
            all(nzchar(output_features)) &&
            !anyDuplicated(output_features)
    )
    expected_features <- source_features[source_features %in% output_features]
    order_matches <- empty_output ||
        identical(output_features, expected_features)
    if (!valid_feature_names || !order_matches) {
        stop(
            "Core output feature names must preserve original order.",
            call. = FALSE
        )
    }

    if (!identical(colnames(data), colnames(source_data))) {
        stop(
            "Core output sample names must match the input order.",
            call. = FALSE
        )
    }

    invisible(data)
}
