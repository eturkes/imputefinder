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

.validate_matrix_data <- function(x) {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("`x` must be an ordinary numeric matrix.", call. = FALSE)
    }

    .validate_axis_names(rownames(x), nrow(x), "feature")
    .validate_axis_names(colnames(x), ncol(x), "sample")

    if (any(is.nan(x))) {
        stop("`x` must not contain NaN values.", call. = FALSE)
    }
    if (any(is.infinite(x))) {
        stop("`x` must not contain Inf or -Inf values.", call. = FALSE)
    }

    invisible(x)
}

.validate_axis_names <- function(x, expected_length, axis) {
    invalid <- expected_length == 0L ||
        is.null(x) ||
        length(x) != expected_length ||
        anyNA(x) ||
        any(!nzchar(x)) ||
        anyDuplicated(x)

    if (invalid) {
        stop(
            sprintf("`x` must have non-empty, unique %s names.", axis),
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
                paste0(
                    "Named `group` must have unique names matching all ",
                    "sample names exactly."
                ),
                call. = FALSE
            )
        }
        group <- group[sample_names]
    }

    .normalise_condition_labels(group, sample_names)
}

.groups_from_design <- function(design, group_col, sample_names) {
    column_index <- if (
        is.character(group_col) &&
            length(group_col) == 1L &&
            !is.na(group_col) &&
            nzchar(group_col)
    ) {
        which(names(design) == group_col)
    } else {
        integer()
    }

    if (length(column_index) != 1L) {
        stop(
            "`group_col` must name exactly one design column.",
            call. = FALSE
        )
    }

    design_names <- rownames(design)
    valid_names <- length(design_names) == length(sample_names) &&
        !anyNA(design_names) &&
        all(nzchar(design_names)) &&
        !anyDuplicated(design_names) &&
        all(sample_names %in% design_names)

    if (!valid_names) {
        stop(
            paste0(
                "Design row names must be unique and match all sample ",
                "names exactly."
            ),
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
