.MISSINGNESS_DESIGN_SCHEMA <- "missingness_design_v1"
.MISSINGNESS_DESIGN_LIFECYCLE <- "experimental"
.MISSINGNESS_DESIGN_FIELDS <- c(
    "spec",
    "sample_data",
    "roles",
    "interactions"
)
.MISSINGNESS_DESIGN_ROLES <- c(
    "condition",
    "nuisance",
    "block",
    "acquisition"
)

#' Declare a Missingness Analysis Design
#'
#' Construct the typed sample design used by ImputeFinder's experimental
#' sidecar workflow. The constructor records declared roles only. It does not
#' fit a model, search for interactions, exclude samples, or decide whether a
#' contrast is estimable.
#'
#' @param sample_data An ordinary data frame with one row per sample and
#'   explicit, non-empty, unique sample row names. Unselected columns are not
#'   retained in the design object.
#' @param condition One column name declaring the biological condition.
#' @param nuisance Optional unique column names declaring technical or other
#'   adjustment terms, such as batch, plate, or run order.
#' @param block An optional single column name declaring subjects, pairs, or
#'   another repeated-measure block.
#' @param acquisition An optional single column name declaring acquisition
#'   mode.
#' @param interactions Optional list of character vectors. Each vector names at
#'   least two distinct declared role columns participating in one explicitly
#'   requested interaction. Ordering is canonicalized; duplicated interaction
#'   sets are rejected.
#'
#' @return An experimental `missingness_design` object containing its
#'   schema/lifecycle specification, declared sample metadata, role map, and
#'   canonical interaction sets. Condition, block, and acquisition identifiers
#'   are stored as character vectors; nuisance-column types are preserved.
#'
#' @details
#' This is an experimental schema: serialized objects identify
#' `missingness_design_v1`, and consumers reject an unknown schema or
#' lifecycle instead of silently upgrading it. Structural support, rank,
#' aliasing, and estimability are assessed by later workflow stages; therefore
#' even a one-sample or one-condition design can be represented here.
#'
#' Construction errors inherit from `imputefinder_design_error` and use
#' the `imputefinder_design_schema_error` or
#' `imputefinder_design_role_error` subclass.
#'
#' @examples
#' sample_data <- data.frame(
#'     condition = c("control", "control", "treated", "treated"),
#'     batch = c("one", "two", "one", "two"),
#'     subject = c("p1", "p2", "p1", "p2"),
#'     row.names = paste0("sample", 1:4)
#' )
#'
#' design <- missingness_design(
#'     sample_data,
#'     condition = "condition",
#'     nuisance = "batch",
#'     block = "subject",
#'     interactions = list(c("condition", "batch"))
#' )
#' design
#'
#' @export
missingness_design <- function(
    sample_data,
    condition,
    nuisance = NULL,
    block = NULL,
    acquisition = NULL,
    interactions = NULL
) {
    .validate_design_sample_data(sample_data)
    roles <- .resolve_missingness_design_roles(
        sample_data,
        condition,
        nuisance,
        block,
        acquisition
    )
    role_columns <- unname(unlist(roles, use.names = FALSE))
    .validate_design_role_values(sample_data, role_columns)
    interactions <- .normalise_design_interactions(
        interactions,
        role_columns
    )

    declared_data <- sample_data[, role_columns, drop = FALSE]
    character_roles <- c(
        roles$condition,
        roles$block,
        roles$acquisition
    )
    for (column in character_roles) {
        declared_data[[column]] <- as.character(declared_data[[column]])
    }

    structure(
        list(
            spec = list(
                schema = .MISSINGNESS_DESIGN_SCHEMA,
                lifecycle = .MISSINGNESS_DESIGN_LIFECYCLE
            ),
            sample_data = declared_data,
            roles = roles,
            interactions = interactions
        ),
        class = "missingness_design"
    )
}

.abort_missingness_design <- function(message, subclass, ...) {
    error <- structure(
        c(
            list(message = message, call = NULL),
            list(...)
        ),
        class = c(
            subclass,
            "imputefinder_design_error",
            "error",
            "condition"
        )
    )
    stop(error)
}

.validate_design_sample_data <- function(sample_data) {
    if (!is.data.frame(sample_data)) {
        .abort_missingness_design(
            "`sample_data` must be an ordinary data frame.",
            "imputefinder_design_schema_error",
            field = "sample_data"
        )
    }

    sample_names <- rownames(sample_data)
    automatic_names <- .row_names_info(sample_data, 1L) < 0L
    valid_sample_names <- nrow(sample_data) > 0L &&
        !automatic_names &&
        length(sample_names) == nrow(sample_data) &&
        !anyNA(sample_names) &&
        all(nzchar(sample_names)) &&
        !anyDuplicated(sample_names)
    if (!valid_sample_names) {
        .abort_missingness_design(
            paste0(
                "`sample_data` must have explicit, non-empty, unique ",
                "sample row names."
            ),
            "imputefinder_design_schema_error",
            field = "row.names"
        )
    }

    column_names <- names(sample_data)
    valid_column_names <- !is.null(column_names) &&
        length(column_names) == ncol(sample_data) &&
        !anyNA(column_names) &&
        all(nzchar(column_names)) &&
        !anyDuplicated(column_names)
    if (!valid_column_names) {
        .abort_missingness_design(
            paste0(
                "`sample_data` must have non-empty, unique column names."
            ),
            "imputefinder_design_schema_error",
            field = "names"
        )
    }

    invisible(sample_data)
}

.valid_design_role_selector <- function(value, required, multiple) {
    valid <- is.character(value) &&
        !anyNA(value) &&
        all(nzchar(value)) &&
        !anyDuplicated(value)
    if (required) {
        return(valid && length(value) == 1L)
    } else if (!multiple) {
        return(valid && length(value) <= 1L)
    }

    valid
}

.design_role_cardinality <- function(required, multiple) {
    if (required) {
        "exactly one"
    } else if (multiple) {
        "zero or more unique"
    } else {
        "zero or one"
    }
}

.validate_design_role_selector <- function(
    value,
    role,
    required,
    multiple
) {
    if (!.valid_design_role_selector(value, required, multiple)) {
        cardinality <- .design_role_cardinality(required, multiple)
        .abort_missingness_design(
            sprintf(
                "`%s` must contain %s non-empty column name%s.",
                role,
                cardinality,
                if (identical(cardinality, "exactly one")) "" else "s"
            ),
            "imputefinder_design_role_error",
            role = role,
            columns = if (is.character(value)) value else character()
        )
    }

    invisible(value)
}

.normalise_design_role <- function(
    value,
    role,
    column_names,
    required = FALSE,
    multiple = FALSE
) {
    if (is.null(value) && !required) {
        return(character())
    }
    .validate_design_role_selector(value, role, required, multiple)

    value <- unname(value)
    absent <- setdiff(value, column_names)
    if (length(absent) > 0L) {
        .abort_missingness_design(
            sprintf(
                "`%s` must name sample metadata columns exactly: %s.",
                role,
                paste(sprintf("`%s`", absent), collapse = ", ")
            ),
            "imputefinder_design_role_error",
            role = role,
            columns = absent
        )
    }

    if (multiple) sort(value, method = "radix") else value
}

.resolve_missingness_design_roles <- function(
    sample_data,
    condition,
    nuisance,
    block,
    acquisition
) {
    column_names <- names(sample_data)
    roles <- list(
        condition = .normalise_design_role(
            condition,
            "condition",
            column_names,
            required = TRUE
        ),
        nuisance = .normalise_design_role(
            nuisance,
            "nuisance",
            column_names,
            multiple = TRUE
        ),
        block = .normalise_design_role(
            block,
            "block",
            column_names
        ),
        acquisition = .normalise_design_role(
            acquisition,
            "acquisition",
            column_names
        )
    )
    role_columns <- unname(unlist(roles, use.names = FALSE))
    duplicated <- unique(role_columns[duplicated(role_columns)])
    if (length(duplicated) > 0L) {
        .abort_missingness_design(
            paste0(
                "Each metadata column must have exactly one design role: ",
                paste(sprintf("`%s`", duplicated), collapse = ", "),
                "."
            ),
            "imputefinder_design_role_error",
            role = "overlap",
            columns = duplicated
        )
    }

    roles
}

.supported_design_role_column <- function(x) {
    is.atomic(x) &&
        is.null(dim(x)) &&
        (
            is.factor(x) ||
                is.character(x) ||
                is.logical(x) ||
                is.integer(x) ||
                is.double(x)
        )
}

.validate_design_role_values <- function(
    sample_data,
    role_columns,
    subclass = "imputefinder_design_role_error"
) {
    for (column in role_columns) {
        values <- sample_data[[column]]
        supported <- .supported_design_role_column(values)
        complete <- supported &&
            !anyNA(values) &&
            all(nzchar(as.character(values)))
        finite <- !supported ||
            !is.numeric(values) ||
            all(is.finite(values))
        if (!supported || !complete || !finite) {
            .abort_missingness_design(
                paste0(
                    "Declared role columns must not contain missing or empty ",
                    "values; they must be atomic vectors and numeric values ",
                    "must be finite: ",
                    sprintf("`%s`.", column)
                ),
                subclass,
                role = "values",
                columns = column
            )
        }
    }

    invisible(sample_data)
}

.design_interaction_key <- function(term) {
    paste0(
        sprintf(
            "%d:%s",
            nchar(term, type = "bytes"),
            term
        ),
        collapse = ""
    )
}

.normalise_design_interaction <- function(term, role_columns, subclass) {
    valid <- is.character(term) &&
        length(term) >= 2L &&
        !anyNA(term) &&
        all(nzchar(term)) &&
        !anyDuplicated(term) &&
        all(term %in% role_columns)
    if (!valid) {
        .abort_missingness_design(
            paste0(
                "Each interaction must contain at least two distinct ",
                "declared role-column names."
            ),
            subclass,
            role = "interactions",
            columns = if (is.character(term)) term else character()
        )
    }

    sort(unname(term), method = "radix")
}

.reject_duplicated_design_interactions <- function(
    normalised,
    keys,
    subclass
) {
    if (!anyDuplicated(keys)) {
        return(invisible(normalised))
    }

    duplicated <- normalised[[which(duplicated(keys))[[1L]]]]
    .abort_missingness_design(
        paste0(
            "`interactions` must not repeat the same role-column set: ",
            paste(sprintf("`%s`", duplicated), collapse = ", "),
            "."
        ),
        subclass,
        role = "interactions",
        columns = duplicated
    )
}

.normalise_design_interactions <- function(
    interactions,
    role_columns,
    subclass = "imputefinder_design_role_error"
) {
    if (is.null(interactions)) {
        return(list())
    }
    if (!is.list(interactions) || is.data.frame(interactions)) {
        .abort_missingness_design(
            "`interactions` must be NULL or a list of role-column sets.",
            subclass,
            role = "interactions",
            columns = character()
        )
    }

    normalised <- lapply(
        unname(interactions),
        function(term) {
            .normalise_design_interaction(term, role_columns, subclass)
        }
    )
    keys <- vapply(normalised, .design_interaction_key, character(1))
    .reject_duplicated_design_interactions(normalised, keys, subclass)

    if (length(keys) == 0L) list() else
        normalised[order(keys, method = "radix")]
}

.missingness_design_spec_value <- function(spec, field) {
    if (!is.list(spec) || !field %in% names(spec)) {
        return(NA_character_)
    }
    value <- spec[[field]]
    if (is.character(value) && length(value) == 1L && !is.na(value)) {
        value
    } else {
        NA_character_
    }
}

.validate_missingness_design_lifecycle <- function(design) {
    if (!inherits(design, "missingness_design")) {
        .abort_missingness_design(
            "`design` must inherit from `missingness_design`.",
            "imputefinder_design_lifecycle_error",
            expected_schema = .MISSINGNESS_DESIGN_SCHEMA,
            actual_schema = NA_character_,
            expected_lifecycle = .MISSINGNESS_DESIGN_LIFECYCLE,
            actual_lifecycle = NA_character_
        )
    }

    spec <- if (is.list(design) && "spec" %in% names(design)) {
        design$spec
    } else {
        NULL
    }
    actual_schema <- .missingness_design_spec_value(spec, "schema")
    actual_lifecycle <- .missingness_design_spec_value(spec, "lifecycle")
    valid_spec <- is.list(spec) &&
        identical(names(spec), c("schema", "lifecycle")) &&
        identical(actual_schema, .MISSINGNESS_DESIGN_SCHEMA) &&
        identical(actual_lifecycle, .MISSINGNESS_DESIGN_LIFECYCLE)
    if (!valid_spec) {
        .abort_missingness_design(
            paste0(
                "Unsupported missingness-design schema or lifecycle; ",
                "reconstruct it with the current `missingness_design()`."
            ),
            "imputefinder_design_lifecycle_error",
            expected_schema = .MISSINGNESS_DESIGN_SCHEMA,
            actual_schema = actual_schema,
            expected_lifecycle = .MISSINGNESS_DESIGN_LIFECYCLE,
            actual_lifecycle = actual_lifecycle
        )
    }

    invisible(design)
}

.missingness_design_role_columns <- function(design) {
    valid_shape <- is.list(design) &&
        identical(names(design), .MISSINGNESS_DESIGN_FIELDS) &&
        identical(class(design), "missingness_design") &&
        is.list(design$roles) &&
        identical(names(design$roles), .MISSINGNESS_DESIGN_ROLES)
    if (!valid_shape) {
        .abort_missingness_design(
            "`design` does not satisfy the missingness-design object schema.",
            "imputefinder_design_schema_error",
            field = "object"
        )
    }

    roles <- design$roles
    valid_roles <- is.character(roles$condition) &&
        length(roles$condition) == 1L &&
        !anyNA(roles$condition) &&
        nzchar(roles$condition) &&
        is.character(roles$nuisance) &&
        !anyNA(roles$nuisance) &&
        all(nzchar(roles$nuisance)) &&
        !anyDuplicated(roles$nuisance) &&
        identical(roles$nuisance, sort(roles$nuisance, method = "radix")) &&
        is.character(roles$block) &&
        length(roles$block) <= 1L &&
        !anyNA(roles$block) &&
        all(nzchar(roles$block)) &&
        is.character(roles$acquisition) &&
        length(roles$acquisition) <= 1L &&
        !anyNA(roles$acquisition) &&
        all(nzchar(roles$acquisition))
    role_columns <- if (valid_roles) {
        unname(unlist(roles, use.names = FALSE))
    } else {
        character()
    }
    valid_roles <- valid_roles &&
        !anyDuplicated(role_columns) &&
        identical(names(design$sample_data), role_columns)
    if (!valid_roles) {
        .abort_missingness_design(
            "Stored roles do not satisfy the missingness-design schema.",
            "imputefinder_design_schema_error",
            field = "roles"
        )
    }

    role_columns
}

.validate_missingness_design_data <- function(design, role_columns) {
    .validate_design_sample_data(design$sample_data)
    roles <- design$roles

    .validate_design_role_values(
        design$sample_data,
        role_columns,
        subclass = "imputefinder_design_schema_error"
    )
    character_roles <- c(roles$condition, roles$block, roles$acquisition)
    if (!all(vapply(
        design$sample_data[character_roles],
        is.character,
        logical(1)
    ))) {
        .abort_missingness_design(
            paste0(
                "Stored condition, block, and acquisition roles must use ",
                "canonical character vectors."
            ),
            "imputefinder_design_schema_error",
            field = "sample_data"
        )
    }

    invisible(design)
}

.validate_missingness_design_interactions <- function(
    design,
    role_columns
) {
    canonical_interactions <- .normalise_design_interactions(
        design$interactions,
        role_columns,
        subclass = "imputefinder_design_schema_error"
    )
    if (!identical(design$interactions, canonical_interactions)) {
        .abort_missingness_design(
            "Stored interactions are not in canonical design order.",
            "imputefinder_design_schema_error",
            field = "interactions"
        )
    }

    invisible(design)
}

.validate_missingness_design <- function(design) {
    .validate_missingness_design_lifecycle(design)
    role_columns <- .missingness_design_role_columns(design)
    .validate_missingness_design_data(design, role_columns)
    .validate_missingness_design_interactions(design, role_columns)

    invisible(design)
}

.validate_design_alignment_names <- function(sample_names) {
    valid_sample_names <- is.character(sample_names) &&
        length(sample_names) > 0L &&
        !anyNA(sample_names) &&
        all(nzchar(sample_names)) &&
        !anyDuplicated(sample_names)
    if (!valid_sample_names) {
        .abort_missingness_design(
            paste0(
                "Analysis sample names must be non-empty and unique before ",
                "design alignment."
            ),
            "imputefinder_design_alignment_error",
            missing_samples = character(),
            unexpected_samples = character()
        )
    }

    unname(sample_names)
}

.missingness_design_alignment_difference <- function(
    design_names,
    sample_names
) {
    list(
        missing_samples = sort(
            setdiff(design_names, sample_names),
            method = "radix"
        ),
        unexpected_samples = sort(
            setdiff(sample_names, design_names),
            method = "radix"
        )
    )
}

.align_missingness_design <- function(design, sample_names) {
    .validate_missingness_design(design)
    sample_names <- .validate_design_alignment_names(sample_names)

    design_names <- rownames(design$sample_data)
    difference <- .missingness_design_alignment_difference(
        design_names,
        sample_names
    )
    exact_set <- length(sample_names) == length(design_names) &&
        length(difference$missing_samples) == 0L &&
        length(difference$unexpected_samples) == 0L
    if (!exact_set) {
        .abort_missingness_design(
            paste0(
                "Design sample row names must match the analysis sample ",
                "names exactly."
            ),
            "imputefinder_design_alignment_error",
            missing_samples = difference$missing_samples,
            unexpected_samples = difference$unexpected_samples
        )
    }

    aligned <- design
    aligned$sample_data <- design$sample_data[
        match(sample_names, design_names),
        ,
        drop = FALSE
    ]
    aligned
}

#' @export
print.missingness_design <- function(x, ...) {
    .validate_missingness_design(x)
    roles <- x$roles
    condition_count <- length(unique(x$sample_data[[roles$condition]]))
    cat("<missingness_design>\n")
    cat(sprintf(
        "Lifecycle: %s (%s)\n",
        x$spec$lifecycle,
        x$spec$schema
    ))
    cat(sprintf(
        "Samples: %d across %d condition%s\n",
        nrow(x$sample_data),
        condition_count,
        if (condition_count == 1L) "" else "s"
    ))
    cat(sprintf("Condition: %s\n", roles$condition))
    cat(sprintf(
        "Nuisance: %s; block: %s; acquisition: %s\n",
        if (length(roles$nuisance)) {
            paste(roles$nuisance, collapse = ", ")
        } else {
            "unavailable"
        },
        if (length(roles$block)) roles$block else "unavailable",
        if (length(roles$acquisition)) roles$acquisition else "unavailable"
    ))
    interaction_text <- if (length(x$interactions)) {
        vapply(x$interactions, paste, collapse = ":", character(1))
    } else {
        "none"
    }
    cat(
        "Interactions: ",
        paste(interaction_text, collapse = ", "),
        "\n",
        sep = ""
    )

    invisible(x)
}
