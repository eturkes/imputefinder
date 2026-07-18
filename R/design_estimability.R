.DESIGN_ESTIMABILITY_SCHEMA <- "design_estimability_v1"
.DESIGN_CONTRAST_SCHEMA <- "design_contrast_estimability_v1"
.DESIGN_ESTIMABILITY_METHODS <- c(
    encoding = "canonical_treatment_contrasts_v1",
    rank = "svd_relative_v1",
    aliasing = "ordered_null_projector_mgs_v1",
    contrast = "row_space_residual_v1",
    units = "declared_block_or_sample_v1"
)
.DESIGN_ALGEBRA_TOLERANCE <- sqrt(.Machine$double.eps)

.abort_design_estimability <- function(message, subclass, ...) {
    .abort_sidecar(message, subclass, ...)
}

.abort_design_contrast <- function(message, ...) {
    .abort_missingness_design(
        message,
        "imputefinder_design_contrast_error",
        ...
    )
}

.canonical_design_text <- function(x) {
    output <- unname(as.character(x))
    text <- Encoding(output) != "bytes"
    output[text] <- enc2utf8(output[text])
    if (anyNA(output)) {
        .abort_design_estimability(
            "Declared design labels must convert losslessly to UTF-8 or bytes.",
            "imputefinder_design_estimability_error",
            field = "design"
        )
    }
    output
}

.design_variable_role_frame <- function(design) {
    roles <- design$roles
    columns <- unname(unlist(roles, use.names = FALSE))
    role <- rep(
        names(roles),
        lengths(roles)
    )
    data.frame(
        column = .canonical_design_text(columns),
        role = unname(role),
        stringsAsFactors = FALSE
    )
}

.design_canonical_data <- function(design) {
    sample <- .canonical_sidecar_names(rownames(design$sample_data))
    ordered <- order(sample, method = "radix")
    data <- design$sample_data[ordered, , drop = FALSE]
    rownames(data) <- sample[ordered]
    data
}

.new_design_variable <- function(values, column, role, sample_names) {
    numeric <- !is.factor(values) &&
        (is.integer(values) || is.double(values))
    if (numeric) {
        basis <- matrix(
            as.double(values),
            ncol = 1L,
            dimnames = list(sample_names, NULL)
        )
        labels <- column
        levels <- character()
        coefficient_level <- NA_character_
        encoding <- "numeric"
        reference <- NA_character_
    } else {
        values <- .canonical_design_text(values)
        levels <- sort(unique(values), method = "radix")
        encoded_levels <- levels[-1L]
        basis <- vapply(
            encoded_levels,
            function(level) as.double(values == level),
            numeric(length(values))
        )
        basis <- matrix(
            basis,
            nrow = length(values),
            ncol = length(encoded_levels),
            dimnames = list(sample_names, NULL)
        )
        labels <- if (length(encoded_levels)) {
            paste0(column, "[", encoded_levels, "]")
        } else {
            character()
        }
        coefficient_level <- encoded_levels
        encoding <- "treatment"
        reference <- levels[[1L]]
    }

    list(
        column = column,
        role = role,
        encoding = encoding,
        reference = reference,
        levels = levels,
        basis = basis,
        labels = labels,
        coefficient_level = coefficient_level
    )
}

.design_variables <- function(design, data) {
    role_frame <- .design_variable_role_frame(design)
    original_columns <- unname(unlist(design$roles, use.names = FALSE))
    variables <- lapply(seq_along(original_columns), function(index) {
        .new_design_variable(
            data[[original_columns[[index]]]],
            role_frame$column[[index]],
            role_frame$role[[index]],
            rownames(data)
        )
    })
    names(variables) <- role_frame$column
    variables
}

.design_variable_frame <- function(variables) {
    data.frame(
        column = vapply(variables, `[[`, character(1L), "column"),
        role = vapply(variables, `[[`, character(1L), "role"),
        encoding = vapply(variables, `[[`, character(1L), "encoding"),
        reference = vapply(
            variables,
            function(variable) variable$reference,
            character(1L)
        ),
        levels = I(lapply(variables, `[[`, "levels")),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

.new_design_term <- function(
    term,
    kind,
    role,
    components,
    basis,
    labels,
    encoding,
    level
) {
    list(
        term = term,
        kind = kind,
        role = role,
        components = components,
        basis = basis,
        labels = labels,
        encoding = encoding,
        level = level
    )
}

.design_main_term <- function(variable) {
    .new_design_term(
        term = variable$column,
        kind = "main",
        role = variable$role,
        components = variable$column,
        basis = variable$basis,
        labels = variable$labels,
        encoding = rep(variable$encoding, ncol(variable$basis)),
        level = variable$coefficient_level
    )
}

.design_interaction_grid <- function(variables) {
    dimensions <- vapply(
        variables,
        function(variable) ncol(variable$basis),
        integer(1L)
    )
    if (any(dimensions == 0L)) {
        return(matrix(integer(), nrow = 0L, ncol = length(variables)))
    }
    grid <- expand.grid(
        lapply(dimensions, seq_len),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    as.matrix(grid)
}

.design_interaction_term <- function(components, variables, sample_names) {
    selected <- variables[components]
    grid <- .design_interaction_grid(selected)
    if (nrow(grid) == 0L) {
        basis <- matrix(
            numeric(),
            nrow = length(sample_names),
            ncol = 0L,
            dimnames = list(sample_names, NULL)
        )
        labels <- level <- character()
    } else {
        basis <- vapply(seq_len(nrow(grid)), function(index) {
            product <- rep(1, length(sample_names))
            for (component_index in seq_along(selected)) {
                product <- product * selected[[component_index]]$basis[
                    , grid[index, component_index]
                ]
            }
            product
        }, numeric(length(sample_names)))
        basis <- matrix(
            basis,
            nrow = length(sample_names),
            ncol = nrow(grid),
            dimnames = list(sample_names, NULL)
        )
        labels <- vapply(seq_len(nrow(grid)), function(index) {
            paste(vapply(seq_along(selected), function(component_index) {
                selected[[component_index]]$labels[
                    grid[index, component_index]
                ]
            }, character(1L)), collapse = ":")
        }, character(1L))
        level <- labels
    }

    .new_design_term(
        term = paste(components, collapse = ":"),
        kind = "interaction",
        role = "interaction",
        components = components,
        basis = basis,
        labels = labels,
        encoding = rep("product", ncol(basis)),
        level = level
    )
}

.design_model_terms <- function(design, variables, sample_names) {
    intercept <- .new_design_term(
        term = "(Intercept)",
        kind = "intercept",
        role = "intercept",
        components = character(),
        basis = matrix(
            1,
            nrow = length(sample_names),
            ncol = 1L,
            dimnames = list(sample_names, NULL)
        ),
        labels = "(Intercept)",
        encoding = "intercept",
        level = NA_character_
    )
    main_columns <- c(
        design$roles$condition,
        design$roles$nuisance,
        design$roles$block
    )
    main_columns <- .canonical_design_text(main_columns)
    main <- lapply(variables[main_columns], .design_main_term)
    interactions <- lapply(design$interactions, function(components) {
        components <- .canonical_design_text(components)
        .design_interaction_term(components, variables, sample_names)
    })

    c(list(intercept), main, interactions)
}

.design_term_frame <- function(terms) {
    count <- vapply(terms, function(term) ncol(term$basis), integer(1L))
    data.frame(
        term_id = sprintf("term_%04d", seq_along(terms)),
        term = vapply(terms, `[[`, character(1L), "term"),
        kind = vapply(terms, `[[`, character(1L), "kind"),
        role = vapply(terms, `[[`, character(1L), "role"),
        coefficient_count = count,
        components = I(lapply(terms, `[[`, "components")),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

.design_column_frame <- function(terms, term_frame) {
    rows <- lapply(seq_along(terms), function(index) {
        term <- terms[[index]]
        count <- ncol(term$basis)
        if (count == 0L) {
            return(NULL)
        }
        data.frame(
            term_id = rep(term_frame$term_id[[index]], count),
            term = rep(term$term, count),
            kind = rep(term$kind, count),
            label = term$labels,
            encoding = term$encoding,
            level = term$level,
            stringsAsFactors = FALSE
        )
    })
    rows <- rows[!vapply(rows, is.null, logical(1L))]
    columns <- do.call(rbind, rows)
    row.names(columns) <- NULL
    columns <- cbind(
        coefficient = sprintf("coef_%04d", seq_len(nrow(columns))),
        columns,
        stringsAsFactors = FALSE
    )
    columns
}

.new_design_model <- function(design) {
    data <- .design_canonical_data(design)
    variables <- .design_variables(design, data)
    terms <- .design_model_terms(design, variables, rownames(data))
    term_frame <- .design_term_frame(terms)
    column_frame <- .design_column_frame(terms, term_frame)
    matrix <- do.call(cbind, lapply(terms, `[[`, "basis"))
    matrix <- matrix(
        as.double(matrix),
        nrow = nrow(data),
        ncol = nrow(column_frame),
        dimnames = list(rownames(data), column_frame$coefficient)
    )

    list(
        matrix = matrix,
        variables = .design_variable_frame(variables),
        terms = term_frame,
        columns = column_frame
    )
}

.design_null_basis <- function(projector, nullity) {
    coefficient_count <- nrow(projector)
    basis <- matrix(numeric(), nrow = coefficient_count, ncol = 0L)
    if (nullity > 0L) {
        for (axis in seq_len(coefficient_count)) {
            candidate <- projector[, axis]
            if (ncol(basis)) {
                for (column in seq_len(ncol(basis))) {
                    candidate <- candidate - basis[, column] *
                        sum(basis[, column] * candidate)
                }
            }
            candidate_norm <- sqrt(sum(candidate^2))
            if (candidate_norm <= .DESIGN_ALGEBRA_TOLERANCE) {
                next
            }
            candidate <- candidate / candidate_norm
            first <- which(
                abs(candidate) > .DESIGN_ALGEBRA_TOLERANCE
            )[[1L]]
            if (candidate[[first]] < 0) {
                candidate <- -candidate
            }
            candidate[
                abs(candidate) <= .Machine$double.eps * coefficient_count
            ] <- 0
            basis <- cbind(basis, candidate)
            if (ncol(basis) == nullity) {
                break
            }
        }
    }
    if (ncol(basis) != nullity) {
        .abort_design_estimability(
            paste0(
                "Canonical null-space construction did not recover its ",
                "SVD nullity."
            ),
            "imputefinder_design_estimability_error",
            field = "aliasing"
        )
    }
    rownames(basis) <- rownames(projector)
    colnames(basis) <- if (nullity) {
        sprintf("alias_%04d", seq_len(nullity))
    } else {
        character()
    }
    basis
}

.empty_affected_coefficients <- function() {
    data.frame(
        coefficient = character(),
        term_id = character(),
        term = character(),
        alias_norm = numeric(),
        stringsAsFactors = FALSE
    )
}

.empty_affected_terms <- function() {
    data.frame(
        term_id = character(),
        term = character(),
        alias_norm = numeric(),
        stringsAsFactors = FALSE
    )
}

.empty_contrast_affected_terms <- function() {
    data.frame(
        term_id = character(),
        term = character(),
        residual_magnitude = numeric(),
        stringsAsFactors = FALSE
    )
}

.design_affected_aliases <- function(alias_norm, model) {
    selected <- which(alias_norm > .DESIGN_ALGEBRA_TOLERANCE)
    if (!length(selected)) {
        return(list(
            coefficients = .empty_affected_coefficients(),
            terms = .empty_affected_terms()
        ))
    }
    columns <- model$columns[selected, , drop = FALSE]
    coefficients <- data.frame(
        coefficient = columns$coefficient,
        term_id = columns$term_id,
        term = columns$term,
        alias_norm = unname(alias_norm[selected]),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    term_selected <- model$terms$term_id %in% unique(columns$term_id)
    term_rows <- model$terms[term_selected, , drop = FALSE]
    term_norm <- vapply(term_rows$term_id, function(term_id) {
        max(coefficients$alias_norm[coefficients$term_id == term_id])
    }, numeric(1L))
    terms <- data.frame(
        term_id = term_rows$term_id,
        term = term_rows$term,
        alias_norm = unname(term_norm),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    list(coefficients = coefficients, terms = terms)
}

.design_svd <- function(x) {
    decomposition <- svd(
        x,
        nu = min(dim(x)),
        nv = min(dim(x))
    )
    threshold <- max(dim(x)) * .Machine$double.eps *
        decomposition$d[[1L]]
    rank <- as.integer(sum(decomposition$d > threshold))
    list(
        decomposition = decomposition,
        threshold = threshold,
        rank = rank,
        nullity = as.integer(ncol(x) - rank)
    )
}

.design_null_projector <- function(decomposition, rank, coefficients) {
    coefficient_count <- length(coefficients)
    row_space <- if (rank) {
        decomposition$v[, seq_len(rank), drop = FALSE]
    } else {
        matrix(numeric(), nrow = coefficient_count, ncol = 0L)
    }
    projector <- diag(coefficient_count) - tcrossprod(row_space)
    projector <- (projector + t(projector)) / 2
    projector[
        abs(projector) <= .Machine$double.eps * coefficient_count
    ] <- 0
    dimnames(projector) <- list(coefficients, coefficients)
    projector
}

.design_rank_record <- function(x, summary) {
    decomposition <- summary$decomposition
    leverage <- if (summary$rank) {
        rowSums(
            decomposition$u[, seq_len(summary$rank), drop = FALSE]^2
        )
    } else {
        rep(0, nrow(x))
    }
    names(leverage) <- rownames(x)
    singular_names <- sprintf("d_%04d", seq_along(decomposition$d))
    singular_values <- stats::setNames(decomposition$d, singular_names)
    scaled <- stats::setNames(
        decomposition$d / decomposition$d[[1L]],
        singular_names
    )
    full_column_rank <- summary$rank == ncol(x)
    condition_number <- if (full_column_rank) {
        decomposition$d[[1L]] / decomposition$d[[summary$rank]]
    } else {
        Inf
    }
    list(
        sample_count = as.integer(nrow(x)),
        coefficient_count = as.integer(ncol(x)),
        rank = summary$rank,
        nullity = summary$nullity,
        threshold = summary$threshold,
        full_column_rank = full_column_rank,
        singular_values = singular_values,
        scaled_singular_values = scaled,
        condition_number = condition_number,
        leverage = leverage
    )
}

.new_design_algebra <- function(model) {
    summary <- .design_svd(model$matrix)
    projector <- .design_null_projector(
        summary$decomposition,
        summary$rank,
        model$columns$coefficient
    )
    null_basis <- .design_null_basis(projector, summary$nullity)
    aliases <- .design_affected_aliases(
        sqrt(pmax(diag(projector), 0)),
        model
    )

    list(
        rank = .design_rank_record(model$matrix, summary),
        aliasing = list(
            null_projector = projector,
            null_basis = null_basis,
            affected_coefficients = aliases$coefficients,
            affected_terms = aliases$terms
        )
    )
}

.design_count_combinations <- function(data, keys, count_name) {
    key_data <- data[keys]
    ordered <- do.call(
        order,
        c(unname(key_data), list(method = "radix"))
    )
    unique_rows <- unique(key_data[ordered, , drop = FALSE])
    counts <- vapply(seq_len(nrow(unique_rows)), function(index) {
        selected <- rep(TRUE, nrow(data))
        for (key in keys) {
            selected <- selected & data[[key]] == unique_rows[[key]][[index]]
        }
        sum(selected)
    }, integer(1L))
    unique_rows[[count_name]] <- counts
    row.names(unique_rows) <- NULL
    unique_rows
}

.design_unit_samples <- function(design) {
    data <- .design_canonical_data(design)
    samples <- rownames(data)
    condition_column <- design$roles$condition
    condition <- .canonical_design_text(data[[condition_column]])
    has_block <- length(design$roles$block) == 1L
    unit <- if (has_block) {
        .canonical_design_text(data[[design$roles$block]])
    } else {
        samples
    }
    list(
        has_block = has_block,
        sample = data.frame(
            sample = samples,
            condition = condition,
            unit = unit,
            weight = rep(1L, length(samples)),
            stringsAsFactors = FALSE
        )
    )
}

.design_unit_tables <- function(sample) {
    unit_condition <- .design_count_combinations(
        sample,
        c("unit", "condition"),
        "sample_count"
    )
    units <- .design_count_combinations(sample, "unit", "sample_count")
    units$condition_count <- vapply(units$unit, function(value) {
        length(unique(sample$condition[sample$unit == value]))
    }, integer(1L))
    conditions <- .design_count_combinations(
        sample,
        "condition",
        "sample_count"
    )
    conditions$unit_count <- vapply(conditions$condition, function(value) {
        length(unique(sample$unit[sample$condition == value]))
    }, integer(1L))
    list(
        unit_condition = unit_condition,
        unit = units,
        condition = conditions
    )
}

.design_attach_unit_counts <- function(sample, tables) {
    unit_condition <- tables$unit_condition
    units <- tables$unit
    sample$unit_sample_count <- units$sample_count[
        match(sample$unit, units$unit)
    ]
    sample$unit_condition_sample_count <- vapply(
        seq_len(nrow(sample)),
        function(index) {
            row <- unit_condition$unit == sample$unit[[index]] &
                unit_condition$condition == sample$condition[[index]]
            unit_condition$sample_count[[which(row)]]
        },
        integer(1L)
    )
    sample
}

.new_design_units <- function(design, model) {
    context <- .design_unit_samples(design)
    tables <- .design_unit_tables(context$sample)
    sample <- .design_attach_unit_counts(context$sample, tables)
    if (!identical(sample$sample, rownames(model$matrix))) {
        .abort_design_estimability(
            "Canonical design rows and unit accounting do not align.",
            "imputefinder_design_estimability_error",
            field = "units"
        )
    }

    list(
        resampling_role = if (context$has_block) "block" else "sample",
        resampling_column = if (context$has_block) {
            .canonical_design_text(design$roles$block)
        } else {
            NA_character_
        },
        independent_unit_count = as.integer(nrow(tables$unit)),
        grouped_technical_siblings = any(
            tables$unit_condition$sample_count > 1L
        ),
        sample = sample,
        unit_condition = tables$unit_condition,
        unit = tables$unit,
        condition = tables$condition
    )
}

.new_design_estimability <- function(design) {
    .validate_missingness_design(design)
    model <- .new_design_model(design)
    algebra <- .new_design_algebra(model)

    list(
        schema = .DESIGN_ESTIMABILITY_SCHEMA,
        methods = .DESIGN_ESTIMABILITY_METHODS,
        model = model,
        rank = algebra$rank,
        aliasing = algebra$aliasing,
        units = .new_design_units(design, model)
    )
}

.design_estimability_schema_value <- function(estimability) {
    if (is.list(estimability) && "schema" %in% names(estimability)) {
        .sidecar_scalar_character(estimability$schema)
    } else {
        NA_character_
    }
}

.validate_design_estimability <- function(estimability, design) {
    actual_schema <- .design_estimability_schema_value(estimability)
    if (!identical(actual_schema, .DESIGN_ESTIMABILITY_SCHEMA)) {
        .abort_design_estimability(
            paste0(
                "Unsupported design-estimability schema; reconstruct the ",
                "analysis with the current `analyze_missingness()`."
            ),
            "imputefinder_design_estimability_lifecycle_error",
            expected_schema = .DESIGN_ESTIMABILITY_SCHEMA,
            actual_schema = actual_schema
        )
    }
    expected <- .new_design_estimability(design)
    valid <- isTRUE(all.equal(
        estimability,
        expected,
        tolerance = .DESIGN_ALGEBRA_TOLERANCE,
        check.attributes = TRUE
    ))
    if (!valid) {
        .abort_design_estimability(
            "Stored design-estimability evidence is malformed or inconsistent.",
            "imputefinder_design_estimability_error",
            field = "design.estimability"
        )
    }
    invisible(estimability)
}

.normalise_coefficient_contrast <- function(core, contrast) {
    valid <- is.numeric(contrast) &&
        is.null(dim(contrast)) &&
        length(contrast) > 0L &&
        !anyNA(contrast) &&
        all(is.finite(contrast)) &&
        !is.null(names(contrast)) &&
        !anyNA(names(contrast)) &&
        all(nzchar(names(contrast))) &&
        !anyDuplicated(names(contrast)) &&
        all(names(contrast) %in% core$model$columns$coefficient)
    if (!valid || all(contrast == 0)) {
        .abort_design_contrast(
            paste0(
                "A coefficient contrast must be finite, nonzero, uniquely ",
                "named, and use model coefficients."
            ),
            contrast_type = "coefficient"
        )
    }
    coefficient <- stats::setNames(
        rep(0, nrow(core$model$columns)),
        core$model$columns$coefficient
    )
    coefficient[names(contrast)] <- as.double(contrast)
    list(
        descriptor = list(
            type = "coefficient",
            weights = coefficient
        ),
        coefficient = coefficient
    )
}

.resolve_design_contrast_variable <- function(core, selector) {
    variables <- core$model$variables
    selected <- which(
        variables$column == selector |
            variables$role == selector
    )
    if (length(selected) != 1L) {
        .abort_design_contrast(
            paste0(
                "A declared contrast must select exactly one categorical ",
                "design variable."
            ),
            selector = selector
        )
    }
    variables[selected, , drop = FALSE]
}

.valid_declared_contrast_container <- function(contrast) {
    is.list(contrast) &&
        !is.data.frame(contrast) &&
        length(contrast) == 1L &&
        !is.null(names(contrast)) &&
        !is.na(names(contrast)) &&
        nzchar(names(contrast))
}

.normalise_declared_weights <- function(variable, weights, selector) {
    levels <- variable$levels[[1L]]
    valid <- identical(variable$encoding, "treatment") &&
        length(levels) >= 2L &&
        is.numeric(weights) &&
        is.null(dim(weights)) &&
        length(weights) >= 2L &&
        !anyNA(weights) &&
        all(is.finite(weights)) &&
        !is.null(names(weights)) &&
        !anyNA(names(weights)) &&
        all(nzchar(names(weights))) &&
        !anyDuplicated(names(weights)) &&
        all(names(weights) %in% levels)
    zero_sum <- valid && abs(sum(weights)) <=
        .DESIGN_ALGEBRA_TOLERANCE * max(1, sum(abs(weights)))
    if (!valid || !zero_sum || sum(weights != 0) < 2L) {
        .abort_design_contrast(
            paste0(
                "Declared categorical weights must be finite, named known ",
                "levels with at least two nonzero sides and sum to zero."
            ),
            selector = selector
        )
    }
    canonical_weights <- stats::setNames(rep(0, length(levels)), levels)
    canonical_weights[names(weights)] <- as.double(weights)
    canonical_weights
}

.declared_contrast_coefficient <- function(
    core,
    variable,
    weights,
    selector
) {
    coefficient <- stats::setNames(
        rep(0, nrow(core$model$columns)),
        core$model$columns$coefficient
    )
    main <- core$model$columns$kind == "main" &
        core$model$columns$term == variable$column
    if (any(main)) {
        coefficient[core$model$columns$coefficient[main]] <-
            weights[core$model$columns$level[main]]
    }
    if (all(coefficient == 0)) {
        .abort_design_contrast(
            "The declared contrast has no encoded coefficient direction.",
            selector = selector
        )
    }
    coefficient
}

.normalise_declared_contrast <- function(core, contrast) {
    if (!.valid_declared_contrast_container(contrast)) {
        .abort_design_contrast(
            "A declared contrast must be a one-entry named list.",
            contrast_type = "declared"
        )
    }
    selector <- names(contrast)[[1L]]
    variable <- .resolve_design_contrast_variable(core, selector)
    canonical_weights <- .normalise_declared_weights(
        variable,
        contrast[[1L]],
        selector
    )
    coefficient <- .declared_contrast_coefficient(
        core,
        variable,
        canonical_weights,
        selector
    )

    list(
        descriptor = list(
            type = "declared_level",
            role = variable$role,
            column = variable$column,
            weights = canonical_weights,
            context = "reference_levels"
        ),
        coefficient = coefficient
    )
}

.contrast_affected_terms <- function(residual, core) {
    affected <- .design_affected_aliases(abs(residual), core$model)$terms
    names(affected)[names(affected) == "alias_norm"] <- "residual_magnitude"
    affected
}

.design_contrast_estimability <- function(core, contrast) {
    actual_schema <- .design_estimability_schema_value(core)
    if (!identical(actual_schema, .DESIGN_ESTIMABILITY_SCHEMA)) {
        .abort_design_contrast(
            paste0(
                "Contrast estimability requires a current ",
                "design-estimability record."
            ),
            actual_schema = actual_schema
        )
    }
    normalised <- if (is.list(contrast)) {
        .normalise_declared_contrast(core, contrast)
    } else {
        .normalise_coefficient_contrast(core, contrast)
    }
    coefficient <- normalised$coefficient
    residual <- as.vector(core$aliasing$null_projector %*% coefficient)
    names(residual) <- names(coefficient)
    contrast_norm <- sqrt(sum(coefficient^2))
    residual_norm <- sqrt(sum(residual^2))
    tolerance <- .DESIGN_ALGEBRA_TOLERANCE * max(1, contrast_norm)
    estimable <- residual_norm <= tolerance

    list(
        schema = .DESIGN_CONTRAST_SCHEMA,
        contrast = normalised$descriptor,
        coefficient = coefficient,
        residual = residual,
        contrast_norm = contrast_norm,
        residual_norm = residual_norm,
        tolerance = tolerance,
        estimable = estimable,
        affected_terms = if (estimable) {
            .empty_contrast_affected_terms()
        } else {
            .contrast_affected_terms(residual, core)
        }
    )
}
