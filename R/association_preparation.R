.ASSOCIATION_PREPARATION_SCHEMA <- "association_preparation_v1"
.ASSOCIATION_RESPONSE_METHOD <- "a_sample_detection_fraction_v3"
.ASSOCIATION_PREPARATION_FIELDS <- c(
    "schema", "protocol", "input_sha256", "identity", "response",
    "hypotheses", "support", "strata"
)
.ASSOCIATION_IDENTITY_FIELDS <- c(
    "dimensions", "feature_names", "sample_names", "original_mask", "design"
)
.ASSOCIATION_RESPONSE_FIELDS <- c(
    "stratum", "acquisition", "sample", "globally_observable_count",
    "detected_count", "detection_fraction"
)
.ASSOCIATION_HYPOTHESIS_FIELDS <- c(
    "hypothesis", "stratum", "coefficient", "label", "term_id", "term",
    "kind", "components", "component_encodings", "role", "eligible",
    "estimable"
)
.ASSOCIATION_SUPPORT_FIELDS <- c(
    "hypothesis", "design", "side_definition", "side_count_low",
    "side_count_high", "complete_block_count", "cell_count_min",
    "numeric_components", "eligible", "code"
)
.ASSOCIATION_STRATUM_FIELDS <- c(
    "stratum", "acquisition", "samples", "design", "core", "response"
)

.association_protocol <- function() {
    c(
        id = .ASSOCIATION_PROTOCOL_ID,
        hash = .ASSOCIATION_PROTOCOL_HASH,
        contract_hash = .ASSOCIATION_CONTRACT_HASH,
        response = .ASSOCIATION_RESPONSE_METHOD,
        encoding = unname(.DESIGN_ESTIMABILITY_METHODS[["encoding"]]),
        rank = unname(.DESIGN_ESTIMABILITY_METHODS[["rank"]]),
        input_hash = .ASSOCIATION_INPUT_HASH
    )
}

.new_association_identity <- function(data, design) {
    .validate_matrix_data(data)
    design <- .align_missingness_design(design, colnames(data))
    list(
        dimensions = c(
            features = as.integer(nrow(data)),
            samples = as.integer(ncol(data))
        ),
        feature_names = .canonical_sidecar_names(rownames(data)),
        sample_names = .canonical_sidecar_names(colnames(data)),
        original_mask = .pack_original_mask(data),
        design = design
    )
}

.association_identity_material <- function(identity) {
    valid <- is.list(identity) &&
        identical(names(identity), .ASSOCIATION_IDENTITY_FIELDS) &&
        is.integer(identity$dimensions) &&
        identical(names(identity$dimensions), c("features", "samples")) &&
        is.character(identity$feature_names) &&
        is.character(identity$sample_names) &&
        inherits(identity$design, "missingness_design")
    if (!valid) {
        .abort_association(
            "Stored association input identity is malformed.",
            "imputefinder_association_preparation_error"
        )
    }
    dimensions <- .validate_mask_dimensions(identity$dimensions)
    .validate_axis_names(
        identity$feature_names,
        dimensions[[1L]],
        "feature",
        "Stored association identity"
    )
    .validate_axis_names(
        identity$sample_names,
        dimensions[[2L]],
        "sample",
        "Stored association identity"
    )
    canonical <- identical(
        identity$feature_names,
        .canonical_sidecar_names(identity$feature_names)
    ) && identical(
        identity$sample_names,
        .canonical_sidecar_names(identity$sample_names)
    )
    .validate_missingness_design(identity$design)
    aligned <- .align_missingness_design(
        identity$design,
        identity$sample_names
    )
    canonical <- canonical && identical(identity$design, aligned)
    if (!canonical) {
        .abort_association(
            "Stored association input identity is not canonical or aligned.",
            "imputefinder_association_preparation_error"
        )
    }
    missing_mask <- .unpack_original_mask(
        identity$original_mask,
        dimensions,
        identity$feature_names,
        identity$sample_names
    )
    data <- matrix(
        0,
        nrow = dimensions[[1L]],
        ncol = dimensions[[2L]],
        dimnames = list(identity$feature_names, identity$sample_names)
    )
    data[missing_mask] <- NA_real_
    list(
        missing_mask = missing_mask,
        input_sha256 = .association_mask_design_sha256(
            missing_mask,
            identity$design$sample_data,
            identity$design$roles,
            identity$design$interactions
        ),
        context = .sentinel_coverage_context(data, identity$design)
    )
}

.validate_association_identity <- function(identity) {
    .association_identity_material(identity)
    invisible(identity)
}

.empty_association_hypotheses <- function() {
    output <- data.frame(
        hypothesis = character(),
        stratum = character(),
        coefficient = character(),
        label = character(),
        term_id = character(),
        term = character(),
        kind = character(),
        role = character(),
        eligible = logical(),
        estimable = logical(),
        stringsAsFactors = FALSE
    )
    output$components <- I(vector("list", 0L))
    output$component_encodings <- I(vector("list", 0L))
    output[.ASSOCIATION_HYPOTHESIS_FIELDS]
}

.empty_association_numeric_support <- function() {
    data.frame(
        component = character(),
        median = double(),
        minimum = double(),
        maximum = double(),
        extrapolates_0_1 = logical(),
        stringsAsFactors = FALSE
    )
}

.empty_association_support <- function() {
    output <- data.frame(
        hypothesis = character(),
        design = character(),
        side_definition = character(),
        side_count_low = integer(),
        side_count_high = integer(),
        complete_block_count = integer(),
        cell_count_min = integer(),
        eligible = logical(),
        code = character(),
        stringsAsFactors = FALSE
    )
    output$numeric_components <- I(vector("list", 0L))
    output[.ASSOCIATION_SUPPORT_FIELDS]
}

.association_bind_rows <- function(rows, empty) {
    rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
    if (!length(rows)) {
        return(empty())
    }
    output <- do.call(rbind, unname(rows))
    row.names(output) <- NULL
    output
}

.association_response <- function(context) {
    globally_observable_count <- as.integer(sum(context$globally_observable))
    if (globally_observable_count == 0L) {
        return(NULL)
    }
    detected_count <- as.integer(colSums(
        context$observed[context$globally_observable, , drop = FALSE]
    ))
    acquisition_column <- context$design$roles$acquisition
    if (length(acquisition_column)) {
        acquisition <- .canonical_design_text(
            context$design$sample_data[[acquisition_column]]
        )
        levels <- sort(unique(acquisition), method = "radix")
    } else {
        acquisition <- rep("undeclared", length(context$sample))
        levels <- "undeclared"
    }

    rows <- lapply(levels, function(level) {
        selected <- acquisition == level
        stratum <- if (length(acquisition_column)) {
            .association_stratum_id(level)
        } else {
            "all"
        }
        data.frame(
            stratum = rep(stratum, sum(selected)),
            acquisition = rep(level, sum(selected)),
            sample = context$sample[selected],
            globally_observable_count = rep(
                globally_observable_count,
                sum(selected)
            ),
            detected_count = detected_count[selected],
            detection_fraction = detected_count[selected] /
                globally_observable_count,
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    })
    output <- do.call(rbind, unname(rows))
    row.names(output) <- NULL
    output[.ASSOCIATION_RESPONSE_FIELDS]
}

.association_stratum_design <- function(design, acquisition) {
    acquisition_column <- design$roles$acquisition
    if (!length(acquisition_column)) {
        return(design)
    }
    values <- .canonical_design_text(
        design$sample_data[[acquisition_column]]
    )
    sample_data <- design$sample_data[values == acquisition, , drop = FALSE]
    interactions <- Filter(
        function(components) !acquisition_column %in% components,
        design$interactions
    )
    missingness_design(
        sample_data,
        condition = design$roles$condition,
        nuisance = design$roles$nuisance,
        block = design$roles$block,
        interactions = interactions
    )
}

.association_component_encodings <- function(components, core) {
    rows <- match(components, core$model$variables$column)
    if (anyNA(rows)) {
        .abort_association(
            "Association hypothesis components do not resolve to variables.",
            "imputefinder_association_preparation_error"
        )
    }
    stats::setNames(core$model$variables$encoding[rows], components)
}

.association_eligible_column <- function(column_index, core) {
    column <- core$model$columns[column_index, , drop = FALSE]
    term_index <- match(column$term_id, core$model$terms$term_id)
    term <- core$model$terms[term_index, , drop = FALSE]
    if (identical(term$kind, "main")) {
        return(term$role %in% c("condition", "nuisance"))
    }
    if (!identical(term$kind, "interaction")) {
        return(FALSE)
    }
    eligible_columns <- c(
        core$model$variables$column[
            core$model$variables$role == "condition"
        ],
        core$model$variables$column[
            core$model$variables$role == "nuisance"
        ]
    )
    all(term$components[[1L]] %in% eligible_columns)
}

.association_hypotheses <- function(core, stratum) {
    selected <- which(vapply(
        seq_len(nrow(core$model$columns)),
        .association_eligible_column,
        logical(1L),
        core = core
    ))
    if (!length(selected)) {
        return(.empty_association_hypotheses())
    }
    columns <- core$model$columns[selected, , drop = FALSE]
    term_index <- match(columns$term_id, core$model$terms$term_id)
    terms <- core$model$terms[term_index, , drop = FALSE]
    components <- unname(terms$components)
    component_encodings <- lapply(
        components,
        .association_component_encodings,
        core = core
    )
    estimable <- vapply(columns$coefficient, function(coefficient) {
        contrast <- stats::setNames(1, coefficient)
        .design_contrast_estimability(core, contrast)$estimable
    }, logical(1L))
    hypothesis <- vapply(seq_len(nrow(columns)), function(index) {
        .association_hypothesis_id(
            stratum,
            columns$coefficient[[index]],
            columns$label[[index]],
            columns$term_id[[index]]
        )
    }, character(1L))
    output <- data.frame(
        hypothesis = hypothesis,
        stratum = rep(stratum, nrow(columns)),
        coefficient = columns$coefficient,
        label = columns$label,
        term_id = columns$term_id,
        term = columns$term,
        kind = columns$kind,
        role = ifelse(
            columns$kind == "interaction",
            "interaction",
            terms$role
        ),
        eligible = rep(TRUE, nrow(columns)),
        estimable = estimable,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    output$components <- I(components)
    output$component_encodings <- I(component_encodings)
    output[.ASSOCIATION_HYPOTHESIS_FIELDS]
}

.association_component_options <- function(variable) {
    if (identical(variable$encoding, "numeric")) {
        return(NA_character_)
    }
    variable$levels[[1L]][-1L]
}

.association_hypothesis_components <- function(core, hypothesis) {
    components <- hypothesis$components[[1L]]
    variable_rows <- match(components, core$model$variables$column)
    variables <- lapply(variable_rows, function(index) {
        core$model$variables[index, , drop = FALSE]
    })
    options <- lapply(variables, .association_component_options)
    if (identical(hypothesis$kind, "main")) {
        column_index <- match(
            hypothesis$coefficient,
            core$model$columns$coefficient
        )
        selected <- if (identical(variables[[1L]]$encoding, "numeric")) {
            1L
        } else {
            match(core$model$columns$level[[column_index]], options[[1L]])
        }
    } else {
        term_columns <- core$model$columns[
            core$model$columns$term_id == hypothesis$term_id,
            ,
            drop = FALSE
        ]
        column_index <- match(hypothesis$coefficient, term_columns$coefficient)
        dimensions <- lengths(options)
        grid <- expand.grid(
            lapply(dimensions, seq_len),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
        selected <- as.integer(grid[column_index, , drop = TRUE])
    }
    output <- lapply(seq_along(components), function(index) {
        variable <- variables[[index]]
        list(
            component = components[[index]],
            encoding = variable$encoding,
            reference = if (identical(variable$encoding, "treatment")) {
                variable$reference
            } else {
                NA_character_
            },
            target = options[[index]][selected[[index]]]
        )
    })
    names(output) <- components
    output
}

.association_position_values <- function(design) {
    data <- .design_canonical_data(design)
    components <- c(design$roles$condition, design$roles$nuisance)
    values <- lapply(components, function(component) {
        value <- data[[component]]
        if (!is.factor(value) && (is.integer(value) || is.double(value))) {
            value <- as.double(value)
            value[value == 0] <- 0
            value
        } else {
            .canonical_design_text(value)
        }
    })
    names(values) <- components
    blocked <- length(design$roles$block) == 1L
    unit <- if (blocked) {
        .canonical_design_text(data[[design$roles$block]])
    } else {
        rownames(data)
    }
    position <- data.frame(
        .sample = rownames(data),
        .unit = unit,
        stringsAsFactors = FALSE
    )
    value_names <- sprintf(".value_%04d", seq_along(values))
    for (index in seq_along(values)) {
        position[[value_names[[index]]]] <- values[[index]]
    }
    if (blocked) {
        keep <- !duplicated(position[c(".unit", value_names)])
        position <- position[keep, , drop = FALSE]
    }
    position_values <- lapply(value_names, function(name) position[[name]])
    names(position_values) <- components
    list(
        blocked = blocked,
        sample = position$.sample,
        unit = position$.unit,
        values = position_values
    )
}

.association_numeric_support <- function(specifications, positions) {
    numeric <- vapply(
        specifications,
        function(specification) identical(specification$encoding, "numeric"),
        logical(1L)
    )
    components <- names(specifications)[numeric]
    if (!length(components)) {
        return(.empty_association_numeric_support())
    }
    minimum <- vapply(components, function(component) {
        min(positions$values[[component]])
    }, numeric(1L))
    maximum <- vapply(components, function(component) {
        max(positions$values[[component]])
    }, numeric(1L))
    median <- vapply(components, function(component) {
        stats::median(positions$values[[component]])
    }, numeric(1L))
    data.frame(
        component = components,
        median = unname(median),
        minimum = unname(minimum),
        maximum = unname(maximum),
        extrapolates_0_1 = unname(minimum > 0 | maximum < 1),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

.association_component_sides <- function(
    specifications,
    positions,
    numeric_support
) {
    output <- lapply(specifications, function(specification) {
        values <- positions$values[[specification$component]]
        if (identical(specification$encoding, "treatment")) {
            side <- rep(NA_integer_, length(values))
            side[values == specification$reference] <- 0L
            side[values == specification$target] <- 1L
            description <- paste0(
                specification$component, "[", specification$reference,
                "]/", specification$component, "[", specification$target,
                "]"
            )
        } else {
            row <- numeric_support$component == specification$component
            median <- numeric_support$median[[which(row)]]
            side <- rep(NA_integer_, length(values))
            side[values < median] <- 0L
            side[values > median] <- 1L
            value <- sprintf("%.17g", median)
            description <- paste0(
                specification$component, "<", value, "/",
                specification$component, ">", value
            )
        }
        list(side = side, description = description)
    })
    names(output) <- names(specifications)
    output
}

.association_unit_side_count <- function(side, value, unit) {
    as.integer(length(unique(unit[!is.na(side) & side == value])))
}

.association_complete_units <- function(sides, unit) {
    low <- unique(unit[!is.na(sides) & sides == 0L])
    high <- unique(unit[!is.na(sides) & sides == 1L])
    as.integer(length(intersect(low, high)))
}

.association_interaction_counts <- function(component_sides, unit) {
    sides <- do.call(cbind, lapply(component_sides, `[[`, "side"))
    sides <- matrix(
        as.integer(sides),
        nrow = length(unit),
        ncol = length(component_sides)
    )
    low <- apply(sides, 2L, .association_unit_side_count, value = 0L, unit = unit)
    high <- apply(sides, 2L, .association_unit_side_count, value = 1L, unit = unit)
    valid <- rowSums(is.na(sides)) == 0L
    cell <- if (any(valid)) {
        apply(sides[valid, , drop = FALSE], 1L, paste0, collapse = "")
    } else {
        character()
    }
    cell_unit <- data.frame(
        unit = unit[valid],
        cell = cell,
        stringsAsFactors = FALSE
    )
    cell_unit <- unique(cell_unit)
    cell_total <- 2^ncol(sides)
    cell_count_min <- if (length(unique(cell_unit$cell)) < cell_total) {
        0L
    } else {
        as.integer(min(table(cell_unit$cell)))
    }
    complete <- if (!nrow(cell_unit)) {
        0L
    } else {
        occupied <- vapply(
            split(cell_unit$cell, cell_unit$unit),
            function(value) length(unique(value)),
            integer(1L)
        )
        as.integer(sum(occupied == cell_total))
    }
    list(
        low = as.integer(min(low)),
        high = as.integer(min(high)),
        cell = cell_count_min,
        complete = complete
    )
}

.association_support_row <- function(design, core, hypothesis, positions) {
    specifications <- .association_hypothesis_components(core, hypothesis)
    numeric_support <- .association_numeric_support(specifications, positions)
    component_sides <- .association_component_sides(
        specifications,
        positions,
        numeric_support
    )
    side_definition <- paste(
        vapply(component_sides, `[[`, character(1L), "description"),
        collapse = " x "
    )
    design_type <- if (positions$blocked) "blocked" else "independent"

    if (identical(hypothesis$kind, "interaction")) {
        counts <- .association_interaction_counts(
            component_sides,
            positions$unit
        )
        complete <- if (positions$blocked) counts$complete else NA_integer_
        eligible <- if (positions$blocked) {
            counts$cell >= 6L && counts$complete >= 6L
        } else {
            counts$cell >= 4L
        }
        code <- if (eligible) NA_character_ else
            "association_low_interaction_support"
        cell_count_min <- counts$cell
    } else {
        side <- component_sides[[1L]]$side
        counts <- list(
            low = .association_unit_side_count(side, 0L, positions$unit),
            high = .association_unit_side_count(side, 1L, positions$unit)
        )
        complete <- if (positions$blocked) {
            .association_complete_units(side, positions$unit)
        } else {
            NA_integer_
        }
        eligible <- if (positions$blocked) {
            counts$low >= 6L && counts$high >= 6L && complete >= 6L
        } else {
            counts$low >= 4L && counts$high >= 4L
        }
        code <- if (eligible) {
            NA_character_
        } else if (positions$blocked) {
            "association_low_block_support"
        } else {
            "association_low_independent_support"
        }
        cell_count_min <- NA_integer_
    }

    output <- data.frame(
        hypothesis = hypothesis$hypothesis,
        design = design_type,
        side_definition = side_definition,
        side_count_low = counts$low,
        side_count_high = counts$high,
        complete_block_count = complete,
        cell_count_min = cell_count_min,
        eligible = eligible,
        code = code,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    output$numeric_components <- I(list(numeric_support))
    output[.ASSOCIATION_SUPPORT_FIELDS]
}

.association_support <- function(design, core, hypotheses) {
    if (!nrow(hypotheses)) {
        return(.empty_association_support())
    }
    positions <- .association_position_values(design)
    rows <- lapply(seq_len(nrow(hypotheses)), function(index) {
        .association_support_row(
            design,
            core,
            hypotheses[index, , drop = FALSE],
            positions
        )
    })
    output <- do.call(rbind, rows)
    row.names(output) <- NULL
    output
}

.association_strata <- function(context, response) {
    acquisitions <- unique(response$acquisition)
    rows <- lapply(acquisitions, function(acquisition) {
        stratum <- unique(response$stratum[response$acquisition == acquisition])
        design <- .association_stratum_design(context$design, acquisition)
        core <- .new_design_estimability(design)
        samples <- rownames(core$model$matrix)
        selected <- match(samples, response$sample)
        values <- response$detection_fraction[selected]
        names(values) <- samples
        list(
            stratum = stratum,
            acquisition = acquisition,
            samples = samples,
            design = design,
            core = core,
            response = values
        )
    })
    names(rows) <- vapply(rows, `[[`, character(1L), "stratum")
    rows
}

.new_association_preparation <- function(data, design) {
    identity <- .new_association_identity(data, design)
    context <- .sentinel_coverage_context(data, identity$design)
    response <- .association_response(context)
    if (is.null(response)) {
        return(.new_unavailable(
            quantity = "association",
            code = "association_no_observable_features",
            message = "No feature is observed in the supplied input.",
            requires = "one globally observable feature"
        ))
    }
    strata <- .association_strata(context, response)
    hypotheses_by_stratum <- lapply(strata, function(stratum) {
        .association_hypotheses(stratum$core, stratum$stratum)
    })
    hypotheses <- .association_bind_rows(
        hypotheses_by_stratum,
        .empty_association_hypotheses
    )
    if (!nrow(hypotheses)) {
        return(.new_unavailable(
            quantity = "association",
            code = "association_no_testable_hypotheses",
            message = "The rebuilt design has no eligible association coefficient.",
            requires = "one eligible condition, nuisance, or interaction coefficient"
        ))
    }
    support_by_stratum <- lapply(seq_along(strata), function(index) {
        .association_support(
            strata[[index]]$design,
            strata[[index]]$core,
            hypotheses_by_stratum[[index]]
        )
    })
    support <- .association_bind_rows(
        support_by_stratum,
        .empty_association_support
    )
    preparation <- structure(
        list(
            schema = .ASSOCIATION_PREPARATION_SCHEMA,
            protocol = .association_protocol(),
            input_sha256 = .association_input_sha256(data, identity$design),
            identity = identity,
            response = response,
            hypotheses = hypotheses,
            support = support,
            strata = strata
        ),
        class = "imputefinder_association_preparation"
    )
    .validate_association_preparation(preparation)
    preparation
}

.valid_association_response <- function(response) {
    valid <- is.data.frame(response) && nrow(response) > 0L &&
        identical(names(response), .ASSOCIATION_RESPONSE_FIELDS) &&
        all(vapply(
            response[c("stratum", "acquisition", "sample")],
            is.character,
            logical(1L)
        )) &&
        is.integer(response$globally_observable_count) &&
        is.integer(response$detected_count) &&
        is.double(response$detection_fraction) && !anyNA(response) &&
        all(nzchar(response$stratum)) && all(nzchar(response$acquisition)) &&
        all(nzchar(response$sample)) &&
        !anyDuplicated(response$sample) &&
        all(response$stratum == "all" |
            grepl("^s_[0-9a-f]{64}$", response$stratum)) &&
        all(response$globally_observable_count > 0L) &&
        length(unique(response$globally_observable_count)) == 1L &&
        all(response$detected_count >= 0L) &&
        all(response$detected_count <= response$globally_observable_count) &&
        identical(
            response$detection_fraction,
            response$detected_count / response$globally_observable_count
        ) &&
        identical(
            do.call(
                order,
                list(response$acquisition, response$sample, method = "radix")
            ),
            seq_len(nrow(response))
        )
    if (!valid) {
        return(FALSE)
    }
    map <- response[
        !duplicated(response$stratum),
        c("stratum", "acquisition"),
        drop = FALSE
    ]
    one_acquisition <- vapply(
        split(response$acquisition, response$stratum),
        function(value) length(unique(value)) == 1L,
        logical(1L)
    )
    undeclared <- identical(map$stratum, "all") &&
        identical(map$acquisition, "undeclared")
    declared <- !any(map$stratum == "all") && identical(
        map$stratum,
        unname(vapply(
            map$acquisition,
            .association_stratum_id,
            character(1L)
        ))
    )
    all(one_acquisition) && (undeclared || declared)
}

.validate_association_stratum <- function(stratum, response) {
    valid <- is.list(stratum) &&
        identical(names(stratum), .ASSOCIATION_STRATUM_FIELDS) &&
        is.character(stratum$stratum) && length(stratum$stratum) == 1L &&
        is.character(stratum$acquisition) &&
        length(stratum$acquisition) == 1L &&
        is.character(stratum$samples) && !anyDuplicated(stratum$samples) &&
        inherits(stratum$design, "missingness_design") &&
        is.double(stratum$response) &&
        identical(names(stratum$response), stratum$samples) &&
        all(is.finite(stratum$response)) &&
        all(stratum$response >= 0 & stratum$response <= 1)
    if (!valid) {
        .abort_association(
            "Stored association stratum is malformed.",
            "imputefinder_association_preparation_error"
        )
    }
    .validate_missingness_design(stratum$design)
    .validate_design_estimability(stratum$core, stratum$design)
    selected <- match(stratum$samples, response$sample)
    expected_samples <- response$sample[response$stratum == stratum$stratum]
    valid <- !anyNA(selected) &&
        identical(stratum$samples, expected_samples) &&
        all(response$stratum[selected] == stratum$stratum) &&
        all(response$acquisition[selected] == stratum$acquisition) &&
        identical(
            response$detection_fraction[selected],
            unname(stratum$response)
        ) &&
        identical(rownames(stratum$core$model$matrix), stratum$samples) &&
        length(stratum$design$roles$acquisition) == 0L
    if (!valid) {
        .abort_association(
            "Stored association stratum does not align with its response.",
            "imputefinder_association_preparation_error"
        )
    }
    invisible(stratum)
}

.validate_association_preparation <- function(preparation) {
    valid <- is.list(preparation) &&
        identical(class(preparation), "imputefinder_association_preparation") &&
        identical(names(preparation), .ASSOCIATION_PREPARATION_FIELDS) &&
        identical(preparation$schema, .ASSOCIATION_PREPARATION_SCHEMA) &&
        identical(preparation$protocol, .association_protocol()) &&
        is.character(preparation$input_sha256) &&
        length(preparation$input_sha256) == 1L &&
        grepl("^[0-9a-f]{64}$", preparation$input_sha256) &&
        .valid_association_response(preparation$response) &&
        is.list(preparation$strata) && length(preparation$strata) > 0L &&
        identical(names(preparation$strata), unique(preparation$response$stratum))
    if (!valid) {
        .abort_association(
            "Stored association preparation has an invalid header.",
            "imputefinder_association_preparation_error"
        )
    }
    identity <- .association_identity_material(preparation$identity)
    expected_response <- .association_response(identity$context)
    valid <- identical(preparation$input_sha256, identity$input_sha256) &&
        identical(preparation$response, expected_response)
    if (!valid) {
        .abort_association(
            "Stored association input identity, hash, and response are detached.",
            "imputefinder_association_preparation_error"
        )
    }
    expected_strata <- .association_strata(identity$context, expected_response)
    valid <- identical(names(preparation$strata), names(expected_strata)) &&
        all(vapply(seq_along(preparation$strata), function(index) {
            stored <- preparation$strata[[index]]
            expected <- expected_strata[[index]]
            identical(names(preparation$strata)[[index]], stored$stratum) &&
                identical(stored$stratum, expected$stratum) &&
                identical(stored$acquisition, expected$acquisition) &&
                identical(stored$samples, expected$samples) &&
                identical(stored$design, expected$design) &&
                identical(stored$response, expected$response)
        }, logical(1L)))
    if (!valid) {
        .abort_association(
            "Stored association strata are detached from the global input.",
            "imputefinder_association_preparation_error"
        )
    }
    invisible(lapply(
        preparation$strata,
        .validate_association_stratum,
        response = preparation$response
    ))
    expected_hypotheses <- .association_bind_rows(
        lapply(preparation$strata, function(stratum) {
            .association_hypotheses(stratum$core, stratum$stratum)
        }),
        .empty_association_hypotheses
    )
    if (!nrow(expected_hypotheses)) {
        .abort_association(
            "An available preparation must contain an eligible hypothesis.",
            "imputefinder_association_preparation_error"
        )
    }
    expected_support <- .association_bind_rows(
        lapply(seq_along(preparation$strata), function(index) {
            stratum <- preparation$strata[[index]]
            selected <- expected_hypotheses$stratum == stratum$stratum
            .association_support(
                stratum$design,
                stratum$core,
                expected_hypotheses[selected, , drop = FALSE]
            )
        }),
        .empty_association_support
    )
    valid <- identical(preparation$hypotheses, expected_hypotheses) &&
        identical(preparation$support, expected_support)
    if (!valid) {
        .abort_association(
            "Stored association hypotheses or support are inconsistent.",
            "imputefinder_association_preparation_error"
        )
    }
    invisible(preparation)
}
