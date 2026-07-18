.SIDECAR_ANALYSIS_FIELDS <- c(
    "classic",
    "spec",
    "design",
    "input",
    "sentinel",
    "stability",
    "provenance"
)
.SIDECAR_DESIGN_FIELDS <- c(
    "declared",
    "estimability",
    "unavailable_roles"
)
.SIDECAR_PROVENANCE_FIELDS <- c(
    "call",
    "seeds",
    "hashes",
    "warnings",
    "failures",
    "assumptions",
    "training_scope"
)
.CLASSIC_RESULT_FIELDS <- c(
    "data",
    "classifications",
    "groups",
    "feature_status",
    "cutoffs",
    "cutoff_diagnostics",
    "profiles",
    "seed_log",
    "groups_by_sample",
    "call"
)
.CLASSIC_FAILURE_FIELDS <- c("stage", "class", "message", "call")

.sidecar_spec_value <- function(spec, field) {
    if (is.list(spec) && field %in% names(spec)) {
        .sidecar_scalar_character(spec[[field]])
    } else {
        NA_character_
    }
}

.validate_sidecar_spec_lifecycle <- function(spec) {
    actual_schema <- .sidecar_spec_value(spec, "schema")
    actual_lifecycle <- .sidecar_spec_value(spec, "lifecycle")
    if (!identical(actual_schema, .SIDECAR_SCHEMA_IDENTIFIER) ||
        !identical(actual_lifecycle, .SIDECAR_LIFECYCLE)) {
        .abort_sidecar(
            paste0(
                "Unsupported analysis schema or lifecycle; reconstruct it ",
                "with the current `analyze_missingness()`."
            ),
            "imputefinder_analysis_lifecycle_error",
            expected_schema = .SIDECAR_SCHEMA_IDENTIFIER,
            actual_schema = actual_schema,
            expected_lifecycle = .SIDECAR_LIFECYCLE,
            actual_lifecycle = actual_lifecycle
        )
    }

    invisible(spec)
}

.valid_sidecar_software <- function(software) {
    is.character(software) &&
        identical(names(software), c("package", "version")) &&
        length(software) == 2L &&
        identical(unname(software[["package"]]), "imputefinder") &&
        !is.na(software[["version"]]) &&
        nzchar(software[["version"]])
}

.validate_sidecar_spec <- function(spec) {
    .validate_sidecar_spec_lifecycle(spec)
    valid <- is.list(spec) &&
        identical(
            names(spec),
            c("schema", "lifecycle", "modules", "scientific", "software")
        ) &&
        identical(spec$modules, .SIDECAR_MODULE_IDENTIFIERS) &&
        identical(spec$scientific, .stable_scientific_identifiers()) &&
        .valid_sidecar_software(spec$software)
    if (!valid) {
        .abort_sidecar(
            "Stored analysis specification is malformed.",
            "imputefinder_analysis_schema_error",
            field = "spec"
        )
    }

    invisible(spec)
}

.design_unavailable_roles <- function(design) {
    optional <- .MISSINGNESS_DESIGN_ROLES[
        .MISSINGNESS_DESIGN_ROLES != "condition"
    ]
    optional[vapply(
        design$roles[optional],
        length,
        integer(1L)
    ) == 0L]
}

.new_sidecar_design_record <- function(design, input) {
    .validate_sidecar_input(input)
    aligned <- .align_missingness_design(design, input$sample_names)

    list(
        declared = aligned,
        estimability = NULL,
        unavailable_roles = .design_unavailable_roles(aligned)
    )
}

.validate_sidecar_design_shape <- function(record) {
    valid <- is.list(record) &&
        identical(names(record), .SIDECAR_DESIGN_FIELDS) &&
        (is.null(record$estimability) || is.list(record$estimability)) &&
        is.character(record$unavailable_roles) &&
        !anyNA(record$unavailable_roles) &&
        !anyDuplicated(record$unavailable_roles)
    if (!valid) {
        .abort_sidecar(
            "Stored design does not satisfy the analysis design schema.",
            "imputefinder_analysis_schema_error",
            field = "design"
        )
    }

    invisible(record)
}

.expected_sidecar_acquisition <- function(design) {
    role <- design$roles$acquisition
    if (length(role) == 0L) {
        return(NULL)
    }
    stats::setNames(
        design$sample_data[[role]],
        rownames(design$sample_data)
    )
}

.validate_sidecar_design <- function(record, input) {
    .validate_sidecar_design_shape(record)
    .validate_missingness_design(record$declared)
    aligned <- .align_missingness_design(
        record$declared,
        input$sample_names
    )
    aligned_names <- .canonical_sidecar_names(
        rownames(record$declared$sample_data)
    )
    valid <- identical(aligned, record$declared) &&
        identical(aligned_names, input$sample_names) &&
        identical(
            record$unavailable_roles,
            .design_unavailable_roles(record$declared)
        ) &&
        identical(
            input$acquisition,
            .expected_sidecar_acquisition(record$declared)
        )
    if (!valid) {
        .abort_sidecar(
            paste0(
                "Stored design must be aligned to input samples and ",
                "declarations."
            ),
            "imputefinder_analysis_schema_error",
            field = "design"
        )
    }

    invisible(record)
}

.new_classic_failure <- function(
    condition,
    stage,
    call = conditionCall(condition)
) {
    valid_stage <- is.character(stage) &&
        length(stage) == 1L &&
        !is.na(stage) &&
        nzchar(stage)
    if (!inherits(condition, "condition") || !valid_stage ||
        !(is.null(call) || is.language(call))) {
        .abort_sidecar(
            "Classic failures require a condition, stage, and portable call.",
            "imputefinder_analysis_schema_error",
            field = "classic"
        )
    }

    structure(
        list(
            stage = unname(stage),
            class = unname(class(condition)),
            message = conditionMessage(condition),
            call = call
        ),
        class = "imputefinder_classic_failure"
    )
}

.validate_classic_failure <- function(failure) {
    valid <- is.list(failure) &&
        identical(class(failure), "imputefinder_classic_failure") &&
        identical(names(failure), .CLASSIC_FAILURE_FIELDS) &&
        is.character(failure$stage) &&
        length(failure$stage) == 1L &&
        !is.na(failure$stage) &&
        nzchar(failure$stage) &&
        is.character(failure$class) &&
        length(failure$class) > 0L &&
        !anyNA(failure$class) &&
        all(nzchar(failure$class)) &&
        is.character(failure$message) &&
        length(failure$message) == 1L &&
        !is.na(failure$message) &&
        (is.null(failure$call) || is.language(failure$call))
    if (!valid) {
        .abort_sidecar(
            "Stored classic failure is malformed.",
            "imputefinder_analysis_schema_error",
            field = "classic"
        )
    }

    invisible(failure)
}

.valid_classic_data <- function(data) {
    (is.matrix(data) && is.numeric(data)) ||
        methods::is(data, "SummarizedExperiment")
}

.validate_classic_result_shape <- function(result) {
    valid <- is.list(result) &&
        identical(class(result), "imputefinder_result") &&
        identical(names(result), .CLASSIC_RESULT_FIELDS) &&
        .valid_classic_data(result$data) &&
        is.data.frame(result$classifications) &&
        is.list(result$groups) &&
        is.data.frame(result$feature_status) &&
        is.numeric(result$cutoffs) &&
        is.list(result$cutoff_diagnostics) &&
        is.list(result$profiles) &&
        is.data.frame(result$seed_log) &&
        is.character(result$groups_by_sample) &&
        (is.null(result$call) || is.language(result$call))
    if (!valid) {
        .abort_sidecar(
            "Stored classic branch is not a complete imputefinder result.",
            "imputefinder_analysis_schema_error",
            field = "classic"
        )
    }

    tryCatch(
        .validate_result_for_presentation(result),
        error = function(error) {
            .abort_sidecar(
                "Stored classic result fails its stable result schema.",
                "imputefinder_analysis_schema_error",
                field = "classic"
            )
        }
    )
    invisible(result)
}

.validate_classic_branch <- function(classic) {
    if (inherits(classic, "imputefinder_classic_failure")) {
        return(.validate_classic_failure(classic))
    }
    .validate_classic_result_shape(classic)
}

.classic_sample_names <- function(classic) {
    if (is.matrix(classic$data)) {
        colnames(classic$data)
    } else {
        colnames(classic$data)
    }
}

.validate_classic_design_alignment <- function(classic, design, input) {
    if (inherits(classic, "imputefinder_classic_failure")) {
        return(invisible(classic))
    }
    condition_role <- design$declared$roles$condition
    expected_groups <- stats::setNames(
        design$declared$sample_data[[condition_role]],
        input$sample_names
    )
    features <- classic$feature_status$feature
    valid <- is.character(features) &&
        identical(
            .canonical_sidecar_names(features),
            input$feature_names
        ) &&
        identical(
            .canonical_sidecar_names(.classic_sample_names(classic)),
            input$sample_names
        ) &&
        identical(classic$groups_by_sample, expected_groups)
    if (!valid) {
        .abort_sidecar(
            "Classic result, design, and input identities do not align.",
            "imputefinder_analysis_schema_error",
            field = "classic"
        )
    }

    invisible(classic)
}

.new_sidecar_provenance <- function(input, classic, call) {
    failures <- if (inherits(classic, "imputefinder_classic_failure")) {
        list(classic = classic)
    } else {
        list()
    }

    list(
        call = call,
        seeds = list(),
        hashes = list(input = input$fingerprint),
        warnings = list(),
        failures = failures,
        assumptions = list(
            scale = input$scale,
            acquisition = if (is.null(input$acquisition)) {
                "unavailable"
            } else {
                "declared"
            }
        ),
        training_scope = list(status = "not_started")
    )
}

.validate_sidecar_provenance <- function(provenance, input, classic) {
    expected_failures <- if (inherits(
        classic,
        "imputefinder_classic_failure"
    )) {
        list(classic = classic)
    } else {
        list()
    }
    expected_acquisition <- if (is.null(input$acquisition)) {
        "unavailable"
    } else {
        "declared"
    }
    valid <- is.list(provenance) &&
        identical(names(provenance), .SIDECAR_PROVENANCE_FIELDS) &&
        (is.null(provenance$call) || is.language(provenance$call)) &&
        is.list(provenance$seeds) &&
        identical(provenance$hashes, list(input = input$fingerprint)) &&
        is.list(provenance$warnings) &&
        identical(provenance$failures, expected_failures) &&
        is.list(provenance$assumptions) &&
        identical(
            names(provenance$assumptions),
            c("scale", "acquisition")
        ) &&
        identical(
            provenance$assumptions$scale,
            .SIDECAR_SCALE_DECLARATION
        ) &&
        identical(
            provenance$assumptions$acquisition,
            expected_acquisition
        ) &&
        is.list(provenance$training_scope) &&
        identical(provenance$training_scope$status, "not_started")
    if (!valid) {
        .abort_sidecar(
            "Stored analysis provenance is malformed.",
            "imputefinder_analysis_schema_error",
            field = "provenance"
        )
    }

    invisible(provenance)
}

.new_imputefinder_analysis <- function(classic, design, input, call = NULL) {
    if (!(is.null(call) || is.language(call))) {
        .abort_sidecar(
            "Analysis provenance call must be language or NULL.",
            "imputefinder_analysis_schema_error",
            field = "call"
        )
    }
    .validate_sidecar_input(input)
    design_record <- .new_sidecar_design_record(design, input)
    .validate_classic_branch(classic)
    .validate_classic_design_alignment(classic, design_record, input)

    analysis <- structure(
        list(
            classic = classic,
            spec = .new_sidecar_spec(),
            design = design_record,
            input = input,
            sentinel = NULL,
            stability = NULL,
            provenance = .new_sidecar_provenance(input, classic, call)
        ),
        class = "imputefinder_analysis"
    )
    .validate_imputefinder_analysis(analysis)
    analysis
}

.validate_imputefinder_analysis_shape <- function(analysis) {
    valid <- is.list(analysis) &&
        identical(class(analysis), "imputefinder_analysis") &&
        identical(names(analysis), .SIDECAR_ANALYSIS_FIELDS) &&
        (is.null(analysis$sentinel) || is.list(analysis$sentinel)) &&
        (is.null(analysis$stability) || is.list(analysis$stability))
    if (!valid) {
        .abort_sidecar(
            "Object does not satisfy the imputefinder-analysis schema.",
            "imputefinder_analysis_schema_error",
            field = "analysis"
        )
    }

    invisible(analysis)
}

.validate_imputefinder_analysis <- function(analysis) {
    .validate_imputefinder_analysis_shape(analysis)
    .validate_sidecar_spec(analysis$spec)
    .validate_sidecar_input(analysis$input)
    .validate_sidecar_design(analysis$design, analysis$input)
    .validate_classic_branch(analysis$classic)
    .validate_classic_design_alignment(
        analysis$classic,
        analysis$design,
        analysis$input
    )
    .validate_sidecar_provenance(
        analysis$provenance,
        analysis$input,
        analysis$classic
    )

    invisible(analysis)
}
