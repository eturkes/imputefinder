.ASSOCIATION_CANDIDATE_EVIDENCE_SCHEMA <- "m13a_candidate_evidence_v1"
.ASSOCIATION_CANDIDATE_EVIDENCE_FIXTURE_SCHEMA <-
    "m13a_candidate_evidence_fixture_v1"
.ASSOCIATION_CANDIDATE_EVIDENCE_FIELDS <- c(
    "schema", "contract_hash", "protocol_hash", "gate_registry_hash",
    "implementation_hash", "effective_manifest_hash",
    "generator_protocol_hash", "artifact_inventory_hash", "candidate_ids",
    "scenario_ids", "replicate_ids", "run_bindings", "outcome_bindings",
    "screening_metrics", "full_gate_results", "ranking_metrics",
    "candidate_order", "selected_candidate", "state", "evidence_hash"
)
.ASSOCIATION_RUN_BINDING_FIELDS <- c(
    "candidate_id", "scenario_id", "replicate", "input_sha256",
    "result_sha256", "result_kind", "elapsed_seconds", "status",
    "failure_code"
)
.ASSOCIATION_OUTCOME_BINDING_FIELDS <- c(
    "candidate_id", "scenario_id", "replicate", "stratum", "hypothesis",
    "status", "code"
)
.ASSOCIATION_SCREENING_METRIC_FIELDS <- c(
    "candidate_id", "metric", "numerator", "denominator", "estimate",
    "licensed_scope", "status", "passed"
)
.ASSOCIATION_FULL_GATE_RESULT_FIELDS <- c(
    "candidate_id", "gate_id", "m12_gate_id", "stage", "metric",
    "estimate", "numerator", "denominator", "lower", "upper",
    "operator", "threshold", "status", "passed", "failure_treatment",
    "evidence_hash", "note"
)
.ASSOCIATION_RANKING_METRIC_FIELDS <- c(
    "candidate_id", "scope_coverage", "power", "absolute_bias",
    "median_runtime"
)
.ASSOCIATION_RANKING_NEAR_TIE <- 0.02
.ASSOCIATION_RANKING_ULP_MULTIPLIER <- 8

.ASSOCIATION_STUDY_RESOLVER_SCHEMA <- "m13a_candidate_study_resolver_v1"
.ASSOCIATION_STUDY_RESOLVER_FIELDS <- c(
    "schema", "contract_hash", "protocol_hash", "generator_protocol_hash",
    "implementation_hash", "mode", "source_sha256", "root",
    "artifact_root", "input", "result"
)
.ASSOCIATION_STUDY_RESOLVER_MODES <- c("fixture", "frozen_protocol")
.ASSOCIATION_STUDY_ARTIFACT_DIRECTORY <- ".agent/m13a-candidate-study"
.ASSOCIATION_STUDY_MANIFEST_ARTIFACT <- "implementation-manifest.rds"
.ASSOCIATION_STUDY_EVIDENCE_ARTIFACT <- "candidate-evidence.rds"
.ASSOCIATION_STUDY_INPUT_SCHEMA <- "m13a_candidate_study_input_v1"
.ASSOCIATION_STUDY_INPUT_FIELDS <- c(
    "schema", "scenario_id", "replicate", "generator_protocol_hash",
    "data", "design", "missing_probability"
)
.ASSOCIATION_STUDY_RESULT_SCHEMA <- "m13a_candidate_study_result_v1"
.ASSOCIATION_STUDY_RESULT_FIELDS <- c(
    "schema", "candidate_id", "scenario_id", "replicate", "input_sha256",
    "result", "elapsed_seconds"
)
.ASSOCIATION_EXECUTION_FAILURE_FIELDS <- c(
    "status", "quantity", "code", "message", "requires"
)
.ASSOCIATION_PANEL_ABSTENTION <- "association_panel_abstention"

.ASSOCIATION_STUDY_SCENARIOS <- c(
    "dda_null_balanced", "dia_null_unequal", "dda_monotone_unequal",
    "dia_monotone_paired", "dda_mixed_outlier", "dia_blockwise_paired",
    "dda_batch_crossed", "dia_batch_partial", "dda_batch_perfect",
    "dia_structural_off", "dda_no_cliff", "dia_no_cliff_low_support",
    "dda_grouped_leakage_trap"
)
.ASSOCIATION_STUDY_REPLICATES <- 1:32
.ASSOCIATION_NULL_SCENARIOS <- c(
    "dda_null_balanced", "dia_null_unequal"
)
.ASSOCIATION_TARGETS <- data.frame(
    scenario_id = c(
        "dda_batch_crossed", "dda_monotone_unequal",
        "dia_batch_partial", "dia_monotone_paired"
    ),
    acquisition = c("DDA", "DDA", "DIA", "DIA"),
    label = c(
        "batch[batch_2]", "condition[D]", "batch[batch_3]",
        "condition[B]"
    ),
    reference = c("batch_1", "A", "batch_1", "A"),
    stringsAsFactors = FALSE
)
.ASSOCIATION_STUDY_MANIFEST <- data.frame(
    scenario_id = .ASSOCIATION_STUDY_SCENARIOS,
    acquisition = c(
        "DDA", "DIA", "DDA", "DIA", "DDA", "DIA", "DDA", "DIA",
        "DDA", "DIA", "DDA", "DIA", "DDA"
    ),
    n_features = as.integer(c(
        3000, 3000, 4000, 4000, 3500, 3500, 4000, 4000, 3000, 3500,
        3000, 1800, 2500
    )),
    condition_sizes = c(
        "A=8;B=8", "A=4;B=7;C=10", "A=4;B=6;C=8;D=10",
        "A=10;B=10", "A=6;B=6;C=6", "A=8;B=8", "A=12;B=12",
        "A=6;B=8;C=10", "A=8;B=8", "A=6;B=6;C=6;D=6",
        "A=8;B=8", "A=3;B=3", "A=12;B=12"
    ),
    nuisance = c(
        "", "", "run_order", "", "run_order", "", "batch", "batch",
        "batch", "", "", "", "technical_replicate"
    ),
    block = c(
        "", "", "", "subject", "", "subject", "", "", "", "", "",
        "", "subject"
    ),
    stringsAsFactors = FALSE
)

.ASSOCIATION_SCREENING_METRICS <- c(
    "common_opportunity_coverage",
    "family_false_flag_clopper_pearson_upper95",
    "maximum_acquisition_false_flag_fraction",
    "candidate_specific_projection_95_interval_coverage",
    "holm_target_detection_fraction",
    "median_absolute_candidate_specific_projection_bias"
)
.ASSOCIATION_SCREENING_OPERATORS <- c(
    family_false_flag_clopper_pearson_upper95 = "<=",
    maximum_acquisition_false_flag_fraction = "<=",
    candidate_specific_projection_95_interval_coverage = ">=",
    holm_target_detection_fraction = ">=",
    median_absolute_candidate_specific_projection_bias = "<="
)
.ASSOCIATION_SCREENING_THRESHOLDS <- c(
    family_false_flag_clopper_pearson_upper95 = 0.15,
    maximum_acquisition_false_flag_fraction = 0.125,
    candidate_specific_projection_95_interval_coverage = 0.90,
    holm_target_detection_fraction = 0.70,
    median_absolute_candidate_specific_projection_bias = 0.03
)
.ASSOCIATION_GATE_BOOTSTRAP_DRAWS <- 9999L

.ASSOCIATION_OPPORTUNITY_AUDIT_FIELDS <- c(
    "candidate_id", "scenario_id", "replicate", "stratum", "hypothesis",
    "available"
)
.ASSOCIATION_NULL_AUDIT_FIELDS <- c(
    "candidate_id", "scenario_id", "replicate", "acquisition",
    "available_p_count", "false_flag"
)
.ASSOCIATION_TARGET_AUDIT_FIELDS <- c(
    "candidate_id", "scenario_id", "replicate", "acquisition",
    "target_label", "outcome_status", "code", "truth_status", "truth",
    "effect", "conf_low", "conf_high", "adjusted_p"
)

.abort_association_evidence <- function(message) {
    .abort_association(
        message,
        "imputefinder_association_evidence_error"
    )
}

.association_study_gate_registry <- function() {
    data.frame(
        gate_id = c(
            "a_assoc_null_upper_v4", "a_assoc_null_stratum_v4",
            "a_assoc_interval_coverage_v4",
            "a_assoc_alternative_power_v4", "a_assoc_effect_bias_v4"
        ),
        m12_gate_id = c(
            "a_assoc_null_upper", "a_assoc_null_stratum",
            "a_assoc_interval_coverage", "a_assoc_alternative_power",
            "a_assoc_effect_bias"
        ),
        registry_version = rep("m13a_gate_registry_v4", 5L),
        metric = c(
            "family_false_flag_clopper_pearson_upper95",
            "maximum_acquisition_false_flag_fraction",
            "fixed_target_95_interval_coverage_fraction",
            "holm_alternative_detection_fraction",
            "median_absolute_effect_bias"
        ),
        operator = c("<=", "<=", ">=", ">=", "<="),
        threshold = c(0.15, 0.125, 0.90, 0.70, 0.03),
        failure_treatment = rep(
            paste0(
                "miss kills or parks the A association panel; mandatory ",
                "design core remains"
            ),
            5L
        ),
        gate_binding_hash = c(
            "4f00acdaa6db9e9035e4b4b7ee5980f7b3beb835474ca1c451bd5fae4cda1c02",
            "de518da9a2aed4bdb1cfdbc50831136723332dfa5497f6f7a1af209eebe00dea",
            "0a1fe11d0760a924af5e91785f1843914b9b7c7a3e95b69282d8b092a045bcc7",
            "9c3ba375943078477f491a4825265e783d8426099437da64c6369070a343418b",
            "63df7d7ef92ee9b6c29f74f1c6f6667d064c5138e33360665b352e11de79de58"
        ),
        stringsAsFactors = FALSE
    )
}

.association_parse_condition_sizes <- function(value) {
    pieces <- strsplit(value, ";", fixed = TRUE)[[1L]]
    fields <- strsplit(pieces, "=", fixed = TRUE)
    stats::setNames(
        as.integer(vapply(fields, `[[`, character(1L), 2L)),
        vapply(fields, `[[`, character(1L), 1L)
    )
}

.association_study_manifest_row <- function(scenario_id) {
    index <- match(scenario_id, .ASSOCIATION_STUDY_MANIFEST$scenario_id)
    if (is.na(index)) {
        .abort_association_evidence("Association study scenario is unknown.")
    }
    .ASSOCIATION_STUDY_MANIFEST[index, , drop = FALSE]
}

.new_association_study_input <- function(
    scenario_id,
    replicate,
    data,
    design,
    missing_probability
) {
    input <- structure(
        list(
            schema = .ASSOCIATION_STUDY_INPUT_SCHEMA,
            scenario_id = scenario_id,
            replicate = replicate,
            generator_protocol_hash = .ASSOCIATION_GENERATOR_PROTOCOL_HASH,
            data = data,
            design = design,
            missing_probability = missing_probability
        ),
        class = "imputefinder_association_study_input"
    )
    .validate_association_study_input(input, scenario_id, replicate)
    input
}

.validate_association_study_input <- function(
    input,
    scenario_id,
    replicate
) {
    valid_key <- .association_scalar_character(scenario_id) &&
        scenario_id %in% .ASSOCIATION_STUDY_SCENARIOS &&
        is.integer(replicate) && length(replicate) == 1L &&
        !is.na(replicate) && replicate %in% .ASSOCIATION_STUDY_REPLICATES
    valid <- valid_key && is.list(input) && identical(
        class(input),
        "imputefinder_association_study_input"
    ) && identical(names(input), .ASSOCIATION_STUDY_INPUT_FIELDS) &&
        identical(input$schema, .ASSOCIATION_STUDY_INPUT_SCHEMA) &&
        identical(input$scenario_id, scenario_id) &&
        identical(input$replicate, replicate) &&
        identical(
            input$generator_protocol_hash,
            .ASSOCIATION_GENERATOR_PROTOCOL_HASH
        ) && is.matrix(input$missing_probability) &&
        is.double(input$missing_probability) &&
        identical(dim(input$missing_probability), dim(input$data)) &&
        identical(dimnames(input$missing_probability), dimnames(input$data)) &&
        !anyNA(input$missing_probability) &&
        all(is.finite(input$missing_probability)) &&
        all(input$missing_probability >= 0 &
            input$missing_probability <= 1)
    if (!valid) {
        .abort_association_evidence(
            "Association study input record is malformed or detached."
        )
    }
    .validate_matrix_data(input$data)
    .validate_missingness_design(input$design)
    aligned <- .align_missingness_design(input$design, colnames(input$data))
    if (!identical(aligned, input$design)) {
        .abort_association_evidence(
            "Association study design is not canonically input-aligned."
        )
    }
    manifest <- .association_study_manifest_row(scenario_id)
    conditions <- .association_parse_condition_sizes(
        manifest$condition_sizes[[1L]]
    )
    condition_column <- input$design$roles$condition
    observed_counts <- table(.canonical_design_text(
        input$design$sample_data[[condition_column]]
    ))
    acquisition_column <- input$design$roles$acquisition
    acquisition <- if (length(acquisition_column)) {
        unique(.canonical_design_text(
            input$design$sample_data[[acquisition_column]]
        ))
    } else {
        character()
    }
    expected_nuisance <- manifest$nuisance[[1L]]
    if (!nzchar(expected_nuisance)) {
        expected_nuisance <- character()
    }
    expected_block <- manifest$block[[1L]]
    if (!nzchar(expected_block)) {
        expected_block <- character()
    }
    structure_valid <- nrow(input$data) == manifest$n_features[[1L]] &&
        ncol(input$data) == sum(conditions) &&
        identical(names(observed_counts), names(conditions)) &&
        identical(as.integer(observed_counts), unname(conditions)) &&
        identical(acquisition, manifest$acquisition[[1L]]) &&
        identical(input$design$roles$condition, "condition") &&
        identical(input$design$roles$nuisance, expected_nuisance) &&
        identical(input$design$roles$block, expected_block) &&
        identical(
            input$design$roles$acquisition,
            "acquisition_mode"
        ) && length(input$design$interactions) == 0L &&
        all(!is.na(input$data[input$missing_probability == 0])) &&
        all(is.na(input$data[input$missing_probability == 1]))
    if (!structure_valid) {
        .abort_association_evidence(
            "Association study input violates its frozen scenario allocation."
        )
    }
    invisible(input)
}

.new_association_execution_failure <- function(code, message) {
    valid <- .association_scalar_character(code) &&
        identical(code, "association_numerical_failure") &&
        .association_scalar_character(message)
    if (!valid) {
        .abort_association_evidence(
            "Association execution failure requires one code and message."
        )
    }
    structure(
        list(
            status = "failed",
            quantity = "association",
            code = code,
            message = message,
            requires = "successful candidate execution"
        ),
        class = "imputefinder_association_execution_failure"
    )
}

.valid_association_execution_failure <- function(result) {
    is.list(result) && identical(
        class(result),
        "imputefinder_association_execution_failure"
    ) && identical(names(result), .ASSOCIATION_EXECUTION_FAILURE_FIELDS) &&
        identical(result$status, "failed") &&
        identical(result$quantity, "association") &&
        identical(result$code, "association_numerical_failure") &&
        .association_scalar_character(result$message) &&
        identical(result$requires, "successful candidate execution")
}

.new_association_study_result <- function(
    candidate_id,
    scenario_id,
    replicate,
    input_sha256,
    result,
    elapsed_seconds
) {
    record <- structure(
        list(
            schema = .ASSOCIATION_STUDY_RESULT_SCHEMA,
            candidate_id = candidate_id,
            scenario_id = scenario_id,
            replicate = replicate,
            input_sha256 = input_sha256,
            result = result,
            elapsed_seconds = elapsed_seconds
        ),
        class = "imputefinder_association_study_result"
    )
    .validate_association_study_result_record(
        record,
        candidate_id,
        scenario_id,
        replicate,
        input_sha256
    )
    record
}

.validate_association_study_result_record <- function(
    record,
    candidate_id,
    scenario_id,
    replicate,
    input_sha256
) {
    valid <- is.list(record) && identical(
        class(record),
        "imputefinder_association_study_result"
    ) && identical(names(record), .ASSOCIATION_STUDY_RESULT_FIELDS) &&
        identical(record$schema, .ASSOCIATION_STUDY_RESULT_SCHEMA) &&
        identical(record$candidate_id, candidate_id) &&
        identical(record$scenario_id, scenario_id) &&
        identical(record$replicate, replicate) &&
        identical(record$input_sha256, input_sha256) &&
        .association_sha256(input_sha256) &&
        .association_double_scalar(record$elapsed_seconds) &&
        record$elapsed_seconds >= 0
    if (!valid) {
        .abort_association_evidence(
            "Association study result record is malformed or detached."
        )
    }
    invisible(record)
}

.new_association_study_resolver <- function(
    input,
    result,
    implementation_hash
) {
    resolver <- structure(
        list(
            schema = .ASSOCIATION_STUDY_RESOLVER_SCHEMA,
            contract_hash = .ASSOCIATION_CONTRACT_HASH,
            protocol_hash = .ASSOCIATION_PROTOCOL_HASH,
            generator_protocol_hash = .ASSOCIATION_GENERATOR_PROTOCOL_HASH,
            implementation_hash = implementation_hash,
            mode = "fixture",
            source_sha256 = NA_character_,
            root = NA_character_,
            artifact_root = NA_character_,
            input = input,
            result = result
        ),
        class = "imputefinder_association_study_resolver"
    )
    .validate_association_study_resolver(resolver, implementation_hash)
    resolver
}

.association_protocol_artifact_layout <- function(include_evidence = FALSE) {
    if (!is.logical(include_evidence) || length(include_evidence) != 1L ||
        is.na(include_evidence)) {
        .abort_association_evidence(
            "Association artifact-layout request is malformed."
        )
    }
    input_directories <- file.path(
        "inputs",
        .ASSOCIATION_STUDY_SCENARIOS
    )
    result_candidate_directories <- file.path(
        "results",
        .ASSOCIATION_CANDIDATES
    )
    result_scenario_directories <- unlist(lapply(
        .ASSOCIATION_CANDIDATES,
        function(candidate_id) file.path(
            "results",
            candidate_id,
            .ASSOCIATION_STUDY_SCENARIOS
        )
    ), use.names = FALSE)
    directories <- c(
        "inputs", input_directories, "results",
        result_candidate_directories, result_scenario_directories
    )
    input_files <- unlist(lapply(
        .ASSOCIATION_STUDY_SCENARIOS,
        function(scenario_id) file.path(
            "inputs",
            scenario_id,
            sprintf(
                "replicate-%03d.rds",
                .ASSOCIATION_STUDY_REPLICATES
            )
        )
    ), use.names = FALSE)
    result_files <- unlist(lapply(
        .ASSOCIATION_CANDIDATES,
        function(candidate_id) unlist(lapply(
            .ASSOCIATION_STUDY_SCENARIOS,
            function(scenario_id) file.path(
                "results",
                candidate_id,
                scenario_id,
                sprintf(
                    "replicate-%03d.rds",
                    .ASSOCIATION_STUDY_REPLICATES
                )
            )
        ), use.names = FALSE)
    ), use.names = FALSE)
    files <- c(
        .ASSOCIATION_STUDY_MANIFEST_ARTIFACT,
        input_files,
        result_files,
        if (include_evidence) .ASSOCIATION_STUDY_EVIDENCE_ARTIFACT
    )
    output <- data.frame(
        relative_path = c(directories, files),
        kind = c(rep("directory", length(directories)), rep("file", length(files))),
        stringsAsFactors = FALSE
    )
    output <- output[
        order(output$relative_path, method = "radix"),
        ,
        drop = FALSE
    ]
    row.names(output) <- NULL
    output
}

.association_read_locked_implementation_manifest <- function(artifact_root) {
    valid <- .association_scalar_character(artifact_root) &&
        dir.exists(artifact_root) && identical(
            artifact_root,
            normalizePath(artifact_root, winslash = "/", mustWork = TRUE)
        )
    path <- file.path(
        artifact_root,
        .ASSOCIATION_STUDY_MANIFEST_ARTIFACT
    )
    valid <- valid && utils::file_test("-f", path) &&
        identical(Sys.readlink(path), "") && identical(
            path,
            normalizePath(path, winslash = "/", mustWork = TRUE)
        )
    if (!valid) {
        .abort_association_evidence(
            "Association locked implementation artifact is absent or unsafe."
        )
    }
    tryCatch(
        readRDS(path),
        error = function(error) {
            .abort_association_evidence(
                "Association locked implementation artifact is unreadable."
            )
        }
    )
}

.association_protocol_artifact_inventory <- function(
    artifact_root,
    implementation_hash,
    complete,
    evidence_policy = c("allow", "absent", "exact"),
    evidence = NULL
) {
    evidence_policy <- match.arg(evidence_policy)
    valid <- .association_scalar_character(artifact_root) &&
        .association_sha256(implementation_hash) &&
        is.logical(complete) && length(complete) == 1L && !is.na(complete) &&
        dir.exists(artifact_root) && identical(
            artifact_root,
            normalizePath(artifact_root, winslash = "/", mustWork = TRUE)
        ) && (!identical(evidence_policy, "exact") || !is.null(evidence))
    if (!valid) {
        .abort_association_evidence(
            "Association artifact-inventory request is malformed."
        )
    }
    actual <- list.files(
        artifact_root,
        all.files = TRUE,
        full.names = FALSE,
        recursive = TRUE,
        include.dirs = TRUE,
        no.. = TRUE
    )
    actual <- sort(actual, method = "radix")
    allowed <- .association_protocol_artifact_layout(TRUE)
    include_evidence <- identical(evidence_policy, "exact")
    required <- .association_protocol_artifact_layout(include_evidence)
    valid <- !anyDuplicated(actual) && all(
        actual %in% allowed$relative_path
    ) && .ASSOCIATION_STUDY_MANIFEST_ARTIFACT %in% actual
    if (complete) {
        valid <- valid && identical(actual, required$relative_path)
    } else if (identical(evidence_policy, "absent")) {
        valid <- valid && !.ASSOCIATION_STUDY_EVIDENCE_ARTIFACT %in% actual
    }
    if (!valid) {
        .abort_association_evidence(
            "Association artifact inventory is incomplete or contains allocation drift."
        )
    }
    paths <- file.path(artifact_root, actual)
    links <- Sys.readlink(paths)
    information <- file.info(paths)
    kinds <- allowed$kind[match(actual, allowed$relative_path)]
    valid <- !anyNA(information$isdir) && !any(nzchar(links)) && all(
        information$isdir == (kinds == "directory")
    ) && all(vapply(seq_along(paths), function(index) {
        if (identical(kinds[[index]], "directory")) {
            dir.exists(paths[[index]])
        } else {
            utils::file_test("-f", paths[[index]])
        }
    }, logical(1L))) && identical(
        paths,
        unname(vapply(
            paths,
            normalizePath,
            character(1L),
            winslash = "/",
            mustWork = TRUE
        ))
    )
    if (!valid) {
        .abort_association_evidence(
            "Association artifact inventory contains an unsafe path or file type."
        )
    }
    locked <- .association_read_locked_implementation_manifest(artifact_root)
    if (!is.list(locked) || !.association_scalar_character(
        locked$manifest_hash
    ) || !identical(locked$manifest_hash, implementation_hash)) {
        .abort_association_evidence(
            "Association artifact inventory is implementation-detached."
        )
    }
    evidence_path <- file.path(
        artifact_root,
        .ASSOCIATION_STUDY_EVIDENCE_ARTIFACT
    )
    if (identical(evidence_policy, "exact")) {
        stored_evidence <- tryCatch(
            readRDS(evidence_path),
            error = function(error) {
                .abort_association_evidence(
                    "Association locked evidence artifact is unreadable."
                )
            }
        )
        if (!identical(stored_evidence, evidence)) {
            .abort_association_evidence(
                "Association locked evidence artifact is detached."
            )
        }
    }
    file_rows <- kinds == "file"
    sha256 <- rep(NA_character_, length(actual))
    sha256[file_rows] <- unname(tools::sha256sum(paths[file_rows]))
    data.frame(
        relative_path = actual,
        kind = kinds,
        sha256 = sha256,
        stringsAsFactors = FALSE
    )
}

.association_protocol_allocation_inventory_hash <- function(inventory) {
    basic <- is.data.frame(inventory) && identical(
        names(inventory),
        c("relative_path", "kind", "sha256")
    ) && is.character(inventory$relative_path) &&
        is.character(inventory$kind) && is.character(inventory$sha256)
    if (!basic) {
        .abort_association_evidence(
            "Association allocation inventory cannot be hashed exactly."
        )
    }
    allocation <- inventory[
        inventory$relative_path != .ASSOCIATION_STUDY_EVIDENCE_ARTIFACT,
        ,
        drop = FALSE
    ]
    row.names(allocation) <- NULL
    expected <- .association_protocol_artifact_layout(FALSE)
    valid <- identical(
        allocation[c("relative_path", "kind")],
        expected
    ) && all(is.na(allocation$sha256[allocation$kind == "directory"])) &&
        all(vapply(
            allocation$sha256[allocation$kind == "file"],
            .association_sha256,
            logical(1L)
        ))
    if (!valid) {
        .abort_association_evidence(
            "Association allocation inventory cannot be hashed exactly."
        )
    }
    descriptor <- data.frame(
        position = 1L,
        field = "schema",
        value = "m13a_protocol_allocation_inventory_v1",
        type = "character",
        stringsAsFactors = FALSE
    )
    .association_evidence_bundle_hash(list(
        descriptor = descriptor,
        allocation = allocation
    ))
}

.new_association_protocol_study_resolver <- function(
    implementation_manifest,
    root = "."
) {
    root <- .association_manifest_root(root)
    expected_artifact_root <- file.path(
        root,
        .ASSOCIATION_STUDY_ARTIFACT_DIRECTORY
    )
    artifact_root <- normalizePath(
        expected_artifact_root,
        winslash = "/",
        mustWork = TRUE
    )
    if (!dir.exists(artifact_root) || !identical(
        artifact_root,
        expected_artifact_root
    ) || !startsWith(artifact_root, paste0(root, "/"))) {
        .abort_association_evidence(
            "Association protocol artifact root is absent or outside the repository."
        )
    }
    locked_manifest <- .association_read_locked_implementation_manifest(
        artifact_root
    )
    if (!identical(locked_manifest, implementation_manifest)) {
        .abort_association_evidence(
            "Association resolver manifest differs from the pre-allocation lock."
        )
    }
    .validate_association_implementation_manifest(
        implementation_manifest,
        root
    )
    resolver <- structure(
        list(
            schema = .ASSOCIATION_STUDY_RESOLVER_SCHEMA,
            contract_hash = .ASSOCIATION_CONTRACT_HASH,
            protocol_hash = .ASSOCIATION_PROTOCOL_HASH,
            generator_protocol_hash = .ASSOCIATION_GENERATOR_PROTOCOL_HASH,
            implementation_hash = implementation_manifest$manifest_hash,
            mode = "frozen_protocol",
            source_sha256 = .association_source_manifest_sha256(
                implementation_manifest$source_files
            ),
            root = root,
            artifact_root = artifact_root,
            input = NULL,
            result = NULL
        ),
        class = "imputefinder_association_study_resolver"
    )
    .validate_association_study_resolver(
        resolver,
        implementation_manifest$manifest_hash
    )
    .association_protocol_artifact_inventory(
        artifact_root,
        implementation_manifest$manifest_hash,
        complete = FALSE,
        evidence_policy = "allow"
    )
    resolver
}

.validate_association_study_resolver <- function(
    resolver,
    implementation_hash
) {
    valid <- is.list(resolver) && identical(
        class(resolver),
        "imputefinder_association_study_resolver"
    ) && identical(names(resolver), .ASSOCIATION_STUDY_RESOLVER_FIELDS) &&
        identical(resolver$schema, .ASSOCIATION_STUDY_RESOLVER_SCHEMA) &&
        identical(resolver$contract_hash, .ASSOCIATION_CONTRACT_HASH) &&
        identical(resolver$protocol_hash, .ASSOCIATION_PROTOCOL_HASH) &&
        identical(
            resolver$generator_protocol_hash,
            .ASSOCIATION_GENERATOR_PROTOCOL_HASH
        ) && identical(resolver$implementation_hash, implementation_hash) &&
        .association_sha256(implementation_hash) &&
        .association_scalar_character(resolver$mode) &&
        resolver$mode %in% .ASSOCIATION_STUDY_RESOLVER_MODES &&
        .association_scalar_character(resolver$mode)
    provenance_valid <- if (valid && identical(
        resolver$mode,
        "frozen_protocol"
    )) {
        .association_sha256(resolver$source_sha256) &&
            .association_scalar_character(resolver$root) &&
            .association_scalar_character(resolver$artifact_root) &&
            dir.exists(resolver$root) && dir.exists(resolver$artifact_root) &&
            identical(
                resolver$root,
                normalizePath(
                    resolver$root,
                    winslash = "/",
                    mustWork = TRUE
                )
            ) && identical(
                resolver$artifact_root,
                normalizePath(
                    file.path(
                        resolver$root,
                        .ASSOCIATION_STUDY_ARTIFACT_DIRECTORY
                    ),
                    winslash = "/",
                    mustWork = TRUE
                )
            ) && startsWith(
                resolver$artifact_root,
                paste0(resolver$root, "/")
            ) && is.null(resolver$input) && is.null(resolver$result)
    } else if (valid) {
        is.na(resolver$source_sha256) &&
            is.na(resolver$root) && is.na(resolver$artifact_root) &&
            is.function(resolver$input) && is.function(resolver$result)
    } else {
        FALSE
    }
    if (!valid || !provenance_valid) {
        .abort_association_evidence(
            "Association study resolver is malformed or implementation-detached."
        )
    }
    invisible(resolver)
}

.association_protocol_artifact_path <- function(
    resolver,
    quantity,
    scenario_id,
    replicate,
    candidate_id = NULL,
    must_exist = TRUE
) {
    .validate_association_study_resolver(
        resolver,
        resolver$implementation_hash
    )
    valid <- identical(resolver$mode, "frozen_protocol") &&
        .association_scalar_character(quantity) &&
        quantity %in% c("input", "result") &&
        scenario_id %in% .ASSOCIATION_STUDY_SCENARIOS &&
        is.integer(replicate) && length(replicate) == 1L &&
        !is.na(replicate) && replicate %in% .ASSOCIATION_STUDY_REPLICATES &&
        is.logical(must_exist) && length(must_exist) == 1L &&
        !is.na(must_exist)
    if (identical(quantity, "result")) {
        valid <- valid && .association_scalar_character(candidate_id) &&
            candidate_id %in% .ASSOCIATION_CANDIDATES
    } else {
        valid <- valid && is.null(candidate_id)
    }
    if (!valid) {
        .abort_association_evidence(
            "Association protocol artifact key is malformed."
        )
    }
    relative <- if (identical(quantity, "input")) {
        file.path(
            "inputs",
            scenario_id,
            sprintf("replicate-%03d.rds", replicate)
        )
    } else {
        file.path(
            "results",
            candidate_id,
            scenario_id,
            sprintf("replicate-%03d.rds", replicate)
        )
    }
    target <- file.path(resolver$artifact_root, relative)
    if (must_exist) {
        information <- file.info(target)
        if (nrow(information) != 1L || is.na(information$isdir) ||
            information$isdir) {
            .abort_association_evidence(
                paste0(
                    "Association protocol ",
                    quantity,
                    " artifact is absent."
                )
            )
        }
        target <- normalizePath(target, winslash = "/", mustWork = TRUE)
    } else {
        ancestor <- dirname(target)
        while (!dir.exists(ancestor) && !identical(
            ancestor,
            resolver$artifact_root
        )) {
            ancestor <- dirname(ancestor)
        }
        ancestor <- normalizePath(
            ancestor,
            winslash = "/",
            mustWork = TRUE
        )
        if (!identical(ancestor, resolver$artifact_root) && !startsWith(
            ancestor,
            paste0(resolver$artifact_root, "/")
        )) {
            .abort_association_evidence(
                "Association artifact directory resolves outside its sealed root."
            )
        }
    }
    if (!startsWith(target, paste0(resolver$artifact_root, "/"))) {
        .abort_association_evidence(
            "Association protocol artifact resolves outside its sealed root."
        )
    }
    target
}

.association_read_protocol_input <- function(
    resolver,
    scenario_id,
    replicate
) {
    path <- .association_protocol_artifact_path(
        resolver,
        "input",
        scenario_id,
        replicate
    )
    tryCatch(
        readRDS(path),
        error = function(error) {
            .abort_association_evidence(
                "Association protocol input artifact is unreadable."
            )
        }
    )
}

.association_read_protocol_result <- function(
    resolver,
    candidate_id,
    scenario_id,
    replicate
) {
    path <- .association_protocol_artifact_path(
        resolver,
        "result",
        scenario_id,
        replicate,
        candidate_id
    )
    tryCatch(
        readRDS(path),
        error = function(error) {
            .abort_association_evidence(
                "Association protocol result artifact is unreadable."
            )
        }
    )
}

.association_protocol_runtime <- function(resolver) {
    .validate_association_study_resolver(
        resolver,
        resolver$implementation_hash
    )
    if (!identical(resolver$mode, "frozen_protocol")) {
        .abort_association_evidence(
            "Association generator runtime requires a protocol resolver."
        )
    }
    previous <- getwd()
    on.exit(setwd(previous), add = TRUE)
    setwd(resolver$root)
    generator <- new.env(parent = baseenv())
    loaded <- tryCatch({
        sys.source(
            "dev/m12-generator-validation.R",
            envir = generator,
            keep.source = FALSE
        )
        TRUE
    }, error = function(error) FALSE)
    valid <- loaded && exists(
        "simulate_m12b_scenario",
        envir = generator,
        inherits = FALSE
    ) && exists(
        "m12b_protocol_hashes",
        envir = generator,
        inherits = FALSE
    )
    protocol_hash <- if (valid) {
        tryCatch(
            unname(generator$m12b_protocol_hashes()[["protocol"]]),
            error = function(error) NA_character_
        )
    } else {
        NA_character_
    }
    if (!valid || !identical(
        protocol_hash,
        .ASSOCIATION_GENERATOR_PROTOCOL_HASH
    )) {
        .abort_association_evidence(
            "Association generator runtime differs from the frozen protocol."
        )
    }
    list(simulate = generator$simulate_m12b_scenario)
}

.association_protocol_input <- function(
    runtime,
    scenario_id,
    replicate
) {
    simulation <- tryCatch(
        runtime$simulate(scenario_id, replicate),
        error = function(error) {
            .abort_association_evidence(
                "Frozen association generator replay failed."
            )
        }
    )
    nuisance <- switch(
        simulation$manifest$nuisance_role,
        none = NULL,
        run_order = "run_order",
        batch_crossed = "batch",
        batch_partial = "batch",
        batch_aliased = "batch",
        technical_replicate = "technical_replicate",
        .abort_association_evidence(
            "Frozen association generator declared an unknown nuisance role."
        )
    )
    block <- if (identical(
        simulation$manifest$block_role,
        "subject"
    )) {
        "subject"
    } else {
        NULL
    }
    columns <- c("condition", nuisance, block, "acquisition_mode")
    metadata <- simulation$design[, columns, drop = FALSE]
    design <- missingness_design(
        metadata,
        condition = "condition",
        nuisance = nuisance,
        block = block,
        acquisition = "acquisition_mode"
    )
    .new_association_study_input(
        scenario_id,
        replicate,
        simulation$data,
        design,
        simulation$missing_probability
    )
}

.association_execute_protocol_candidate <- function(
    candidate_id,
    preparation
) {
    if (inherits(preparation, "imputefinder_unavailable")) {
        return(preparation)
    }
    runner <- switch(
        candidate_id,
        a_fraction_ols_hc3_cr2 = .run_association_ols_hc3_cr2,
        a_fraction_freedman_lane = .run_association_freedman_lane,
        a_fraction_quasibinomial = .run_association_quasibinomial,
        NULL
    )
    if (is.null(runner)) {
        .abort_association_evidence(
            "Association protocol candidate is unknown."
        )
    }
    tryCatch(
        runner(preparation),
        error = function(error) .new_association_execution_failure(
            "association_numerical_failure",
            "Candidate engine raised an error under sealed execution."
        )
    )
}

.association_serialized_sha256 <- function(value) {
    unname(tools::sha256sum(bytes = serialize(
        value,
        connection = NULL,
        version = 3L
    )))
}

.association_study_input_payload_sha256 <- function(input) {
    .association_serialized_sha256(list(
        data = input$data,
        design = input$design,
        missing_probability = input$missing_probability
    ))
}

.association_study_input_sha256 <- function(input, payload_sha256 = NULL) {
    if (is.null(payload_sha256)) {
        payload_sha256 <- .association_study_input_payload_sha256(input)
    }
    valid <- is.list(input) &&
        identical(input$schema, .ASSOCIATION_STUDY_INPUT_SCHEMA) &&
        .association_scalar_character(input$scenario_id) &&
        is.integer(input$replicate) && length(input$replicate) == 1L &&
        !is.na(input$replicate) && .association_sha256(payload_sha256)
    if (!valid) {
        .abort_association_evidence(
            "Association study input cannot be hashed exactly."
        )
    }
    bytes <- charToRaw(enc2utf8(paste(
        "imputefinder:association-study-input:v1",
        .ASSOCIATION_GENERATOR_PROTOCOL_HASH,
        input$scenario_id,
        input$replicate,
        payload_sha256,
        sep = "\t"
    )))
    unname(tools::sha256sum(bytes = bytes))
}

.empty_association_run_bindings <- function() {
    data.frame(
        candidate_id = character(),
        scenario_id = character(),
        replicate = integer(),
        input_sha256 = character(),
        result_sha256 = character(),
        result_kind = character(),
        elapsed_seconds = double(),
        status = character(),
        failure_code = character(),
        stringsAsFactors = FALSE
    )[.ASSOCIATION_RUN_BINDING_FIELDS]
}

.empty_association_outcome_bindings <- function() {
    data.frame(
        candidate_id = character(),
        scenario_id = character(),
        replicate = integer(),
        stratum = character(),
        hypothesis = character(),
        status = character(),
        code = character(),
        stringsAsFactors = FALSE
    )[.ASSOCIATION_OUTCOME_BINDING_FIELDS]
}

.empty_association_opportunity_audit <- function() {
    data.frame(
        candidate_id = character(),
        scenario_id = character(),
        replicate = integer(),
        stratum = character(),
        hypothesis = character(),
        available = logical(),
        stringsAsFactors = FALSE
    )[.ASSOCIATION_OPPORTUNITY_AUDIT_FIELDS]
}

.empty_association_null_audit <- function() {
    data.frame(
        candidate_id = character(),
        scenario_id = character(),
        replicate = integer(),
        acquisition = character(),
        available_p_count = integer(),
        false_flag = logical(),
        stringsAsFactors = FALSE
    )[.ASSOCIATION_NULL_AUDIT_FIELDS]
}

.empty_association_target_audit <- function() {
    data.frame(
        candidate_id = character(),
        scenario_id = character(),
        replicate = integer(),
        acquisition = character(),
        target_label = character(),
        outcome_status = character(),
        code = character(),
        truth_status = character(),
        truth = double(),
        effect = double(),
        conf_low = double(),
        conf_high = double(),
        adjusted_p = double(),
        stringsAsFactors = FALSE
    )[.ASSOCIATION_TARGET_AUDIT_FIELDS]
}

.association_bind_audit_rows <- function(rows, empty) {
    rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
    if (!length(rows)) {
        return(empty())
    }
    output <- do.call(rbind, unname(rows))
    row.names(output) <- NULL
    output
}

.association_oracle_response <- function(input, preparation) {
    if (!inherits(
        preparation,
        "imputefinder_association_preparation"
    )) {
        return(NULL)
    }
    feature_order <- match(
        preparation$identity$feature_names,
        rownames(input$data)
    )
    sample_order <- match(
        preparation$identity$sample_names,
        colnames(input$data)
    )
    if (anyNA(feature_order) || anyNA(sample_order)) {
        return(NULL)
    }
    data <- input$data[feature_order, sample_order, drop = FALSE]
    probability <- input$missing_probability[
        feature_order,
        sample_order,
        drop = FALSE
    ]
    globally_observable <- rowSums(!is.na(data)) > 0L
    count <- as.integer(sum(globally_observable))
    if (count < 1L) {
        return(NULL)
    }
    response <- colMeans(
        1 - probability[globally_observable, , drop = FALSE]
    )
    names(response) <- colnames(data)
    if (any(!is.finite(response)) || any(response < 0 | response > 1)) {
        return(NULL)
    }
    list(count = count, response = response)
}

.association_target_hypothesis <- function(preparation, target_label) {
    selected <- which(preparation$hypotheses$label == target_label)
    if (length(selected) != 1L) {
        return(NULL)
    }
    index <- selected[[1L]]
    stratum <- preparation$strata[[
        preparation$hypotheses$stratum[[index]]
    ]]
    reference_index <- match(target_label, .ASSOCIATION_TARGETS$label)
    component <- preparation$hypotheses$components[[index]]
    variable_index <- if (length(component) == 1L) {
        match(component, stratum$core$model$variables$column)
    } else {
        NA_integer_
    }
    reference_valid <- !is.na(reference_index) && !is.na(variable_index) &&
        identical(preparation$hypotheses$kind[[index]], "main") &&
        identical(
            stratum$core$model$variables$encoding[[variable_index]],
            "treatment"
        ) && identical(
            stratum$core$model$variables$reference[[variable_index]],
            .ASSOCIATION_TARGETS$reference[[reference_index]]
        )
    residual_df <- length(stratum$samples) - stratum$core$rank$rank
    valid <- reference_valid && preparation$hypotheses$eligible[[index]] &&
        preparation$hypotheses$estimable[[index]] &&
        preparation$support$eligible[[index]] && residual_df >= 3L
    if (!valid) {
        return(NULL)
    }
    list(
        index = index,
        hypothesis = preparation$hypotheses[index, , drop = FALSE],
        stratum = stratum
    )
}

.validate_association_null_oracle <- function(oracle) {
    if (is.null(oracle)) {
        return(invisible(NULL))
    }
    reference <- oracle$response[[1L]]
    valid <- all(vapply(oracle$response, function(value) {
        .association_numeric_close(value, reference)
    }, logical(1L)))
    if (!valid) {
        .abort_association_evidence(
            "Association null scenario violates its zero plug-in oracle."
        )
    }
    invisible(oracle)
}

.association_ols_oracle_effect <- function(stratum, response, coefficient) {
    samples <- rownames(stratum$core$model$matrix)
    selected <- match(samples, names(response))
    column <- match(coefficient, colnames(stratum$core$model$matrix))
    if (anyNA(selected) || is.na(column)) {
        return(NA_real_)
    }
    x <- stratum$core$model$matrix
    y <- unname(response[selected])
    decomposition <- tryCatch(
        .design_svd(x),
        error = function(error) NULL
    )
    rank <- stratum$core$rank$rank
    if (is.null(decomposition) || decomposition$rank != rank || rank < 1L ||
        any(!is.finite(y))) {
        return(NA_real_)
    }
    u <- decomposition$decomposition$u[, seq_len(rank), drop = FALSE]
    v <- decomposition$decomposition$v[, seq_len(rank), drop = FALSE]
    singular <- decomposition$decomposition$d[seq_len(rank)]
    coordinates <- as.vector(crossprod(u, y)) / singular
    coefficient_values <- as.vector(v %*% coordinates)
    effect <- coefficient_values[[column]]
    if (is.finite(effect)) as.double(effect) else NA_real_
}

.association_quasibinomial_oracle_effect <- function(
    stratum,
    response,
    total,
    hypothesis
) {
    x <- stratum$core$model$matrix
    samples <- rownames(x)
    selected <- match(samples, names(response))
    column <- match(hypothesis$coefficient, colnames(x))
    if (anyNA(selected) || is.na(column)) {
        return(NA_real_)
    }
    y <- unname(response[selected])
    decomposition <- tryCatch(
        .design_svd(x),
        error = function(error) NULL
    )
    rank <- stratum$core$rank$rank
    if (is.null(decomposition) || decomposition$rank != rank || rank < 1L ||
        any(!is.finite(y)) || any(y < 0 | y > 1)) {
        return(NA_real_)
    }
    u <- decomposition$decomposition$u[, seq_len(rank), drop = FALSE]
    basis <- decomposition$decomposition$v[, seq_len(rank), drop = FALSE]
    singular <- decomposition$decomposition$d[seq_len(rank)]
    z <- sweep(u, 2L, singular, `*`)
    fit <- tryCatch(
        .association_quasibinomial_glm_fit(
            z,
            as.double(total * y),
            as.double(total)
        ),
        error = function(error) NULL
    )
    valid <- !is.null(fit) && isTRUE(fit$converged) &&
        identical(fit$boundary, FALSE) &&
        length(fit$coefficients) == rank &&
        all(is.finite(fit$coefficients)) &&
        length(fit$fitted.values) == nrow(x) &&
        all(is.finite(fit$fitted.values))
    if (!valid) {
        return(NA_real_)
    }
    probability_floor <- sqrt(.Machine$double.eps)
    if (any(fit$fitted.values < probability_floor |
        fit$fitted.values > 1 - probability_floor)) {
        return(NA_real_)
    }
    weights <- as.double(total * fit$fitted.values *
        (1 - fit$fitted.values))
    weighted <- tryCatch(
        .design_svd(sweep(z, 1L, sqrt(weights), `*`)),
        error = function(error) NULL
    )
    if (is.null(weighted) || weighted$rank != rank) {
        return(NA_real_)
    }
    axis <- numeric(ncol(x))
    axis[[column]] <- 1
    contrast <- as.vector(crossprod(basis, axis))
    observed_column <- x[, column]
    z0 <- z - tcrossprod(observed_column, contrast)
    z1 <- sweep(z0, 2L, contrast, `+`)
    mu0 <- stats::plogis(as.vector(z0 %*% fit$coefficients))
    mu1 <- stats::plogis(as.vector(z1 %*% fit$coefficients))
    effect <- mean(mu1 - mu0)
    if (is.finite(effect)) as.double(effect) else NA_real_
}

.association_target_truth <- function(
    candidate_id,
    preparation,
    oracle,
    target_label
) {
    if (is.null(oracle) || !inherits(
        preparation,
        "imputefinder_association_preparation"
    )) {
        return(NA_real_)
    }
    target <- .association_target_hypothesis(preparation, target_label)
    if (is.null(target)) {
        return(NA_real_)
    }
    if (identical(candidate_id, .ASSOCIATION_QUASIBINOMIAL_CANDIDATE)) {
        .association_quasibinomial_oracle_effect(
            target$stratum,
            oracle$response,
            oracle$count,
            target$hypothesis
        )
    } else {
        .association_ols_oracle_effect(
            target$stratum,
            oracle$response,
            target$hypothesis$coefficient[[1L]]
        )
    }
}

.association_common_opportunities <- function(preparation) {
    if (!inherits(
        preparation,
        "imputefinder_association_preparation"
    )) {
        return(integer())
    }
    residual_ok <- vapply(
        preparation$hypotheses$stratum,
        function(stratum_id) {
            stratum <- preparation$strata[[stratum_id]]
            length(stratum$samples) - stratum$core$rank$rank >= 3L
        },
        logical(1L)
    )
    which(
        preparation$hypotheses$eligible &
            preparation$hypotheses$estimable &
            preparation$support$eligible & residual_ok
    )
}

.association_artifact_outcome_bindings <- function(
    artifact,
    candidate_id,
    scenario_id,
    replicate
) {
    available <- vapply(artifact$outcomes, function(outcome) {
        identical(class(outcome), "imputefinder_association")
    }, logical(1L))
    code <- vapply(seq_along(artifact$outcomes), function(index) {
        if (available[[index]]) NA_character_ else
            artifact$outcomes[[index]]$code
    }, character(1L))
    data.frame(
        candidate_id = rep(candidate_id, length(artifact$outcomes)),
        scenario_id = rep(scenario_id, length(artifact$outcomes)),
        replicate = rep(replicate, length(artifact$outcomes)),
        stratum = artifact$hypotheses$stratum,
        hypothesis = artifact$hypotheses$hypothesis,
        status = ifelse(available, "available", "unavailable"),
        code = code,
        stringsAsFactors = FALSE
    )[.ASSOCIATION_OUTCOME_BINDING_FIELDS]
}

.association_panel_abstention_binding <- function(
    result,
    candidate_id,
    scenario_id,
    replicate
) {
    data.frame(
        candidate_id = candidate_id,
        scenario_id = scenario_id,
        replicate = replicate,
        stratum = NA_character_,
        hypothesis = .ASSOCIATION_PANEL_ABSTENTION,
        status = "unavailable",
        code = result$code,
        stringsAsFactors = FALSE
    )[.ASSOCIATION_OUTCOME_BINDING_FIELDS]
}

.association_execution_failure_bindings <- function(
    preparation,
    candidate_id,
    scenario_id,
    replicate
) {
    if (!inherits(
        preparation,
        "imputefinder_association_preparation"
    )) {
        .abort_association_evidence(
            "Association execution failure lacks a hypothesis disposition set."
        )
    }
    data.frame(
        candidate_id = rep(candidate_id, nrow(preparation$hypotheses)),
        scenario_id = rep(scenario_id, nrow(preparation$hypotheses)),
        replicate = rep(replicate, nrow(preparation$hypotheses)),
        stratum = preparation$hypotheses$stratum,
        hypothesis = preparation$hypotheses$hypothesis,
        status = rep("unavailable", nrow(preparation$hypotheses)),
        code = rep(
            "association_numerical_failure",
            nrow(preparation$hypotheses)
        ),
        stringsAsFactors = FALSE
    )[.ASSOCIATION_OUTCOME_BINDING_FIELDS]
}

.association_opportunity_rows <- function(
    preparation,
    artifact,
    candidate_id,
    scenario_id,
    replicate
) {
    selected <- .association_common_opportunities(preparation)
    if (!length(selected)) {
        return(.empty_association_opportunity_audit())
    }
    available <- rep(FALSE, length(selected))
    if (!is.null(artifact)) {
        available <- vapply(artifact$outcomes[selected], function(outcome) {
            identical(class(outcome), "imputefinder_association")
        }, logical(1L))
    }
    data.frame(
        candidate_id = rep(candidate_id, length(selected)),
        scenario_id = rep(scenario_id, length(selected)),
        replicate = rep(replicate, length(selected)),
        stratum = preparation$hypotheses$stratum[selected],
        hypothesis = preparation$hypotheses$hypothesis[selected],
        available = unname(available),
        stringsAsFactors = FALSE
    )[.ASSOCIATION_OPPORTUNITY_AUDIT_FIELDS]
}

.association_null_row <- function(
    artifact,
    candidate_id,
    scenario_id,
    replicate,
    acquisition
) {
    available <- if (is.null(artifact)) {
        logical()
    } else {
        vapply(artifact$outcomes, function(outcome) {
            identical(class(outcome), "imputefinder_association")
        }, logical(1L))
    }
    flagged <- if (any(available)) {
        any(vapply(
            artifact$outcomes[available],
            `[[`,
            logical(1L),
            "flag"
        ))
    } else {
        FALSE
    }
    data.frame(
        candidate_id = candidate_id,
        scenario_id = scenario_id,
        replicate = replicate,
        acquisition = acquisition,
        available_p_count = as.integer(sum(available)),
        false_flag = flagged,
        stringsAsFactors = FALSE
    )[.ASSOCIATION_NULL_AUDIT_FIELDS]
}

.association_target_row <- function(
    artifact,
    preparation,
    oracle,
    candidate_id,
    scenario_id,
    replicate,
    acquisition,
    target_label,
    fallback_code
) {
    truth <- .association_target_truth(
        candidate_id,
        preparation,
        oracle,
        target_label
    )
    target <- if (is.null(artifact)) integer() else
        which(artifact$hypotheses$label == target_label)
    outcome <- if (length(target) == 1L) {
        artifact$outcomes[[target[[1L]]]]
    } else {
        NULL
    }
    available <- !is.null(outcome) &&
        identical(class(outcome), "imputefinder_association")
    code <- if (available) {
        NA_character_
    } else if (!is.null(outcome) && inherits(
        outcome,
        "imputefinder_unavailable"
    )) {
        outcome$code
    } else {
        fallback_code
    }
    values <- if (available) {
        c(
            effect = outcome$effect,
            conf_low = outcome$conf_low,
            conf_high = outcome$conf_high,
            adjusted_p = outcome$adjusted_p
        )
    } else {
        c(
            effect = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
            adjusted_p = NA_real_
        )
    }
    data.frame(
        candidate_id = candidate_id,
        scenario_id = scenario_id,
        replicate = replicate,
        acquisition = acquisition,
        target_label = target_label,
        outcome_status = if (available) "available" else "unavailable",
        code = code,
        truth_status = if (is.finite(truth) && truth != 0) {
            "available"
        } else {
            "failed"
        },
        truth = as.double(truth),
        effect = as.double(values[["effect"]]),
        conf_low = as.double(values[["conf_low"]]),
        conf_high = as.double(values[["conf_high"]]),
        adjusted_p = as.double(values[["adjusted_p"]]),
        stringsAsFactors = FALSE
    )[.ASSOCIATION_TARGET_AUDIT_FIELDS]
}

.association_resolve_run <- function(
    record,
    preparation,
    oracle,
    candidate_id,
    scenario_id,
    replicate,
    acquisition,
    target
) {
    result <- record$result
    artifact <- NULL
    if (inherits(
        result,
        "imputefinder_association_candidate_artifact"
    )) {
        if (!inherits(
            preparation,
            "imputefinder_association_preparation"
        ) || !identical(result$candidate, candidate_id)) {
            .abort_association_evidence(
                "Association candidate artifact is detached from its run."
            )
        }
        .validate_association_candidate_artifact(result, preparation)
        artifact <- result
        result_kind <- "panel"
        status <- "measured"
        failure_code <- NA_character_
        result_sha256 <- .association_serialized_sha256(result)
        outcomes <- .association_artifact_outcome_bindings(
            result,
            candidate_id,
            scenario_id,
            replicate
        )
    } else if (inherits(result, "imputefinder_unavailable")) {
        if (!identical(result, preparation)) {
            .abort_association_evidence(
                "Association panel abstention is detached from preparation."
            )
        }
        result_kind <- "unavailable"
        status <- "failed"
        failure_code <- result$code
        result_sha256 <- NA_character_
        outcomes <- .association_panel_abstention_binding(
            result,
            candidate_id,
            scenario_id,
            replicate
        )
    } else if (.valid_association_execution_failure(result)) {
        result_kind <- "execution_failure"
        status <- "failed"
        failure_code <- result$code
        result_sha256 <- NA_character_
        outcomes <- .association_execution_failure_bindings(
            preparation,
            candidate_id,
            scenario_id,
            replicate
        )
    } else {
        .abort_association_evidence(
            "Association study result has no exact permitted disposition."
        )
    }
    run <- data.frame(
        candidate_id = candidate_id,
        scenario_id = scenario_id,
        replicate = replicate,
        input_sha256 = record$input_sha256,
        result_sha256 = result_sha256,
        result_kind = result_kind,
        elapsed_seconds = record$elapsed_seconds,
        status = status,
        failure_code = failure_code,
        stringsAsFactors = FALSE
    )[.ASSOCIATION_RUN_BINDING_FIELDS]
    opportunity <- .association_opportunity_rows(
        preparation,
        artifact,
        candidate_id,
        scenario_id,
        replicate
    )
    null <- if (scenario_id %in% .ASSOCIATION_NULL_SCENARIOS) {
        .association_null_row(
            artifact,
            candidate_id,
            scenario_id,
            replicate,
            acquisition
        )
    } else {
        .empty_association_null_audit()
    }
    target_row <- if (!is.null(target)) {
        .association_target_row(
            artifact,
            preparation,
            oracle,
            candidate_id,
            scenario_id,
            replicate,
            acquisition,
            target$label[[1L]],
            if (identical(result_kind, "execution_failure")) {
                "association_numerical_failure"
            } else if (identical(result_kind, "unavailable")) {
                failure_code
            } else {
                "association_target_absent"
            }
        )
    } else {
        .empty_association_target_audit()
    }
    list(
        run = run,
        outcomes = outcomes,
        opportunity = opportunity,
        null = null,
        target = target_row
    )
}

.association_run_grid <- function() {
    rows <- lapply(.ASSOCIATION_CANDIDATES, function(candidate_id) {
        do.call(rbind, lapply(
            .ASSOCIATION_STUDY_SCENARIOS,
            function(scenario_id) {
                data.frame(
                    candidate_id = candidate_id,
                    scenario_id = scenario_id,
                    replicate = .ASSOCIATION_STUDY_REPLICATES,
                    stringsAsFactors = FALSE
                )
            }
        ))
    })
    output <- do.call(rbind, rows)
    row.names(output) <- NULL
    output
}

.association_null_grid <- function() {
    grid <- .association_run_grid()
    grid <- grid[
        grid$scenario_id %in% .ASSOCIATION_NULL_SCENARIOS,
        ,
        drop = FALSE
    ]
    manifest <- .ASSOCIATION_STUDY_MANIFEST
    grid$acquisition <- manifest$acquisition[
        match(grid$scenario_id, manifest$scenario_id)
    ]
    row.names(grid) <- NULL
    grid
}

.association_target_grid <- function() {
    rows <- lapply(.ASSOCIATION_CANDIDATES, function(candidate_id) {
        do.call(rbind, lapply(seq_len(nrow(.ASSOCIATION_TARGETS)), function(
            index
        ) {
            target <- .ASSOCIATION_TARGETS[index, , drop = FALSE]
            data.frame(
                candidate_id = candidate_id,
                scenario_id = target$scenario_id,
                replicate = .ASSOCIATION_STUDY_REPLICATES,
                acquisition = target$acquisition,
                target_label = target$label,
                stringsAsFactors = FALSE
            )
        }))
    })
    output <- do.call(rbind, rows)
    row.names(output) <- NULL
    output
}

.association_order_index <- function(frame, fields) {
    candidate <- if ("candidate_id" %in% fields) {
        match(frame$candidate_id, .ASSOCIATION_CANDIDATES)
    } else {
        NULL
    }
    scenario <- if ("scenario_id" %in% fields) {
        match(frame$scenario_id, .ASSOCIATION_STUDY_SCENARIOS)
    } else {
        NULL
    }
    target <- if ("target_label" %in% fields) {
        match(frame$target_label, .ASSOCIATION_TARGETS$label)
    } else {
        NULL
    }
    values <- list()
    for (field in fields) {
        values[[length(values) + 1L]] <- switch(
            field,
            candidate_id = candidate,
            scenario_id = scenario,
            target_label = target,
            frame[[field]]
        )
    }
    do.call(order, c(values, list(na.last = TRUE, method = "radix")))
}

.association_canonicalize_audit <- function(frame, fields) {
    if (!nrow(frame)) {
        row.names(frame) <- NULL
        return(frame)
    }
    frame <- frame[
        .association_order_index(frame, fields),
        ,
        drop = FALSE
    ]
    row.names(frame) <- NULL
    frame
}

.validate_association_run_bindings <- function(run_bindings) {
    grid <- .association_run_grid()
    valid <- is.data.frame(run_bindings) &&
        identical(names(run_bindings), .ASSOCIATION_RUN_BINDING_FIELDS) &&
        nrow(run_bindings) == nrow(grid) &&
        all(vapply(run_bindings[c(
            "candidate_id", "scenario_id", "input_sha256", "result_sha256",
            "result_kind", "status", "failure_code"
        )], is.character, logical(1L))) &&
        is.integer(run_bindings$replicate) &&
        is.double(run_bindings$elapsed_seconds) &&
        identical(
            run_bindings[c("candidate_id", "scenario_id", "replicate")],
            grid
        ) && all(vapply(
            run_bindings$input_sha256,
            .association_sha256,
            logical(1L)
        )) && all(is.finite(run_bindings$elapsed_seconds)) &&
        all(run_bindings$elapsed_seconds >= 0) &&
        all(run_bindings$result_kind %in% c(
            "panel", "unavailable", "execution_failure"
        )) && all(run_bindings$status %in% c("measured", "failed"))
    panel <- if (is.data.frame(run_bindings)) {
        run_bindings$result_kind == "panel"
    } else {
        logical()
    }
    valid <- valid && all(run_bindings$status[panel] == "measured") &&
        all(!is.na(run_bindings$result_sha256[panel])) &&
        all(vapply(
            run_bindings$result_sha256[panel],
            .association_sha256,
            logical(1L)
        )) && all(is.na(run_bindings$failure_code[panel])) &&
        all(run_bindings$status[!panel] == "failed") &&
        all(is.na(run_bindings$result_sha256[!panel])) &&
        all(!is.na(run_bindings$failure_code[!panel])) &&
        all(nzchar(run_bindings$failure_code[!panel]))
    execution_failure <- run_bindings$result_kind == "execution_failure"
    panel_unavailable <- run_bindings$result_kind == "unavailable"
    valid <- valid && all(
        run_bindings$failure_code[execution_failure] ==
            "association_numerical_failure"
    ) && all(run_bindings$failure_code[panel_unavailable] %in% c(
        "association_no_observable_features",
        "association_no_testable_hypotheses"
    ))
    if (!valid) {
        .abort_association_evidence(
            "Association run bindings are malformed or incomplete."
        )
    }
    invisible(run_bindings)
}

.validate_association_outcome_bindings <- function(outcome_bindings) {
    valid <- is.data.frame(outcome_bindings) &&
        identical(
            names(outcome_bindings),
            .ASSOCIATION_OUTCOME_BINDING_FIELDS
    ) &&
        all(vapply(outcome_bindings[c(
            "candidate_id", "scenario_id", "stratum", "hypothesis",
            "status", "code"
        )], is.character, logical(1L))) &&
        is.integer(outcome_bindings$replicate) &&
        all(outcome_bindings$candidate_id %in% .ASSOCIATION_CANDIDATES) &&
        all(outcome_bindings$scenario_id %in%
            .ASSOCIATION_STUDY_SCENARIOS) &&
        all(outcome_bindings$replicate %in%
            .ASSOCIATION_STUDY_REPLICATES) &&
        all(outcome_bindings$status %in% c("available", "unavailable")) &&
        !anyDuplicated(outcome_bindings[c(
            "candidate_id", "scenario_id", "replicate", "stratum",
            "hypothesis"
        )])
    abstention <- if (is.data.frame(outcome_bindings)) {
        outcome_bindings$hypothesis == .ASSOCIATION_PANEL_ABSTENTION
    } else {
        logical()
    }
    allowed_codes <- c(
        .ASSOCIATION_UNAVAILABLE_CODES,
        "association_no_observable_features",
        "association_no_testable_hypotheses"
    )
    valid <- valid && all(is.na(outcome_bindings$code[
        outcome_bindings$status == "available"
    ])) && all(
        !is.na(outcome_bindings$code[
            outcome_bindings$status == "unavailable"
        ])
    ) && all(outcome_bindings$code[
        outcome_bindings$status == "unavailable"
    ] %in% allowed_codes) && all(is.na(outcome_bindings$stratum[abstention])) &&
        all(outcome_bindings$status[abstention] == "unavailable") &&
        all(!is.na(outcome_bindings$stratum[!abstention])) &&
        all(grepl("^a_[0-9a-f]{64}$", outcome_bindings$hypothesis[!abstention]))
    canonical <- if (valid) {
        .association_canonicalize_audit(
            outcome_bindings,
            c(
                "candidate_id", "scenario_id", "replicate", "stratum",
                "hypothesis"
            )
        )
    } else {
        NULL
    }
    if (!valid || !identical(outcome_bindings, canonical)) {
        .abort_association_evidence(
            "Association outcome bindings are malformed or noncanonical."
        )
    }
    invisible(outcome_bindings)
}

.validate_association_opportunity_audit <- function(opportunity) {
    valid <- is.data.frame(opportunity) && identical(
        names(opportunity),
        .ASSOCIATION_OPPORTUNITY_AUDIT_FIELDS
    ) &&
        all(vapply(opportunity[c(
            "candidate_id", "scenario_id", "stratum", "hypothesis"
        )], is.character, logical(1L))) &&
        is.integer(opportunity$replicate) &&
        is.logical(opportunity$available) && !anyNA(opportunity) &&
        all(opportunity$candidate_id %in% .ASSOCIATION_CANDIDATES) &&
        all(opportunity$scenario_id %in% .ASSOCIATION_STUDY_SCENARIOS) &&
        all(opportunity$replicate %in% .ASSOCIATION_STUDY_REPLICATES) &&
        all(nzchar(opportunity$stratum)) &&
        all(grepl("^a_[0-9a-f]{64}$", opportunity$hypothesis)) &&
        !anyDuplicated(opportunity[c(
            "candidate_id", "scenario_id", "replicate", "stratum",
            "hypothesis"
        )])
    canonical <- if (valid) {
        .association_canonicalize_audit(
            opportunity,
            c(
                "candidate_id", "scenario_id", "replicate", "stratum",
                "hypothesis"
            )
        )
    } else {
        NULL
    }
    if (valid) {
        keys <- lapply(.ASSOCIATION_CANDIDATES, function(candidate_id) {
            selected <- opportunity$candidate_id == candidate_id
            output <- opportunity[
                selected,
                c(
                    "scenario_id", "replicate", "stratum", "hypothesis"
                ),
                drop = FALSE
            ]
            row.names(output) <- NULL
            output
        })
        common <- all(vapply(keys[-1L], identical, logical(1L), keys[[1L]]))
    } else {
        common <- FALSE
    }
    if (!valid || !identical(opportunity, canonical) || !common) {
        .abort_association_evidence(
            paste0(
                "Association opportunity audit lacks one canonical common ",
                "candidate-agnostic denominator."
            )
        )
    }
    invisible(opportunity)
}

.validate_association_null_audit <- function(null) {
    grid <- .association_null_grid()
    valid <- is.data.frame(null) && identical(
        names(null),
        .ASSOCIATION_NULL_AUDIT_FIELDS
    ) && nrow(null) == nrow(grid) &&
        all(vapply(null[c(
            "candidate_id", "scenario_id", "acquisition"
        )], is.character, logical(1L))) &&
        is.integer(null$replicate) && is.integer(null$available_p_count) &&
        is.logical(null$false_flag) && !anyNA(null) &&
        identical(
            null[c(
                "candidate_id", "scenario_id", "replicate", "acquisition"
            )],
            grid
        ) && all(null$available_p_count >= 0L) &&
        all(!null$false_flag | null$available_p_count > 0L)
    if (!valid) {
        .abort_association_evidence(
            "Association null-family audit is malformed or incomplete."
        )
    }
    invisible(null)
}

.validate_association_target_audit <- function(target) {
    grid <- .association_target_grid()
    valid <- is.data.frame(target) && identical(
        names(target),
        .ASSOCIATION_TARGET_AUDIT_FIELDS
    ) && nrow(target) == nrow(grid) &&
        all(vapply(target[c(
            "candidate_id", "scenario_id", "acquisition", "target_label",
            "outcome_status", "code", "truth_status"
        )], is.character, logical(1L))) &&
        is.integer(target$replicate) && all(vapply(target[c(
            "truth", "effect", "conf_low", "conf_high", "adjusted_p"
        )], is.double, logical(1L))) && identical(
            target[c(
                "candidate_id", "scenario_id", "replicate", "acquisition",
                "target_label"
            )],
            grid
        ) && all(target$outcome_status %in% c("available", "unavailable")) &&
        all(target$truth_status %in% c("available", "failed"))
    outcome_available <- if (is.data.frame(target)) {
        target$outcome_status == "available"
    } else {
        logical()
    }
    truth_available <- if (is.data.frame(target)) {
        target$truth_status == "available"
    } else {
        logical()
    }
    outcome_values <- c("effect", "conf_low", "conf_high", "adjusted_p")
    valid <- valid && all(is.na(target$code[outcome_available])) &&
        all(!is.na(target$code[!outcome_available])) &&
        all(nzchar(target$code[!outcome_available])) &&
        all(target$code[!outcome_available] %in% c(
            .ASSOCIATION_UNAVAILABLE_CODES,
            "association_no_observable_features",
            "association_no_testable_hypotheses",
            "association_target_absent"
        )) &&
        all(is.finite(as.matrix(target[outcome_available, outcome_values]))) &&
        all(is.na(as.matrix(target[!outcome_available, outcome_values]))) &&
        all(target$conf_low[outcome_available] <=
            target$conf_high[outcome_available]) &&
        all(target$adjusted_p[outcome_available] >= 0 &
            target$adjusted_p[outcome_available] <= 1) &&
        all(is.finite(target$truth[truth_available])) &&
        all(target$truth[truth_available] != 0) &&
        all(is.na(target$truth[!truth_available]))
    quasi <- target$candidate_id == .ASSOCIATION_QUASIBINOMIAL_CANDIDATE
    valid <- valid && all(abs(target$truth[quasi & truth_available]) <= 1) &&
        all(abs(target$effect[quasi & outcome_available]) <= 1)
    if (!valid) {
        .abort_association_evidence(
            "Association fixed-target audit is malformed or off scale."
        )
    }
    invisible(target)
}

.association_compare_threshold <- function(value, operator, threshold) {
    switch(
        operator,
        "<=" = value <= threshold,
        ">=" = value >= threshold,
        "==" = value == threshold,
        "<" = value < threshold,
        ">" = value > threshold,
        FALSE
    )
}

.association_clopper_pearson_upper <- function(numerator, denominator) {
    valid <- length(numerator) == 1L && length(denominator) == 1L &&
        is.finite(numerator) && is.finite(denominator) && denominator > 0 &&
        numerator >= 0 && numerator <= denominator &&
        numerator == floor(numerator) && denominator == floor(denominator)
    if (!valid) {
        .abort_association_evidence(
            "Clopper-Pearson counts are malformed."
        )
    }
    if (numerator == denominator) {
        return(1.0)
    }
    as.double(stats::qbeta(
        0.95,
        numerator + 1,
        denominator - numerator
    ))
}

.association_exact_binomial_interval <- function(numerator, denominator) {
    interval <- stats::binom.test(
        as.integer(numerator),
        as.integer(denominator),
        conf.level = 0.95
    )$conf.int
    as.double(interval)
}

.association_screening_scope <- function(candidate_id, metric) {
    if (identical(metric, "common_opportunity_coverage")) {
        return("all_candidate_agnostic_opportunities")
    }
    if (metric %in% names(.ASSOCIATION_SCREENING_OPERATORS)[1:2]) {
        return(paste0(
            paste(.ASSOCIATION_NULL_SCENARIOS, collapse = ";"),
            ":replicates_1-32"
        ))
    }
    scenarios <- .ASSOCIATION_TARGETS$scenario_id
    if (identical(candidate_id, .ASSOCIATION_QUASIBINOMIAL_CANDIDATE)) {
        scenarios <- setdiff(scenarios, "dia_monotone_paired")
    }
    paste0(paste(scenarios, collapse = ";"), ":replicates_1-32")
}

.association_screening_metric_row <- function(
    candidate_id,
    metric,
    numerator = NA_real_,
    denominator = NA_real_,
    estimate = NA_real_,
    status = "failed",
    passed = FALSE
) {
    data.frame(
        candidate_id = candidate_id,
        metric = metric,
        numerator = as.double(numerator),
        denominator = as.double(denominator),
        estimate = as.double(estimate),
        licensed_scope = .association_screening_scope(
            candidate_id,
            metric
        ),
        status = status,
        passed = passed,
        stringsAsFactors = FALSE
    )[.ASSOCIATION_SCREENING_METRIC_FIELDS]
}

.association_screening_measured_row <- function(
    candidate_id,
    metric,
    numerator,
    denominator,
    estimate
) {
    operator <- .ASSOCIATION_SCREENING_OPERATORS[[metric]]
    threshold <- .ASSOCIATION_SCREENING_THRESHOLDS[[metric]]
    .association_screening_metric_row(
        candidate_id,
        metric,
        numerator,
        denominator,
        estimate,
        "measured",
        .association_compare_threshold(estimate, operator, threshold)
    )
}

.association_candidate_screening <- function(
    candidate_id,
    run_bindings,
    opportunity,
    null,
    target
) {
    candidate_runs <- run_bindings$candidate_id == candidate_id
    candidate_opportunity <- opportunity$candidate_id == candidate_id
    denominator <- sum(candidate_opportunity)
    coverage <- if (denominator > 0L) {
        .association_screening_metric_row(
            candidate_id,
            "common_opportunity_coverage",
            sum(opportunity$available[candidate_opportunity]),
            denominator,
            sum(opportunity$available[candidate_opportunity]) / denominator,
            "measured",
            TRUE
        )
    } else {
        .association_screening_metric_row(
            candidate_id,
            "common_opportunity_coverage"
        )
    }
    failure <- any(
        run_bindings$result_kind[candidate_runs] == "execution_failure"
    )
    candidate_null <- null[null$candidate_id == candidate_id, , drop = FALSE]
    null_complete <- all(candidate_null$available_p_count > 0L)
    if (!failure && null_complete) {
        false_count <- sum(candidate_null$false_flag)
        null_upper <- .association_screening_measured_row(
            candidate_id,
            "family_false_flag_clopper_pearson_upper95",
            false_count,
            nrow(candidate_null),
            .association_clopper_pearson_upper(
                false_count,
                nrow(candidate_null)
            )
        )
        by_acquisition <- split(
            candidate_null$false_flag,
            candidate_null$acquisition
        )
        fractions <- vapply(by_acquisition, mean, numeric(1L))
        maximum <- max(fractions)
        maximum_count <- max(vapply(by_acquisition, sum, integer(1L)))
        stratum <- .association_screening_measured_row(
            candidate_id,
            "maximum_acquisition_false_flag_fraction",
            maximum_count,
            unique(lengths(by_acquisition)),
            maximum
        )
    } else {
        null_upper <- .association_screening_metric_row(
            candidate_id,
            "family_false_flag_clopper_pearson_upper95"
        )
        stratum <- .association_screening_metric_row(
            candidate_id,
            "maximum_acquisition_false_flag_fraction"
        )
    }
    licensed_scenarios <- .ASSOCIATION_TARGETS$scenario_id
    if (identical(candidate_id, .ASSOCIATION_QUASIBINOMIAL_CANDIDATE)) {
        licensed_scenarios <- setdiff(
            licensed_scenarios,
            "dia_monotone_paired"
        )
    }
    candidate_target <- target[
        target$candidate_id == candidate_id &
            target$scenario_id %in% licensed_scenarios,
        ,
        drop = FALSE
    ]
    truth_complete <- all(candidate_target$truth_status == "available")
    available <- candidate_target$outcome_status == "available"
    if (!failure && truth_complete) {
        covered <- available & candidate_target$conf_low <=
            candidate_target$truth & candidate_target$truth <=
            candidate_target$conf_high
        detected <- available & candidate_target$adjusted_p <= 0.05
        interval <- .association_screening_measured_row(
            candidate_id,
            "candidate_specific_projection_95_interval_coverage",
            sum(covered),
            nrow(candidate_target),
            mean(covered)
        )
        power <- .association_screening_measured_row(
            candidate_id,
            "holm_target_detection_fraction",
            sum(detected),
            nrow(candidate_target),
            mean(detected)
        )
    } else {
        interval <- .association_screening_metric_row(
            candidate_id,
            "candidate_specific_projection_95_interval_coverage"
        )
        power <- .association_screening_metric_row(
            candidate_id,
            "holm_target_detection_fraction"
        )
    }
    if (!failure && truth_complete && all(available)) {
        absolute <- abs(candidate_target$effect - candidate_target$truth)
        bias <- .association_screening_measured_row(
            candidate_id,
            "median_absolute_candidate_specific_projection_bias",
            sum(absolute),
            length(absolute),
            stats::median(absolute)
        )
    } else {
        bias <- .association_screening_metric_row(
            candidate_id,
            "median_absolute_candidate_specific_projection_bias"
        )
    }
    do.call(rbind, list(
        coverage, null_upper, stratum, interval, power, bias
    ))
}

.association_screening_metrics <- function(
    run_bindings,
    opportunity,
    null,
    target
) {
    rows <- lapply(.ASSOCIATION_CANDIDATES, function(candidate_id) {
        .association_candidate_screening(
            candidate_id,
            run_bindings,
            opportunity,
            null,
            target
        )
    })
    output <- do.call(rbind, rows)
    row.names(output) <- NULL
    output[.ASSOCIATION_SCREENING_METRIC_FIELDS]
}

.validate_association_screening_metrics <- function(screening) {
    expected <- do.call(rbind, lapply(.ASSOCIATION_CANDIDATES, function(
        candidate_id
    ) {
        data.frame(
            candidate_id = candidate_id,
            metric = .ASSOCIATION_SCREENING_METRICS,
            stringsAsFactors = FALSE
        )
    }))
    row.names(expected) <- NULL
    valid <- is.data.frame(screening) && identical(
        names(screening),
        .ASSOCIATION_SCREENING_METRIC_FIELDS
    ) && nrow(screening) == nrow(expected) && identical(
        screening[c("candidate_id", "metric")],
        expected
    ) && all(vapply(screening[c(
        "numerator", "denominator", "estimate"
    )], is.double, logical(1L))) && is.logical(screening$passed) &&
        !anyNA(screening[c(
            "candidate_id", "metric", "licensed_scope", "status", "passed"
        )]) && all(screening$status %in% c("measured", "failed")) &&
        identical(
            screening$licensed_scope,
            unname(mapply(
                .association_screening_scope,
                screening$candidate_id,
                screening$metric,
                USE.NAMES = FALSE
            ))
        )
    measured <- if (is.data.frame(screening)) {
        screening$status == "measured"
    } else {
        logical()
    }
    numeric_fields <- c("numerator", "denominator", "estimate")
    valid <- valid && all(is.finite(as.matrix(screening[
        measured,
        numeric_fields,
        drop = FALSE
    ]))) && all(screening$numerator[measured] >= 0) &&
        all(screening$denominator[measured] > 0) &&
        all(is.na(as.matrix(screening[
            !measured,
            numeric_fields,
            drop = FALSE
        ]))) && all(!screening$passed[!measured])
    fractions <- measured & screening$metric !=
        "median_absolute_candidate_specific_projection_bias"
    valid <- valid && all(screening$estimate[fractions] >= 0 &
        screening$estimate[fractions] <= 1)
    expected_pass <- vapply(seq_len(nrow(screening)), function(index) {
        if (!measured[[index]]) {
            return(FALSE)
        }
        metric <- screening$metric[[index]]
        if (identical(metric, "common_opportunity_coverage")) {
            return(TRUE)
        }
        .association_compare_threshold(
            screening$estimate[[index]],
            .ASSOCIATION_SCREENING_OPERATORS[[metric]],
            .ASSOCIATION_SCREENING_THRESHOLDS[[metric]]
        )
    }, logical(1L))
    if (!valid || !identical(screening$passed, expected_pass)) {
        .abort_association_evidence(
            "Association screening metrics are malformed or off scale."
        )
    }
    invisible(screening)
}

.association_bootstrap_seed <- function(candidate_id, gate_id) {
    registry <- .association_study_gate_registry()
    valid <- .association_scalar_character(candidate_id) &&
        candidate_id %in% .ASSOCIATION_CANDIDATES &&
        .association_scalar_character(gate_id) &&
        gate_id %in% registry$gate_id[3:5]
    if (!valid) {
        .abort_association_evidence(
            "Association gate-bootstrap identity is malformed."
        )
    }
    digest <- unname(tools::sha256sum(bytes = charToRaw(enc2utf8(paste(
        .ASSOCIATION_PROTOCOL_ID,
        candidate_id,
        gate_id,
        "scenario_stratified_whole_replicate_bootstrap_v1",
        sep = "\t"
    )))))
    as.integer(strtoi(substr(digest, 1L, 7L), 16L) + 1L)
}

.association_bootstrap_interval <- function(
    values,
    scenarios,
    replicates,
    candidate_id,
    gate_id,
    statistic
) {
    expected_scenarios <- rep(
        .ASSOCIATION_TARGETS$scenario_id,
        each = length(.ASSOCIATION_STUDY_REPLICATES)
    )
    expected_replicates <- rep(
        .ASSOCIATION_STUDY_REPLICATES,
        times = nrow(.ASSOCIATION_TARGETS)
    )
    valid <- is.double(values) && length(values) == 128L &&
        all(is.finite(values)) && is.character(scenarios) &&
        identical(scenarios, expected_scenarios) &&
        is.integer(replicates) && identical(replicates, expected_replicates) &&
        .association_scalar_character(statistic) &&
        statistic %in% c("mean", "median")
    if (!valid) {
        .abort_association_evidence(
            "Association gate bootstrap input is malformed."
        )
    }
    groups <- lapply(.ASSOCIATION_TARGETS$scenario_id, function(scenario_id) {
        which(scenarios == scenario_id)
    })
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    on.exit(
        .association_restore_rng(old_kind, had_seed, old_seed),
        add = TRUE
    )
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(.association_bootstrap_seed(candidate_id, gate_id))
    estimates <- double(.ASSOCIATION_GATE_BOOTSTRAP_DRAWS)
    fun <- if (identical(statistic, "mean")) mean else stats::median
    for (draw in seq_len(.ASSOCIATION_GATE_BOOTSTRAP_DRAWS)) {
        selected <- unlist(lapply(groups, function(group) {
            group[sample.int(length(group), length(group), replace = TRUE)]
        }), use.names = FALSE)
        estimates[[draw]] <- fun(values[selected])
    }
    interval <- as.double(stats::quantile(
        estimates,
        c(0.025, 0.975),
        names = FALSE,
        type = 8L
    ))
    estimate <- fun(values)
    .association_expand_bootstrap_interval(interval, estimate)
}

.association_expand_bootstrap_interval <- function(interval, estimate) {
    valid <- is.double(interval) && length(interval) == 2L &&
        all(is.finite(interval)) && interval[[1L]] <= interval[[2L]] &&
        is.double(estimate) && length(estimate) == 1L &&
        is.finite(estimate)
    if (!valid) {
        .abort_association_evidence(
            "Association bootstrap interval expansion is malformed."
        )
    }
    c(min(interval[[1L]], estimate), max(interval[[2L]], estimate))
}

.association_gate_evidence_hash <- function(
    candidate_id,
    gate_id,
    evidence
) {
    if (!is.data.frame(evidence)) {
        .abort_association_evidence(
            "Association gate evidence must be one exact audit frame."
        )
    }
    descriptor <- data.frame(
        position = 1:3,
        field = c("candidate_id", "gate_id", "protocol_hash"),
        value = c(candidate_id, gate_id, .ASSOCIATION_PROTOCOL_HASH),
        type = "character",
        stringsAsFactors = FALSE
    )
    keyed <- data.frame(
        row_id = sprintf("row_%05d", seq_len(nrow(evidence))),
        evidence,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    .association_evidence_bundle_hash(list(
        descriptor = descriptor,
        evidence = keyed
    ))
}

.association_gate_raw_fields <- c(
    "gate_id", "registry_version", "metric", "estimate", "numerator",
    "denominator", "lower", "upper", "status", "gate_binding_hash",
    "evidence_hash", "note"
)

.association_gate_raw_row <- function(
    registry,
    evidence_hash,
    estimate = NA_real_,
    numerator = NA_real_,
    denominator = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    status = "failed",
    note
) {
    data.frame(
        gate_id = registry$gate_id,
        registry_version = registry$registry_version,
        metric = registry$metric,
        estimate = as.double(estimate),
        numerator = as.double(numerator),
        denominator = as.double(denominator),
        lower = as.double(lower),
        upper = as.double(upper),
        status = status,
        gate_binding_hash = registry$gate_binding_hash,
        evidence_hash = evidence_hash,
        note = note,
        stringsAsFactors = FALSE
    )[.association_gate_raw_fields]
}

.association_evaluate_gate_results <- function(results) {
    registry <- .association_study_gate_registry()
    valid <- is.data.frame(results) && identical(
        names(results),
        .association_gate_raw_fields
    ) && nrow(results) == nrow(registry) &&
        identical(results$gate_id, registry$gate_id) &&
        identical(results$registry_version, registry$registry_version) &&
        identical(results$metric, registry$metric) &&
        identical(results$gate_binding_hash, registry$gate_binding_hash) &&
        all(vapply(results[c(
            "estimate", "numerator", "denominator", "lower", "upper"
        )], is.double, logical(1L))) &&
        all(results$status %in% c("measured", "failed", "unavailable")) &&
        all(vapply(results$evidence_hash, .association_sha256, logical(1L))) &&
        is.character(results$note) && !anyNA(results$note) &&
        all(nzchar(results$note))
    measured <- if (is.data.frame(results)) {
        results$status == "measured"
    } else {
        logical()
    }
    numeric_fields <- c(
        "estimate", "numerator", "denominator", "lower", "upper"
    )
    valid <- valid && all(is.finite(as.matrix(
        results[measured, numeric_fields, drop = FALSE]
    ))) && all(results$denominator[measured] > 0) &&
        all(results$lower[measured] <= results$estimate[measured]) &&
        all(results$estimate[measured] <= results$upper[measured]) &&
        all(is.na(as.matrix(
            results[!measured, numeric_fields, drop = FALSE]
        )))
    if (!valid) {
        .abort_association_evidence(
            "Association selected-winner gate rows are malformed."
        )
    }
    passed <- measured & mapply(
        .association_compare_threshold,
        results$estimate,
        registry$operator,
        registry$threshold,
        USE.NAMES = FALSE
    )
    data.frame(
        gate_id = registry$gate_id,
        m12_gate_id = registry$m12_gate_id,
        stage = rep("selected_winner_synthetic", nrow(registry)),
        metric = registry$metric,
        estimate = results$estimate,
        operator = registry$operator,
        threshold = registry$threshold,
        status = results$status,
        passed = passed,
        failure_treatment = registry$failure_treatment,
        evidence_hash = results$evidence_hash,
        stringsAsFactors = FALSE
    )
}

.association_candidate_full_gates <- function(
    candidate_id,
    null,
    target
) {
    registry <- .association_study_gate_registry()
    candidate_null <- null[null$candidate_id == candidate_id, , drop = FALSE]
    candidate_target <- target[
        target$candidate_id == candidate_id,
        ,
        drop = FALSE
    ]
    null_hash <- vapply(registry$gate_id[1:2], function(gate_id) {
        .association_gate_evidence_hash(
            candidate_id,
            gate_id,
            candidate_null
        )
    }, character(1L))
    target_hash <- vapply(registry$gate_id[3:5], function(gate_id) {
        .association_gate_evidence_hash(
            candidate_id,
            gate_id,
            candidate_target
        )
    }, character(1L))
    null_complete <- all(candidate_null$available_p_count > 0L)
    if (null_complete) {
        false_count <- sum(candidate_null$false_flag)
        upper <- .association_clopper_pearson_upper(
            false_count,
            nrow(candidate_null)
        )
        null_upper <- .association_gate_raw_row(
            registry[1L, , drop = FALSE],
            null_hash[[1L]],
            upper,
            false_count,
            nrow(candidate_null),
            0,
            upper,
            "measured",
            "Pooled one-sided 95% Clopper-Pearson upper bound."
        )
        groups <- split(
            seq_len(nrow(candidate_null)),
            candidate_null$acquisition
        )
        counts <- vapply(groups, function(index) {
            sum(candidate_null$false_flag[index])
        }, integer(1L))
        denominators <- lengths(groups)
        fractions <- counts / denominators
        worst <- order(-fractions, names(fractions), method = "radix")[[1L]]
        interval <- .association_exact_binomial_interval(
            counts[[worst]],
            denominators[[worst]]
        )
        null_stratum <- .association_gate_raw_row(
            registry[2L, , drop = FALSE],
            null_hash[[2L]],
            fractions[[worst]],
            counts[[worst]],
            denominators[[worst]],
            interval[[1L]],
            interval[[2L]],
            "measured",
            paste0(
                "Maximum acquisition fraction; exact 95% interval for ",
                names(fractions)[[worst]], "."
            )
        )
    } else {
        note <- paste0(
            "Failed: at least one null family has zero available raw p-values."
        )
        null_upper <- .association_gate_raw_row(
            registry[1L, , drop = FALSE],
            null_hash[[1L]],
            note = note
        )
        null_stratum <- .association_gate_raw_row(
            registry[2L, , drop = FALSE],
            null_hash[[2L]],
            note = note
        )
    }
    complete <- candidate_target$outcome_status == "available" &
        candidate_target$truth_status == "available"
    target_complete <- all(complete)
    if (target_complete) {
        covered <- as.double(
            candidate_target$conf_low <= candidate_target$truth &
                candidate_target$truth <= candidate_target$conf_high
        )
        detected <- as.double(candidate_target$adjusted_p <= 0.05)
        absolute <- abs(candidate_target$effect - candidate_target$truth)
        coverage_estimate <- mean(covered)
        power_estimate <- mean(detected)
        bias_estimate <- stats::median(absolute)
        coverage_interval <- .association_bootstrap_interval(
            covered,
            candidate_target$scenario_id,
            candidate_target$replicate,
            candidate_id,
            registry$gate_id[[3L]],
            "mean"
        )
        power_interval <- .association_bootstrap_interval(
            detected,
            candidate_target$scenario_id,
            candidate_target$replicate,
            candidate_id,
            registry$gate_id[[4L]],
            "mean"
        )
        bias_interval <- .association_bootstrap_interval(
            as.double(absolute),
            candidate_target$scenario_id,
            candidate_target$replicate,
            candidate_id,
            registry$gate_id[[5L]],
            "median"
        )
        coverage <- .association_gate_raw_row(
            registry[3L, , drop = FALSE],
            target_hash[[1L]],
            coverage_estimate,
            sum(covered),
            length(covered),
            coverage_interval[[1L]],
            coverage_interval[[2L]],
            "measured",
            paste0(
                "9999-draw scenario-stratified whole-generator-replicate ",
                "bootstrap 95% interval."
            )
        )
        power <- .association_gate_raw_row(
            registry[4L, , drop = FALSE],
            target_hash[[2L]],
            power_estimate,
            sum(detected),
            length(detected),
            power_interval[[1L]],
            power_interval[[2L]],
            "measured",
            paste0(
                "9999-draw scenario-stratified whole-generator-replicate ",
                "bootstrap 95% interval."
            )
        )
        bias <- .association_gate_raw_row(
            registry[5L, , drop = FALSE],
            target_hash[[3L]],
            bias_estimate,
            sum(absolute),
            length(absolute),
            bias_interval[[1L]],
            bias_interval[[2L]],
            "measured",
            paste0(
                "9999-draw scenario-stratified whole-generator-replicate ",
                "percentile 95% interval."
            )
        )
    } else {
        note <- paste0(
            "Failed: all 128 fixed targets require finite truth and outcome."
        )
        coverage <- .association_gate_raw_row(
            registry[3L, , drop = FALSE],
            target_hash[[1L]],
            note = note
        )
        power <- .association_gate_raw_row(
            registry[4L, , drop = FALSE],
            target_hash[[2L]],
            note = note
        )
        bias <- .association_gate_raw_row(
            registry[5L, , drop = FALSE],
            target_hash[[3L]],
            note = note
        )
    }
    raw <- do.call(rbind, list(
        null_upper, null_stratum, coverage, power, bias
    ))
    row.names(raw) <- NULL
    output <- .association_evaluate_gate_results(raw)
    output$candidate_id <- candidate_id
    output$numerator <- raw$numerator
    output$denominator <- raw$denominator
    output$lower <- raw$lower
    output$upper <- raw$upper
    output$note <- raw$note
    output <- output[c(
        "candidate_id",
        setdiff(names(output), "candidate_id")
    )]
    output[.ASSOCIATION_FULL_GATE_RESULT_FIELDS]
}

.association_ranking_near <- function(value, anchor) {
    valid <- is.double(value) && is.double(anchor) &&
        length(anchor) == 1L && is.finite(anchor) &&
        length(value) > 0L && all(is.finite(value))
    if (!valid) {
        .abort_association_evidence(
            "Association ranking comparison is malformed."
        )
    }
    slack <- .ASSOCIATION_RANKING_ULP_MULTIPLIER *
        .Machine$double.eps * pmax(
            1,
            abs(value),
            abs(anchor),
            .ASSOCIATION_RANKING_NEAR_TIE
        )
    abs(value - anchor) <= .ASSOCIATION_RANKING_NEAR_TIE + slack
}

.association_rank_candidates <- function(metrics) {
    valid <- is.data.frame(metrics) && identical(
        names(metrics),
        .ASSOCIATION_RANKING_METRIC_FIELDS
    ) && nrow(metrics) > 0L && !anyDuplicated(metrics$candidate_id) &&
        all(metrics$candidate_id %in% .ASSOCIATION_CANDIDATES) &&
        all(vapply(metrics[-1L], is.double, logical(1L))) && !anyNA(metrics) &&
        all(is.finite(as.matrix(metrics[-1L]))) &&
        all(metrics$scope_coverage >= 0 & metrics$scope_coverage <= 1) &&
        all(metrics$power >= 0 & metrics$power <= 1) &&
        all(metrics$absolute_bias >= 0) &&
        all(metrics$median_runtime >= 0)
    if (!valid) {
        .abort_association_evidence(
            "Association ranking metrics are malformed."
        )
    }
    simplicity <- stats::setNames(
        seq_along(.ASSOCIATION_CANDIDATES),
        .ASSOCIATION_CANDIDATES
    )
    metrics$simplicity_rank <- unname(simplicity[metrics$candidate_id])
    output <- character()
    remaining <- metrics
    while (nrow(remaining)) {
        anchor_order <- order(
            -remaining$scope_coverage,
            remaining$candidate_id,
            method = "radix"
        )
        anchor <- remaining[anchor_order[[1L]], , drop = FALSE]
        near <- .association_ranking_near(
            remaining$scope_coverage,
            anchor$scope_coverage
        ) & .association_ranking_near(
            remaining$power,
            anchor$power
        ) & .association_ranking_near(
            remaining$absolute_bias,
            anchor$absolute_bias
        )
        cluster <- remaining[near, , drop = FALSE]
        chosen_order <- order(
            cluster$simplicity_rank,
            cluster$median_runtime,
            cluster$candidate_id,
            method = "radix"
        )
        chosen <- cluster$candidate_id[[chosen_order[[1L]]]]
        output <- c(output, chosen)
        remaining <- remaining[
            remaining$candidate_id != chosen,
            ,
            drop = FALSE
        ]
    }
    output
}

.validate_association_full_gate_results <- function(full_gates) {
    registry <- .association_study_gate_registry()
    valid <- is.data.frame(full_gates) && identical(
        names(full_gates),
        .ASSOCIATION_FULL_GATE_RESULT_FIELDS
    ) && nrow(full_gates) %% nrow(registry) == 0L &&
        all(vapply(full_gates[c(
            "candidate_id", "gate_id", "m12_gate_id", "stage", "metric",
            "operator", "status", "failure_treatment", "evidence_hash",
            "note"
        )], is.character, logical(1L))) &&
        all(vapply(full_gates[c(
            "estimate", "numerator", "denominator", "lower", "upper",
            "threshold"
        )], is.double, logical(1L))) && is.logical(full_gates$passed) &&
        !anyNA(full_gates[c(
            "candidate_id", "gate_id", "m12_gate_id", "stage", "metric",
            "operator", "threshold", "status", "passed",
            "failure_treatment", "evidence_hash", "note"
        )]) && all(full_gates$candidate_id %in% .ASSOCIATION_CANDIDATES) &&
        !anyDuplicated(full_gates[c("candidate_id", "gate_id")])
    candidates <- if (is.data.frame(full_gates)) {
        unique(full_gates$candidate_id)
    } else {
        character()
    }
    valid <- valid && identical(
        candidates,
        .ASSOCIATION_CANDIDATES[
            .ASSOCIATION_CANDIDATES %in% candidates
        ]
    )
    if (valid && length(candidates)) {
        expected <- do.call(rbind, lapply(candidates, function(candidate_id) {
            data.frame(
                candidate_id = candidate_id,
                gate_id = registry$gate_id,
                m12_gate_id = registry$m12_gate_id,
                stage = "selected_winner_synthetic",
                metric = registry$metric,
                operator = registry$operator,
                threshold = registry$threshold,
                failure_treatment = registry$failure_treatment,
                stringsAsFactors = FALSE
            )
        }))
        row.names(expected) <- NULL
        comparison <- c(
            "candidate_id", "gate_id", "m12_gate_id", "stage", "metric",
            "operator", "threshold", "failure_treatment"
        )
        valid <- identical(full_gates[comparison], expected[comparison])
    }
    measured <- if (is.data.frame(full_gates)) {
        full_gates$status == "measured"
    } else {
        logical()
    }
    raw_numeric <- c(
        "estimate", "numerator", "denominator", "lower", "upper"
    )
    valid <- valid && all(full_gates$status %in% c(
        "measured", "failed", "unavailable"
    )) && all(is.finite(as.matrix(full_gates[
        measured,
        raw_numeric,
        drop = FALSE
    ]))) && all(is.na(as.matrix(full_gates[
        !measured,
        raw_numeric,
        drop = FALSE
    ]))) && all(full_gates$numerator[measured] >= 0) &&
        all(full_gates$denominator[measured] > 0) &&
        all(full_gates$lower[measured] <= full_gates$estimate[measured]) &&
        all(full_gates$estimate[measured] <= full_gates$upper[measured]) &&
        all(!full_gates$passed[!measured]) && all(vapply(
            full_gates$evidence_hash,
            .association_sha256,
            logical(1L)
        )) && all(nzchar(full_gates$note))
    if (any(measured)) {
        bias <- full_gates$gate_id == "a_assoc_effect_bias_v4"
        valid <- valid && all(full_gates$estimate[measured & !bias] >= 0 &
            full_gates$estimate[measured & !bias] <= 1) &&
            all(full_gates$estimate[measured & bias] >= 0)
    }
    expected_pass <- if (nrow(full_gates)) {
        measured & mapply(
            .association_compare_threshold,
            full_gates$estimate,
            full_gates$operator,
            full_gates$threshold,
            USE.NAMES = FALSE
        )
    } else {
        logical()
    }
    if (!valid || !identical(full_gates$passed, expected_pass)) {
        .abort_association_evidence(
            "Association full-gate results are malformed or detached."
        )
    }
    invisible(full_gates)
}

.validate_association_audit_consistency <- function(
    run_bindings,
    opportunity,
    null,
    target
) {
    key <- function(frame) {
        paste(
            frame$candidate_id,
            frame$scenario_id,
            frame$replicate,
            sep = "\r"
        )
    }
    run_keys <- key(run_bindings)
    nonpanel <- run_bindings$result_kind != "panel"
    invalid_opportunity <- if (nrow(opportunity)) {
        rows <- match(key(opportunity), run_keys)
        anyNA(rows) || any(opportunity$available & nonpanel[rows])
    } else {
        FALSE
    }
    null_rows <- match(key(null), run_keys)
    target_rows <- match(key(target), run_keys)
    invalid_null <- anyNA(null_rows) || any(
        nonpanel[null_rows] &
            (null$available_p_count != 0L | null$false_flag)
    )
    invalid_target <- anyNA(target_rows) || any(
        nonpanel[target_rows] & target$outcome_status != "unavailable"
    )
    if (invalid_opportunity || invalid_null || invalid_target) {
        .abort_association_evidence(
            "Association audits contradict their run dispositions."
        )
    }
    invisible(TRUE)
}

.association_candidate_metrics <- function(
    run_bindings,
    opportunity,
    null,
    target
) {
    .validate_association_run_bindings(run_bindings)
    .validate_association_opportunity_audit(opportunity)
    .validate_association_null_audit(null)
    .validate_association_target_audit(target)
    .validate_association_audit_consistency(
        run_bindings,
        opportunity,
        null,
        target
    )
    screening <- .association_screening_metrics(
        run_bindings,
        opportunity,
        null,
        target
    )
    .validate_association_screening_metrics(screening)
    threshold_rows <- screening$metric != "common_opportunity_coverage"
    qualified <- vapply(.ASSOCIATION_CANDIDATES, function(candidate_id) {
        rows <- screening$candidate_id == candidate_id
        scope <- rows & !threshold_rows
        thresholds <- rows & threshold_rows
        all(screening$status[scope] == "measured") &&
            all(screening$status[thresholds] == "measured") &&
            all(screening$passed[thresholds])
    }, logical(1L))
    gate_rows <- lapply(.ASSOCIATION_CANDIDATES[qualified], function(
        candidate_id
    ) {
        .association_candidate_full_gates(candidate_id, null, target)
    })
    full_gates <- if (length(gate_rows)) {
        output <- do.call(rbind, gate_rows)
        row.names(output) <- NULL
        output[.ASSOCIATION_FULL_GATE_RESULT_FIELDS]
    } else {
        data.frame(
            candidate_id = character(), gate_id = character(),
            m12_gate_id = character(), stage = character(),
            metric = character(), estimate = double(),
            numerator = double(), denominator = double(), lower = double(),
            upper = double(),
            operator = character(), threshold = double(),
            status = character(), passed = logical(),
            failure_treatment = character(), evidence_hash = character(),
            note = character(),
            stringsAsFactors = FALSE
        )[.ASSOCIATION_FULL_GATE_RESULT_FIELDS]
    }
    passers <- vapply(.ASSOCIATION_CANDIDATES, function(candidate_id) {
        rows <- full_gates$candidate_id == candidate_id
        sum(rows) == 5L && all(full_gates$status[rows] == "measured") &&
            all(full_gates$passed[rows])
    }, logical(1L))
    .validate_association_full_gate_results(full_gates)
    ranking_rows <- lapply(.ASSOCIATION_CANDIDATES[passers], function(
        candidate_id
    ) {
        scope <- screening$candidate_id == candidate_id &
            screening$metric == "common_opportunity_coverage"
        power <- full_gates$candidate_id == candidate_id &
            full_gates$gate_id == "a_assoc_alternative_power_v4"
        bias <- full_gates$candidate_id == candidate_id &
            full_gates$gate_id == "a_assoc_effect_bias_v4"
        runtime <- run_bindings$elapsed_seconds[
            run_bindings$candidate_id == candidate_id
        ]
        data.frame(
            candidate_id = candidate_id,
            scope_coverage = as.double(screening$estimate[scope]),
            power = as.double(full_gates$estimate[power]),
            absolute_bias = as.double(full_gates$estimate[bias]),
            median_runtime = as.double(stats::median(runtime)),
            stringsAsFactors = FALSE
        )
    })
    ranking <- if (length(ranking_rows)) {
        output <- do.call(rbind, ranking_rows)
        row.names(output) <- NULL
        output[.ASSOCIATION_RANKING_METRIC_FIELDS]
    } else {
        data.frame(
            candidate_id = character(), scope_coverage = double(),
            power = double(), absolute_bias = double(),
            median_runtime = double(), stringsAsFactors = FALSE
        )[.ASSOCIATION_RANKING_METRIC_FIELDS]
    }
    candidate_order <- if (nrow(ranking)) {
        .association_rank_candidates(ranking)
    } else {
        character()
    }
    list(
        screening_metrics = screening,
        full_gate_results = full_gates,
        ranking_metrics = ranking,
        candidate_order = candidate_order,
        selected_candidate = if (length(candidate_order)) {
            candidate_order[[1L]]
        } else {
            NA_character_
        },
        state = if (length(candidate_order)) "winner_locked" else "no_winner"
    )
}

.association_call_study_resolver <- function(fun, ..., quantity) {
    tryCatch(
        fun(...),
        error = function(error) {
            .abort_association_evidence(paste0(
                "Association study resolver could not read ",
                quantity,
                ": ",
                conditionMessage(error)
            ))
        }
    )
}

.association_candidate_evidence_components <- function(evidence) {
    scalar_names <- c(
        "schema", "contract_hash", "protocol_hash", "gate_registry_hash",
        "implementation_hash", "effective_manifest_hash",
        "generator_protocol_hash", "artifact_inventory_hash",
        "selected_candidate", "state"
    )
    scalar <- data.frame(
        position = match(
            scalar_names,
            .ASSOCIATION_CANDIDATE_EVIDENCE_FIELDS
        ),
        field = scalar_names,
        value = unname(vapply(scalar_names, function(field) {
            evidence[[field]]
        }, character(1L))),
        type = "character",
        stringsAsFactors = FALSE
    )
    vector_frame <- function(value, type) {
        data.frame(
            position = seq_along(value),
            value = value,
            type = rep(type, length(value)),
            stringsAsFactors = FALSE
        )
    }
    nested_frame <- function(value) {
        data.frame(
            position = seq_len(nrow(value)),
            value,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }
    list(
        scalar_fields = scalar,
        candidate_ids = vector_frame(evidence$candidate_ids, "character"),
        scenario_ids = vector_frame(evidence$scenario_ids, "character"),
        replicate_ids = vector_frame(evidence$replicate_ids, "integer"),
        run_bindings = nested_frame(evidence$run_bindings),
        outcome_bindings = nested_frame(evidence$outcome_bindings),
        screening_metrics = nested_frame(evidence$screening_metrics),
        full_gate_results = nested_frame(evidence$full_gate_results),
        ranking_metrics = nested_frame(evidence$ranking_metrics),
        candidate_order = vector_frame(
            evidence$candidate_order,
            "character"
        )
    )
}

.association_candidate_evidence_hash <- function(evidence) {
    .association_evidence_bundle_hash(
        .association_candidate_evidence_components(evidence)
    )
}

.resolve_association_candidate_evidence <- function(
    resolver,
    artifact_inventory_hash = NA_character_
) {
    .validate_association_study_resolver(
        resolver,
        resolver$implementation_hash
    )
    production <- identical(resolver$mode, "frozen_protocol")
    inventory_valid <- if (production) {
        .association_sha256(artifact_inventory_hash)
    } else {
        is.character(artifact_inventory_hash) &&
            length(artifact_inventory_hash) == 1L &&
            is.na(artifact_inventory_hash)
    }
    if (!inventory_valid) {
        .abort_association_evidence(
            "Association evidence inventory binding is malformed."
        )
    }
    run_rows <- outcome_rows <- opportunity_rows <- null_rows <- target_rows <-
        list()
    runtime <- if (identical(resolver$mode, "frozen_protocol")) {
        .association_protocol_runtime(resolver)
    } else {
        NULL
    }
    position <- 0L
    for (scenario_id in .ASSOCIATION_STUDY_SCENARIOS) {
        manifest <- .association_study_manifest_row(scenario_id)
        target_index <- match(scenario_id, .ASSOCIATION_TARGETS$scenario_id)
        target <- if (is.na(target_index)) NULL else
            .ASSOCIATION_TARGETS[target_index, , drop = FALSE]
        previous_input <- NULL
        previous_material <- NULL
        for (replicate in .ASSOCIATION_STUDY_REPLICATES) {
            input <- if (identical(resolver$mode, "frozen_protocol")) {
                .association_read_protocol_input(
                    resolver,
                    scenario_id,
                    replicate
                )
            } else {
                .association_call_study_resolver(
                    resolver$input,
                    scenario_id,
                    replicate,
                    quantity = "one fixture input"
                )
            }
            if (identical(resolver$mode, "frozen_protocol")) {
                regenerated <- .association_protocol_input(
                    runtime,
                    scenario_id,
                    replicate
                )
                .validate_association_study_input(
                    regenerated,
                    scenario_id,
                    replicate
                )
                if (!identical(input, regenerated)) {
                    .abort_association_evidence(
                        paste0(
                            "Stored association input differs from exact ",
                            "frozen-generator replay."
                        )
                    )
                }
            }
            repeated <- FALSE
            if (!is.null(previous_input)) {
                expected <- previous_input
                expected$replicate <- replicate
                repeated <- identical(input, expected)
            }
            if (repeated) {
                input_sha256 <- .association_study_input_sha256(
                    input,
                    previous_material$payload_sha256
                )
                preparation <- previous_material$preparation
                oracle <- previous_material$oracle
            } else {
                .validate_association_study_input(
                    input,
                    scenario_id,
                    replicate
                )
                payload_sha256 <-
                    .association_study_input_payload_sha256(input)
                input_sha256 <- .association_study_input_sha256(
                    input,
                    payload_sha256
                )
                preparation <- .new_association_preparation(
                    input$data,
                    input$design
                )
                oracle <- .association_oracle_response(input, preparation)
                if (scenario_id %in% .ASSOCIATION_NULL_SCENARIOS) {
                    .validate_association_null_oracle(oracle)
                }
                previous_material <- list(
                    payload_sha256 = payload_sha256,
                    preparation = preparation,
                    oracle = oracle
                )
            }
            previous_input <- input
            for (candidate_id in .ASSOCIATION_CANDIDATES) {
                record <- if (identical(
                    resolver$mode,
                    "frozen_protocol"
                )) {
                    .association_read_protocol_result(
                        resolver,
                        candidate_id,
                        scenario_id,
                        replicate
                    )
                } else {
                    .association_call_study_resolver(
                        resolver$result,
                        candidate_id,
                        scenario_id,
                        replicate,
                        quantity = "one fixture candidate result"
                    )
                }
                .validate_association_study_result_record(
                    record,
                    candidate_id,
                    scenario_id,
                    replicate,
                    input_sha256
                )
                if (identical(resolver$mode, "frozen_protocol")) {
                    regenerated_result <-
                        .association_execute_protocol_candidate(
                            candidate_id,
                            preparation
                        )
                    if (!identical(record$result, regenerated_result)) {
                        .abort_association_evidence(
                            paste0(
                                "Stored association result differs from ",
                                "sealed candidate-engine replay."
                            )
                        )
                    }
                }
                resolved <- .association_resolve_run(
                    record,
                    preparation,
                    oracle,
                    candidate_id,
                    scenario_id,
                    replicate,
                    manifest$acquisition[[1L]],
                    target
                )
                position <- position + 1L
                run_rows[[position]] <- resolved$run
                outcome_rows[[position]] <- resolved$outcomes
                opportunity_rows[[position]] <- resolved$opportunity
                null_rows[[position]] <- resolved$null
                target_rows[[position]] <- resolved$target
            }
        }
    }
    run_bindings <- .association_bind_audit_rows(
        run_rows,
        .empty_association_run_bindings
    )
    run_bindings <- .association_canonicalize_audit(
        run_bindings,
        c("candidate_id", "scenario_id", "replicate")
    )
    outcome_bindings <- .association_bind_audit_rows(
        outcome_rows,
        .empty_association_outcome_bindings
    )
    outcome_bindings <- .association_canonicalize_audit(
        outcome_bindings,
        c(
            "candidate_id", "scenario_id", "replicate", "stratum",
            "hypothesis"
        )
    )
    opportunity <- .association_bind_audit_rows(
        opportunity_rows,
        .empty_association_opportunity_audit
    )
    opportunity <- .association_canonicalize_audit(
        opportunity,
        c(
            "candidate_id", "scenario_id", "replicate", "stratum",
            "hypothesis"
        )
    )
    null <- .association_bind_audit_rows(
        null_rows,
        .empty_association_null_audit
    )
    null <- .association_canonicalize_audit(
        null,
        c("candidate_id", "scenario_id", "replicate")
    )
    target_audit <- .association_bind_audit_rows(
        target_rows,
        .empty_association_target_audit
    )
    target_audit <- .association_canonicalize_audit(
        target_audit,
        c("candidate_id", "target_label", "replicate")
    )
    .validate_association_run_bindings(run_bindings)
    .validate_association_outcome_bindings(outcome_bindings)
    metrics <- .association_candidate_metrics(
        run_bindings,
        opportunity,
        null,
        target_audit
    )
    evidence_schema <- if (identical(resolver$mode, "frozen_protocol")) {
        .ASSOCIATION_CANDIDATE_EVIDENCE_SCHEMA
    } else {
        .ASSOCIATION_CANDIDATE_EVIDENCE_FIXTURE_SCHEMA
    }
    evidence_class <- if (identical(resolver$mode, "frozen_protocol")) {
        "imputefinder_association_candidate_evidence"
    } else {
        "imputefinder_association_candidate_evidence_fixture"
    }
    evidence <- structure(
        list(
            schema = evidence_schema,
            contract_hash = .ASSOCIATION_CONTRACT_HASH,
            protocol_hash = .ASSOCIATION_PROTOCOL_HASH,
            gate_registry_hash = .ASSOCIATION_GATE_REGISTRY_HASH,
            implementation_hash = resolver$implementation_hash,
            effective_manifest_hash = .ASSOCIATION_EFFECTIVE_MANIFEST_HASH,
            generator_protocol_hash = .ASSOCIATION_GENERATOR_PROTOCOL_HASH,
            artifact_inventory_hash = artifact_inventory_hash,
            candidate_ids = .ASSOCIATION_CANDIDATES,
            scenario_ids = .ASSOCIATION_STUDY_SCENARIOS,
            replicate_ids = .ASSOCIATION_STUDY_REPLICATES,
            run_bindings = run_bindings,
            outcome_bindings = outcome_bindings,
            screening_metrics = metrics$screening_metrics,
            full_gate_results = metrics$full_gate_results,
            ranking_metrics = metrics$ranking_metrics,
            candidate_order = metrics$candidate_order,
            selected_candidate = metrics$selected_candidate,
            state = metrics$state,
            evidence_hash = NA_character_
        ),
        class = evidence_class
    )
    evidence$evidence_hash <- .association_candidate_evidence_hash(evidence)
    evidence
}

.validate_association_candidate_evidence_replay <- function(
    evidence,
    resolver,
    implementation_hash
) {
    .validate_association_study_resolver(resolver, implementation_hash)
    if (!identical(resolver$mode, "fixture")) {
        .abort_association_evidence(
            "Association fixture replay validator rejects production evidence."
        )
    }
    expected <- .resolve_association_candidate_evidence(resolver)
    valid <- is.list(evidence) && identical(
        class(evidence),
        "imputefinder_association_candidate_evidence_fixture"
    ) && identical(names(evidence), .ASSOCIATION_CANDIDATE_EVIDENCE_FIELDS) &&
        identical(evidence, expected) &&
        identical(
            evidence$evidence_hash,
            .association_candidate_evidence_hash(evidence)
        )
    if (!valid) {
        .abort_association_evidence(
            paste0(
                "Association candidate evidence differs from complete ",
                "semantic resolver replay."
            )
        )
    }
    invisible(evidence)
}

.authorize_association_candidate_evidence <- function(
    resolver,
    implementation_manifest,
    root = ".",
    existing_evidence = NULL
) {
    .validate_association_implementation_manifest(
        implementation_manifest,
        root
    )
    .validate_association_study_resolver(
        resolver,
        implementation_manifest$manifest_hash
    )
    canonical_root <- .association_manifest_root(root)
    valid <- identical(resolver$mode, "frozen_protocol") && identical(
        resolver$implementation_hash,
        implementation_manifest$manifest_hash
    ) && identical(
        resolver$root,
        canonical_root
    ) && identical(
        resolver$source_sha256,
        .association_source_manifest_sha256(
            implementation_manifest$source_files
        )
    )
    if (!valid) {
        .abort_association_evidence(
            paste0(
                "Association evidence authorization requires the exact ",
                "locked implementation and protocol-replay resolver."
            )
        )
    }
    evidence_policy <- if (is.null(existing_evidence)) "absent" else "exact"
    inventory_before <- .association_protocol_artifact_inventory(
        resolver$artifact_root,
        implementation_manifest$manifest_hash,
        complete = TRUE,
        evidence_policy = evidence_policy,
        evidence = existing_evidence
    )
    artifact_inventory_hash <-
        .association_protocol_allocation_inventory_hash(inventory_before)
    locked_before <- .association_read_locked_implementation_manifest(
        resolver$artifact_root
    )
    if (!identical(locked_before, implementation_manifest)) {
        .abort_association_evidence(
            "Association production replay differs from its on-disk implementation lock."
        )
    }
    evidence <- .resolve_association_candidate_evidence(
        resolver,
        artifact_inventory_hash
    )
    .validate_association_implementation_manifest(
        implementation_manifest,
        canonical_root
    )
    valid_after <- identical(
        resolver$source_sha256,
        .association_source_manifest_sha256(
            implementation_manifest$source_files
        )
    )
    inventory_after <- .association_protocol_artifact_inventory(
        resolver$artifact_root,
        implementation_manifest$manifest_hash,
        complete = TRUE,
        evidence_policy = evidence_policy,
        evidence = existing_evidence
    )
    locked_after <- .association_read_locked_implementation_manifest(
        resolver$artifact_root
    )
    if (!valid_after || !identical(locked_after, implementation_manifest) ||
        !identical(inventory_after, inventory_before) || !identical(
            .association_protocol_allocation_inventory_hash(inventory_after),
            artifact_inventory_hash
        )) {
        .abort_association_evidence(
            "Association implementation or artifacts changed during evidence replay."
        )
    }
    evidence
}

.validate_association_candidate_evidence <- function(
    evidence,
    resolver,
    implementation_manifest,
    root = "."
) {
    expected <- .authorize_association_candidate_evidence(
        resolver,
        implementation_manifest,
        root,
        existing_evidence = evidence
    )
    valid <- identical(
        class(evidence),
        "imputefinder_association_candidate_evidence"
    ) && identical(
        evidence$schema,
        .ASSOCIATION_CANDIDATE_EVIDENCE_SCHEMA
    ) && identical(evidence, expected) && identical(
        evidence$evidence_hash,
        .association_candidate_evidence_hash(evidence)
    )
    if (!valid) {
        .abort_association_evidence(
            paste0(
                "Association candidate evidence differs from its exact ",
                "authorized production replay."
            )
        )
    }
    invisible(evidence)
}
