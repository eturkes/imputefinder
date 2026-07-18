#!/usr/bin/env Rscript

# M13a mandatory algebraic-core development gate. Uses only frozen synthetic
# replicates 33-64; A association candidates and external artifacts stay closed.
# Run from repository root after installing the current package:
# Rscript --vanilla dev/m13-design-core-validation.R --verify

.m13_contract <- local({
    environment <- new.env(parent = baseenv())
    sys.source("dev/m12-validation-contract.R", envir = environment)
    environment
})

.m13_generator <- local({
    environment <- new.env(parent = globalenv())
    sys.source("dev/m12-generator-validation.R", envir = environment)
    environment
})

.M13_PROTOCOL_ID <- "m13_design_core_validation_v1"
.M13_DEVELOPMENT_REPLICATES <- 33:64
.M13_MANDATORY_GATES <- c(
    "a_rank_deficient_sensitivity",
    "a_rank_full_specificity",
    "a_alias_exact_sensitivity",
    "a_alias_false_positive",
    "a_contrast_nonestimable_rejection",
    "a_contrast_estimable_retention",
    "a_block_unit_accounting",
    "a_unequal_design_accounting",
    "a_named_order_invariance",
    "a_side_effect_global",
    "a_side_effect_classic"
)

.M13_EXPECTED_AUDIT_HASH <-
    "37528216590dcfff3338501ef61b766697ea63c14b7324912fcb00f84d58790e"
.M13_EXPECTED_GATE_RESULT_HASH <-
    "cb72a69301100287cee569ae7a0d709f42edc10e3b10a62dbff3a86ab88aec6f"

.m13_fail <- function(...) {
    stop(..., call. = FALSE)
}

.m13_hash_object <- function(value) {
    unname(tools::sha256sum(bytes = serialize(value, NULL, version = 3L)))
}

.m13_package_functions <- function() {
    if (!requireNamespace("imputefinder", quietly = TRUE)) {
        .m13_fail("Install the current imputefinder package before M13a validation.")
    }
    namespace <- asNamespace("imputefinder")
    list(
        design = getExportedValue("imputefinder", "missingness_design"),
        analyze = getExportedValue("imputefinder", "analyze_missingness"),
        classify = getExportedValue("imputefinder", "classify_missingness"),
        core = get(".new_design_estimability", envir = namespace),
        contrast = get(".design_contrast_estimability", envir = namespace),
        validate = get(".validate_design_estimability", envir = namespace)
    )
}

.m13_validate_evidence_boundary <- function() {
    contract <- .m13_contract$m12_validation_contract()
    generator_hash <- unname(
        .m13_generator$m12b_protocol_hashes()[["protocol"]]
    )
    expected_generator_hash <- contract$generator_protocol$protocol_hash[[1L]]
    closed <- all(contract$artifact_manifest$open_state == "metadata_only") &&
        all(is.na(contract$artifact_manifest$local_sha256)) &&
        all(contract$protocol_registry$result_state == "frozen_unrun")
    if (!identical(generator_hash, expected_generator_hash) || !closed) {
        .m13_fail(
            "M13a requires the frozen generator plus sealed candidate and external-result boundary."
        )
    }
    invisible(contract)
}

.m13_nuisance_columns <- function(manifest) {
    switch(
        manifest$nuisance_role,
        none = NULL,
        run_order = "run_order",
        batch_crossed = "batch",
        batch_partial = "batch",
        batch_aliased = "batch",
        technical_replicate = "technical_replicate",
        .m13_fail("Unknown M13a nuisance role: ", manifest$nuisance_role)
    )
}

.m13_design_columns <- function(simulation) {
    nuisance <- .m13_nuisance_columns(simulation$manifest)
    block <- if (identical(simulation$manifest$block_role, "subject")) {
        "subject"
    } else {
        NULL
    }
    c("condition", nuisance, block, "acquisition_mode")
}

.m13_typed_design <- function(simulation, functions, permuted = FALSE) {
    columns <- .m13_design_columns(simulation)
    metadata <- simulation$design[, columns, drop = FALSE]
    if (permuted) {
        metadata <- metadata[
            rev(rownames(metadata)),
            rev(names(metadata)),
            drop = FALSE
        ]
        for (column in columns) {
            values <- metadata[[column]]
            if (is.character(values)) {
                metadata[[column]] <- factor(
                    values,
                    levels = rev(sort(unique(values), method = "radix"))
                )
            } else if (is.numeric(values)) {
                metadata[[column]] <- if (is.integer(values)) {
                    as.double(values)
                } else {
                    as.integer(values)
                }
            }
        }
    }
    functions$design(
        metadata,
        condition = "condition",
        nuisance = .m13_nuisance_columns(simulation$manifest),
        block = if (identical(
            simulation$manifest$block_role,
            "subject"
        )) "subject" else NULL,
        acquisition = "acquisition_mode"
    )
}

.m13_condition_contrast <- function(simulation) {
    levels <- sort(unique(simulation$design$condition), method = "radix")
    weights <- stats::setNames(c(-1, 1), levels[c(1L, length(levels))])
    list(condition = weights)
}

.m13_snapshot <- function() {
    has_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
    list(
        rng_kind = RNGkind(),
        has_seed = has_seed,
        seed = if (has_seed) {
            get(".Random.seed", globalenv(), inherits = FALSE)
        } else {
            NULL
        },
        options = options(),
        directory = getwd()
    )
}

.m13_snapshot_identical <- function(before, after) {
    identical(before, after)
}

.m13_classic_without_call <- function(result) {
    result$call <- NULL
    result
}

.m13_classic_pair <- function(simulation, design, functions) {
    group <- stats::setNames(
        simulation$design$condition,
        colnames(simulation$data)
    )
    before <- functions$classify(
        simulation$data,
        group = group,
        cutoffs = simulation$manual_cutoffs,
        seed = simulation$rescue_seed
    )
    analysis <- functions$analyze(
        simulation$data,
        design,
        cutoffs = simulation$manual_cutoffs,
        seed = simulation$rescue_seed,
        modules = "sentinel"
    )
    after <- analysis$classic
    before_canonical <- .m13_classic_without_call(before)
    after_canonical <- .m13_classic_without_call(after)
    list(
        drift = !identical(before_canonical, after_canonical),
        classic_hash = .m13_hash_object(before_canonical),
        analysis = analysis
    )
}

.m13_block_match <- function(core, simulation) {
    expected_role <- if (identical(
        simulation$manifest$block_role,
        "subject"
    )) "block" else "sample"
    expected <- stats::setNames(
        simulation$design$resampling_unit,
        rownames(simulation$design)
    )
    sample <- core$units$sample
    expected <- unname(expected[sample$sample])
    expected_grouped <- any(duplicated(data.frame(
        unit = simulation$design$resampling_unit,
        condition = simulation$design$condition
    )))
    identical(core$units$resampling_role, expected_role) &&
        identical(sample$unit, expected) &&
        identical(
            core$units$independent_unit_count,
            as.integer(length(unique(expected)))
        ) &&
        identical(
            core$units$grouped_technical_siblings,
            expected_grouped
        )
}

.m13_unequal_match <- function(core, simulation) {
    expected <- table(simulation$design$condition)
    condition <- core$units$condition
    expected_count <- as.integer(expected[condition$condition])
    identical(
        rownames(core$model$matrix),
        sort(rownames(simulation$design), method = "radix")
    ) &&
        identical(core$units$sample$weight, rep(1L, nrow(simulation$design))) &&
        identical(condition$sample_count, expected_count)
}

.m13_alias_detected <- function(core) {
    affected <- core$aliasing$affected_terms$term
    !core$rank$full_column_rank &&
        all(c("condition", "batch") %in% affected)
}

.m13_run_instance <- function(scenario_id, replicate, functions) {
    simulation <- .m13_generator$simulate_m12b_scenario(
        scenario_id,
        as.integer(replicate)
    )
    design <- .m13_typed_design(simulation, functions)
    original_data <- simulation$data
    original_design <- design
    caller <- .m13_snapshot()

    core <- functions$core(design)
    functions$validate(core, design)
    contrast <- functions$contrast(
        core,
        .m13_condition_contrast(simulation)
    )
    permuted_design <- .m13_typed_design(
        simulation,
        functions,
        permuted = TRUE
    )
    permuted_core <- functions$core(permuted_design)
    classic <- .m13_classic_pair(simulation, design, functions)

    caller_after <- .m13_snapshot()
    global_mutation <- !.m13_snapshot_identical(caller, caller_after) ||
        !identical(simulation$data, original_data) ||
        !identical(design, original_design)
    rank_oracle_full <- simulation$design_truth$model_rank ==
        simulation$design_truth$model_columns

    data.frame(
        scenario_id = scenario_id,
        replicate = as.integer(replicate),
        acquisition = simulation$manifest$acquisition,
        rank_oracle_full = rank_oracle_full,
        rank_full = core$rank$full_column_rank,
        rank_oracle_match = identical(
            core$rank$full_column_rank,
            rank_oracle_full
        ),
        alias_detected = .m13_alias_detected(core),
        contrast_estimable = contrast$estimable,
        block_unit_match = .m13_block_match(core, simulation),
        unequal_design_match = .m13_unequal_match(core, simulation),
        order_invariant = identical(core, permuted_core),
        global_mutation = global_mutation,
        classic_drift = classic$drift,
        core_sha256 = .m13_hash_object(core),
        classic_sha256 = classic$classic_hash,
        stringsAsFactors = FALSE
    )
}

m13_design_core_audit <- function() {
    functions <- .m13_package_functions()
    scenarios <- .m13_generator$m12b_scenario_parameters()$scenario_id
    grid <- expand.grid(
        scenario_id = scenarios,
        replicate = .M13_DEVELOPMENT_REPLICATES,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    grid <- grid[order(
        match(grid$scenario_id, scenarios),
        grid$replicate
    ), , drop = FALSE]
    rows <- lapply(seq_len(nrow(grid)), function(index) {
        .m13_run_instance(
            grid$scenario_id[[index]],
            as.integer(grid$replicate[[index]]),
            functions
        )
    })
    audit <- do.call(rbind, rows)
    row.names(audit) <- NULL
    audit
}

.m13_normative_classic_drift <- function(functions) {
    x <- rbind(
        on_off = c(NA, NA, NA, NA, 15, 16, 15, 17),
        mar_both = c(14, 15, NA, 16, 14, NA, 15, 16),
        sparse_mar = c(20, NA, NA, NA, 20, 21, 22, 23),
        sparse_mnar = c(8, NA, NA, NA, 14, 15, 16, 17),
        all_mnar = c(8, NA, NA, NA, 9, NA, NA, NA),
        globally_absent = rep(NA, 8),
        complete_low = c(8, 8, 8, 8, 9, 9, 9, 9)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    group <- rep(c("A", "B"), each = 4L)
    cutoffs <- c(A = 12, B = 12)
    design <- functions$design(
        data.frame(condition = group, row.names = colnames(x)),
        condition = "condition"
    )
    before <- functions$classify(x, group = group, cutoffs = cutoffs)
    after <- functions$analyze(
        x,
        design,
        cutoffs = cutoffs,
        modules = "sentinel"
    )$classic
    !identical(
        .m13_classic_without_call(before),
        .m13_classic_without_call(after)
    )
}

.m13_scenarios <- function(...) {
    sort(unique(c(...)), method = "radix")
}

.m13_gate_definition <- function(audit, normative_drift) {
    all_rows <- rep(TRUE, nrow(audit))
    perfect <- audit$scenario_id == "dda_batch_perfect"
    full <- audit$scenario_id %in% .m13_scenarios(
        "dda_batch_crossed", "dia_batch_partial", "dia_monotone_paired"
    )
    alias_nonperfect <- audit$scenario_id %in% .m13_scenarios(
        "dda_batch_crossed", "dia_batch_partial"
    )
    paired <- audit$scenario_id %in% .m13_scenarios(
        "dda_grouped_leakage_trap", "dia_blockwise_paired",
        "dia_monotone_paired"
    )
    unequal <- audit$scenario_id %in% .m13_scenarios(
        "dda_monotone_unequal", "dia_batch_partial", "dia_null_unequal"
    )
    list(
        a_rank_deficient_sensitivity = list(
            selected = perfect,
            outcome = !audit$rank_full,
            kind = "fraction"
        ),
        a_rank_full_specificity = list(
            selected = full,
            outcome = audit$rank_full,
            kind = "fraction"
        ),
        a_alias_exact_sensitivity = list(
            selected = perfect,
            outcome = audit$alias_detected,
            kind = "fraction"
        ),
        a_alias_false_positive = list(
            selected = alias_nonperfect,
            outcome = audit$alias_detected,
            kind = "fraction"
        ),
        a_contrast_nonestimable_rejection = list(
            selected = perfect,
            outcome = !audit$contrast_estimable,
            kind = "fraction"
        ),
        a_contrast_estimable_retention = list(
            selected = full,
            outcome = audit$contrast_estimable,
            kind = "fraction"
        ),
        a_block_unit_accounting = list(
            selected = paired,
            outcome = audit$block_unit_match,
            kind = "fraction"
        ),
        a_unequal_design_accounting = list(
            selected = unequal,
            outcome = audit$unequal_design_match,
            kind = "fraction"
        ),
        a_named_order_invariance = list(
            selected = all_rows,
            outcome = audit$order_invariant,
            kind = "fraction"
        ),
        a_side_effect_global = list(
            selected = all_rows,
            outcome = audit$global_mutation,
            kind = "count"
        ),
        a_side_effect_classic = list(
            selected = all_rows,
            outcome = c(audit$classic_drift, normative_drift),
            kind = "count",
            denominator = nrow(audit) + 1L
        )
    )
}

.m13_gate_evidence_hash <- function(gate_id, audit, definition) {
    selected <- audit[definition$selected, , drop = FALSE]
    outcome <- definition$outcome
    if (length(outcome) == nrow(audit) + 1L) {
        selected <- audit
    }
    .m13_hash_object(list(
        protocol_id = .M13_PROTOCOL_ID,
        gate_id = gate_id,
        instances = selected,
        outcome = outcome
    ))
}

m13_gate_results <- function(audit, normative_drift) {
    contract <- .m13_contract$m12_validation_contract()
    registry <- contract$gate_registry
    registry <- registry[
        match(.M13_MANDATORY_GATES, registry$gate_id),
        ,
        drop = FALSE
    ]
    definitions <- .m13_gate_definition(audit, normative_drift)
    rows <- lapply(.M13_MANDATORY_GATES, function(gate_id) {
        definition <- definitions[[gate_id]]
        selected_outcome <- if (length(definition$outcome) == nrow(audit)) {
            definition$outcome[definition$selected]
        } else {
            definition$outcome
        }
        denominator <- if (!is.null(definition$denominator)) {
            definition$denominator
        } else {
            sum(definition$selected)
        }
        numerator <- sum(selected_outcome)
        estimate <- if (identical(definition$kind, "fraction")) {
            numerator / denominator
        } else {
            numerator
        }
        row <- registry[registry$gate_id == gate_id, , drop = FALSE]
        if (!identical(as.integer(denominator), row$sample_size)) {
            .m13_fail("M13a denominator drift for ", gate_id)
        }
        data.frame(
            gate_id = gate_id,
            registry_version = row$registry_version,
            metric = row$metric,
            estimate = as.numeric(estimate),
            numerator = as.numeric(numerator),
            denominator = as.numeric(denominator),
            lower = NA_real_,
            upper = NA_real_,
            status = "measured",
            gate_binding_hash = unname(
                .m13_contract$m12_gate_binding_hashes()[gate_id]
            ),
            evidence_hash = .m13_gate_evidence_hash(
                gate_id,
                audit,
                definition
            ),
            note = "complete frozen development enumeration; association candidates unopened",
            stringsAsFactors = FALSE
        )
    })
    results <- do.call(rbind, rows)
    row.names(results) <- NULL
    results
}

m13_validation_bundle <- function() {
    .m13_validate_evidence_boundary()
    audit <- m13_design_core_audit()
    normative_drift <- .m13_normative_classic_drift(
        .m13_package_functions()
    )
    results <- m13_gate_results(audit, normative_drift)
    report <- .m13_contract$m12_evaluate_gate_results(
        results,
        gate_ids = .M13_MANDATORY_GATES
    )
    list(
        audit = audit,
        normative_classic_drift = normative_drift,
        gate_results = results,
        gate_report = report,
        audit_hash = .m13_hash_object(audit),
        gate_result_hash = .m13_hash_object(results)
    )
}

.m13_main <- function(args) {
    if (!identical(args, "--verify")) {
        .m13_fail(
            "usage: Rscript --vanilla dev/m13-design-core-validation.R --verify"
        )
    }
    bundle <- m13_validation_bundle()
    if (!all(bundle$gate_report$passed)) {
        .m13_fail(
            "M13a mandatory gate failures: ",
            paste(
                bundle$gate_report$gate_id[!bundle$gate_report$passed],
                collapse = ", "
            )
        )
    }
    frozen <- !anyNA(c(
        .M13_EXPECTED_AUDIT_HASH,
        .M13_EXPECTED_GATE_RESULT_HASH
    ))
    if (frozen && !identical(
        c(bundle$audit_hash, bundle$gate_result_hash),
        c(.M13_EXPECTED_AUDIT_HASH, .M13_EXPECTED_GATE_RESULT_HASH)
    )) {
        .m13_fail("M13a frozen evidence hash mismatch.")
    }

    cat("protocol_id: ", .M13_PROTOCOL_ID, "\n", sep = "")
    cat("development_replicates: 33-64; instances=", nrow(bundle$audit), "\n", sep = "")
    cat("association_state: frozen_unrun; external_artifacts: metadata_only\n")
    cat("audit_sha256: ", bundle$audit_hash, "\n", sep = "")
    cat("gate_results_sha256: ", bundle$gate_result_hash, "\n", sep = "")
    print(bundle$gate_report[, c(
        "gate_id", "estimate", "operator", "threshold", "passed"
    )], row.names = FALSE)
    if (!frozen) {
        cat("hash_state: bootstrap_unfrozen\n")
    }
    invisible(TRUE)
}

if (sys.nframe() == 0L) {
    .m13_main(commandArgs(trailingOnly = TRUE))
}
