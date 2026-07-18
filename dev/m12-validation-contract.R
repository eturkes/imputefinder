#!/usr/bin/env Rscript

# M12c validation contract. This file freezes scope, generator/data roles,
# A-C candidate protocols, and numeric gates; it neither downloads nor reads
# external result artifacts.
# Run from the repository root:
# Rscript --vanilla dev/m12-validation-contract.R --verify

.M12_CONTRACT_VERSION <- "m12_validation_contract_v3"
.M12_CHECKED_ON <- "2026-07-18"

.m12_schema_catalog <- function() {
    rows <- list(
        generator_manifest = c(
            scenario_id = "character:id,unique",
            generator_version = "character:id",
            evidence_tier = "character:tier1|tier2",
            acquisition = "character:DDA|DIA",
            acquisition_profile = "character:id",
            n_features = "integer:positive",
            n_conditions = "integer:2..4",
            condition_sizes = "character:named-positive-counts",
            sampling_design = "character:independent|paired",
            block_role = "character:none|subject",
            nuisance_role = "character:declared-role",
            confounding = "character:none|partial|perfect",
            missingness_patterns = "character:token-set",
            intensity_stress = "character:token-set",
            evidence_support = "character:adequate|low",
            truth_targets = "character:token-set",
            negative_controls = "character:token-set|none",
            required_tracks = "character:A/B/C-token-set",
            seed_stream = "character:id",
            implementation_state = "character:implemented_verified"
        ),
        generator_protocol = c(
            protocol_id = "character:id,unique",
            generator_version = "character:id",
            seed_stream = "character:id",
            development_replicates = "integer:positive",
            implementation_file = "character:repository-path",
            protocol_hash = "character:sha256",
            generator_audit_hash = "character:sha256",
            stable_v1_audit_hash = "character:sha256",
            implementation_state = "character:implemented_verified"
        ),
        claim_inventory = c(
            claim_id = "character:id,unique",
            track = "character:A|B|C",
            component = "character:id",
            claim = "character:nonempty",
            evidence_tiers = "character:token-set",
            required_strata = "character:token-set",
            disposition = "character:mandatory|independent"
        ),
        protocol_registry = c(
            protocol_id = "character:id,unique",
            track = "character:A|B|C",
            purpose = "character:nonempty",
            implementation_file = "character:repository-path",
            protocol_hash = "character:sha256",
            result_state = "character:frozen_unrun"
        ),
        internal_evidence = c(
            data_id = "character:id,unique",
            evidence_tier = "character:tier0",
            role = "character:compatibility",
            scope = "character:nonempty",
            unit = "character:nonempty",
            sample_size = "integer:positive",
            source_path = "character:repository-path",
            content_sha256 = "character:sha256",
            state = "character:frozen"
        ),
        gate_registry = c(
            gate_id = "character:id,unique",
            registry_version = "character:id",
            claim_id = "character:claim-foreign-key",
            track = "character:A|B|C",
            metric = "character:id",
            unit = "character:nonempty",
            comparator = "character:nonempty",
            strata = "character:nonempty",
            sample_size = "integer:positive",
            uncertainty_method = "character:nonempty",
            threshold_operator = "character:<=|>=|==|<|>",
            threshold_value = "numeric:finite",
            threshold_scale = "character:nonempty",
            failure_treatment = "character:nonempty",
            data_ids = "character:generator/dataset-token-set",
            protocol_id = "character:id",
            data_hash = "character:sha256",
            protocol_hash = "character:sha256"
        ),
        gate_results = c(
            gate_id = "character:gate-foreign-key,unique",
            registry_version = "character:id",
            metric = "character:id",
            estimate = "numeric:finite-if-measured-else-NA",
            numerator = "numeric:finite-if-measured-else-NA",
            denominator = "numeric:finite-positive-if-measured-else-NA",
            lower = "numeric:finite-or-NA-if-measured;NA-otherwise",
            upper = "numeric:finite-or-NA-if-measured;NA-otherwise",
            status = "character:measured|failed|unavailable",
            gate_binding_hash = "character:sha256",
            evidence_hash = "character:sha256",
            note = "character:nonempty"
        ),
        data_cards = c(
            dataset_id = "character:id,unique",
            card_version = "character:id",
            title = "character:nonempty",
            accession = "character:nonempty",
            repository = "character:nonempty",
            metadata_url = "character:https",
            evidence_tiers = "character:token-set",
            acquisition = "character:DDA/DIA-token-set",
            deployment_target = "character:nonempty",
            experimental_unit = "character:nonempty",
            related_sample_rule = "character:nonempty",
            n_runs = "integer:positive",
            n_conditions = "integer:2..5-candidate",
            condition_layout = "character:nonempty",
            nuisance_roles = "character:token-set|none",
            block_roles = "character:token-set|none",
            truth_scope = "character:nonempty",
            truth_limit = "character:nonempty",
            metadata_scope = "character:nonempty",
            licence_id = "character:nonempty",
            licence_url = "character:https",
            licence_status = "character:verified_open|terms_reviewed",
            exclusions = "character:nonempty",
            split_feasibility = "character:whole_family_only|grouped_unit_only",
            split_unit = "character:nonempty",
            related_dataset_ids = "character:token-set|none",
            role = "character:unassigned|training|development|confirmation",
            split_hash = "character:sha256-or-NA",
            checked_on = "character:YYYY-MM-DD"
        ),
        data_roles = c(
            dataset_id = "character:dataset-foreign-key,unique",
            protocol_id = "character:id",
            role = "character:development|confirmation",
            included_artifacts = "character:artifact-token-set",
            included_conditions = "character:condition-token-set",
            excluded_conditions = "character:condition-token-set|none",
            grouping_rule = "character:nonempty",
            analysis_unit = "character:nonempty",
            linked_tracks = "character:A/B/C-token-set",
            linked_claims = "character:claim-token-set",
            opening_stage = "character:nonempty",
            licence_snapshot = "character:nonempty"
        ),
        artifact_manifest = c(
            artifact_id = "character:id,unique",
            dataset_id = "character:dataset-foreign-key",
            acquisition = "character:DDA|DIA",
            file_name = "character:nonempty",
            source_url = "character:https",
            file_bytes = "numeric:positive-integer",
            checksum_algorithm = "character:sha1",
            upstream_checksum = "character:sha1",
            checksum_source = "character:nonempty",
            content_scope = "character:nonempty",
            role = "character:unassigned|training|development|confirmation",
            open_state = "character:metadata_only|opened",
            local_sha256 = "character:sha256-or-NA"
        ),
        candidate_exclusions = c(
            candidate_id = "character:id,unique",
            accession = "character:nonempty",
            reason_code = "character:id",
            reason = "character:nonempty",
            reconsider_when = "character:nonempty",
            metadata_url = "character:https"
        ),
        source_manifest = c(
            source_id = "character:id,unique",
            scope = "character:nonempty",
            url = "character:https",
            checked_on = "character:YYYY-MM-DD",
            note = "character:nonempty"
        )
    )
    do.call(
        rbind,
        lapply(names(rows), function(schema_name) {
            fields <- rows[[schema_name]]
            data.frame(
                schema = schema_name,
                position = seq_along(fields),
                field = names(fields),
                rule = unname(fields),
                stringsAsFactors = FALSE,
                row.names = NULL
            )
        })
    )
}

.m12_generator_manifest <- function() {
    data.frame(
        scenario_id = c(
            "dda_null_balanced", "dia_null_unequal",
            "dda_monotone_unequal", "dia_monotone_paired",
            "dda_mixed_outlier", "dia_blockwise_paired",
            "dda_batch_crossed", "dia_batch_partial",
            "dda_batch_perfect", "dia_structural_off",
            "dda_no_cliff", "dia_no_cliff_low_support",
            "dda_grouped_leakage_trap"
        ),
        generator_version = rep("m12_generator_v1", 13L),
        evidence_tier = c(
            "tier1", "tier1", "tier1", "tier1", "tier2", "tier2",
            "tier1", "tier1", "tier1", "tier2", "tier1", "tier2",
            "tier2"
        ),
        acquisition = c(
            "DDA", "DIA", "DDA", "DIA", "DDA", "DIA", "DDA",
            "DIA", "DDA", "DIA", "DDA", "DIA", "DDA"
        ),
        acquisition_profile = paste0(
            "synthetic_",
            tolower(c(
                "DDA", "DIA", "DDA", "DIA", "DDA", "DIA", "DDA",
                "DIA", "DDA", "DIA", "DDA", "DIA", "DDA"
            )),
            "_v1"
        ),
        n_features = as.integer(c(
            3000, 3000, 4000, 4000, 3500, 3500, 4000,
            4000, 3000, 3500, 3000, 1800, 2500
        )),
        n_conditions = as.integer(c(2, 3, 4, 2, 3, 2, 2, 3, 2, 4, 2, 2, 2)),
        condition_sizes = c(
            "A=8;B=8", "A=4;B=7;C=10", "A=4;B=6;C=8;D=10",
            "A=10;B=10", "A=6;B=6;C=6", "A=8;B=8",
            "A=12;B=12", "A=6;B=8;C=10", "A=8;B=8",
            "A=6;B=6;C=6;D=6", "A=8;B=8", "A=3;B=3",
            "A=12;B=12"
        ),
        sampling_design = c(
            "independent", "independent", "independent", "paired",
            "independent", "paired", "independent", "independent",
            "independent", "independent", "independent", "independent",
            "paired"
        ),
        block_role = c(
            "none", "none", "none", "subject", "none", "subject",
            "none", "none", "none", "none", "none", "none", "subject"
        ),
        nuisance_role = c(
            "none", "none", "run_order", "none", "run_order", "none",
            "batch_crossed", "batch_partial", "batch_aliased", "none",
            "none", "none", "technical_replicate"
        ),
        confounding = c(
            "none", "none", "none", "none", "none", "none", "none",
            "partial", "perfect", "none", "none", "none", "none"
        ),
        missingness_patterns = c(
            "intensity_independent", "intensity_independent",
            "monotone_abundance", "monotone_abundance",
            "intensity_independent;mixed;monotone_abundance", "blockwise",
            "batch_associated", "batch_associated;mixed",
            "batch_associated", "structural_off_compatible",
            "no_cliff", "intensity_independent;no_cliff",
            "blockwise;intensity_independent"
        ),
        intensity_stress = c(
            "none", "heteroscedastic", "differing_support;heteroscedastic",
            "none", "heteroscedastic;outliers", "differing_support",
            "none", "differing_support", "none", "differing_support",
            "none", "heteroscedastic", "none"
        ),
        evidence_support = c(
            "adequate", "adequate", "adequate", "adequate", "adequate",
            "adequate", "adequate", "adequate", "adequate", "adequate",
            "low", "low", "adequate"
        ),
        truth_targets = c(
            "null_association;state;retention",
            "null_association;state;retention",
            "cutoff;detection_contrast;state;retention",
            "detection_contrast;paired_effect;state",
            "detection_contrast;state;retention",
            "block_effect;detection_contrast;state",
            "batch_association;estimability;state",
            "batch_association;estimability;state",
            "aliasing;nonestimability",
            "detection_contrast;structural_off_compatibility;state",
            "cutoff_failure;weak_identification",
            "abstention;cutoff_failure;weak_identification",
            "grouped_split;resampling_unit"
        ),
        negative_controls = c(
            "association_null", "association_null", "none", "none", "none",
            "none", "none", "none", "confounded_biology", "causal_wording",
            "flat_profile", "low_support", "leakage_trap"
        ),
        required_tracks = c(
            "A;B;C", "A;B;C", "B;C", "A;B;C", "B;C", "A;B;C",
            "A;B;C", "A;B;C", "A;C", "B;C", "B", "A;B;C", "A;B;C"
        ),
        seed_stream = rep("m12_sim_seed_stream_v1", 13L),
        implementation_state = rep("implemented_verified", 13L),
        stringsAsFactors = FALSE
    )
}

.m12_generator_protocol <- function() {
    data.frame(
        protocol_id = "m12_generator_protocol_v1",
        generator_version = "m12_generator_v1",
        seed_stream = "m12_sim_seed_stream_v1",
        development_replicates = 64L,
        implementation_file = "dev/m12-generator-validation.R",
        protocol_hash =
            "cdea1bf874152e63fba08c49e390da1dd18690105e43a2dc7a03fc6866d0d080",
        generator_audit_hash =
            "4d216966c8cc54b1b23685d8801c7c65b2afcc5879dd15430d57a88f765e8a00",
        stable_v1_audit_hash =
            "8463f4306fe9b6192336487e4ca365c8c7e3525dc06aeca5b66b1a99b1a10250",
        implementation_state = "implemented_verified",
        stringsAsFactors = FALSE
    )
}

.m12_claim_inventory <- function() {
    rows <- list(
        c("A_design_rank", "A", "design_core", "Constructed rank deficiency is detected exactly.", "tier1", "DDA;DIA", "mandatory"),
        c("A_design_alias", "A", "design_core", "Perfect condition-nuisance aliasing is detected exactly.", "tier1", "DDA;DIA", "mandatory"),
        c("A_contrast_estimability", "A", "design_core", "Non-estimable contrasts are rejected and estimable contrasts retained.", "tier1", "DDA;DIA", "mandatory"),
        c("A_block_accounting", "A", "design_core", "Paired blocks and unequal replication are accounted for at their declared units.", "tier1;tier2", "DDA;DIA", "mandatory"),
        c("A_order_invariance", "A", "design_core", "Label-preserving order and re-encoding leave named outputs invariant.", "tier1", "DDA;DIA", "mandatory"),
        c("A_side_effects", "A", "design_core", "Inputs, classic branch, RNG, options, and global state remain unchanged.", "tier0;tier1", "DDA;DIA", "mandatory"),
        c("A_assoc_null", "A", "association_panel", "Association flags satisfy the frozen null error gate.", "tier1;tier2", "DDA;DIA", "independent"),
        c("A_assoc_public", "A", "association_panel", "A frozen multibatch public case meets coverage and wording gates.", "tier3;tier4", "DDA;DIA", "independent"),
        c("A_assoc_calibration", "A", "association_panel", "Declared alternatives meet uncertainty and multiplicity calibration gates.", "tier1;tier2", "DDA;DIA", "independent"),
        c("B_cutoff_coverage", "B", "estimator_uncertainty", "Cutoff ranges meet coverage where the generating model licenses coverage.", "tier1;tier2", "DDA;DIA", "independent"),
        c("B_false_confidence", "B", "estimator_uncertainty", "False-confidence frequency stays within the frozen null bound.", "tier1;tier2", "DDA;DIA", "independent"),
        c("B_no_cliff", "B", "estimator_uncertainty", "Flat and no-cliff controls report weak identification or failure.", "tier1;tier2", "DDA;DIA", "independent"),
        c("B_sparse_abstention", "B", "sampling_stability", "Unsupported sparse panels return structured unavailable independently.", "tier1;tier2", "DDA;DIA", "independent"),
        c("B_reproducibility", "B", "all_families", "Instance/input order and seed streams reproduce canonical named output.", "tier1", "DDA;DIA", "independent"),
        c("B_side_effects", "B", "all_families", "Perturbations leave caller RNG and global state unchanged.", "tier0;tier1", "DDA;DIA", "independent"),
        c("B_stable_overhead", "B", "execution", "Stable default calls perform zero sidecar work and retain v1 gates.", "tier0", "stable", "independent"),
        c("B_external_behavior", "B", "all_families", "Controlled public matrices produce complete stratified panels without an inferential cutoff-coverage claim.", "tier3", "DDA;DIA", "independent"),
        c("B_display_calibration", "B", "communication", "Optional stability labels pass frozen held-out risk-coverage gates or remain disabled.", "tier1;tier2", "DDA;DIA", "independent"),
        c("C_null_error", "C", "detection_abundance", "Correlated-feature null rejection and FDR satisfy separate component bounds.", "tier1;tier2", "DDA;DIA", "independent"),
        c("C_false_sign", "C", "detection_abundance", "Component-specific false-sign behavior satisfies frozen bounds.", "tier1;tier2", "DDA;DIA", "independent"),
        c("C_effect_bias", "C", "detection_abundance", "Detection and observed-abundance effect bias satisfy separate gates.", "tier1;tier2;tier3", "DDA;DIA", "independent"),
        c("C_interval_coverage", "C", "detection_abundance", "Intervals meet estimand-specific coverage gates where licensed.", "tier1;tier2;tier3", "DDA;DIA", "independent"),
        c("C_alternative_utility", "C", "detection_abundance", "Frozen alternatives meet component-specific power and precision-recall gates.", "tier1;tier2;tier3", "DDA;DIA", "independent"),
        c("C_edge_cases", "C", "design_support", "Paired, unequal, separated, low-support, and confounded cases follow frozen outcomes.", "tier1;tier2", "DDA;DIA", "independent"),
        c("C_reference_recovery", "C", "detection_abundance", "Controlled reference detection and recovery meet experimental gates.", "tier3", "DDA;DIA", "independent"),
        c("C_order_invariance", "C", "execution", "Untouched development data and named-order invariance gates pass.", "tier1;tier3", "DDA;DIA", "independent"),
        c("C_wording", "C", "communication", "Output separates differential detection, observed abundance, and causal absence.", "tier1;tier3;tier4", "DDA;DIA", "independent")
    )
    out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(out) <- c(
        "claim_id", "track", "component", "claim", "evidence_tiers",
        "required_strata", "disposition"
    )
    out$claim_id <- tolower(out$claim_id)
    out
}

.m12_empty_gate_registry <- function() {
    data.frame(
        gate_id = character(), registry_version = character(),
        claim_id = character(), track = character(), metric = character(),
        unit = character(), comparator = character(), strata = character(),
        sample_size = integer(), uncertainty_method = character(),
        threshold_operator = character(), threshold_value = numeric(),
        threshold_scale = character(), failure_treatment = character(),
        data_ids = character(), protocol_id = character(),
        data_hash = character(), protocol_hash = character(),
        stringsAsFactors = FALSE
    )
}

.m12_protocol_registry <- function() {
    data.frame(
        protocol_id = c(
            "m12_a_candidate_protocol_v1",
            "m12_b_perturbation_protocol_v1",
            "m12_c_candidate_protocol_v1"
        ),
        track = c("A", "B", "C"),
        purpose = c(
            "algebraic design core plus association candidate comparison",
            "separated perturbation families cutoff policy streams and display calibration",
            "detection plus conditional-observed-abundance estimator comparison"
        ),
        implementation_file = rep("dev/m12-candidate-protocol.R", 3L),
        protocol_hash = c(
            "e20bdcc8e040eed65457480be7ae2ce2542761165c3d27161b953f86f39e5edc",
            "57ffa1220c740dcb69508d3ad1b9beb05e7013c3a92af528293211e41a13f256",
            "58b6e3c6d7f35f5561549a3a52420e49b0f3eb296752353d1bb3c9d433cf9d96"
        ),
        result_state = rep("frozen_unrun", 3L),
        stringsAsFactors = FALSE
    )
}

.m12_internal_evidence <- function() {
    data.frame(
        data_id = c(
            "v1_normative_fixture",
            "v1_scientific_oracle",
            "v1_performance_10000x50"
        ),
        evidence_tier = rep("tier0", 3L),
        role = rep("compatibility", 3L),
        scope = c(
            "released-shape normative matrix/result fixture",
            "frozen M7 scientific protocol and routine oracle",
            "10,000 feature x 50 sample stable-call performance fixture"
        ),
        unit = c("fixture", "scenario", "benchmark repetition"),
        sample_size = c(1L, 13L, 7L),
        source_path = c(
            "tests/testthat/helper-fixtures.R",
            "dev/scientific-validation.R",
            "dev/performance-validation.R"
        ),
        content_sha256 = c(
            "856b62f23c41140f07ac7965f4afe62f7dc23e25036a7769306f8bdd4a731d4f",
            "5b6b2f861d867378fa950a7f79857711b1a8056f29b557505194090214d5dc47",
            "74c53ccaf0f0c400693ac76607c575f686012c5b3980ae4d758a5dc4893fc2c2"
        ),
        state = rep("frozen", 3L),
        stringsAsFactors = FALSE
    )
}

.m12_ids <- function(...) {
    paste(sort(unique(c(...)), method = "radix"), collapse = ";")
}

.m12_gate_specs <- function() {
    scenario_ids <- .m12_generator_manifest()$scenario_id
    all_sim <- .m12_ids(scenario_ids)
    full_design <- .m12_ids(
        "dda_batch_crossed", "dia_batch_partial", "dia_monotone_paired"
    )
    paired <- .m12_ids(
        "dda_grouped_leakage_trap", "dia_blockwise_paired",
        "dia_monotone_paired"
    )
    unequal <- .m12_ids(
        "dda_monotone_unequal", "dia_batch_partial", "dia_null_unequal"
    )
    nulls <- .m12_ids("dda_null_balanced", "dia_null_unequal")
    no_cliff <- .m12_ids("dda_no_cliff", "dia_no_cliff_low_support")
    assoc_alt <- .m12_ids(
        "dda_batch_crossed", "dda_monotone_unequal",
        "dia_batch_partial", "dia_monotone_paired"
    )
    c_detection_alt <- .m12_ids(
        "dda_mixed_outlier", "dda_monotone_unequal", "dia_batch_partial",
        "dia_monotone_paired", "dia_structural_off"
    )
    c_abundance_alt <- .m12_ids(
        "dda_grouped_leakage_trap", "dda_mixed_outlier",
        "dda_monotone_unequal", "dia_batch_partial",
        "dia_blockwise_paired", "dia_monotone_paired", "dia_structural_off"
    )
    g <- function(
        gate_id,
        claim_id,
        metric,
        unit,
        comparator,
        strata,
        sample_size,
        uncertainty_method,
        threshold_operator,
        threshold_value,
        threshold_scale,
        failure_treatment,
        data_ids
    ) {
        data.frame(
            gate_id = gate_id,
            claim_id = claim_id,
            metric = metric,
            unit = unit,
            comparator = comparator,
            strata = strata,
            sample_size = as.integer(sample_size),
            uncertainty_method = uncertainty_method,
            threshold_operator = threshold_operator,
            threshold_value = as.numeric(threshold_value),
            threshold_scale = threshold_scale,
            failure_treatment = failure_treatment,
            data_ids = data_ids,
            stringsAsFactors = FALSE
        )
    }
    exact <- "complete enumeration of frozen development-test scenario-replicates"
    mandatory <- "any miss blocks the mandatory design core and therefore M13-M15"
    kill_a <- "miss kills or parks the A association panel; mandatory design core remains"
    kill_b <- "miss kills or parks the affected B claim without weakening another panel"
    disable_b_display <- "miss disables stable/fragile labels; continuous B panels remain eligible and thresholds stay frozen"
    kill_c <- "miss kills or parks the affected C component/track; threshold stays frozen"
    rows <- list(
        g("a_rank_deficient_sensitivity", "a_design_rank", "rank_deficiency_detection_fraction", "scenario-replicate", "constructed SVD rank oracle", "DDA perfect confounding", 32L, exact, "==", 1, "fraction", mandatory, "dda_batch_perfect"),
        g("a_rank_full_specificity", "a_design_rank", "full_rank_retention_fraction", "scenario-replicate", "constructed full-rank SVD oracle", "DDA crossed; DIA partial; DIA paired", 96L, exact, "==", 1, "fraction", mandatory, full_design),
        g("a_alias_exact_sensitivity", "a_design_alias", "perfect_alias_detection_fraction", "scenario-replicate", "condition=batch constructed alias", "DDA perfect confounding", 32L, exact, "==", 1, "fraction", mandatory, "dda_batch_perfect"),
        g("a_alias_false_positive", "a_design_alias", "false_exact_alias_fraction", "scenario-replicate", "crossed/partial nonalias oracle", "DDA crossed; DIA partial", 64L, exact, "==", 0, "fraction", mandatory, .m12_ids("dda_batch_crossed", "dia_batch_partial")),
        g("a_contrast_nonestimable_rejection", "a_contrast_estimability", "nonestimable_rejection_fraction", "scenario-replicate", "row-space projection oracle", "DDA perfect confounding", 32L, exact, "==", 1, "fraction", mandatory, "dda_batch_perfect"),
        g("a_contrast_estimable_retention", "a_contrast_estimability", "estimable_retention_fraction", "scenario-replicate", "row-space projection oracle", "crossed; partial; paired", 96L, exact, "==", 1, "fraction", mandatory, full_design),
        g("a_block_unit_accounting", "a_block_accounting", "block_resampling_unit_match_fraction", "scenario-replicate", "generator subject/technical-sibling oracle", "DDA/DIA paired", 96L, exact, "==", 1, "fraction", mandatory, paired),
        g("a_unequal_design_accounting", "a_block_accounting", "design_row_and_weight_match_fraction", "scenario-replicate", "generator named unequal-count oracle", "DDA/DIA unequal", 96L, exact, "==", 1, "fraction", mandatory, unequal),
        g("a_named_order_invariance", "a_order_invariance", "canonical_named_hash_match_fraction", "scenario-replicate", "baseline versus row/column/label/re-encoding permutations", "all 13 DDA/DIA scenarios", 416L, exact, "==", 1, "fraction", mandatory, all_sim),
        g("a_side_effect_global", "a_side_effects", "caller_state_mutation_count", "scenario-replicate", "input/RNG/options/working-directory snapshots", "all 13 DDA/DIA scenarios", 416L, exact, "==", 0, "count", mandatory, all_sim),
        g("a_side_effect_classic", "a_side_effects", "classic_exact_hash_drift_count", "scenario-replicate", "stable result/failure before versus after sentinel", "all scenarios plus normative fixture", 417L, exact, "==", 0, "count", mandatory, .m12_ids(scenario_ids, "v1_normative_fixture")),
        g("a_assoc_null_upper", "a_assoc_null", "family_false_flag_clopper_pearson_upper95", "scenario-replicate association family", "Holm-adjusted alpha=0.05", "pooled DDA+DIA nulls; strata also reported", 64L, "one-sided exact binomial 95 percent bound", "<=", 0.15, "probability", kill_a, nulls),
        g("a_assoc_null_stratum", "a_assoc_null", "maximum_acquisition_false_flag_fraction", "scenario-replicate association family", "Holm-adjusted alpha=0.05", "DDA and DIA separately n=32", 64L, "exact binomial intervals reported per acquisition", "<=", 0.125, "fraction", kill_a, nulls),
        g("a_public_development_coverage", "a_assoc_public", "eligible_declared_term_completion_fraction", "mouse", "all algebraically estimable condition/batch terms return effect+interval+raw/adjusted evidence", "HarmonizR DDA four-batch development", 25L, "whole-family descriptive audit; no cross-dataset inference", ">=", 0.95, "fraction", kill_a, "harmonizr_mouse_pxd027467"),
        g("a_public_development_wording", "a_assoc_public", "causal_or_correction_wording_violation_count", "wording rubric item", "associated/aliased/cannot-separate vocabulary", "HarmonizR development", 12L, "exact adversarial wording rubric", "==", 0, "count", kill_a, "harmonizr_mouse_pxd027467"),
        g("a_public_confirmation_coverage", "a_assoc_public", "eligible_declared_term_completion_fraction", "biological-replicate family", "all estimable condition/instrument/acquisition terms return complete evidence", "MultiPro DDA+DIA sealed confirmation", 6L, "whole-family confirmation audit after synchronized opening", ">=", 0.95, "fraction", kill_a, "multipro_hcc_pxd041391"),
        g("a_public_confirmation_wording", "a_assoc_public", "causal_or_correction_wording_violation_count", "wording rubric item", "associated/aliased/cannot-separate vocabulary", "MultiPro sealed confirmation", 12L, "exact adversarial wording rubric", "==", 0, "count", kill_a, "multipro_hcc_pxd041391"),
        g("a_assoc_interval_coverage", "a_assoc_calibration", "median_replicate_95_interval_coverage", "scenario-replicate", "generator expected sample-detection-fraction contrast", "DDA/DIA condition/batch alternatives", 128L, "feature-module-aware scenario bootstrap interval", ">=", 0.90, "fraction", kill_a, assoc_alt),
        g("a_assoc_alternative_power", "a_assoc_calibration", "holm_alternative_detection_fraction", "scenario-replicate declared term", "nonzero generator expected detection-fraction term", "DDA/DIA condition/batch alternatives", 128L, "scenario-replicate clustered 95 percent interval", ">=", 0.70, "fraction", kill_a, assoc_alt),
        g("a_assoc_effect_bias", "a_assoc_calibration", "median_absolute_effect_bias", "scenario-replicate declared term", "generator expected detection-fraction contrast", "DDA/DIA condition/batch alternatives", 128L, "scenario-replicate percentile interval", "<=", 0.03, "probability points", kill_a, assoc_alt),

        g("b_cutoff_range_coverage", "b_cutoff_coverage", "type8_95_range_boundary_coverage", "scenario-replicate condition", "single-monotone generating midpoint", "DDA unequal; DIA paired", 64L, "scenario-replicate clustered 95 percent interval", ">=", 0.90, "fraction", kill_b, .m12_ids("dda_monotone_unequal", "dia_monotone_paired")),
        g("b_cutoff_median_error", "b_cutoff_coverage", "median_absolute_cutoff_error", "scenario-replicate condition", "single-monotone generating midpoint", "DDA unequal; DIA paired", 64L, "scenario-replicate percentile 95 percent interval", "<=", 0.35, "log2 intensity", kill_b, .m12_ids("dda_monotone_unequal", "dia_monotone_paired")),
        g("b_false_confidence_upper", "b_false_confidence", "false_confidence_clopper_pearson_upper95", "scenario-replicate", "frozen success/range-width event", "DDA/DIA no-cliff", 64L, "one-sided exact binomial 95 percent bound", "<=", 0.10, "probability", kill_b, no_cliff),
        g("b_no_cliff_abstention", "b_no_cliff", "weak_identification_or_failure_fraction", "scenario-replicate condition", "flat/no-cliff generator control", "pooled DDA+DIA", 64L, "exact binomial interval", ">=", 0.95, "fraction", kill_b, no_cliff),
        g("b_no_cliff_stratum", "b_no_cliff", "minimum_acquisition_abstention_fraction", "scenario-replicate condition", "flat/no-cliff generator control", "DDA and DIA separately n=32", 64L, "exact binomial intervals by acquisition", ">=", 0.90, "fraction", kill_b, no_cliff),
        g("b_sparse_correct_abstention", "b_sparse_abstention", "sampling_panel_unavailable_match_fraction", "scenario-replicate", "frozen low-support rule/code", "DIA low support", 32L, exact, "==", 1, "fraction", kill_b, "dia_no_cliff_low_support"),
        g("b_sparse_panel_independence", "b_sparse_abstention", "independent_panel_completion_fraction", "scenario-replicate panel", "unsupported sampling leaves other eligible panels measured", "DIA low support", 32L, exact, ">=", 0.95, "fraction", kill_b, "dia_no_cliff_low_support"),
        g("b_reproducibility_exact", "b_reproducibility", "canonical_perturbation_hash_match_fraction", "scenario-replicate", "instance/input/draw order and seed-stream permutations", "all 13 DDA/DIA scenarios", 416L, exact, "==", 1, "fraction", kill_b, all_sim),
        g("b_side_effects_exact", "b_side_effects", "caller_state_mutation_count", "scenario-replicate", "input/RNG/options/working-directory snapshots", "all 13 DDA/DIA scenarios", 416L, exact, "==", 0, "count", kill_b, all_sim),
        g("b_stable_zero_work", "b_stable_overhead", "stable_default_sidecar_invocation_count", "normative call", "instrumented stable classify_missingness default", "v1 normative fixture", 1L, "exact call-path counter", "==", 0, "count", kill_b, "v1_normative_fixture"),
        g("b_stable_runtime", "b_stable_overhead", "median_runtime_ratio", "benchmark repetition", "stable default after/before frozen M10 baseline", "10000x50 performance fixture", 7L, "paired repetition ratios with full distribution reported", "<=", 1.05, "ratio", kill_b, "v1_performance_10000x50"),
        g("b_external_development_schema", "b_external_behavior", "required_panel_field_completion_fraction", "DDA/DIA run", "all supported B continuous fields/failures/provenance present", "MS-DAP linked development family", 18L, "whole-family descriptive audit", "==", 1, "fraction", kill_b, "msdap_hy_pxd036134"),
        g("b_external_development_wording", "b_external_behavior", "inferential_cutoff_coverage_wording_violation_count", "wording rubric item", "public ranges labelled finite-dataset descriptive", "MS-DAP development", 12L, "exact adversarial wording rubric", "==", 0, "count", kill_b, "msdap_hy_pxd036134"),
        g("b_external_confirmation_schema", "b_external_behavior", "required_panel_field_completion_fraction", "selected technical run", "all supported B continuous fields/failures/provenance present", "UPS1 DDA + LFQbench DIA sealed confirmation", 18L, "synchronized whole-family confirmation audit", "==", 1, "fraction", kill_b, .m12_ids("lfqbench_pxd002952", "ups1_yeast_pxd002099")),
        g("b_external_confirmation_wording", "b_external_behavior", "inferential_cutoff_coverage_wording_violation_count", "wording rubric item", "public ranges labelled finite-dataset descriptive", "UPS1/LFQbench sealed confirmation", 12L, "exact adversarial wording rubric", "==", 0, "count", kill_b, .m12_ids("lfqbench_pxd002952", "ups1_yeast_pxd002099")),
        g("b_display_heldout_risk", "b_display_calibration", "maximum_family_wilson_upper95_instability_risk", "scenario-replicate", "held-out draws 500-999 after threshold learned on 1-499", "sampling and estimator; DDA/DIA", 416L, "Wilson one-sided 95 percent by family", "<=", 0.15, "probability", disable_b_display, all_sim),
        g("b_display_heldout_coverage", "b_display_calibration", "minimum_family_label_coverage", "supported feature-condition", "frozen risk-qualified threshold", "sampling and estimator; DDA/DIA", 416L, "feature modules clustered within scenario-replicate", ">=", 0.10, "fraction", disable_b_display, all_sim),
        g("b_display_no_cliff", "b_display_calibration", "false_stable_label_clopper_pearson_upper95", "scenario-replicate", "no-cliff negative control", "DDA/DIA no-cliff", 64L, "one-sided exact binomial 95 percent bound", "<=", 0.10, "probability", disable_b_display, no_cliff),

        g("c_null_fdr_upper", "c_null_error", "empirical_fdr_upper95", "scenario-replicate detection family", "BY q=0.05; all-null FDR equals family false-flag probability", "pooled correlated-feature DDA+DIA nulls", 64L, "one-sided exact binomial 95 percent bound", "<=", 0.10, "probability", kill_c, nulls),
        g("c_null_fdr_stratum", "c_null_error", "maximum_acquisition_empirical_fdr", "scenario-replicate detection family", "BY q=0.05", "DDA and DIA separately n=32", 64L, "exact binomial intervals by acquisition", "<=", 0.10, "fraction", kill_c, nulls),
        g("c_abundance_null_fdr_upper", "c_null_error", "empirical_abundance_fdr_upper95", "scenario-replicate observed-abundance family", "BY q=0.05; all-null FDR equals family false-flag probability", "pooled correlated-feature DDA+DIA nulls", 64L, "one-sided exact binomial 95 percent bound", "<=", 0.10, "probability", kill_c, nulls),
        g("c_abundance_null_fdr_stratum", "c_null_error", "maximum_acquisition_abundance_empirical_fdr", "scenario-replicate observed-abundance family", "BY q=0.05", "DDA and DIA separately n=32", 64L, "exact binomial intervals by acquisition", "<=", 0.10, "fraction", kill_c, nulls),
        g("c_false_sign_upper", "c_false_sign", "false_sign_proportion_upper95", "scenario-replicate discovery family", "nonzero generator standardized detection-effect sign", "five detection-contrast DDA/DIA alternatives", 160L, "feature-module clustered percentile 95 percent bound", "<=", 0.10, "fraction", kill_c, c_detection_alt),
        g("c_abundance_false_sign_upper", "c_false_sign", "abundance_false_sign_proportion_upper95", "scenario-replicate discovery family", "nonzero generator conditional-observed-abundance effect sign", "seven abundance-contrast DDA/DIA alternatives", 224L, "feature-module clustered percentile 95 percent bound", "<=", 0.10, "fraction", kill_c, c_abundance_alt),
        g("c_detection_effect_bias", "c_effect_bias", "median_absolute_detection_risk_difference_bias", "scenario-replicate true detection-alternative feature", "nonzero generator standardized marginal detection difference", "five detection-contrast DDA/DIA alternatives", 160L, "feature-module clustered percentile 95 percent interval", "<=", 0.075, "probability points", kill_c, c_detection_alt),
        g("c_abundance_effect_bias", "c_effect_bias", "median_absolute_observed_abundance_bias", "scenario-replicate true abundance-alternative feature", "nonzero generator conditional-observed-abundance contrast", "seven abundance-contrast DDA/DIA alternatives", 224L, "feature-module clustered percentile 95 percent interval", "<=", 0.25, "log2 intensity", kill_c, c_abundance_alt),
        g("c_detection_effect_rmse", "c_effect_bias", "detection_risk_difference_rmse", "scenario-replicate true detection-alternative feature", "nonzero generator standardized marginal detection difference", "five detection-contrast DDA/DIA alternatives", 160L, "feature-module clustered percentile 95 percent interval", "<=", 0.15, "probability points", kill_c, c_detection_alt),
        g("c_detection_interval_coverage", "c_interval_coverage", "median_replicate_detection_interval_coverage", "scenario-replicate eligible feature", "generator standardized marginal detection difference", "five detection-contrast DDA/DIA alternatives", 160L, "feature-module clustered percentile 95 percent interval", ">=", 0.90, "fraction", kill_c, c_detection_alt),
        g("c_abundance_interval_coverage", "c_interval_coverage", "median_replicate_abundance_interval_coverage", "scenario-replicate eligible feature", "generator conditional-observed-abundance contrast", "seven abundance-contrast DDA/DIA alternatives", 224L, "feature-module clustered percentile 95 percent interval", ">=", 0.90, "fraction", kill_c, c_abundance_alt),
        g("c_alternative_power", "c_alternative_utility", "by_detection_power", "scenario-replicate true detection-alternative feature", "BY q=0.05 and correct nonzero generator detection-effect sign", "five detection-contrast DDA/DIA alternatives", 160L, "feature-module clustered percentile 95 percent interval", ">=", 0.50, "fraction", kill_c, c_detection_alt),
        g("c_alternative_precision", "c_alternative_utility", "detection_discovery_precision", "scenario-replicate discovery", "nonzero generator standardized detection effect", "five detection-contrast DDA/DIA alternatives", 160L, "feature-module clustered percentile 95 percent interval", ">=", 0.90, "fraction", kill_c, c_detection_alt),
        g("c_abundance_alternative_power", "c_alternative_utility", "by_abundance_power", "scenario-replicate true abundance-alternative feature", "BY q=0.05 and correct nonzero generator conditional-abundance sign", "seven abundance-contrast DDA/DIA alternatives", 224L, "feature-module clustered percentile 95 percent interval", ">=", 0.50, "fraction", kill_c, c_abundance_alt),
        g("c_abundance_alternative_precision", "c_alternative_utility", "abundance_discovery_precision", "scenario-replicate discovery", "nonzero generator conditional-observed-abundance effect", "seven abundance-contrast DDA/DIA alternatives", 224L, "feature-module clustered percentile 95 percent interval", ">=", 0.90, "fraction", kill_c, c_abundance_alt),
        g("c_edge_paired", "c_edge_cases", "paired_supported_completion_fraction", "scenario-replicate eligible feature", "frozen paired support + finite/structured outcome", "three paired DDA/DIA scenarios", 96L, exact, ">=", 0.90, "fraction", kill_c, paired),
        g("c_edge_unequal", "c_edge_cases", "unequal_supported_completion_fraction", "scenario-replicate eligible feature", "frozen unequal support + finite/structured outcome", "three unequal DDA/DIA scenarios", 96L, exact, ">=", 0.90, "fraction", kill_c, unequal),
        g("c_edge_separation", "c_edge_cases", "separated_finite_or_unavailable_fraction", "scenario-replicate eligible feature", "finite selected estimator or named separation code", "DIA structural-off-compatible", 32L, exact, "==", 1, "fraction", kill_c, "dia_structural_off"),
        g("c_edge_low_support", "c_edge_cases", "low_support_code_match_fraction", "scenario-replicate", "frozen component-specific support/code", "DIA low support", 32L, exact, "==", 1, "fraction", kill_c, "dia_no_cliff_low_support"),
        g("c_edge_confounding", "c_edge_cases", "nonestimable_contrast_rejection_fraction", "scenario-replicate", "A row-space oracle before fitting", "DDA perfect confounding", 32L, exact, "==", 1, "fraction", kill_c, "dda_batch_perfect"),
        g("c_reference_development_direction", "c_reference_recovery", "known_direction_agreement_fraction", "acquisition x yeast-level contrast", "12.5/15.625/18.75 ng known ordering", "MS-DAP linked DDA+DIA development", 18L, "technical-family bootstrap; no biological generalization", ">=", 0.90, "fraction", kill_c, "msdap_hy_pxd036134"),
        g("c_reference_development_ratio", "c_reference_recovery", "median_absolute_log2_ratio_error", "acquisition x yeast-level contrast", "known MS-DAP mixture ratios", "MS-DAP linked DDA+DIA development", 18L, "technical-family percentile interval", "<=", 0.50, "log2 ratio", kill_c, "msdap_hy_pxd036134"),
        g("c_reference_confirmation_ups1", "c_reference_recovery", "median_absolute_log2_slope_error", "selected UPS1 concentration contrast", "2/4/10/25 fmol known DDA series", "UPS1 sealed confirmation", 12L, "technical-family percentile interval", "<=", 0.35, "log2 ratio", kill_c, "ups1_yeast_pxd002099"),
        g("c_reference_confirmation_lfqbench", "c_reference_recovery", "median_absolute_species_ratio_error", "species x A/B contrast", "HYE124 human 1:1 yeast 2:1 E_coli 1:4", "LFQbench DIA sealed confirmation", 6L, "technical-family percentile interval", "<=", 0.50, "log2 ratio", kill_c, "lfqbench_pxd002952"),
        g("c_reference_confirmation_direction", "c_reference_recovery", "known_direction_agreement_fraction", "controlled species/concentration contrast", "UPS1 and HYE124 known direction", "DDA+DIA sealed confirmation", 18L, "whole-family descriptive interval", ">=", 0.90, "fraction", kill_c, .m12_ids("lfqbench_pxd002952", "ups1_yeast_pxd002099")),
        g("c_order_invariance_synthetic", "c_order_invariance", "canonical_named_result_hash_match_fraction", "scenario-replicate contrast", "feature/sample/condition/term order permutations", "all 13 DDA/DIA scenarios", 416L, exact, "==", 1, "fraction", kill_c, all_sim),
        g("c_order_invariance_external", "c_order_invariance", "canonical_named_result_hash_match_fraction", "development technical run", "original versus named-order permutations", "MS-DAP linked DDA+DIA development", 18L, "exact whole-family comparison", "==", 1, "fraction", kill_c, "msdap_hy_pxd036134"),
        g("c_wording_rubric", "c_wording", "estimand_wording_violation_count", "wording rubric item", "detection/conditional-abundance/causal-absence separation", "synthetic plus development external outputs", 18L, "exact adversarial wording rubric", "==", 0, "count", kill_c, .m12_ids("dda_batch_perfect", "dia_structural_off", "msdap_hy_pxd036134"))
    )
    out <- do.call(rbind, rows)
    row.names(out) <- NULL
    out
}

.m12_evidence_hash <- function(
    data_ids,
    generator,
    generator_protocol,
    internal,
    cards,
    roles,
    artifacts
) {
    ids <- sort(.m12_tokens(data_ids), method = "radix")
    pieces <- vapply(ids, function(id) {
        if (id %in% generator$scenario_id) {
            frames <- list(
                generator[generator$scenario_id == id, , drop = FALSE],
                generator_protocol
            )
        } else if (id %in% internal$data_id) {
            frames <- list(internal[internal$data_id == id, , drop = FALSE])
        } else if (id %in% cards$dataset_id) {
            frames <- list(
                cards[cards$dataset_id == id, , drop = FALSE],
                roles[roles$dataset_id == id, , drop = FALSE],
                artifacts[artifacts$dataset_id == id, , drop = FALSE]
            )
        } else {
            .m12_fail("unknown evidence ID while hashing: ", id)
        }
        hashes <- vapply(frames, .m12_hash_frame, character(1L))
        unname(tools::sha256sum(bytes = charToRaw(enc2utf8(paste(
            seq_along(hashes), hashes, sep = "\t", collapse = "\n"
        )))))
    }, character(1L))
    unname(tools::sha256sum(bytes = charToRaw(enc2utf8(paste(
        ids, pieces, sep = "\t", collapse = "\n"
    )))))
}

.m12_gate_registry <- function(
    claims,
    protocols,
    generator,
    generator_protocol,
    internal,
    cards,
    roles,
    artifacts
) {
    specs <- .m12_gate_specs()
    track <- claims$track[match(specs$claim_id, claims$claim_id)]
    protocol_index <- match(track, protocols$track)
    out <- data.frame(
        gate_id = specs$gate_id,
        registry_version = rep("m12_gate_registry_v1", nrow(specs)),
        claim_id = specs$claim_id,
        track = track,
        metric = specs$metric,
        unit = specs$unit,
        comparator = specs$comparator,
        strata = specs$strata,
        sample_size = specs$sample_size,
        uncertainty_method = specs$uncertainty_method,
        threshold_operator = specs$threshold_operator,
        threshold_value = specs$threshold_value,
        threshold_scale = specs$threshold_scale,
        failure_treatment = specs$failure_treatment,
        data_ids = specs$data_ids,
        protocol_id = protocols$protocol_id[protocol_index],
        data_hash = vapply(specs$data_ids, .m12_evidence_hash, character(1L),
            generator = generator,
            generator_protocol = generator_protocol,
            internal = internal,
            cards = cards,
            roles = roles,
            artifacts = artifacts
        ),
        protocol_hash = protocols$protocol_hash[protocol_index],
        stringsAsFactors = FALSE
    )
    row.names(out) <- NULL
    out
}

.m12_data_roles <- function() {
    data.frame(
        dataset_id = c(
            "harmonizr_mouse_pxd027467", "multipro_hcc_pxd041391",
            "msdap_hy_pxd036134", "ups1_yeast_pxd002099",
            "lfqbench_pxd002952"
        ),
        protocol_id = rep("m12_data_role_protocol_v1", 5L),
        role = c(
            "development", "confirmation", "development",
            "confirmation", "confirmation"
        ),
        included_artifacts = c(
            paste(c(
                "harmonizr_mouse_metadata", "harmonizr_mouse_ffpe_2018",
                "harmonizr_mouse_ff_2018", "harmonizr_mouse_ffpe_2020",
                "harmonizr_mouse_ff_2020"
            ), collapse = ";"),
            "multipro_hcc_dda_fragpipe;multipro_hcc_dia_diann",
            "msdap_hy_dda_maxquant;msdap_hy_dia_diann",
            "ups1_yeast_nonnormalized",
            "lfqbench_hye124_6600_64var"
        ),
        included_conditions = c(
            "tumor;control", "HCC1806;HS578T",
            "yeast_12.5ng;yeast_15.625ng;yeast_18.75ng",
            "UPS1_2fmol;UPS1_4fmol;UPS1_10fmol;UPS1_25fmol",
            "HYE124_A;HYE124_B"
        ),
        excluded_conditions = c(
            "none", "related_PXD041421_family", "none", "UPS1_50fmol",
            "none"
        ),
        grouping_rule = c(
            "whole four-batch family; animals and batch siblings never split",
            "whole HCC1806/HS578T DDA+DIA family; biological, instrument, and technical siblings synchronized",
            "whole matched DDA+DIA concentration family; mixtures and technical replicates synchronized",
            "whole protein artifact under one role; selected concentration triplicates synchronized",
            "whole inventoried instrument/software family; master-mixture derivatives synchronized"
        ),
        analysis_unit = c(
            "mouse within declared preservation/timepoint batch",
            "cell-culture biological replicate",
            "prepared concentration mixture",
            "prepared concentration mixture",
            "master species mixture"
        ),
        linked_tracks = c("A", "A", "B;C", "B;C", "B;C"),
        linked_claims = c(
            "a_assoc_public",
            "a_assoc_public",
            paste(c(
                "b_external_behavior", "c_effect_bias",
                "c_interval_coverage", "c_alternative_utility",
                "c_reference_recovery"
            ), collapse = ";"),
            paste(c(
                "b_external_behavior", "c_effect_bias",
                "c_interval_coverage", "c_alternative_utility",
                "c_reference_recovery"
            ), collapse = ";"),
            paste(c(
                "b_external_behavior", "c_effect_bias",
                "c_interval_coverage", "c_alternative_utility",
                "c_reference_recovery"
            ), collapse = ";")
        ),
        opening_stage = c(
            "after linked numeric gates and A candidate protocol hashes freeze",
            "M15P only after A association becomes a confirmation candidate",
            "after linked numeric gates and B/C candidate protocol hashes freeze",
            "M15P only after linked B/C tracks become confirmation candidates",
            "M15P only after linked B/C tracks become confirmation candidates"
        ),
        licence_snapshot = c(
            rep("PRIDE project licence=CC0; checked 2026-07-18", 4L),
            "EMBL-EBI terms revised 2024-02-05; owner-rights caveat retained; checked 2026-07-18"
        ),
        stringsAsFactors = FALSE
    )
}

.m12_role_hashes <- function(roles = .m12_data_roles()) {
    stats::setNames(
        vapply(seq_len(nrow(roles)), function(index) {
            .m12_hash_frame(roles[index, , drop = FALSE])
        }, character(1L)),
        roles$dataset_id
    )
}

.m12_data_cards <- function() {
    role_hashes <- .m12_role_hashes()
    data.frame(
        dataset_id = c(
            "harmonizr_mouse_pxd027467", "multipro_hcc_pxd041391",
            "msdap_hy_pxd036134", "ups1_yeast_pxd002099",
            "lfqbench_pxd002952"
        ),
        card_version = rep("m12_data_card_v2", 5L),
        title = c(
            "HarmonizR mouse tumour/control four-batch DDA study",
            "MultiPro HCC1806/HS578T deliberate-batch benchmark",
            "MS-DAP HeLa/yeast matched DDA/DIA spike-in",
            "UPS1-in-yeast DDA concentration series",
            "LFQbench HYE124 DIA benchmark"
        ),
        accession = c(
            "PXD027467", "PXD041391", "PXD036134", "PXD002099",
            "PXD002952"
        ),
        repository = rep("PRIDE Archive", 5L),
        metadata_url = paste0(
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/",
            c(
                "PXD027467", "PXD041391", "PXD036134", "PXD002099",
                "PXD002952"
            )
        ),
        evidence_tiers = c("tier4", "tier3;tier4", "tier3", "tier3", "tier3"),
        acquisition = c("DDA", "DDA;DIA", "DDA;DIA", "DDA", "DIA"),
        deployment_target = c(
            "two-condition bulk LFQ DDA across four declared preservation/timepoint batches",
            "balanced two-class bulk LFQ with declared machine, biological, and technical replication",
            "three-condition bulk LFQ technical benchmark with matched acquisition modes",
            "four-level bulk LFQ DDA concentration benchmark using a pre-role <=4-level subset",
            "two-condition bulk LFQ DIA mixed-species benchmark"
        ),
        experimental_unit = c(
            "mouse",
            "cell-culture biological replicate",
            "prepared concentration mixture; available repeats are technical",
            "prepared concentration mixture; available repeats are technical",
            "master species mixture; available repeats are technical"
        ),
        related_sample_rule = c(
            "keep the complete four-batch family together; use mouse as the biological unit within batch",
            "keep every injection, machine, and acquisition output from one biological replicate together",
            "keep DDA/DIA runs and technical replicates from one mixture together; treat the whole family as one role",
            "keep triplicates and every selected concentration together; treat the whole artifact as one role",
            "keep instruments, software outputs, and triplicates derived from each A/B master mixture together"
        ),
        n_runs = as.integer(c(25, 72, 18, 15, 6)),
        n_conditions = as.integer(c(2, 2, 3, 4, 2)),
        condition_layout = c(
            "tumour/control across four DDA batches: n=6,5,7,7 animals with preservation and analysis-time variation",
            "2 cell lines x 3 biological replicates x 3 technical replicates x 2 machines x 2 acquisition modes",
            "3 yeast spike levels x 3 technical replicates x 2 acquisition modes",
            "selected 2/4/10/25 fmol UPS1 levels x 3 technical replicates; 50 fmol excluded before role assignment",
            "2 HYE mixtures x 3 technical replicates for the inventoried 6600/64-variable-window artifact"
        ),
        nuisance_roles = c(
            "preservation;analysis_timepoint;batch", "instrument;technical_replicate",
            "acquisition_mode", "run_order", "instrument;software"
        ),
        block_roles = c(
            "mouse", "biological_replicate", "mixture", "concentration",
            "master_mixture"
        ),
        truth_scope = c(
            "declared tumour/control phenotype and four preservation/timepoint batches",
            "declared cell class, deliberate machine batch, and replicate hierarchy",
            "fixed human background plus known 12.5/15.625/18.75 ng yeast inputs",
            "fixed yeast background plus selected known 2/4/10/25 fmol UPS1 inputs",
            "known species composition and A:B abundance ratios"
        ),
        truth_limit = c(
            "natural missing cells and causal mechanisms are unlabeled; phenotype and technical effects remain associational",
            "individual natural missing cells and causal mechanisms are unlabeled; cell-line abundance is not a spike-in truth",
            "technical repeats do not identify biological-unit generalization or cell-level missingness mechanisms",
            "technical repeats do not identify biological-unit generalization; excluded 50 fmol level receives no claim",
            "technical repeats and shared master mixtures preclude biological-unit or software-independent confirmation claims"
        ),
        metadata_scope = c(
            "primary paper plus PRIDE project/file metadata; result archives unopened",
            "PRIDE project protocol plus file metadata; result archives unopened",
            "PRIDE project protocol plus file metadata; result archives unopened",
            "PRIDE project protocol plus file metadata; CSV unopened",
            "PRIDE project protocol plus file metadata; result archive unopened"
        ),
        licence_id = c(
            "CC0-1.0", "CC0-1.0", "CC0-1.0", "CC0-1.0",
            "EMBL-EBI-terms"
        ),
        licence_url = c(
            rep("https://creativecommons.org/publicdomain/zero/1.0/", 4L),
            "https://www.ebi.ac.uk/about/terms-of-use/"
        ),
        licence_status = c(
            rep("verified_open", 4L),
            "terms_reviewed"
        ),
        exclusions = c(
            "protein-level DDA outputs only; whole family is development-only; no causal missing-cell label",
            "protein-level outputs only; no peptide-hierarchy claim; PXD041421 is related and cannot supply independent confirmation",
            "protein-level outputs only; no random run split; DDA and DIA are linked artifacts",
            "exact <=4-level subset = 2/4/10/25 fmol; exclude 50 fmol; no peptide-level claim",
            "one software/instrument artifact only after role freeze; no causal label for naturally missing cells"
        ),
        split_feasibility = c(
            "whole_family_only", "grouped_unit_only", "whole_family_only",
            "whole_family_only", "whole_family_only"
        ),
        split_unit = c(
            "entire PXD027467 four-batch mouse family",
            "biological replicate across all technical siblings",
            "entire PXD036134 matched-acquisition family",
            "entire 2/4/10/25 fmol selected concentration artifact",
            "entire PXD002952 family of master-mixture derivatives"
        ),
        related_dataset_ids = c(
            "none", "PXD041421", "none", "none", "PXD020529"
        ),
        role = c(
            "development", "confirmation", "development", "confirmation",
            "confirmation"
        ),
        split_hash = unname(role_hashes[c(
            "harmonizr_mouse_pxd027467", "multipro_hcc_pxd041391",
            "msdap_hy_pxd036134", "ups1_yeast_pxd002099",
            "lfqbench_pxd002952"
        )]),
        checked_on = rep(.M12_CHECKED_ON, 5L),
        stringsAsFactors = FALSE
    )
}

.m12_artifact_manifest <- function() {
    accessions <- c(
        rep("PXD027467", 5L),
        "PXD041391", "PXD041391", "PXD036134", "PXD036134",
        "PXD002099", "PXD002952"
    )
    files <- c(
        "Mouse_Metadata.xlsx",
        "Mouse_FFPE_2018_SEARCH.zip",
        "Mouse_FF_2018_SEARCH.zip",
        "Mouse_FFPE_2020_SEARCH.zip",
        "Mouse_FF_2020_SEARCH.zip",
        "HCC1806_HS578T_DDA_FragPipe_Output.zip",
        "HCC1806_HS578T_DIA_DIANN_Output.zip",
        "timstofpro2_30SPD_two-proteome_spike-in_series_DDA.zip",
        "timstofpro2_30SPD_two-proteome_spike-in_series_DIA.zip",
        "YEAST-Data-NonNormalized.csv",
        "lfqbench_analysis_HYE124_TTOF6600_64var.zip"
    )
    data.frame(
        artifact_id = c(
            "harmonizr_mouse_metadata", "harmonizr_mouse_ffpe_2018",
            "harmonizr_mouse_ff_2018", "harmonizr_mouse_ffpe_2020",
            "harmonizr_mouse_ff_2020",
            "multipro_hcc_dda_fragpipe", "multipro_hcc_dia_diann",
            "msdap_hy_dda_maxquant", "msdap_hy_dia_diann",
            "ups1_yeast_nonnormalized", "lfqbench_hye124_6600_64var"
        ),
        dataset_id = c(
            rep("harmonizr_mouse_pxd027467", 5L),
            "multipro_hcc_pxd041391", "multipro_hcc_pxd041391",
            "msdap_hy_pxd036134", "msdap_hy_pxd036134",
            "ups1_yeast_pxd002099", "lfqbench_pxd002952"
        ),
        acquisition = c(
            rep("DDA", 5L), "DDA", "DIA", "DDA", "DIA", "DDA", "DIA"
        ),
        file_name = files,
        source_url = paste0(
            "https://www.ebi.ac.uk/pride/ws/archive-file-downloader/files/s3/",
            accessions,
            "/",
            files
        ),
        file_bytes = c(
            12361, 97768338, 168611708, 184092139, 304630271,
            2942612390, 3083089287, 101116244, 127676636, 273919, 263968636
        ),
        checksum_algorithm = rep("sha1", 11L),
        upstream_checksum = c(
            "8392a710eeeff2b801336697d2134ec74e85ea01",
            "d6788a617342ffc676dcc7a81ce371e40e67e3b5",
            "634314ff9fa1bf5eb4c53aa22fac91ecf863d4c4",
            "8379e4c92ce2bcdabcc6f8d264efe38681d1a604",
            "3fc9ee403270c9e6b147ee8ac279c3469f83d71d",
            "4e372285a6ee8c8301ee3f1e7ff9bc548e0ee4bd",
            "d027d5f89b2eea04f6c90b8189760856f70f039b",
            "9113637ddfdcd9c97d10ab0358bc918e948eee19",
            "f03cb375f06ba243f73d02580881ee44a9ccd734",
            "90e76ee3fa6238ba5325aac7d82bd2b8b6bcafb0",
            "a4420dd9f5c5833371b5623cebd6742d36888241"
        ),
        checksum_source = rep(
            paste0("PRIDE Archive v3 file record, checked ", .M12_CHECKED_ON),
            11L
        ),
        content_scope = c(
            "mouse sample and four-batch design metadata workbook",
            "2018 FFPE DDA search and protein-quantification output archive",
            "2018 fresh-frozen DDA search and protein-quantification output archive",
            "2020 FFPE DDA search and protein-quantification output archive",
            "2020 fresh-frozen DDA search and protein-quantification output archive",
            "FragPipe DDA search and protein-quantification output archive",
            "DIA-NN DIA search and protein-quantification output archive",
            "MaxQuant DDA processed spike-in output archive",
            "DIA-NN DIA processed spike-in output archive",
            "non-normalized protein quantitative CSV",
            "LFQbench analysis archive for HYE124 TripleTOF 6600 64-variable-window runs"
        ),
        role = c(
            rep("development", 5L), rep("confirmation", 2L),
            rep("development", 2L), rep("confirmation", 2L)
        ),
        open_state = rep("metadata_only", 11L),
        local_sha256 = rep(NA_character_, 11L),
        stringsAsFactors = FALSE
    )
}

.m12_candidate_exclusions <- function() {
    data.frame(
        candidate_id = c(
            "multipro_a549_pxd041421", "lfq_multiplatform_pxd028735",
            "ra_swath_btz898", "ionstar_pxd003881"
        ),
        accession = c("PXD041421", "PXD028735", "doi:10.1093/bioinformatics/btz898", "PXD003881"),
        reason_code = c(
            "related_family", "artifact_identity_unresolved",
            "artifact_and_licence_unresolved", "processed_matrix_unavailable"
        ),
        reason = c(
            "same MultiPro laboratory/protocol family as PXD041391 and inventoried processed archives lacked upstream checksums",
            "more than 700 cross-platform runs; no minimum processed protein artifact and checksummed subset was identifiable without opening data",
            "valuable blocked longitudinal DIA case, but no stable deidentified protein artifact/accession and data licence were identified",
            "strong controlled DDA truth, but PRIDE exposes raw/mzXML/mzIdentML rather than a compact deposited protein matrix"
        ),
        reconsider_when = c(
            "a non-independent exploratory role is needed and exact artifacts receive checksums",
            "platform-shift evidence is promoted and an exact <=4-condition artifact/protocol is frozen",
            "a stable deidentified protein matrix, subject/batch metadata, and reuse terms are verified",
            "a versioned protein-level derivation with immutable preprocessing and checksum is frozen"
        ),
        metadata_url = c(
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD041421",
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD028735",
            "https://doi.org/10.1093/bioinformatics/btz898",
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD003881"
        ),
        stringsAsFactors = FALSE
    )
}

.m12_source_manifest <- function() {
    data_sources <- data.frame(
        source_id = c(
            "pride_api", "pride_checksum", "ebi_terms", "harmonizr_paper",
            "harmonizr_record", "multipro_paper",
            "multipro_record", "msdap_paper", "msdap_record", "ups1_record",
            "ups1_paper", "lfqbench_paper", "lfqbench_record"
        ),
        scope = c(
            "repository metadata interface", "repository checksum algorithm",
            "reuse terms", "HarmonizR mouse batch design",
            "HarmonizR licence/files", "MultiPro design", "MultiPro licence/files",
            "MS-DAP design", "MS-DAP licence/files", "UPS1 design/licence/files",
            "UPS1 scientific truth", "LFQbench design/truth",
            "LFQbench licence/files"
        ),
        url = c(
            "https://www.ebi.ac.uk/pride/markdownpage/prideapi",
            "https://github.com/PRIDE-Archive/pride-checksum",
            "https://www.ebi.ac.uk/about/terms-of-use/",
            "https://doi.org/10.1038/s41467-022-31007-x",
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD027467",
            "https://doi.org/10.1038/s41597-023-02779-8",
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD041391",
            "https://doi.org/10.1021/acs.jproteome.2c00513",
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD036134",
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD002099",
            "https://doi.org/10.1021/acs.jproteome.5b00183",
            "https://doi.org/10.1038/nbt.3685",
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD002952"
        ),
        checked_on = rep(.M12_CHECKED_ON, 13L),
        note = c(
            "official project/file metadata API", "official SHA-1 checksum generator",
            "repository adds no owner-independent restrictions; attribution expected",
            "two-condition mouse DDA study across four preservation/timepoint batches",
            "CC0 record and four processed-batch archive identities",
            "balanced DDA/DIA deliberate-batch replicate design",
            "CC0 record and processed archive identities",
            "matched DDA/DIA three-level spike-in design",
            "CC0 record and processed archive identities",
            "CC0 five-level UPS1 DDA design and compact CSV identity",
            "controlled UPS1/background interpretation",
            "mixed-species DIA reference ratios and replicate design",
            "EBI-terms record and checksummed analysis archive"
        ),
        stringsAsFactors = FALSE
    )
    method_sources <- data.frame(
        source_id = c(
            "r_p_adjust", "r_quantile", "brglm2_current",
            "bias_reduction_2020", "survival_clogit", "geepack_current",
            "gee_small_sample_md", "limma_current", "limma_manual",
            "permutation_glm_2014", "cluster_bootstrap_2013",
            "cluster_cr2_2018",
            "selective_risk_coverage_2010"
        ),
        scope = c(
            "Holm/BY multiplicity", "sample quantile definition",
            "finite bias-reduced binomial implementation",
            "adjusted-score bias reduction", "conditional logistic likelihood",
            "clustered binary GEE implementation", "small-cluster GEE covariance",
            "omics linear-model implementation", "missing/blocked linear models",
            "nuisance-aware permutation", "cluster bootstrap",
            "small-sample cluster-robust fixed-effect inference",
            "risk-coverage display calibration"
        ),
        url = c(
            "https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html",
            "https://stat.ethz.ch/R-manual/R-devel/library/stats/html/quantile.html",
            "https://cran.r-project.org/package=brglm2",
            "https://doi.org/10.1007/s11222-019-09860-6",
            "https://stat.ethz.ch/R-manual/R-devel/library/survival/html/clogit.html",
            "https://cran.r-project.org/package=geepack",
            "https://doi.org/10.1111/j.0006-341X.2001.00126.x",
            "https://bioconductor.org/packages/release/bioc/html/limma.html",
            "https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf",
            "https://doi.org/10.1016/j.neuroimage.2014.01.060",
            "https://doi.org/10.1016/j.jmva.2012.10.006",
            "https://doi.org/10.1080/07350015.2016.1247004",
            "https://www.jmlr.org/papers/v11/el-yaniv10a.html"
        ),
        checked_on = rep(.M12_CHECKED_ON, 13L),
        note = c(
            "Holm strong FWER and BY arbitrary-dependence FDR definitions",
            "type 8 approximately median-unbiased and Hyndman-Fan recommended",
            "mean bias reduction returns finite full-rank logistic estimates",
            "adjusted-score framework and frequency properties",
            "exact conditional likelihood with explicit large-tie limits",
            "current GEE solver; standard sandwich remains asymptotic",
            "small-sample robust covariance can inflate type-I error; MD correction",
            "current maintained Bioconductor linear-model engine",
            "lmFit accepts missing values and declared block/correlation structures",
            "permutation validity depends on exchangeability/restricted design",
            "resampling subjects preserves within-subject dependence",
            "CR2 with Satterthwaite tests corrects small-cluster fixed-effect inference",
            "abstention labels require explicit risk-coverage tradeoff"
        ),
        stringsAsFactors = FALSE
    )
    rbind(data_sources, method_sources)
}

m12_validation_contract <- function() {
    generator <- .m12_generator_manifest()
    generator_protocol <- .m12_generator_protocol()
    claims <- .m12_claim_inventory()
    protocols <- .m12_protocol_registry()
    internal <- .m12_internal_evidence()
    roles <- .m12_data_roles()
    cards <- .m12_data_cards()
    artifacts <- .m12_artifact_manifest()
    gates <- .m12_gate_registry(
        claims, protocols, generator, generator_protocol, internal,
        cards, roles, artifacts
    )
    list(
        schema_catalog = .m12_schema_catalog(),
        generator_manifest = generator,
        generator_protocol = generator_protocol,
        claim_inventory = claims,
        protocol_registry = protocols,
        internal_evidence = internal,
        gate_registry = gates,
        data_cards = cards,
        data_roles = roles,
        artifact_manifest = artifacts,
        candidate_exclusions = .m12_candidate_exclusions(),
        source_manifest = .m12_source_manifest()
    )
}

.m12_fail <- function(...) {
    stop(..., call. = FALSE)
}

.m12_nonempty <- function(x) {
    is.character(x) && !anyNA(x) && all(nzchar(x))
}

.m12_id <- function(x) {
    .m12_nonempty(x) && all(grepl("^[a-z0-9][a-z0-9_]*$", x))
}

.m12_sha <- function(x, n) {
    is.character(x) && !anyNA(x) && all(grepl(
        paste0("^[0-9a-f]{", n, "}$"),
        x
    ))
}

.m12_tokens <- function(x) {
    unique(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
}

.m12_require_columns <- function(x, schema, catalog) {
    expected <- catalog$field[catalog$schema == schema]
    if (!identical(names(x), expected)) {
        .m12_fail(schema, " columns/order differ from schema")
    }
    if (any(vapply(x, is.factor, logical(1L)))) {
        .m12_fail(schema, " must not contain factors")
    }
}

.m12_validate_generator <- function(x, catalog) {
    .m12_require_columns(x, "generator_manifest", catalog)
    if (!.m12_id(x$scenario_id) || anyDuplicated(x$scenario_id)) {
        .m12_fail("generator scenario IDs must be unique canonical IDs")
    }
    if (!.m12_id(x$generator_version) || !.m12_id(x$acquisition_profile) ||
        !.m12_id(x$seed_stream)) {
        .m12_fail("generator version/profile/seed IDs must be canonical")
    }
    if (!all(x$evidence_tier %in% c("tier1", "tier2")) ||
        !all(x$acquisition %in% c("DDA", "DIA")) ||
        !all(x$sampling_design %in% c("independent", "paired")) ||
        !all(x$block_role %in% c("none", "subject")) ||
        !all(x$confounding %in% c("none", "partial", "perfect")) ||
        !all(x$evidence_support %in% c("adequate", "low")) ||
        !all(x$implementation_state == "implemented_verified")) {
        .m12_fail("generator manifest contains an invalid vocabulary value")
    }
    if (!is.integer(x$n_features) || any(x$n_features <= 0L) ||
        !is.integer(x$n_conditions) || any(x$n_conditions < 2L) ||
        any(x$n_conditions > 4L)) {
        .m12_fail("generator dimensions leave Section 9.1 scope")
    }
    parsed_sizes <- strsplit(x$condition_sizes, ";", fixed = TRUE)
    size_ok <- vapply(seq_along(parsed_sizes), function(index) {
        pieces <- strsplit(parsed_sizes[[index]], "=", fixed = TRUE)
        length(pieces) == x$n_conditions[[index]] &&
            all(lengths(pieces) == 2L) &&
            !anyDuplicated(vapply(pieces, `[[`, character(1L), 1L)) &&
            all(grepl("^[1-9][0-9]*$", vapply(
                pieces,
                `[[`,
                character(1L),
                2L
            )))
    }, logical(1L))
    if (!all(size_ok)) {
        .m12_fail("condition_sizes must match n_conditions")
    }
    patterns <- unique(unlist(lapply(x$missingness_patterns, .m12_tokens)))
    required_patterns <- c(
        "intensity_independent", "monotone_abundance", "mixed", "blockwise",
        "batch_associated", "structural_off_compatible", "no_cliff"
    )
    if (!all(required_patterns %in% patterns)) {
        .m12_fail("generator manifest does not cover every Section 9.1 pattern")
    }
    stresses <- unique(unlist(lapply(x$intensity_stress, .m12_tokens)))
    required_stress <- c("outliers", "heteroscedastic", "differing_support")
    if (!all(required_stress %in% stresses) ||
        !all(c("DDA", "DIA") %in% x$acquisition) ||
        !all(c("independent", "paired") %in% x$sampling_design) ||
        !all(c("partial", "perfect") %in% x$confounding)) {
        .m12_fail("generator manifest does not cover required design/stress axes")
    }
    invisible(TRUE)
}

.m12_validate_generator_protocol <- function(x, generator, catalog) {
    .m12_require_columns(x, "generator_protocol", catalog)
    if (nrow(x) != 1L || !.m12_id(x$protocol_id) ||
        !.m12_id(x$generator_version) || !.m12_id(x$seed_stream) ||
        !identical(x$generator_version, unique(generator$generator_version)) ||
        !identical(x$seed_stream, unique(generator$seed_stream)) ||
        !is.integer(x$development_replicates) ||
        x$development_replicates <= 0L ||
        !identical(x$implementation_file, "dev/m12-generator-validation.R") ||
        !.m12_sha(x$protocol_hash, 64L) ||
        !.m12_sha(x$generator_audit_hash, 64L) ||
        !.m12_sha(x$stable_v1_audit_hash, 64L) ||
        !identical(x$implementation_state, "implemented_verified")) {
        .m12_fail("generator protocol linkage is incomplete or malformed")
    }
    invisible(TRUE)
}

.m12_validate_claims <- function(x, catalog) {
    .m12_require_columns(x, "claim_inventory", catalog)
    if (!.m12_id(x$claim_id) || anyDuplicated(x$claim_id) ||
        !all(x$track %in% c("A", "B", "C")) ||
        !.m12_id(x$component) || !.m12_nonempty(x$claim) ||
        !.m12_nonempty(x$evidence_tiers) ||
        !.m12_nonempty(x$required_strata) ||
        !all(x$disposition %in% c("mandatory", "independent"))) {
        .m12_fail("claim inventory violates its schema")
    }
    if (!all(c("A", "B", "C") %in% x$track)) {
        .m12_fail("claim inventory must cover A-C")
    }
    invisible(TRUE)
}

.m12_validate_protocols <- function(x, catalog) {
    .m12_require_columns(x, "protocol_registry", catalog)
    if (nrow(x) != 3L || !.m12_id(x$protocol_id) ||
        anyDuplicated(x$protocol_id) || !identical(x$track, c("A", "B", "C")) ||
        !.m12_nonempty(x$purpose) ||
        !all(x$implementation_file == "dev/m12-candidate-protocol.R") ||
        !.m12_sha(x$protocol_hash, 64L) ||
        !all(x$result_state == "frozen_unrun")) {
        .m12_fail("candidate protocol registry is malformed or results are open")
    }
    invisible(TRUE)
}

.m12_validate_internal <- function(x, catalog) {
    .m12_require_columns(x, "internal_evidence", catalog)
    if (!.m12_id(x$data_id) || anyDuplicated(x$data_id) ||
        !all(x$evidence_tier == "tier0") ||
        !all(x$role == "compatibility") || !.m12_nonempty(x$scope) ||
        !.m12_nonempty(x$unit) || !is.integer(x$sample_size) ||
        any(x$sample_size <= 0L) || !.m12_nonempty(x$source_path) ||
        !.m12_sha(x$content_sha256, 64L) || !all(x$state == "frozen")) {
        .m12_fail("internal evidence registry is malformed")
    }
    if (any(!file.exists(x$source_path))) {
        .m12_fail("internal evidence source path is absent")
    }
    live <- unname(tools::sha256sum(x$source_path))
    if (!identical(live, x$content_sha256)) {
        .m12_fail("internal evidence content hash drifted")
    }
    invisible(TRUE)
}

.m12_validate_gates <- function(
    x,
    claims,
    protocols,
    generator,
    generator_protocol,
    internal,
    cards,
    roles,
    artifacts,
    catalog
) {
    .m12_require_columns(x, "gate_registry", catalog)
    if (!nrow(x)) {
        .m12_fail("numeric gate registry must be populated before candidates run")
    }
    character_fields <- setdiff(names(x), c("sample_size", "threshold_value"))
    if (!all(vapply(x[character_fields], .m12_nonempty, logical(1L))) ||
        !.m12_id(x$gate_id) || anyDuplicated(x$gate_id) ||
        !.m12_id(x$registry_version) || !.m12_id(x$metric) ||
        !.m12_id(x$protocol_id) || !all(x$track %in% c("A", "B", "C")) ||
        !all(x$claim_id %in% claims$claim_id) ||
        !all(x$track == claims$track[match(x$claim_id, claims$claim_id)]) ||
        !is.integer(x$sample_size) || any(x$sample_size <= 0L) ||
        !is.numeric(x$threshold_value) || any(!is.finite(x$threshold_value)) ||
        !all(x$threshold_operator %in% c("<=", ">=", "==", "<", ">")) ||
        !.m12_sha(x$data_hash, 64L) || !.m12_sha(x$protocol_hash, 64L)) {
        .m12_fail("gate registry accepts only complete, numeric, hash-frozen rows")
    }
    if (!setequal(unique(x$claim_id), claims$claim_id)) {
        .m12_fail("every retained A-C claim needs at least one numeric gate")
    }
    protocol_index <- match(x$protocol_id, protocols$protocol_id)
    if (anyNA(protocol_index) ||
        any(x$track != protocols$track[protocol_index]) ||
        any(x$protocol_hash != protocols$protocol_hash[protocol_index])) {
        .m12_fail("gate rows do not match their frozen track protocol")
    }
    if (any(grepl("tbd|pending|provisional", x$failure_treatment,
                  ignore.case = TRUE))) {
        .m12_fail("gate failure treatment must be final and explicit")
    }
    token_lists <- lapply(x$data_ids, .m12_tokens)
    canonical <- vapply(token_lists, function(ids) {
        identical(ids, sort(unique(ids), method = "radix"))
    }, logical(1L))
    if (!all(canonical)) {
        .m12_fail("gate data IDs must be unique canonical sorted token sets")
    }
    data_ids <- c(generator$scenario_id, internal$data_id, cards$dataset_id)
    referenced <- unique(unlist(token_lists))
    if (!all(referenced %in% data_ids)) {
        .m12_fail("gate registry references unknown generator/dataset IDs")
    }
    expected_data_hash <- unname(vapply(x$data_ids, .m12_evidence_hash, character(1L),
        generator = generator,
        generator_protocol = generator_protocol,
        internal = internal,
        cards = cards,
        roles = roles,
        artifacts = artifacts
    ))
    if (!identical(x$data_hash, expected_data_hash)) {
        .m12_fail("gate data hashes do not match frozen evidence bindings")
    }
    invisible(TRUE)
}

.m12_validate_cards <- function(x, catalog) {
    .m12_require_columns(x, "data_cards", catalog)
    chars <- setdiff(names(x), c("n_runs", "n_conditions", "split_hash"))
    if (!all(vapply(x[chars], .m12_nonempty, logical(1L))) ||
        !.m12_id(x$dataset_id) || anyDuplicated(x$dataset_id) ||
        !.m12_id(x$card_version) ||
        !all(grepl("^https://", x$metadata_url)) ||
        !is.integer(x$n_runs) || any(x$n_runs <= 0L) ||
        !is.integer(x$n_conditions) || any(x$n_conditions < 2L) ||
        any(x$n_conditions > 5L) ||
        !all(x$licence_status %in% c("verified_open", "terms_reviewed")) ||
        !all(x$split_feasibility %in% c("whole_family_only", "grouped_unit_only")) ||
        !all(x$role %in% c("unassigned", "training", "development", "confirmation"))) {
        .m12_fail("data cards violate field vocabularies/types")
    }
    assigned <- x$role != "unassigned"
    if (any(!assigned & !is.na(x$split_hash)) ||
        any(assigned & !vapply(x$split_hash, function(value) {
            !is.na(value) && grepl("^[0-9a-f]{64}$", value)
        }, logical(1L)))) {
        .m12_fail("split hashes must be absent before role assignment and frozen after")
    }
    out_of_scope <- x$n_conditions > 4L
    if (any(out_of_scope & !grepl("<=4", x$exclusions[out_of_scope], fixed = TRUE))) {
        .m12_fail("candidate artifacts above four conditions need an explicit subset rail")
    }
    invisible(TRUE)
}

.m12_validate_roles <- function(x, cards, artifacts, claims, catalog) {
    .m12_require_columns(x, "data_roles", catalog)
    if (nrow(x) != nrow(cards) || !identical(x$dataset_id, cards$dataset_id) ||
        !.m12_id(x$dataset_id) || anyDuplicated(x$dataset_id) ||
        !.m12_id(x$protocol_id) ||
        !all(x$role %in% c("development", "confirmation")) ||
        !all(x$role == cards$role) ||
        !all(vapply(x, .m12_nonempty, logical(1L)))) {
        .m12_fail("data-role protocol violates identity/assignment rails")
    }
    artifact_references <- unlist(lapply(x$included_artifacts, .m12_tokens))
    artifact_ids <- unique(artifact_references)
    claim_ids <- unique(unlist(lapply(x$linked_claims, .m12_tokens)))
    tracks <- unique(unlist(lapply(x$linked_tracks, .m12_tokens)))
    if (anyDuplicated(artifact_references) ||
        length(artifact_references) != nrow(artifacts) ||
        !setequal(artifact_ids, artifacts$artifact_id) ||
        !all(claim_ids %in% claims$claim_id) ||
        !all(tracks %in% c("A", "B", "C"))) {
        .m12_fail("data-role protocol references unknown/omitted artifacts or claims")
    }
    for (index in seq_len(nrow(x))) {
        row_artifacts <- .m12_tokens(x$included_artifacts[[index]])
        row_claims <- .m12_tokens(x$linked_claims[[index]])
        row_tracks <- .m12_tokens(x$linked_tracks[[index]])
        artifact_rows <- artifacts$artifact_id %in% row_artifacts
        if (!all(artifacts$dataset_id[artifact_rows] == x$dataset_id[[index]]) ||
            !all(artifacts$role[artifact_rows] == x$role[[index]]) ||
            !all(claims$track[match(row_claims, claims$claim_id)] %in% row_tracks) ||
            !identical(
                unname(cards$split_hash[[index]]),
                unname(.m12_hash_frame(x[index, , drop = FALSE]))
            )) {
            .m12_fail("data-role artifact/split hash linkage is inconsistent")
        }
    }
    ups1 <- x[x$dataset_id == "ups1_yeast_pxd002099", , drop = FALSE]
    if (nrow(ups1) != 1L ||
        !identical(
            ups1$included_conditions,
            "UPS1_2fmol;UPS1_4fmol;UPS1_10fmol;UPS1_25fmol"
        ) || !identical(ups1$excluded_conditions, "UPS1_50fmol")) {
        .m12_fail("UPS1 role must freeze the exact four-level subset")
    }
    role_track <- paste(
        rep(x$role, lengths(lapply(x$linked_tracks, .m12_tokens))),
        unlist(lapply(x$linked_tracks, .m12_tokens)),
        sep = "\r"
    )
    required_role_track <- as.vector(outer(
        c("development", "confirmation"),
        c("A", "B", "C"),
        paste,
        sep = "\r"
    ))
    if (!all(required_role_track %in% role_track)) {
        .m12_fail("data roles must cover development + confirmation for A-C")
    }
    invisible(TRUE)
}

.m12_validate_artifacts <- function(x, cards, catalog) {
    .m12_require_columns(x, "artifact_manifest", catalog)
    chars <- setdiff(names(x), c("file_bytes", "local_sha256"))
    if (!all(vapply(x[chars], .m12_nonempty, logical(1L))) ||
        !.m12_id(x$artifact_id) || anyDuplicated(x$artifact_id) ||
        !all(x$dataset_id %in% cards$dataset_id) ||
        !all(x$acquisition %in% c("DDA", "DIA")) ||
        !all(grepl("^https://", x$source_url)) ||
        !is.numeric(x$file_bytes) || any(!is.finite(x$file_bytes)) ||
        any(x$file_bytes <= 0) || any(x$file_bytes != floor(x$file_bytes)) ||
        !all(x$checksum_algorithm == "sha1") ||
        !.m12_sha(x$upstream_checksum, 40L) ||
        !all(x$role %in% c("unassigned", "training", "development", "confirmation")) ||
        !all(x$open_state %in% c("metadata_only", "opened"))) {
        .m12_fail("artifact manifest violates identity/type vocabularies")
    }
    card_role <- cards$role[match(x$dataset_id, cards$dataset_id)]
    if (any(x$role != card_role)) {
        .m12_fail("artifact roles must equal their data-card roles")
    }
    opened <- x$open_state == "opened"
    local_ok <- vapply(x$local_sha256, function(value) {
        !is.na(value) && grepl("^[0-9a-f]{64}$", value)
    }, logical(1L))
    if (any(opened & (x$role == "unassigned" | !local_ok)) ||
        any(!opened & !is.na(x$local_sha256))) {
        .m12_fail("opening requires assigned role + local SHA-256; metadata-only has none")
    }
    card_acquisition <- strsplit(cards$acquisition, ";", fixed = TRUE)
    acquisition_ok <- vapply(seq_len(nrow(x)), function(index) {
        x$acquisition[[index]] %in% card_acquisition[[match(
            x$dataset_id[[index]],
            cards$dataset_id
        )]]
    }, logical(1L))
    if (!all(acquisition_ok) || !all(cards$dataset_id %in% x$dataset_id)) {
        .m12_fail("artifact/card acquisition coverage is inconsistent")
    }
    invisible(TRUE)
}

.m12_validate_simple <- function(x, schema, catalog) {
    .m12_require_columns(x, schema, catalog)
    id_field <- if (schema == "candidate_exclusions") "candidate_id" else "source_id"
    if (!.m12_id(x[[id_field]]) || anyDuplicated(x[[id_field]]) ||
        !all(vapply(x, .m12_nonempty, logical(1L)))) {
        .m12_fail(schema, " violates its nonempty/unique-ID schema")
    }
    url_field <- if (schema == "candidate_exclusions") "metadata_url" else "url"
    if (!all(grepl("^https://", x[[url_field]]))) {
        .m12_fail(schema, " URLs must use HTTPS")
    }
    invisible(TRUE)
}

m12_validate_contract <- function(contract = m12_validation_contract()) {
    expected_names <- c(
        "schema_catalog", "generator_manifest", "generator_protocol",
        "claim_inventory", "protocol_registry", "internal_evidence",
        "gate_registry", "data_cards", "data_roles", "artifact_manifest",
        "candidate_exclusions", "source_manifest"
    )
    if (!identical(names(contract), expected_names)) {
        .m12_fail("contract components/order differ from v3")
    }
    catalog <- contract$schema_catalog
    if (!identical(
        names(catalog),
        c("schema", "position", "field", "rule")
    ) || !.m12_nonempty(catalog$schema) || !is.integer(catalog$position) ||
        !.m12_nonempty(catalog$field) || !.m12_nonempty(catalog$rule) ||
        anyDuplicated(paste(catalog$schema, catalog$field, sep = "\r"))) {
        .m12_fail("schema catalog is malformed")
    }
    .m12_validate_generator(contract$generator_manifest, catalog)
    .m12_validate_generator_protocol(
        contract$generator_protocol,
        contract$generator_manifest,
        catalog
    )
    .m12_validate_claims(contract$claim_inventory, catalog)
    .m12_validate_protocols(contract$protocol_registry, catalog)
    .m12_validate_internal(contract$internal_evidence, catalog)
    .m12_validate_cards(contract$data_cards, catalog)
    .m12_validate_artifacts(contract$artifact_manifest, contract$data_cards, catalog)
    .m12_validate_roles(
        contract$data_roles,
        contract$data_cards,
        contract$artifact_manifest,
        contract$claim_inventory,
        catalog
    )
    .m12_validate_gates(
        contract$gate_registry,
        contract$claim_inventory,
        contract$protocol_registry,
        contract$generator_manifest,
        contract$generator_protocol,
        contract$internal_evidence,
        contract$data_cards,
        contract$data_roles,
        contract$artifact_manifest,
        catalog
    )
    .m12_validate_simple(contract$candidate_exclusions, "candidate_exclusions", catalog)
    .m12_validate_simple(contract$source_manifest, "source_manifest", catalog)
    if (any(contract$data_cards$role == "unassigned") ||
        any(contract$artifact_manifest$role == "unassigned") ||
        !all(contract$artifact_manifest$open_state == "metadata_only") ||
        any(!is.na(contract$artifact_manifest$local_sha256)) ||
        !all(contract$protocol_registry$result_state == "frozen_unrun") ||
        !setequal(unique(contract$gate_registry$claim_id),
                  contract$claim_inventory$claim_id)) {
        .m12_fail("M12c must be gate-complete, metadata-only, and result-free")
    }
    invisible(TRUE)
}

m12_gate_binding_hashes <- function(
    gate_registry = m12_validation_contract()$gate_registry
) {
    frozen_registry <- m12_validation_contract()$gate_registry
    if (!identical(gate_registry, frozen_registry)) {
        .m12_fail("gate bindings require the canonical frozen registry")
    }
    stats::setNames(
        vapply(seq_len(nrow(gate_registry)), function(index) {
            .m12_hash_frame(gate_registry[index, , drop = FALSE])
        }, character(1L)),
        gate_registry$gate_id
    )
}

m12_evaluate_gate_results <- function(
    results,
    gate_registry = m12_validation_contract()$gate_registry,
    gate_ids = gate_registry$gate_id
) {
    frozen_registry <- m12_validation_contract()$gate_registry
    if (!identical(gate_registry, frozen_registry)) {
        .m12_fail("gate evaluation requires the canonical frozen registry")
    }
    catalog <- .m12_schema_catalog()
    .m12_require_columns(results, "gate_results", catalog)
    if (!length(gate_ids) || !.m12_id(gate_ids) || anyDuplicated(gate_ids) ||
        !all(gate_ids %in% gate_registry$gate_id)) {
        .m12_fail("requested gate IDs are unknown or duplicated")
    }
    registry <- gate_registry[match(gate_ids, gate_registry$gate_id), , drop = FALSE]
    if (nrow(results) != length(gate_ids) ||
        anyDuplicated(results$gate_id) || !setequal(results$gate_id, gate_ids)) {
        .m12_fail("gate results must contain every requested gate exactly once")
    }
    results <- results[match(gate_ids, results$gate_id), , drop = FALSE]
    row.names(results) <- NULL
    numeric_required <- c("estimate", "numerator", "denominator")
    measured <- results$status == "measured"
    numeric_values <- as.matrix(results[numeric_required])
    if (!all(vapply(results[numeric_required], is.numeric, logical(1L))) ||
        !all(results$status %in% c("measured", "failed", "unavailable")) ||
        any(!is.finite(numeric_values[measured, , drop = FALSE])) ||
        any(results$denominator[measured] <= 0) ||
        any(!is.na(numeric_values[!measured, , drop = FALSE])) ||
        !is.numeric(results$lower) ||
        !is.numeric(results$upper) ||
        any(!is.na(results$lower[measured]) &
            !is.finite(results$lower[measured])) ||
        any(!is.na(results$upper[measured]) &
            !is.finite(results$upper[measured])) ||
        any(!is.na(results$lower[measured]) &
            results$lower[measured] > results$estimate[measured]) ||
        any(!is.na(results$upper[measured]) &
            results$upper[measured] < results$estimate[measured]) ||
        any(!is.na(results$lower[!measured])) ||
        any(!is.na(results$upper[!measured])) ||
        !.m12_sha(results$gate_binding_hash, 64L) ||
        !.m12_sha(results$evidence_hash, 64L) || !.m12_nonempty(results$note)) {
        .m12_fail("gate result values/status/bounds violate the result schema")
    }
    bindings <- unname(m12_gate_binding_hashes(gate_registry)[gate_ids])
    if (any(results$registry_version != registry$registry_version) ||
        any(results$metric != registry$metric) ||
        any(results$gate_binding_hash != bindings)) {
        .m12_fail("gate results are detached from their frozen registry rows")
    }
    compare <- function(value, operator, threshold) {
        switch(operator,
            "<=" = value <= threshold,
            ">=" = value >= threshold,
            "==" = value == threshold,
            "<" = value < threshold,
            ">" = value > threshold,
            .m12_fail("unknown threshold operator: ", operator)
        )
    }
    numeric_pass <- rep(FALSE, nrow(results))
    numeric_pass[measured] <- mapply(
        compare,
        results$estimate[measured],
        registry$threshold_operator[measured],
        registry$threshold_value[measured],
        USE.NAMES = FALSE
    )
    passed <- measured & numeric_pass
    data.frame(
        gate_id = registry$gate_id,
        claim_id = registry$claim_id,
        track = registry$track,
        metric = registry$metric,
        estimate = results$estimate,
        operator = registry$threshold_operator,
        threshold = registry$threshold_value,
        status = results$status,
        passed = passed,
        failure_treatment = registry$failure_treatment,
        evidence_hash = results$evidence_hash,
        stringsAsFactors = FALSE
    )
}

.m12_escape <- function(x) {
    x <- enc2utf8(x)
    x <- gsub("%", "%25", x, fixed = TRUE)
    x <- gsub("\t", "%09", x, fixed = TRUE)
    x <- gsub("\r", "%0D", x, fixed = TRUE)
    gsub("\n", "%0A", x, fixed = TRUE)
}

.m12_canonical_value <- function(x) {
    if (is.integer(x)) {
        out <- as.character(x)
    } else if (is.numeric(x)) {
        out <- vapply(x, function(value) {
            if (is.na(value)) {
                return("<NA>")
            }
            format(value, digits = 17L, scientific = FALSE, trim = TRUE)
        }, character(1L))
    } else if (is.logical(x)) {
        out <- ifelse(is.na(x), "<NA>", ifelse(x, "TRUE", "FALSE"))
    } else {
        out <- x
    }
    out[is.na(out)] <- "<NA>"
    .m12_escape(out)
}

.m12_hash_frame <- function(x) {
    if (nrow(x)) {
        key <- x[[1L]]
        order_index <- order(key, method = "radix")
        x <- x[order_index, , drop = FALSE]
    }
    header <- paste(.m12_escape(names(x)), collapse = "\t")
    rows <- if (nrow(x)) {
        vapply(seq_len(nrow(x)), function(index) {
            paste(vapply(x, function(column) {
                .m12_canonical_value(column[[index]])
            }, character(1L)), collapse = "\t")
        }, character(1L))
    } else {
        character()
    }
    canonical <- paste(c(header, rows), collapse = "\n")
    unname(tools::sha256sum(bytes = charToRaw(enc2utf8(canonical))))
}

m12_contract_hashes <- function(contract = m12_validation_contract()) {
    m12_validate_contract(contract)
    component_hashes <- vapply(contract, .m12_hash_frame, character(1L))
    aggregate <- paste(
        names(component_hashes),
        component_hashes,
        sep = "\t",
        collapse = "\n"
    )
    c(
        component_hashes,
        contract = unname(tools::sha256sum(
            bytes = charToRaw(enc2utf8(aggregate))
        ))
    )
}

.M12_EXPECTED_HASHES <- c(
    schema_catalog = "3c117b427cda433a58548855d25ebe7b6bed0f235cbae7a12aab69a3b2e1689c",
    generator_manifest = "9d1381833cd59e8775bd99e730d5c783b98c5a82afff6d52276e83a46a8be76a",
    generator_protocol = "9f2166bc7a31112d5529d81adcd68cae70b99d12374f3d77bc45844d3d16bd7a",
    claim_inventory = "c31032238677acc5e57478057852481f29dd24bc3f6fb88d629653be5ba8d296",
    protocol_registry = "61f7c1d9c42e80e1dafb13693e7a1818c7835ee9fd219e8a953167fe0c84ab6d",
    internal_evidence = "7c2deb16729dc06720bea8ce14b1ed4024dcc51f836c03e00b3ccc7696464a1a",
    gate_registry = "12e1fb5f9886b4ca6a4accb18909b3aa40cf33ac220800d970d79d456d9fcd87",
    data_cards = "0ba7fa86286afd4cc72dfea6eda4b0bc1324553229c2a3b7c412d41c777d684c",
    data_roles = "500d08b89b44c3adde9782176ae1065bf20c2eb135b89e666c8623ec2c83665c",
    artifact_manifest = "7c3267020d484ad9bc1bf4b8ac442cf0916d554e9e467e5c696ab2a20fb1c652",
    candidate_exclusions = "7ef14dd357171ffaad2e3793dfd1e038d1cc4108d56d052a2f751e4061e71af8",
    source_manifest = "e413eb42bb5f173dde1f994664b2dd9b081eee304323cea3cf8ac390ac2651b4",
    contract = "2cd5ecd8bf5763da1c5d9e1d8994207e701bcea9a7efcf4ae539dfd9b52d7431"
)

.m12_expect_error <- function(code) {
    inherits(try(force(code), silent = TRUE), "try-error")
}

.m12_gate_evaluator_self_tests <- function(contract) {
    registry <- contract$gate_registry
    estimate <- registry$threshold_value
    estimate[registry$threshold_operator == "<"] <-
        estimate[registry$threshold_operator == "<"] - 1e-8
    estimate[registry$threshold_operator == ">"] <-
        estimate[registry$threshold_operator == ">"] + 1e-8
    results <- data.frame(
        gate_id = registry$gate_id,
        registry_version = registry$registry_version,
        metric = registry$metric,
        estimate = estimate,
        numerator = estimate,
        denominator = rep(1, nrow(registry)),
        lower = rep(NA_real_, nrow(registry)),
        upper = rep(NA_real_, nrow(registry)),
        status = rep("measured", nrow(registry)),
        gate_binding_hash = unname(m12_gate_binding_hashes(registry)),
        evidence_hash = rep(paste(rep("0", 64L), collapse = ""), nrow(registry)),
        note = rep("synthetic evaluator boundary self-test", nrow(registry)),
        stringsAsFactors = FALSE
    )
    baseline <- m12_evaluate_gate_results(results, registry)
    subset_ids <- registry$gate_id[1:2]
    subset_report <- m12_evaluate_gate_results(
        results[results$gate_id %in% subset_ids, , drop = FALSE],
        registry,
        subset_ids
    )
    omitted <- results[-1L, , drop = FALSE]
    detached <- results
    detached$gate_binding_hash[[1L]] <- paste0(
        "1", substr(detached$gate_binding_hash[[1L]], 2L, 64L)
    )
    unavailable <- results
    unavailable$status[[1L]] <- "unavailable"
    unavailable[1L, c(
        "estimate", "numerator", "denominator", "lower", "upper"
    )] <- NA_real_
    unavailable_report <- m12_evaluate_gate_results(unavailable, registry)
    fabricated_unavailable <- results
    fabricated_unavailable$status[[1L]] <- "unavailable"
    wrong_side <- results
    first <- 1L
    wrong_side$estimate[[first]] <- switch(
        registry$threshold_operator[[first]],
        "<=" = registry$threshold_value[[first]] + 1,
        ">=" = registry$threshold_value[[first]] - 1,
        "==" = registry$threshold_value[[first]] + 1,
        "<" = registry$threshold_value[[first]],
        ">" = registry$threshold_value[[first]]
    )
    wrong_report <- m12_evaluate_gate_results(wrong_side, registry)
    altered_registry <- registry
    altered_registry$threshold_value[[1L]] <-
        altered_registry$threshold_value[[1L]] + 1
    altered_results <- results
    altered_results$estimate[[1L]] <- altered_registry$threshold_value[[1L]]
    altered_results$numerator[[1L]] <- altered_results$estimate[[1L]]
    altered_results$gate_binding_hash <- vapply(
        seq_len(nrow(altered_registry)),
        function(index) {
            .m12_hash_frame(altered_registry[index, , drop = FALSE])
        },
        character(1L)
    )
    c(
        boundary_pass = all(baseline$passed),
        staged_subset_pass = all(subset_report$passed),
        empty_subset_rejected = .m12_expect_error(
            m12_evaluate_gate_results(
                results[FALSE, , drop = FALSE],
                registry,
                character()
            )
        ),
        omission_rejected = .m12_expect_error(
            m12_evaluate_gate_results(omitted, registry)
        ),
        binding_rejected = .m12_expect_error(
            m12_evaluate_gate_results(detached, registry)
        ),
        unavailable_fails = !unavailable_report$passed[[1L]],
        unavailable_endpoint_rejected = .m12_expect_error(
            m12_evaluate_gate_results(fabricated_unavailable, registry)
        ),
        wrong_side_fails = !wrong_report$passed[[1L]],
        altered_registry_rejected = .m12_expect_error(
            m12_evaluate_gate_results(altered_results, altered_registry)
        )
    )
}

.m12_self_tests <- function(contract) {
    out_of_scope <- contract
    out_of_scope$generator_manifest$n_conditions[[1L]] <- 5L
    missing_split <- contract
    missing_split$data_cards$split_hash[[1L]] <- NA_character_
    role_mismatch <- contract
    role_mismatch$artifact_manifest$role[[1L]] <- "confirmation"
    opened <- contract
    opened$artifact_manifest$open_state[[1L]] <- "opened"
    wrong_checksum <- contract
    wrong_checksum$artifact_manifest$upstream_checksum[[1L]] <- "bad"
    wrong_subset <- contract
    ups1 <- wrong_subset$data_roles$dataset_id == "ups1_yeast_pxd002099"
    wrong_subset$data_roles$included_conditions[ups1] <-
        "UPS1_2fmol;UPS1_4fmol;UPS1_10fmol;UPS1_50fmol"
    missing_claim_gate <- contract
    missing_claim_gate$gate_registry <- missing_claim_gate$gate_registry[
        missing_claim_gate$gate_registry$claim_id != "a_design_rank",
        ,
        drop = FALSE
    ]
    protocol_drift <- contract
    protocol_drift$gate_registry$protocol_hash[[1L]] <- paste0(
        "0", substr(protocol_drift$gate_registry$protocol_hash[[1L]], 2L, 64L)
    )
    data_drift <- contract
    data_drift$gate_registry$data_hash[[1L]] <- paste0(
        "0", substr(data_drift$gate_registry$data_hash[[1L]], 2L, 64L)
    )
    nonnumeric_gate <- contract
    nonnumeric_gate$gate_registry$threshold_value[[1L]] <- Inf
    results_open <- contract
    results_open$protocol_registry$result_state[[1L]] <- "results_opened"
    contract_tests <- c(
        generator_scope = .m12_expect_error(m12_validate_contract(out_of_scope)),
        role_requires_split = .m12_expect_error(m12_validate_contract(missing_split)),
        role_artifact_agreement = .m12_expect_error(m12_validate_contract(role_mismatch)),
        opening_requires_role_hash = .m12_expect_error(m12_validate_contract(opened)),
        checksum_shape = .m12_expect_error(m12_validate_contract(wrong_checksum)),
        exact_subset = .m12_expect_error(m12_validate_contract(wrong_subset)),
        every_claim_gated = .m12_expect_error(m12_validate_contract(missing_claim_gate)),
        protocol_hash_link = .m12_expect_error(m12_validate_contract(protocol_drift)),
        data_hash_link = .m12_expect_error(m12_validate_contract(data_drift)),
        numeric_threshold = .m12_expect_error(m12_validate_contract(nonnumeric_gate)),
        results_sealed = .m12_expect_error(m12_validate_contract(results_open))
    )
    c(contract_tests, .m12_gate_evaluator_self_tests(contract))
}

.m12_main <- function(args) {
    if (!identical(args, "--verify")) {
        .m12_fail("usage: Rscript --vanilla dev/m12-validation-contract.R --verify")
    }
    contract <- m12_validation_contract()
    hashes <- m12_contract_hashes(contract)
    self_tests <- .m12_self_tests(contract)
    if (!all(self_tests)) {
        .m12_fail(
            "M12 schema self-test failures: ",
            paste(names(self_tests)[!self_tests], collapse = ", ")
        )
    }
    frozen <- !anyNA(.M12_EXPECTED_HASHES)
    if (frozen && !identical(hashes, .M12_EXPECTED_HASHES)) {
        mismatch <- names(hashes)[hashes != .M12_EXPECTED_HASHES]
        .m12_fail("M12 frozen hash mismatch: ", paste(mismatch, collapse = ", "))
    }
    cat("contract_version: ", .M12_CONTRACT_VERSION, "\n", sep = "")
    cat("state: metadata_only; roles=assigned; protocols=frozen_unrun; gates=", nrow(contract$gate_registry), "\n", sep = "")
    cat("self_tests: ", paste(names(self_tests), self_tests, sep = "=", collapse = "; "), "\n", sep = "")
    cat(paste(names(hashes), hashes, sep = ": "), sep = "\n")
    cat("\n")
    if (!frozen) {
        cat("hash_state: bootstrap_unfrozen\n")
    }
    invisible(TRUE)
}

if (sys.nframe() == 0L) {
    .m12_main(commandArgs(trailingOnly = TRUE))
}
