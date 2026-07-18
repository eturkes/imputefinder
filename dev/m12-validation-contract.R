#!/usr/bin/env Rscript

# M12b validation contract. This file freezes scope, generator linkage, and
# whole-family evidence roles; it neither downloads nor reads result artifacts.
# Run from the repository root:
# Rscript --vanilla dev/m12-validation-contract.R --verify

.M12_CONTRACT_VERSION <- "m12_validation_contract_v2"
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
        c("C_null_error", "C", "detection", "Correlated-feature null rejection and FDR satisfy frozen bounds.", "tier1;tier2", "DDA;DIA", "independent"),
        c("C_false_sign", "C", "detection", "False-sign behavior satisfies its frozen bound.", "tier1;tier2", "DDA;DIA", "independent"),
        c("C_effect_bias", "C", "detection_abundance", "Detection and observed-abundance effect bias satisfy separate gates.", "tier1;tier2;tier3", "DDA;DIA", "independent"),
        c("C_interval_coverage", "C", "detection_abundance", "Intervals meet estimand-specific coverage gates where licensed.", "tier1;tier2;tier3", "DDA;DIA", "independent"),
        c("C_alternative_utility", "C", "detection", "Frozen alternatives meet power and precision-recall gates.", "tier1;tier2;tier3", "DDA;DIA", "independent"),
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
                "b_cutoff_coverage", "b_false_confidence", "c_effect_bias",
                "c_interval_coverage", "c_alternative_utility",
                "c_reference_recovery"
            ), collapse = ";"),
            paste(c(
                "b_cutoff_coverage", "b_false_confidence", "c_effect_bias",
                "c_interval_coverage", "c_alternative_utility",
                "c_reference_recovery"
            ), collapse = ";"),
            paste(c(
                "b_cutoff_coverage", "b_false_confidence", "c_effect_bias",
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
    data.frame(
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
}

m12_validation_contract <- function() {
    list(
        schema_catalog = .m12_schema_catalog(),
        generator_manifest = .m12_generator_manifest(),
        generator_protocol = .m12_generator_protocol(),
        claim_inventory = .m12_claim_inventory(),
        gate_registry = .m12_empty_gate_registry(),
        data_cards = .m12_data_cards(),
        data_roles = .m12_data_roles(),
        artifact_manifest = .m12_artifact_manifest(),
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

.m12_validate_gates <- function(x, claims, data_ids, catalog) {
    .m12_require_columns(x, "gate_registry", catalog)
    if (!nrow(x)) {
        return(invisible(TRUE))
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
    referenced <- unique(unlist(lapply(x$data_ids, .m12_tokens)))
    if (!all(referenced %in% data_ids)) {
        .m12_fail("gate registry references unknown generator/dataset IDs")
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
        "claim_inventory", "gate_registry", "data_cards", "data_roles",
        "artifact_manifest",
        "candidate_exclusions", "source_manifest"
    )
    if (!identical(names(contract), expected_names)) {
        .m12_fail("contract components/order differ from v1")
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
        c(contract$generator_manifest$scenario_id, contract$data_cards$dataset_id),
        catalog
    )
    .m12_validate_simple(contract$candidate_exclusions, "candidate_exclusions", catalog)
    .m12_validate_simple(contract$source_manifest, "source_manifest", catalog)
    if (any(contract$data_cards$role == "unassigned") ||
        any(contract$artifact_manifest$role == "unassigned") ||
        !all(contract$artifact_manifest$open_state == "metadata_only") ||
        any(!is.na(contract$artifact_manifest$local_sha256)) ||
        nrow(contract$gate_registry) != 0L) {
        .m12_fail("M12b must be role-assigned, metadata-only, and result-free")
    }
    invisible(TRUE)
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
    schema_catalog = "24e4aa3c93f5b0d87ff5ed4f92e70b08d0e258359a57b7c88fd1bd931c8b18f3",
    generator_manifest = "9d1381833cd59e8775bd99e730d5c783b98c5a82afff6d52276e83a46a8be76a",
    generator_protocol = "9f2166bc7a31112d5529d81adcd68cae70b99d12374f3d77bc45844d3d16bd7a",
    claim_inventory = "c6265115e6e2abb9ed30ef52a93d1605427de6869beed70f642a358bb9e3ce40",
    gate_registry = "0f10cd353879c95a851be3925a3c2f9570179cdf4592134215cccfdc58d80305",
    data_cards = "ea8605fe7f0b5c0765cb13eecb9851fee53ea527038984a44b8839d20d576188",
    data_roles = "650b6492d2f2a8d5542bc7097602c67f1f68db71d706a18e0bfd4dcfb4637aa7",
    artifact_manifest = "7c3267020d484ad9bc1bf4b8ac442cf0916d554e9e467e5c696ab2a20fb1c652",
    candidate_exclusions = "7ef14dd357171ffaad2e3793dfd1e038d1cc4108d56d052a2f751e4061e71af8",
    source_manifest = "b56898ee3ce5cd8ad2f2656a0b1cfb94792d8f1fb7f51b54472555279d8033ed",
    contract = "bb5624fecd9cd52980bd901c1b75eac0d58dcb1869d95f92014e0d9e086a53ed"
)

.m12_expect_error <- function(code) {
    inherits(try(force(code), silent = TRUE), "try-error")
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
    c(
        generator_scope = .m12_expect_error(m12_validate_contract(out_of_scope)),
        role_requires_split = .m12_expect_error(m12_validate_contract(missing_split)),
        role_artifact_agreement = .m12_expect_error(m12_validate_contract(role_mismatch)),
        opening_requires_role_hash = .m12_expect_error(m12_validate_contract(opened)),
        checksum_shape = .m12_expect_error(m12_validate_contract(wrong_checksum)),
        exact_subset = .m12_expect_error(m12_validate_contract(wrong_subset))
    )
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
    cat("state: metadata_only; roles=assigned; gates=0\n")
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
