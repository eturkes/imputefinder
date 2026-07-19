#!/usr/bin/env Rscript

# M12c A-C candidate/perturbation protocol. This freezes comparisons, support,
# multiplicity, calibration, and streams; it runs no candidate and reads no
# external result artifact. Run from the repository root:
# Rscript --vanilla dev/m12-candidate-protocol.R --verify

.m12c_contract <- local({
    path <- "dev/m12-validation-contract.R"
    if (!file.exists(path)) {
        stop("Run the M12c protocol from the repository root.", call. = FALSE)
    }
    environment <- new.env(parent = baseenv())
    sys.source(path, envir = environment)
    environment
})

.M12C_VERSION <- "m12_candidate_protocol_v2"
.M12C_CHECKED_ON <- "2026-07-19"
.M12C_CANDIDATE_REPLICATES <- 1:32
.M12C_DEVELOPMENT_REPLICATES <- 33:64
.M12C_BOOTSTRAP_DRAWS <- 999L
.M12C_PERMUTATION_DRAWS <- 9999L
.M12C_PROTOCOL_IDS <- c(
    A = "m12_a_candidate_protocol_v2",
    B = "m12_b_perturbation_protocol_v2",
    C = "m12_c_candidate_protocol_v2"
)

.m12c_fail <- function(...) {
    stop(..., call. = FALSE)
}

.m12c_nonempty <- function(x) {
    is.character(x) && !anyNA(x) && all(nzchar(x))
}

.m12c_id <- function(x) {
    .m12c_nonempty(x) && all(grepl("^[a-z0-9][a-z0-9_]*$", x))
}

.m12c_tokens <- function(x) {
    unique(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
}

m12c_evaluation_allocation <- function() {
    scenario_ids <- .m12c_contract$.m12_generator_manifest()$scenario_id
    data.frame(
        allocation_id = paste0(scenario_ids, "_allocation_v1"),
        scenario_id = scenario_ids,
        candidate_replicates = rep("1-32", length(scenario_ids)),
        candidate_count = rep(32L, length(scenario_ids)),
        development_replicates = rep("33-64", length(scenario_ids)),
        development_count = rep(32L, length(scenario_ids)),
        candidate_role = rep("method_selection_and_calibration", length(scenario_ids)),
        development_role = rep("frozen_untouched_test", length(scenario_ids)),
        result_state = rep("unrun", length(scenario_ids)),
        stringsAsFactors = FALSE
    )
}

m12c_descriptors <- function() {
    data.frame(
        protocol_id = c(
            unname(.M12C_PROTOCOL_IDS)
        ),
        track = c("A", "B", "C"),
        objective = c(
            "freeze algebraic design core and association candidate comparison",
            "freeze separated perturbation families cutoff policy streams and display calibration",
            "freeze detection and conditional-observed-abundance estimator comparison"
        ),
        input_evidence = c(
            "original pre-rescue mask plus typed design",
            "original pre-rescue evidence plus stable rules_v1 path",
            "fingerprint-matched original mask/intensities plus typed design and contrast"
        ),
        allocation = rep("replicates 1-32 candidate; 33-64 frozen development", 3L),
        implementation_file = rep("dev/m12-candidate-protocol.R", 3L),
        checked_on = rep(.M12C_CHECKED_ON, 3L),
        result_state = rep("frozen_unrun", 3L),
        stringsAsFactors = FALSE
    )
}

m12c_a_design_core <- function() {
    data.frame(
        core_id = "a_design_core_svd_v1",
        model_terms = "condition plus declared nuisance main effects plus only explicit interactions; block fixed effects where declared",
        encoding = "canonical sorted labels; intercept plus treatment contrasts; named term-to-column map retained",
        rank = "SVD rank; d_j > max(nrow,ncol)*double_epsilon*d_1",
        alias = "form null projector I-V_r V_r^T; scan canonical coefficient axes, modified-Gram-Schmidt projected axes, accept norm>sqrt(double_epsilon), unit-normalize, make first coefficient above tolerance positive; report basis, exact aliases, affected terms",
        estimability = "contrast row-space residual L2 <= sqrt(double_epsilon)*max(1,contrast_L2)",
        near_confounding = "report scaled singular values condition number overlap cells and leverage; no automatic correction",
        block = "subject/pair is one resampling unit and a fixed design term; technical siblings remain inside their biological unit",
        failure = "rank or estimability failure returns named nonestimable evidence before statistical association",
        state = "frozen_unrun",
        stringsAsFactors = FALSE
    )
}

m12c_a_candidates <- function() {
    data.frame(
        candidate_id = c(
            "a_fraction_ols_hc3_cr2",
            "a_fraction_freedman_lane",
            "a_fraction_quasibinomial"
        ),
        target = rep("declared_term_association_with_sample_detection_fraction", 3L),
        design_scope = c(
            "independent HC3 or paired/repeated CR2; estimable declared terms",
            "independent or paired/repeated only when frozen transformations preserve declared exchangeability and variance groups; estimable declared terms",
            "independent units only; empirical aggregate candidate with no feature-level peptide-effect analogue; estimable declared terms"
        ),
        response = rep("detected globally-observable proteins / globally-observable proteins per sample before rescue", 3L),
        fit = c(
            "ordinary least squares on detection fraction",
            "same OLS effect; reduced fitted values plus transformed reduced-model residuals; full-model scalar robust Wald statistic recomputed per transform",
            "one aggregate detected/undetected count per sample; quasibinomial logit with one estimated dispersion"
        ),
        effect = c(
            "design-adjusted percentage-point contrast",
            "design-adjusted percentage-point contrast",
            "standardized marginal probability contrast plus log-odds coefficient"
        ),
        uncertainty = c(
            "independent HC3 t interval with residual df; paired/repeated full-design identity-working-model CR2 t interval, unavailable when computed Satterthwaite df<5",
            "same robust OLS interval and reference-df floors as OLS candidate; two-sided W*>=Wobs; exact exceed/total including identity when allowable transformations<=100000, else 9999 unique keyed draws and p=(1+exceed)/(1+9999)",
            "dispersion-scaled covariance; t critical value with residual degrees of freedom"
        ),
        randomization = c(
            "none",
            "independent: permute whole-unit reduced residuals within the intersection of identical nuisance, exchangeability, and variance strata; paired/repeated: use only null-invariant condition-position swaps within complete blocks, moving technical siblings together; required crossing of a declared group => unavailable",
            "none"
        ),
        multiplicity = rep("Holm 0.05 within acquisition x declared association panel", 3L),
        separation = c("not applicable", "not applicable", "boundary/nonconvergence => structured unavailable"),
        reference_implementation = c(
            "base R matrix algebra; explicit HC3 and full-design identity-working-model CR2/Satterthwaite covariance",
            "base R matrix algebra; HC3/CR2 robust Wald statistic and frozen exchangeability/variance groups",
            "stats::glm family=quasibinomial"
        ),
        simplicity_rank = 1:3,
        state = rep("candidate_unrun", 3L),
        stringsAsFactors = FALSE
    )
}

m12c_a_support <- function() {
    data.frame(
        support_id = c(
            "a_association_independent",
            "a_association_paired"
        ),
        design = c("independent", "paired_or_repeated"),
        min_units_per_contrast_side = c(4L, 6L),
        min_complete_blocks = c(0L, 6L),
        residual_df_floor = c(3L, 3L),
        min_satterthwaite_df = c(NA_integer_, 5L),
        min_allowable_permutations = c(20L, 20L),
        extra_rule = c(
            "whole biological units only; every permutation stratum must preserve the declared nuisance layout",
            "complete blocks for the contrasted levels; technical siblings move together; CR2 uses the full design with identity working model and requires computed Satterthwaite df>=5"
        ),
        failure_code = c(
            "association_low_independent_support",
            "association_low_block_support"
        ),
        satterthwaite_failure_code = c(
            "not_applicable",
            "association_low_reference_df"
        ),
        stringsAsFactors = FALSE
    )
}

m12c_b_perturbations <- function() {
    data.frame(
        perturbation_id = c(
            "b_sampling_unit_bootstrap",
            "b_sampling_leave_one_unit_out",
            "b_estimator_feature_bootstrap",
            "b_assumption_policy_panel"
        ),
        family = c("sampling", "sampling", "estimator", "assumption"),
        draw_count = c(.M12C_BOOTSTRAP_DRAWS, 0L, .M12C_BOOTSTRAP_DRAWS, 6L),
        count_rule = c(
            "exactly 999",
            "exhaustive one draw per unique resampling unit",
            "exactly 999",
            "six named deterministic policies"
        ),
        unit = c(
            "declared subject/block; otherwise biological resampling unit",
            "declared subject/block; otherwise biological resampling unit",
            "declared synthetic feature module; otherwise individual protein descriptive only",
            "original named feature-condition evidence"
        ),
        replacement = c("with replacement", "without replacement", "with replacement", "none"),
        stratification = c(
            "independent units within condition; paired blocks jointly across conditions; technical siblings inseparable",
            "whole unit removed jointly; condition support rechecked",
            "synthetic modules within acquisition/scenario; public proteins by acquisition with no inferential coverage claim",
            "none"
        ),
        cutoff_policy = c(
            "fix successful baseline cutoff to isolate sampling",
            "fix successful baseline cutoff to isolate influence",
            "re-estimate automatic cutoff per weighted draw",
            "named fixed/re-estimated/sweep policies"
        ),
        weighting = c(
            "integer multiplicity weights equivalent to repeated units; original names never duplicated",
            "0/1 unit weights",
            "integer multiplicity weights equivalent to repeated feature clusters",
            "unit weights"
        ),
        failure = c(
            "unsupported draw recorded unavailable and remains in unconditional denominator",
            "unsupported omission recorded unavailable for that unit",
            "each cutoff failure retained by class; no fabricated endpoint",
            "unsupported policy retained as unavailable independently"
        ),
        stringsAsFactors = FALSE
    )
}

m12c_b_cutoff_policy <- function() {
    data.frame(
        policy_id = "b_cutoff_range_policy_v1",
        successful_draw_floor = 900L,
        successful_fraction_floor = 0.90,
        range_probability_low = 0.025,
        range_probability_high = 0.975,
        quantile_type = 8L,
        coverage_truth = "synthetic single-monotone generating detection midpoint only",
        public_scope = "empirical finite-dataset range only; never inferential coverage",
        false_confidence = "no-cliff draw family has success_fraction>=0.90 and range_width<=0.50 log2 without weak-identification flag",
        manual_cutoff = "estimator panel unavailable; sampling and named non-cutoff policy sensitivity remain",
        denominator = "all 999 draws for failure; successful draws only for range, with counts always reported",
        state = "frozen_unrun",
        stringsAsFactors = FALSE
    )
}

m12c_b_policy_scenarios <- function() {
    data.frame(
        policy_id = c(
            "b_policy_v1_baseline",
            "b_policy_v1_reestimated",
            "b_policy_v1_cutoff_sweep",
            "b_policy_v1_median",
            "b_policy_v1_trimmed",
            "b_policy_v1_half_or_more"
        ),
        cutoff = c(
            "fixed baseline", "re-estimated automatic",
            "type-8 q025;q25;q50;q75;q975 successful estimator cutoffs",
            "fixed baseline", "fixed baseline", "fixed baseline"
        ),
        summary = c(
            "arithmetic observed mean", "arithmetic observed mean",
            "arithmetic observed mean", "observed median",
            "20-percent observed trimmed mean", "arithmetic observed mean"
        ),
        majority = c(
            rep("strictly more than half", 5L),
            "at least half"
        ),
        interpretation = rep("sensitivity only; never changes rules_v1", 6L),
        state = rep("frozen_unrun", 6L),
        stringsAsFactors = FALSE
    )
}

m12c_b_display_calibration <- function() {
    data.frame(
        family = c("sampling", "estimator", "assumption"),
        score = c(
            "minimum of baseline-state frequency and retention agreement using draws 1-499",
            "baseline-state agreement with failures counted disagreement using draws 1-499",
            "continuous named-policy agreement only"
        ),
        validation_event = c(
            "draws 500-999 disagreement fraction exceeds 0.10",
            "draws 500-999 disagreement fraction exceeds 0.10",
            "not calibrated: six deterministic assumptions lack an independent draw law"
        ),
        threshold_grid = c(
            "0.80;0.85;0.90;0.95;0.975;0.99",
            "0.80;0.85;0.90;0.95;0.975;0.99",
            "none"
        ),
        selection = c(
            rep("replicates 1-32: lowest threshold with Wilson upper-95 risk<=0.10 and Wilson lower-95 coverage>=0.20", 2L),
            "continuous_only"
        ),
        development_gate = c(
            rep("replicates 33-64: Wilson upper-95 risk<=0.15 and point coverage>=0.10", 2L),
            "not_applicable"
        ),
        display = c(
            rep("stable/fragile only if calibration+development gates pass; unavailable stays distinct; continuous values retained", 2L),
            "no stable/fragile label; continuous scenario matrix retained"
        ),
        state = rep("frozen_unrun", 3L),
        stringsAsFactors = FALSE
    )
}

m12c_seed_protocol <- function() {
    data.frame(
        protocol_id = unname(.M12C_PROTOCOL_IDS),
        stream_id = c(
            "m12c_a_permutation_seed_stream_v2",
            "m12c_b_perturbation_seed_stream_v2",
            "m12c_c_resampling_seed_stream_v2"
        ),
        purpose = c(
            "A Monte Carlo restricted permutations",
            "B sampling and estimator bootstrap draws",
            "C paired-detection whole-block bootstrap intervals"
        ),
        key = rep("protocol_id|input_sha256|perturbation_id|canonical_instance_id|draw_id|nonce", 3L),
        digest = rep("SHA-256 UTF-8 with literal | separators over pipe-free canonical-token fields", 3L),
        integer_map = rep("first seven hex digits plus one => 1..268435456", 3L),
        collision = rep("canonical-order open addressing by nonce within one protocol+input+perturbation manifest; all stream instances/draws enter one call", 3L),
        rng_kind = rep("Mersenne-Twister/Inversion/Rejection in local restored scope", 3L),
        instance_order = rep("radix-sorted acquisition/condition/feature/unit IDs, then integer draw ID", 3L),
        state = rep("frozen_unrun", 3L),
        stringsAsFactors = FALSE
    )
}

m12c_c_estimands <- function() {
    data.frame(
        estimand_id = c(
            "c_detection_standardized_risk_difference_v1",
            "c_observed_abundance_contrast_v1"
        ),
        component = c("detection", "observed_abundance"),
        response = c(
            "D_ij=1 iff original x_ij finite before rescue",
            "original finite log2 x_ij only; no rescue or invented abundance"
        ),
        estimand = c(
            "mean_i[p_i(condition=contrast_right)-p_i(condition=contrast_left)] after setting each eligible empirical declared nuisance/block row to each contrasted level; log-odds contrast secondary",
            "design-adjusted mean log2 contrast conditional on detection"
        ),
        synthetic_truth = c(
            "same empirical-row counterfactual standardization applied to generator probabilities; alternative iff the standardized difference is nonzero",
            "within each contrasted level sum(p_ij*y_ij)/sum(p_ij) over eligible empirical design cells, then take the declared contrast; alternative iff nonzero"
        ),
        interval = rep("two-sided 95 percent; candidate-specific method retained", 2L),
        limitation = c(
            "model-conditional association; not structural absence or causal mechanism",
            "selected observed distribution on supplied unnormalised scale; not unconditional abundance"
        ),
        combined_decision = rep("none", 2L),
        stringsAsFactors = FALSE
    )
}

m12c_c_candidates <- function() {
    data.frame(
        candidate_id = c(
            "c_detection_glm_hc3",
            "c_detection_mean_br",
            "c_detection_gee_md",
            "c_detection_conditional_exact_aux",
            "c_abundance_ols_hc3_cr2",
            "c_abundance_limma_ebayes",
            "c_abundance_limma_robust"
        ),
        component = c(rep("detection", 4L), rep("observed_abundance", 3L)),
        design_scope = c(
            "independent; fixed-block paired when full rank after support filter",
            "independent; fixed-block paired when full rank after support filter",
            "paired/repeated or independent clusters with at least 10 clusters",
            "paired blocks with discordance; auxiliary odds-only comparator",
            "independent HC3 or fixed-block paired/repeated CR2",
            "independent or fixed-block paired",
            "independent or fixed-block paired"
        ),
        fit = c(
            "binomial-logit maximum likelihood",
            "mean bias-reduced binomial-logit adjusted scores",
            "marginal logit GEE with independence/exchangeable working correlation",
            "exact conditional logistic likelihood stratified by block",
            "feature-wise least squares on observed cells",
            "feature-wise linear model plus empirical-Bayes residual variance",
            "feature-wise linear model plus robust empirical-Bayes residual variance"
        ),
        uncertainty = c(
            "delta-method standardized-risk interval: HC3 t with residual df independent; 999 whole-block bootstrap percentile interval paired/repeated",
            "delta-method standardized-risk interval from adjusted-score covariance independent; 999 whole-block bootstrap percentile interval paired/repeated",
            "delta-method standardized-risk interval from Mancl-DeRouen corrected sandwich; t df=clusters-rank",
            "conditional profile-likelihood interval; no marginal risk-difference interval",
            "HC3 t interval with residual df independent; full-design identity-working-model CR2 t interval paired/repeated, unavailable when computed Satterthwaite df<5",
            "moderated t interval with missing values retained",
            "robust moderated t interval with missing values retained"
        ),
        separation = c(
            "detected before fit; separated fit unavailable",
            "finite full-rank estimates required; profile failure unavailable",
            "nonconvergence/infinite interval unavailable",
            "one-direction discordance may yield unbounded interval and is reported",
            rep("not applicable", 3L)
        ),
        reference_implementation = c(
            "base R stats plus explicit covariance",
            "brglm2 mean bias reduction reference; package choice deferred to passing study",
            "geepack score reference; exchangeable repeated/paired and independence single-observation clusters; explicit Mancl-DeRouen covariance",
            "survival exact conditional-likelihood reference",
            "base R matrix algebra; explicit HC3 and full-design identity-working-model CR2/Satterthwaite covariance",
            "Bioconductor limma lmFit/eBayes",
            "Bioconductor limma lmFit/eBayes robust"
        ),
        release_eligible = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
        simplicity_rank = c(1L, 2L, 3L, 4L, 1L, 2L, 3L),
        state = rep("candidate_unrun", 7L),
        stringsAsFactors = FALSE
    )
}

m12c_c_resampling <- function() {
    data.frame(
        resampling_id = "c_detection_paired_block_bootstrap_v2",
        candidate_ids = "c_detection_glm_hc3;c_detection_mean_br",
        draw_count = .M12C_BOOTSTRAP_DRAWS,
        unit = "complete declared subject/pair block",
        replacement = "with replacement via integer block weights; original names remain unique",
        stratification = "all contrasted condition positions and technical siblings stay inside the sampled block",
        estimand_grid = "refit each weighted draw; standardize fitted risks over the original empirical declared design grid",
        interval = "type-8 0.025/0.975 percentiles over successful draws",
        failure = "retain every failure class; interval unavailable unless at least 900 of 999 draws succeed",
        seed_stream = "m12c_c_resampling_seed_stream_v2",
        state = "frozen_unrun",
        stringsAsFactors = FALSE
    )
}

m12c_c_support <- function() {
    data.frame(
        support_id = c(
            "c_detection_independent",
            "c_detection_paired",
            "c_detection_gee",
            "c_abundance_independent",
            "c_abundance_paired",
            "c_all_components"
        ),
        component = c(rep("detection", 3L), rep("observed_abundance", 2L), "all"),
        design = c("independent", "paired", "clustered", "independent", "paired", "all"),
        min_units_per_contrast_side = c(4L, 6L, 10L, 3L, 6L, 1L),
        min_informative_events = c(2L, 3L, 2L, 3L, 6L, 1L),
        residual_df_floor = c(3L, 3L, 4L, 3L, 3L, 1L),
        min_satterthwaite_df = c(
            NA_integer_, NA_integer_, NA_integer_, NA_integer_, 5L, NA_integer_
        ),
        extra_rule = c(
            "at least two detected and two undetected overall",
            "at least three discordant complete blocks for contrasted levels",
            "at least ten independent clusters and both outcomes overall",
            "at least three observed units on each contrast side",
            "at least six independent blocks observed on both contrast sides; full-design identity-working-model CR2 requires computed Satterthwaite df>=5",
            "A contrast estimable; batch-exclusive condition contrast always unavailable"
        ),
        failure_code = c(
            "detection_low_support", "detection_low_discordance",
            "detection_too_few_clusters", "abundance_low_support",
            "abundance_low_paired_support", "contrast_nonestimable"
        ),
        satterthwaite_failure_code = c(
            rep("not_applicable", 4L), "abundance_low_reference_df",
            "not_applicable"
        ),
        stringsAsFactors = FALSE
    )
}

m12c_multiplicity <- function() {
    data.frame(
        family_id = c(
            "a_declared_terms_holm_v1",
            "c_detection_by_v1",
            "c_abundance_by_v1"
        ),
        track = c("A", "C", "C"),
        component = c("association", "detection", "observed_abundance"),
        family = c(
            "one acquisition x one analysis x all estimable declared association terms/interactions",
            "one requested contrast x detection component x all original nonglobal features",
            "one requested contrast x observed-abundance component x all original nonglobal features"
        ),
        method = c("Holm", "Benjamini-Yekutieli", "Benjamini-Yekutieli"),
        level = c(0.05, 0.05, 0.05),
        unavailable = c(
            "not tested and absent from term family",
            "p=1 for adjustment count; output remains unavailable",
            "p=1 for adjustment count; output remains unavailable"
        ),
        raw_evidence = rep("retain raw p interval effect support and adjusted value", 3L),
        stringsAsFactors = FALSE
    )
}

m12c_candidate_selection <- function() {
    data.frame(
        selection_id = c(
            "a_association_selection_v1",
            "c_detection_independent_selection_v1",
            "c_detection_paired_selection_v1",
            "c_abundance_selection_v1"
        ),
        track = c("A", "C", "C", "C"),
        candidate_ids = c(
            "a_fraction_ols_hc3_cr2;a_fraction_freedman_lane;a_fraction_quasibinomial",
            "c_detection_glm_hc3;c_detection_mean_br",
            "c_detection_glm_hc3;c_detection_mean_br;c_detection_gee_md;c_detection_conditional_exact_aux",
            "c_abundance_ols_hc3_cr2;c_abundance_limma_ebayes;c_abundance_limma_robust"
        ),
        hard_rejection = c(
            "any null/calibration/support/finite-output gate miss or causal wording violation",
            "any null/bias/coverage/support/finite-output gate miss",
            "any null/bias/coverage/support/finite-output gate miss; auxiliary conditional candidate cannot win",
            "any null/false-sign/bias/coverage/power/support/finite-output gate miss"
        ),
        ranking = c(
            rep("passing scope coverage descending; simplicity rank ascending when coverage/power within 0.02 and bias within 0.02; median runtime final tie-break", 3L),
            "passing scope coverage descending; simplicity rank ascending when coverage within 0.02 and bias within 0.10 log2; median runtime final tie-break"
        ),
        candidate_data = rep("only simulation replicates 1-32", 4L),
        development_data = rep("winner locked before untouched replicates 33-64", 4L),
        no_winner = rep("track association/estimator panel killed or parked; no gate weakening", 4L),
        state = rep("frozen_unrun", 4L),
        stringsAsFactors = FALSE
    )
}

m12c_excluded_candidates <- function() {
    data.frame(
        candidate_id = c(
            "a_cell_level_feature_fixed_logit",
            "c_detection_beta_binomial",
            "c_detection_quasibinomial_feature",
            "c_detection_bayesian_hierarchical"
        ),
        tracks = c("A", "C", "C", "C"),
        reason = c(
            "thousands of feature intercepts plus two-way dependence add complexity without a distinct global sentinel estimand",
            "grouping binary cells erases sample-level nuisance/block structure and does not identify the requested feature contrast",
            "one Bernoulli response per sample supplies no feature-level dispersion replication beyond declared clustering",
            "no validation-backed prior/transfer regime and thousands-feature computation add prior sensitivity beyond first-release need"
        ),
        reconsider_when = c(
            "global sample-fraction candidates all fail with a diagnosed estimand gap",
            "a grouped repeated-binomial acquisition regime is promoted",
            "replicated binomial denominators are available per feature/sample unit",
            "a prior is externally calibrated and a simpler candidate fails"
        ),
        decision = rep("excluded_before_results", 4L),
        stringsAsFactors = FALSE
    )
}

m12c_method_sources <- function() {
    data.frame(
        source_id = c(
            "r_p_adjust", "r_quantile", "brglm2_current", "bias_reduction_2020",
            "survival_clogit", "geepack_current", "gee_small_sample_md",
            "limma_current", "limma_manual", "permutation_glm_2014",
            "permutation_robust_wald_2019", "cluster_bootstrap_2013",
            "cluster_cr2_2018", "cluster_cr2_corrigendum_2023",
            "msqrob_hurdle_2020",
            "selective_risk_coverage_2010"
        ),
        tracks = c(
            "A;C", "B", "C", "C", "C", "C", "C", "C", "C", "A",
            "A", "B;C", "A;C", "A;C", "A", "B"
        ),
        scope = c(
            "Holm/BY definitions and dependency guarantees",
            "named sample-quantile algorithms; type 8 median-unbiased recommendation",
            "current finite mean-bias-reduced binomial reference implementation",
            "adjusted-score mean/median bias-reduction properties",
            "exact conditional logistic likelihood and computational limits",
            "current clustered categorical GEE reference implementation",
            "small-cluster sandwich bias and Mancl-DeRouen correction",
            "current omics linear-model implementation",
            "missing-value and blocked linear-model behavior",
            "nuisance-aware restricted permutation conditions",
            "heteroscedastic robust-W asymptotics and small-sample warnings",
            "cluster/subject bootstrap preserves within-cluster structure",
            "small-sample CR2/Satterthwaite inference with arbitrary fixed effects",
            "corrected absorbed-fixed-effect shortcut limited to OLS with identity working model",
            "protein-within-peptide quasibinomial precedent; not validation of the global all-protein denominator",
            "risk-coverage framing for abstaining displays"
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
            "https://doi.org/10.1016/j.neuroimage.2019.116030",
            "https://doi.org/10.1016/j.jmva.2012.09.003",
            "https://doi.org/10.1080/07350015.2016.1247004",
            "https://doi.org/10.1080/07350015.2023.2174123",
            "https://doi.org/10.1021/acs.analchem.9b04375",
            "https://www.jmlr.org/papers/v11/el-yaniv10a.html"
        ),
        checked_on = rep(.M12C_CHECKED_ON, 16L),
        stringsAsFactors = FALSE
    )
}

m12c_protocol <- function() {
    list(
        descriptors = m12c_descriptors(),
        evaluation_allocation = m12c_evaluation_allocation(),
        a_design_core = m12c_a_design_core(),
        a_candidates = m12c_a_candidates(),
        a_support = m12c_a_support(),
        b_perturbations = m12c_b_perturbations(),
        b_cutoff_policy = m12c_b_cutoff_policy(),
        b_policy_scenarios = m12c_b_policy_scenarios(),
        b_display_calibration = m12c_b_display_calibration(),
        seed_protocol = m12c_seed_protocol(),
        c_estimands = m12c_c_estimands(),
        c_candidates = m12c_c_candidates(),
        c_resampling = m12c_c_resampling(),
        c_support = m12c_c_support(),
        multiplicity = m12c_multiplicity(),
        candidate_selection = m12c_candidate_selection(),
        excluded_candidates = m12c_excluded_candidates(),
        method_sources = m12c_method_sources()
    )
}

m12c_seed_manifest <- function(
    protocol_id,
    input_sha256,
    perturbation_id,
    instance_ids,
    draw_ids
) {
    if (!is.character(protocol_id) || length(protocol_id) != 1L ||
        !protocol_id %in% unname(.M12C_PROTOCOL_IDS) ||
        !is.character(input_sha256) || length(input_sha256) != 1L ||
        !grepl("^[0-9a-f]{64}$", input_sha256) ||
        !is.character(perturbation_id) || length(perturbation_id) != 1L ||
        !grepl("^[a-z0-9][a-z0-9_]*$", perturbation_id) ||
        !length(instance_ids) || !.m12c_nonempty(instance_ids) ||
        anyDuplicated(instance_ids) ||
        any(grepl("|", instance_ids, fixed = TRUE)) ||
        any(grepl("\r", instance_ids, fixed = TRUE)) ||
        any(grepl("\n", instance_ids, fixed = TRUE)) ||
        !is.integer(draw_ids) || !length(draw_ids) || anyNA(draw_ids) ||
        any(draw_ids <= 0L) ||
        anyDuplicated(draw_ids)) {
        .m12c_fail("M12c seed manifest input is malformed.")
    }
    grid <- expand.grid(
        instance_id = sort(instance_ids, method = "radix"),
        draw_id = sort(draw_ids, method = "radix"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    used <- new.env(hash = TRUE, parent = emptyenv())
    seeds <- integer(nrow(grid))
    nonces <- integer(nrow(grid))
    for (index in seq_len(nrow(grid))) {
        nonce <- 0L
        repeat {
            key <- paste(
                protocol_id, input_sha256, perturbation_id,
                grid$instance_id[[index]], grid$draw_id[[index]], nonce,
                sep = "|"
            )
            digest <- unname(tools::sha256sum(
                bytes = charToRaw(enc2utf8(key))
            ))
            seed <- as.integer(strtoi(substr(digest, 1L, 7L), 16L) + 1L)
            seed_key <- as.character(seed)
            if (!exists(seed_key, envir = used, inherits = FALSE)) {
                assign(seed_key, TRUE, envir = used)
                break
            }
            nonce <- nonce + 1L
        }
        seeds[[index]] <- seed
        nonces[[index]] <- nonce
    }
    data.frame(
        protocol_id = protocol_id,
        perturbation_id = perturbation_id,
        instance_id = grid$instance_id,
        draw_id = as.integer(grid$draw_id),
        seed = seeds,
        nonce = nonces,
        stringsAsFactors = FALSE
    )
}

.m12c_validate_columns <- function(x, expected, label) {
    if (!is.data.frame(x) || !identical(names(x), expected) ||
        any(vapply(x, is.factor, logical(1L)))) {
        .m12c_fail(label, " schema/order is malformed.")
    }
}

m12c_validate_protocol <- function(protocol = m12c_protocol()) {
    expected_components <- c(
        "descriptors", "evaluation_allocation", "a_design_core",
        "a_candidates", "a_support", "b_perturbations", "b_cutoff_policy",
        "b_policy_scenarios", "b_display_calibration", "seed_protocol",
        "c_estimands", "c_candidates", "c_resampling", "c_support", "multiplicity",
        "candidate_selection", "excluded_candidates", "method_sources"
    )
    if (!identical(names(protocol), expected_components)) {
        .m12c_fail("M12c protocol components/order differ from v2.")
    }
    if (any(!vapply(protocol, is.data.frame, logical(1L))) ||
        any(vapply(protocol, function(x) any(vapply(x, is.factor, logical(1L))),
                   logical(1L)))) {
        .m12c_fail("M12c components must be factor-free data frames.")
    }

    descriptors <- protocol$descriptors
    .m12c_validate_columns(descriptors, c(
        "protocol_id", "track", "objective", "input_evidence", "allocation",
        "implementation_file", "checked_on", "result_state"
    ), "descriptor")
    if (!identical(descriptors$track, c("A", "B", "C")) ||
        !.m12c_id(descriptors$protocol_id) || anyDuplicated(descriptors$protocol_id) ||
        !all(descriptors$result_state == "frozen_unrun")) {
        .m12c_fail("M12c descriptor identities/state are malformed.")
    }

    allocation <- protocol$evaluation_allocation
    .m12c_validate_columns(allocation, c(
        "allocation_id", "scenario_id", "candidate_replicates",
        "candidate_count", "development_replicates", "development_count",
        "candidate_role", "development_role", "result_state"
    ), "allocation")
    scenarios <- .m12c_contract$.m12_generator_manifest()$scenario_id
    if (!identical(allocation$scenario_id, scenarios) ||
        anyDuplicated(allocation$allocation_id) ||
        !all(allocation$candidate_count == 32L) ||
        !all(allocation$development_count == 32L) ||
        !all(allocation$candidate_replicates == "1-32") ||
        !all(allocation$development_replicates == "33-64") ||
        !all(allocation$result_state == "unrun")) {
        .m12c_fail("M12c allocation must be disjoint 32/32 and unrun.")
    }

    if (nrow(protocol$a_design_core) != 1L || nrow(protocol$a_candidates) != 3L ||
        !.m12c_id(protocol$a_candidates$candidate_id) ||
        anyDuplicated(protocol$a_candidates$candidate_id) ||
        !identical(protocol$a_candidates$simplicity_rank, 1:3) ||
        !all(protocol$a_candidates$state == "candidate_unrun") ||
        nrow(protocol$a_support) != 2L ||
        any(protocol$a_support$min_units_per_contrast_side < 4L) ||
        any(protocol$a_support$residual_df_floor < 3L) ||
        !identical(
            protocol$a_support$min_satterthwaite_df,
            c(NA_integer_, 5L)
        ) ||
        !identical(
            protocol$a_support$satterthwaite_failure_code,
            c("not_applicable", "association_low_reference_df")
        ) ||
        any(protocol$a_support$min_allowable_permutations < 20L)) {
        .m12c_fail("M12c A core/candidate comparison is malformed.")
    }

    b <- protocol$b_perturbations
    if (!identical(b$family, c("sampling", "sampling", "estimator", "assumption")) ||
        !identical(b$draw_count, c(999L, 0L, 999L, 6L)) ||
        !identical(b$replacement, c(
            "with replacement", "without replacement", "with replacement", "none"
        ))) {
        .m12c_fail("M12c B perturbation units/counts/replacement drifted.")
    }
    cutoff <- protocol$b_cutoff_policy
    if (nrow(cutoff) != 1L || cutoff$successful_draw_floor != 900L ||
        cutoff$successful_fraction_floor != 0.90 || cutoff$quantile_type != 8L ||
        cutoff$range_probability_low != 0.025 ||
        cutoff$range_probability_high != 0.975 ||
        nrow(protocol$b_policy_scenarios) != 6L ||
        !identical(protocol$b_display_calibration$family,
                   c("sampling", "estimator", "assumption"))) {
        .m12c_fail("M12c B cutoff/policy/display calibration drifted.")
    }

    candidates <- protocol$c_candidates
    if (nrow(protocol$c_estimands) != 2L || nrow(candidates) != 7L ||
        !.m12c_id(candidates$candidate_id) || anyDuplicated(candidates$candidate_id) ||
        !identical(sum(candidates$release_eligible), 6L) ||
        nrow(protocol$c_resampling) != 1L ||
        protocol$c_resampling$draw_count != 999L ||
        protocol$c_resampling$replacement !=
            "with replacement via integer block weights; original names remain unique" ||
        nrow(protocol$c_support) != 6L ||
        any(protocol$c_support$min_units_per_contrast_side <= 0L) ||
        any(protocol$c_support$min_informative_events <= 0L) ||
        any(protocol$c_support$residual_df_floor <= 0L) ||
        !identical(
            protocol$c_support$min_satterthwaite_df,
            c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, 5L,
              NA_integer_)
        ) ||
        !identical(
            protocol$c_support$satterthwaite_failure_code,
            c(rep("not_applicable", 4L), "abundance_low_reference_df",
              "not_applicable")
        )) {
        .m12c_fail("M12c C estimand/candidate/support comparison is malformed.")
    }

    multiplicity <- protocol$multiplicity
    if (!identical(multiplicity$method,
                   c("Holm", "Benjamini-Yekutieli", "Benjamini-Yekutieli")) ||
        !all(multiplicity$level == 0.05) ||
        nrow(protocol$candidate_selection) != 4L ||
        !all(protocol$candidate_selection$state == "frozen_unrun") ||
        !all(protocol$excluded_candidates$decision == "excluded_before_results")) {
        .m12c_fail("M12c multiplicity/selection/exclusion rules drifted.")
    }
    sources <- protocol$method_sources
    if (!.m12c_id(sources$source_id) || anyDuplicated(sources$source_id) ||
        !all(grepl("^https://", sources$url)) ||
        !all(sources$checked_on == .M12C_CHECKED_ON)) {
        .m12c_fail("M12c primary method-source manifest is malformed.")
    }

    seed_protocol <- protocol$seed_protocol
    if (!identical(seed_protocol$protocol_id, unname(.M12C_PROTOCOL_IDS)) ||
        !all(seed_protocol$key ==
            "protocol_id|input_sha256|perturbation_id|canonical_instance_id|draw_id|nonce") ||
        !all(seed_protocol$state == "frozen_unrun")) {
        .m12c_fail("M12c seed stream identities/key/state drifted.")
    }

    seed_a <- m12c_seed_manifest(
        unname(.M12C_PROTOCOL_IDS[["B"]]),
        paste(rep("0", 64L), collapse = ""),
        "b_sampling_unit_bootstrap",
        c("condition_b", "condition_a"),
        as.integer(c(3, 1, 2))
    )
    seed_b <- m12c_seed_manifest(
        unname(.M12C_PROTOCOL_IDS[["B"]]),
        paste(rep("0", 64L), collapse = ""),
        "b_sampling_unit_bootstrap",
        c("condition_a", "condition_b"),
        1:3
    )
    if (!identical(seed_a, seed_b) || anyDuplicated(seed_a$seed) ||
        any(seed_a$seed <= 0L)) {
        .m12c_fail("M12c seed stream is not canonical and collision-safe.")
    }
    invisible(TRUE)
}

.m12c_hash_frames <- function(frames) {
    vapply(frames, .m12c_contract$.m12_hash_frame, character(1L))
}

.m12c_hash_bundle <- function(frames) {
    hashes <- .m12c_hash_frames(frames)
    aggregate <- paste(names(hashes), hashes, sep = "\t", collapse = "\n")
    unname(tools::sha256sum(bytes = charToRaw(enc2utf8(aggregate))))
}

m12c_protocol_hashes <- function(protocol = m12c_protocol()) {
    m12c_validate_protocol(protocol)
    components <- .m12c_hash_frames(protocol)
    descriptor <- function(track) {
        protocol$descriptors[protocol$descriptors$track == track, , drop = FALSE]
    }
    sources <- function(track) {
        keep <- vapply(protocol$method_sources$tracks, function(value) {
            track %in% .m12c_tokens(value)
        }, logical(1L))
        protocol$method_sources[keep, , drop = FALSE]
    }
    excluded <- function(track) {
        keep <- vapply(protocol$excluded_candidates$tracks, function(value) {
            track %in% .m12c_tokens(value)
        }, logical(1L))
        protocol$excluded_candidates[keep, , drop = FALSE]
    }
    a_frames <- list(
        descriptor = descriptor("A"),
        allocation = protocol$evaluation_allocation,
        design_core = protocol$a_design_core,
        candidates = protocol$a_candidates,
        support = protocol$a_support,
        seed = protocol$seed_protocol[protocol$seed_protocol$protocol_id == .M12C_PROTOCOL_IDS[["A"]], , drop = FALSE],
        multiplicity = protocol$multiplicity[protocol$multiplicity$track == "A", , drop = FALSE],
        selection = protocol$candidate_selection[protocol$candidate_selection$track == "A", , drop = FALSE],
        excluded = excluded("A"),
        sources = sources("A")
    )
    b_frames <- list(
        descriptor = descriptor("B"),
        allocation = protocol$evaluation_allocation,
        perturbations = protocol$b_perturbations,
        cutoff = protocol$b_cutoff_policy,
        policies = protocol$b_policy_scenarios,
        display = protocol$b_display_calibration,
        seed = protocol$seed_protocol[protocol$seed_protocol$protocol_id == .M12C_PROTOCOL_IDS[["B"]], , drop = FALSE],
        sources = sources("B")
    )
    c_frames <- list(
        descriptor = descriptor("C"),
        allocation = protocol$evaluation_allocation,
        estimands = protocol$c_estimands,
        candidates = protocol$c_candidates,
        resampling = protocol$c_resampling,
        support = protocol$c_support,
        seed = protocol$seed_protocol[protocol$seed_protocol$protocol_id == .M12C_PROTOCOL_IDS[["C"]], , drop = FALSE],
        multiplicity = protocol$multiplicity[protocol$multiplicity$track == "C", , drop = FALSE],
        selection = protocol$candidate_selection[protocol$candidate_selection$track == "C", , drop = FALSE],
        excluded = excluded("C"),
        sources = sources("C")
    )
    track_hashes <- c(
        a_protocol = .m12c_hash_bundle(a_frames),
        b_protocol = .m12c_hash_bundle(b_frames),
        c_protocol = .m12c_hash_bundle(c_frames)
    )
    bundle <- paste(names(track_hashes), track_hashes, sep = "\t", collapse = "\n")
    c(
        components,
        track_hashes,
        protocol_bundle = unname(tools::sha256sum(
            bytes = charToRaw(enc2utf8(bundle))
        ))
    )
}

.M12C_EXPECTED_HASHES <- c(
    descriptors = "28b4fb982a7d48747ace1e701a2250b7d1663e9a79948940cfc9d111bf896203",
    evaluation_allocation = "67eab76bd0e2d9507f4f558d9b823e33fc87ced797645a2848697540a1222754",
    a_design_core = "06899b13aff59d6bd727fd747611bc0a5c901ccedc1d6f376aee441e55b6f6c6",
    a_candidates = "ea6da31838d22bf079d812864c0ecfd204bc431ff7e3fd4ce97eea7efb124c9f",
    a_support = "2b1b65f3726798d04ec466fc8ef8f6b3e9ce390cf3ea822276cdff5b1a791a4e",
    b_perturbations = "0ce3cb0b8d781aa829b42e022e201ce83c572532885dbde5cce90ac055d65c19",
    b_cutoff_policy = "9079d59a76896874e469982f9308fd01f933ca967390a4cd2300ac52ebe2a052",
    b_policy_scenarios = "9e94fb14ffd9a4f1a4521468cc8f0821aaf03070b33a40b7d06da9d42763b04b",
    b_display_calibration = "0d7985d628f85cde7ded6a10bf721b70e7f0b64e099d3e79dbe116c775002034",
    seed_protocol = "efd82566f20f2fd3e035943e5c2854b383974351fe64fc3a2b3a9726aa33d574",
    c_estimands = "8049c2bb6bc135b1bd34627b16bbe262f89223443eb56b76a09440539e0c7585",
    c_candidates = "0ffb487bce24814fdebf570b83fea7d11c15d21dc226c3d2cdcba90cf033c7c0",
    c_resampling = "1759fff01d569b22f08657861785da9bcd0ec358dfa1f1e0c307bccae11c8629",
    c_support = "81c8726aa1cf68dd5af07af551a9d0b8c750cb5ced15074b09359be28067d26c",
    multiplicity = "75482e90cb35d5491bea00573bedcc09fcb0797aba852e4aebeacbb611e1e7c8",
    candidate_selection = "aa11b90dc459cd38baa6539e0f989b482394e71bba4c18e133a60b0c90b2c01a",
    excluded_candidates = "27dabc8a910dee3fc2e513a3a14041a812670abe6db90520d3f4d7207102085b",
    method_sources = "4f677fd08b19268a547395e1870878a24e00c58fb11b01880f43f8d8f51b56cf",
    a_protocol = "ca0cf8dbbf082446d9116ce280f0acf2c6517d35ec6137edb7b415585ce92683",
    b_protocol = "20912c4dfadeebb0c1da7100cb54dfa27f202ea39c3e0d991fb0ee38bab383c2",
    c_protocol = "f28b01a4cf5373807b4c71241a3563bc64cf55aa8d446cb62e242c9ecae76b04",
    protocol_bundle = "3bc4ef531d6c84475114befc5c153e50680459f1d8125fdbabac5fd840a1de65"
)

.m12c_expect_error <- function(code) {
    inherits(try(force(code), silent = TRUE), "try-error")
}

.m12c_self_tests <- function(protocol) {
    overlap <- protocol
    overlap$evaluation_allocation$development_replicates[[1L]] <- "32-63"
    wrong_replacement <- protocol
    wrong_replacement$b_perturbations$replacement[[1L]] <- "without replacement"
    weak_range <- protocol
    weak_range$b_cutoff_policy$successful_draw_floor <- 899L
    wrong_multiplicity <- protocol
    wrong_multiplicity$multiplicity$method[[2L]] <- "BH"
    opened_results <- protocol
    opened_results$descriptors$result_state[[1L]] <- "results_opened"
    wrong_seed_key <- protocol
    wrong_seed_key$seed_protocol$key[[2L]] <- "version|input|draw"
    wrong_c_resampling <- protocol
    wrong_c_resampling$c_resampling$draw_count <- 998L
    weak_reference_df <- protocol
    weak_reference_df$a_support$min_satterthwaite_df[[2L]] <- 4L
    empty_seed_manifest <- .m12c_expect_error(m12c_seed_manifest(
        unname(.M12C_PROTOCOL_IDS[["B"]]),
        paste(rep("0", 64L), collapse = ""),
        "b_sampling_unit_bootstrap",
        character(),
        integer()
    ))
    c(
        allocation_seal = .m12c_expect_error(m12c_validate_protocol(overlap)),
        replacement_rail = .m12c_expect_error(m12c_validate_protocol(wrong_replacement)),
        cutoff_floor = .m12c_expect_error(m12c_validate_protocol(weak_range)),
        dependence_control = .m12c_expect_error(m12c_validate_protocol(wrong_multiplicity)),
        results_sealed = .m12c_expect_error(m12c_validate_protocol(opened_results)),
        seed_key_seal = .m12c_expect_error(m12c_validate_protocol(wrong_seed_key)),
        c_resampling_seal = .m12c_expect_error(
            m12c_validate_protocol(wrong_c_resampling)
        ),
        satterthwaite_floor = .m12c_expect_error(
            m12c_validate_protocol(weak_reference_df)
        ),
        empty_seed_manifest = empty_seed_manifest
    )
}

.m12c_main <- function(args) {
    if (!identical(args, "--verify")) {
        .m12c_fail("usage: Rscript --vanilla dev/m12-candidate-protocol.R --verify")
    }
    protocol <- m12c_protocol()
    hashes <- m12c_protocol_hashes(protocol)
    registry <- .m12c_contract$.m12_protocol_registry()
    linked <- c(
        a_protocol = registry$protocol_hash[registry$track == "A"],
        b_protocol = registry$protocol_hash[registry$track == "B"],
        c_protocol = registry$protocol_hash[registry$track == "C"]
    )
    if (!identical(hashes[names(linked)], linked)) {
        .m12c_fail("M12c track hashes do not match the validation contract.")
    }
    self_tests <- .m12c_self_tests(protocol)
    if (!all(self_tests)) {
        .m12c_fail(
            "M12c protocol self-test failures: ",
            paste(names(self_tests)[!self_tests], collapse = ", ")
        )
    }
    frozen <- !anyNA(.M12C_EXPECTED_HASHES)
    if (frozen && !identical(hashes, .M12C_EXPECTED_HASHES)) {
        mismatch <- names(hashes)[hashes != .M12C_EXPECTED_HASHES]
        .m12c_fail("M12c frozen hash mismatch: ", paste(mismatch, collapse = ", "))
    }
    cat("protocol_version: ", .M12C_VERSION, "\n", sep = "")
    cat("state: frozen_unrun; external_artifacts=metadata_only\n")
    cat("allocation: candidate=1-32; development=33-64\n")
    cat("self_tests: ", paste(names(self_tests), self_tests, sep = "=", collapse = "; "), "\n", sep = "")
    cat(paste(names(hashes), hashes, sep = ": "), sep = "\n")
    cat("\n")
    if (!frozen) {
        cat("hash_state: bootstrap_unfrozen\n")
    }
    invisible(TRUE)
}

if (sys.nframe() == 0L) {
    .m12c_main(commandArgs(trailingOnly = TRUE))
}
