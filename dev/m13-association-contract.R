#!/usr/bin/env Rscript

# M13c pre-result association addendum. This supersedes only Track A details
# that M12 v2 left ambiguous. It opens no result artifact and runs no candidate.
# Run from the repository root:
# Rscript --vanilla dev/m13-association-contract.R --verify

.m13a_source <- function(path) {
    if (!file.exists(path)) {
        stop("Run the M13c contract from the repository root.", call. = FALSE)
    }
    environment <- new.env(parent = baseenv())
    sys.source(path, envir = environment)
    environment
}

.m13a_m12 <- .m13a_source("dev/m12-validation-contract.R")
.m13a_candidates <- .m13a_source("dev/m12-candidate-protocol.R")
.m13a_generator <- .m13a_source("dev/m12-generator-validation.R")

.M13A_VERSION <- "m13_association_contract_v4"
.M13A_PROTOCOL_ID <- "m13_a_association_protocol_v4"
.M13A_CHECKED_ON <- "2026-07-19"
.M13A_CANDIDATE_REPLICATES <- 1:32
.M13A_DEVELOPMENT_REPLICATES <- 33:64
.M13A_EXPECTED_UPSTREAM <- c(
    m12_contract = "45d1cda936b21a42d6735d900917734398eecb3ff1c136cab486cf90ee2b21e5",
    m12_gate_registry = "70d101b890de687be8973c9ff263caa454b61505f99470e5ec6061cb344db92f",
    m12_a_protocol = "ca0cf8dbbf082446d9116ce280f0acf2c6517d35ec6137edb7b415585ce92683",
    m12_protocol_bundle = "3bc4ef531d6c84475114befc5c153e50680459f1d8125fdbabac5fd840a1de65",
    generator_protocol = "cdea1bf874152e63fba08c49e390da1dd18690105e43a2dc7a03fc6866d0d080",
    generator_audit = "4d216966c8cc54b1b23685d8801c7c65b2afcc5879dd15430d57a88f765e8a00"
)
.M13A_EXPECTED_HASHES <- c(
    descriptor = "2e4785f5276ae64c23c9fd6fb42be12555d70a861957dbcf979cfe18b5df2de6",
    upstream_bindings = "5dc11e3dd5e43e1b638d8faf377836ad3b39b59319e0bfe3098d273891277bdd",
    normative_digest = "23c8288d81e08d7b87fd1d2d0f85cb0c4319edca680c276d3334a2ef94a74782",
    normative_source_digest = "ee49e77725a7bfc3e85a094496193e73869af4be4d1f9bceee8010bb55765a0b",
    implementation_bindings = "8ad1d822c6f626c45b9d55631a824eb5c1b165a54dc0ef32e3ebe51c4cb497f9",
    effective_manifest = "5dcf4f6408aed4ac7f9ccf07f4502805d34d2ea513ac68ac1c8e1259f23372f6",
    response = "baeb69d8ce3164629e8a50c5290549291d29b2a1c2dac155391be9da67ef57eb",
    strata = "4b7df275274bd127eed9a402695b2cfc9196f54a4c59271cde866aed17ceec92",
    hypotheses = "6a395c7dea901373cda91e22b719fa4dfad19ca2dcfaa9405d694dd9e8905e5c",
    support = "df04c882ea0928dac9e7f13ca813f2789091953cbf3fc237d857f6160df7d0a8",
    unit_positions = "6a69ec1ad4a3e4d85cbcdae4a1809e0d39728b12a3dc003e9d7335131f15d701",
    candidate_engines = "790608cd7875c28f2b5c2eae1445b9da1cdfe6aecfc4d39d68496f511a6835bf",
    numerics = "997a7c8d6c0eb3a20b3b841bb9d9688f586ca26ec6e96d1a3f3ba4ee61a56b85",
    permutation = "62ed7f671c9ee871eb6c9311cf6fba7eef0afe4d7ad384021c6b824c83d3430b",
    permutation_counts = "8a605d537ec00e71c2cdb0a6e8b27126775bb3d70824cd43e1e555e54b772816",
    permutation_maps = "fd4a0288087ca7ff0f402b9b0bc8d1a657810ddf32e4bce37d8160ae4f174203",
    seed_protocol = "d384a4f446acba3f4e36b3f3b662e6910402c43c43fc51909f476d0899140ebe",
    input_hash = "35871cc878f941af0b299b7cc680842df35b760b85c10b849f5168692c6fa64c",
    multiplicity = "1e6d84fe9da0b0623aa97f0bd0f7cba88051e8e703798d901badc06379e7d0ae",
    result_schema = "1b924dd61b42f24d1480da585fe4166cb2d2fe7ef90be83f99a9f1b74834a5cf",
    result_records = "1cfbdf3fae1be20dead577bf961d3e2091b8edbca36f66479c7a5b57a024c21e",
    lifecycle = "d7b358c0835b2fef76892e654d7477a66a8d745445c53b82339c88862d303040",
    unavailable_codes = "75ff8476e75dde29007d7b879409df4f157a844ddf8ee09e2abad30627bdc021",
    truth = "821ad25d39770f46d84ea48b2a76ecae2ec1641068bebbb29e4d2709432946cc",
    truth_methods = "13dabc75ae7f8a780831586b35e7f7900bd727b412260f74efcdf5731607b261",
    opportunities = "21eb96391a0efa21f410390f2e77a67fd174fe165ea2631cac969ecf78ca4e6b",
    evaluation_allocation = "f2fdcbc9116fb76c9dbd7a28cba8b8d30e685fc97f785547037ea8b0b941cc38",
    selection = "f0e99d3f86baf2a8993491932f9a1f5378d629986e592a3c8b39d2fe39cc20e2",
    screening_allocation = "d861cba428cc76eac35cb3bd6648f672985e10afc6c74779428f1f5895e8da09",
    screening_thresholds = "c8ca1c22cba1e7da8bdd2898d97370efa508eeccc93f41129182cf4a1d4b20aa",
    gate_statistics = "7deed5afc046265b327229f6fc8756cc0ee2525d9ef2da4d3fc7d367fba5b10c",
    gate_registry = "65f01b20e74b216399e17433df393f32e9bda189f73992ab419069cf7d5c343b",
    completed_gate_lineage = "9d7f5c30e28cd437a1dde001a1715e6e143dfbe8b1f399fb7d7f3029b1198a6e",
    structural_gates = "c29ac1fe23bb78f10640fb65a1756d38d4b18279dd967eb422f08f8024ec1db8",
    execution_seal = "db611f42f94c2823fffc9dae039ab26c2d1df97debae6a3fac99529e321879bf",
    evidence_hash_protocol = "0d850d7d013d199bc1f6ba0d707f7b3390d26b82a2cbe10fe85039f8f4ec7df3",
    implementation_manifest_schema = "8af3ebaaff8172968896500731747f0cbebe716bd4258143d271717f3aaf678a",
    candidate_evidence_schema = "1e2d6e09eac4754cd8b39ac90d858188edb09a93c8422468018050f079708205",
    release_boundary = "5e25356ea2e3cb33b31439bf78ac01e254e99435fe0e12fc8be716c6bf1348e0",
    contract = "422954c88ed7683f58ac296cb4d9e785d61c3f80a276914235d4b285bbc16547"
)

# BEGIN M13A NORMATIVE SOURCE

.m13a_fail <- function(...) {
    stop(..., call. = FALSE)
}

.m13a_tokens <- function(...) {
    paste(sort(unique(c(...)), method = "radix"), collapse = ";")
}

.m13a_live_upstream <- function() {
    m12 <- .m13a_m12$m12_contract_hashes()
    candidates <- .m13a_candidates$m12c_protocol_hashes()
    generator <- .m13a_generator$m12b_protocol_hashes()
    c(
        m12_contract = unname(m12[["contract"]]),
        m12_gate_registry = unname(m12[["gate_registry"]]),
        m12_a_protocol = unname(candidates[["a_protocol"]]),
        m12_protocol_bundle = unname(candidates[["protocol_bundle"]]),
        generator_protocol = unname(generator[["protocol"]]),
        generator_audit = unname(
            .m13a_m12$.m12_generator_protocol()$generator_audit_hash
        )
    )
}

m13a_descriptor <- function() {
    data.frame(
        protocol_id = .M13A_PROTOCOL_ID,
        version = .M13A_VERSION,
        track = "A",
        role = "pre_result_superseding_addendum",
        supersedes = "ambiguous Track A details in m12_a_candidate_protocol_v2",
        preserves = "M12 thresholds;candidates;generator identities;external data roles;sealed evidence",
        gate_version = "nine association/public rows superseded by m13a_gate_registry_v4; eleven completed core rows retain exact M12 bindings",
        result_state = "frozen_unrun",
        external_state = "metadata_only",
        checked_on = .M13A_CHECKED_ON,
        implementation_file = "dev/m13-association-contract.R",
        detail_file = "dev/m13-association-contract.md",
        stringsAsFactors = FALSE
    )
}

m13a_upstream_bindings <- function() {
    data.frame(
        binding = names(.M13A_EXPECTED_UPSTREAM),
        sha256 = unname(.M13A_EXPECTED_UPSTREAM),
        drift_rule = "reject_before_candidate_computation",
        stringsAsFactors = FALSE
    )
}

m13a_implementation_bindings <- function() {
    data.frame(
        component = c(
            "design_schema", "design_encoding", "design_rank",
            "design_aliasing", "design_contrast", "design_units",
            "input_hash", "cr2_residual_projector"
        ),
        identifier = c(
            "design_estimability_v1", "canonical_treatment_contrasts_v1",
            "svd_relative_v1", "ordered_null_projector_mgs_v1",
            "row_space_residual_v1", "declared_block_or_sample_v1",
            "association_mask_design_be_v1",
            "complete_qr_residual_projector_v1"
        ),
        drift_rule = "new version plus explicit migration; never silent reinterpretation",
        stringsAsFactors = FALSE
    )
}

m13a_response <- function() {
    data.frame(
        response_id = "a_sample_detection_fraction_v3",
        source = "original finite mask before rescue",
        denominator = "features finite in at least one sample of the complete analysis input",
        denominator_scope = "fixed before acquisition stratification and shared by all strata",
        numerator = "denominator features finite in the named sample",
        scale = "probability fraction; percentage points only for presentation",
        zero_denominator = "association_no_observable_features",
        intensity_use = "none",
        stringsAsFactors = FALSE
    )
}

m13a_strata <- function() {
    data.frame(
        stratum_id = c("declared_acquisition", "no_acquisition_role"),
        selector = c(
            "one canonical sorted acquisition level per stratum",
            "one synthetic stratum named all"
        ),
        model_rebuild = c(
            "subset samples then remove acquisition role/interactions and rebuild encoding/SVD/units",
            "rebuild encoding/SVD/units on all samples"
        ),
        family = rep("one independent Holm family per rebuilt stratum", 2L),
        identity = c(
            "s_ plus full SHA-256 of length-prefixed acquisition label",
            "literal all; acquisition display label literal undeclared"
        ),
        forbidden = rep("subsetting a globally built design core", 2L),
        stringsAsFactors = FALSE
    )
}

m13a_stratum_id <- function(acquisition = NULL) {
    if (is.null(acquisition)) {
        return("all")
    }
    if (!.m13a_scalar_character(acquisition)) {
        .m13a_fail("M13c acquisition stratum labels must be nonempty scalars.")
    }
    bytes <- c(
        charToRaw("imputefinder:association-stratum:v1"), as.raw(0L),
        .m13a_encode_text(acquisition)
    )
    paste0("s_", unname(tools::sha256sum(bytes = bytes)))
}

m13a_hypotheses <- function() {
    data.frame(
        hypothesis_class = c(
            "condition_main", "nuisance_main", "declared_interaction",
            "intercept", "block_main", "block_interaction",
            "acquisition_term"
        ),
        family_role = c(
            "eligible", "eligible", "eligible", "adjustment_only",
            "adjustment_only", "adjustment_only", "stratifier_only"
        ),
        rule = c(
            "each canonical nonreference treatment column or numeric column is one 1-df coefficient hypothesis",
            "each canonical nonreference treatment column or numeric column is one 1-df coefficient hypothesis",
            "each canonical product column is one 1-df hypothesis only when every component role is condition or nuisance",
            "retained in fit and absent from association family",
            "retained in fit and absent from association family",
            "retained only if declared and absent from association family",
            "removed after defining strata; any containing interaction removed"
        ),
        effect_direction = c(
            "encoded level minus reference or plus one numeric encoded unit",
            "encoded level minus reference or plus one numeric encoded unit",
            "plus one product-column unit conditional on all retained design columns",
            "not_reported", "not_reported", "not_reported", "not_reported"
        ),
        estimability = c(
            rep("raw coefficient axis must pass frozen row-space residual check", 3L),
            rep("not_applicable", 4L)
        ),
        stringsAsFactors = FALSE
    )
}

m13a_support <- function() {
    data.frame(
        support_id = c(
            "independent_categorical", "independent_numeric",
            "blocked_categorical", "blocked_numeric", "interaction_cells"
        ),
        unit = c(
            "unique biological sample units",
            "unique biological sample units",
            "unique declared blocks",
            "unique declared blocks",
            "Cartesian component-side cells in biological units"
        ),
        side_rule = c(
            "canonical reference versus encoded target",
            "below versus above unit-level median; exact ties excluded",
            "blocks containing both canonical reference and encoded target",
            "blocks containing unit-positions below and above the global unit-position median; exact ties excluded",
            "Cartesian product of every categorical reference/target and every numeric below/above-median component side"
        ),
        floor = c(4L, 4L, 6L, 6L, 4L),
        blocked_floor = c(NA_integer_, NA_integer_, 6L, 6L, 6L),
        extra = c(
            "each side independently",
            "median computed over independent units; record range and whether encoded 0-to-1 leaves range",
            "at least six complete blocks and six blocks per side",
            "at least six complete blocks and six blocks per side; record range and 0-to-1 extrapolation",
            "every 2^k relevant cell meets 4 independent units or 6 complete blocks; stored low/high counts are the minima over component marginal sides; numeric support is stored once per numeric component"
        ),
        failure_code = c(
            "association_low_independent_support",
            "association_low_independent_support",
            "association_low_block_support",
            "association_low_block_support",
            "association_low_interaction_support"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_unit_positions <- function() {
    data.frame(
        design = c("independent", "blocked", "all"),
        position_key = c(
            "canonical sample ID",
            "block ID plus canonical values of every condition and nuisance column",
            "acquisition omitted because the model is already stratified"
        ),
        sibling_rule = c(
            "one row per sample",
            "rows with identical complete position key are technical siblings and collapse once; their tested numeric values must be bit-identical after integer-to-double normalization",
            "a duplicated key with discordant declared values is association_incompatible_unit_positions"
        ),
        weighting = c(
            "each sample once",
            "each unique unit-position once for medians; each block once for side/cell counts regardless of repeated positions",
            "exact median ties belong to neither side"
        ),
        interaction = c(
            "each unit contributes once to its Cartesian cell",
            "a block supports a cell when at least one unique position occupies it; repeated positions never increase the block count",
            "all 2^k cells are materialized in canonical binary order"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_candidate_engines <- function() {
    data.frame(
        candidate_id = c(
            "a_fraction_ols_hc3_cr2",
            "a_fraction_freedman_lane",
            "a_fraction_quasibinomial"
        ),
        licensed_scope = c(
            "independent HC3; declared-block full-design CR2",
            "independent HC3 or declared-block CR2 plus compatible restricted transformations",
            "independent units only"
        ),
        estimand = c(
            "OLS encoded-column coefficient on response fraction",
            "same OLS coefficient and robust interval; restricted-permutation p",
            "empirical standardized encoded-column 0-to-1 probability difference; log-odds coefficient secondary"
        ),
        reference_df = c(
            "HC3 residual n-r; CR2 scalar Satterthwaite >=5",
            "same robust interval floor; permutation p has no parametric df",
            "residual n-r for log-odds test and delta interval"
        ),
        p_value = c(
            "two-sided robust coefficient t",
            "two-sided robust W permutation exceedance",
            "two-sided log-odds coefficient t"
        ),
        numerical_abstention = c(
            "nonpositive robust variance; invalid leverage; CR2 spectral failure",
            "base robust failure; transformed-fit failure; incompatible or fewer than 20 transforms",
            "boundary counts; separation risk; nonconvergence; weighted rank loss; invalid dispersion/covariance"
        ),
        simplicity_rank = 1:3,
        stringsAsFactors = FALSE
    )
}

m13a_numerics <- function() {
    data.frame(
        component = c(
            "rank_coordinates", "hc3_leverage", "covariance_psd", "cr2",
            "null_invariance", "wald_exceedance", "quasibinomial_irls",
            "quasibinomial_probability", "quasibinomial_delta",
            "global_floors"
        ),
        frozen_rule = c(
            "X=U_r D_r V_r'; Z=X V_r; K=(Z'Z)^-1; coefficient direction c=V_r' e_j",
            "tau(A)=max(dim(A))*double_epsilon*max(1,max_abs(A)); clamp leverage into [0,1] only within tau(H); require 1-h>tau(H); HC3 uses K Z' diag((e/(1-h))^2) Z K",
            "one tau(A) rule for every matrix; symmetrize covariance; eigenvalue below -tau(V) fails; values in [-tau(V),tau(V)] become zero; scalar contrast variance must exceed tau(V)",
            "full design including block fixed effects; identity working model; M=Q_perp Q_perp' from symmetrized complete LAPACK QR orthogonal complement of U_r; symmetric MP inverse-square-root correction using tau(B_g); scalar Satterthwaite",
            "require max_abs(P Z0-Z0)<=sqrt(double_epsilon)*max(1,max_abs(Z0)) plus exact label preservation",
            "after finite W computation compare literal IEEE-754 Wstar>=Wobs; no tie jitter or discarded transform",
            "stats::glm.fit in rank coordinates with quasibinomial; start/etastart NULL; family mustart; glm.control(epsilon=1e-10,maxit=100,trace=FALSE); final weighted SVD rank r",
            "every grouped count strictly between 0 and G; every fitted mu in [sqrt(double_epsilon),1-sqrt(double_epsilon)]; Pearson dispersion finite and strictly positive; V_Q then passes the common covariance PSD/scalar rules",
            "g=mean(mu1*(1-mu1)*z1-mu0*(1-mu0)*z0); SE=sqrt(g' V_Q g); no interval clamp",
            "residual df >=3; CR2 Satterthwaite df >=5"
        ),
        tolerance = c(
            "existing SVD relative threshold and row-space tolerance",
            "common entrywise tau(A)",
            "common entrywise tau(A)",
            "common entrywise tau(B_g)",
            "sqrt(double_epsilon) relative max norm",
            "literal >= on finite doubles",
            "IRLS deviance relative change <=1e-10 within 100 iterations; frozen weighted SVD threshold",
            "sqrt(double_epsilon) probability boundary; exact finite phi>0; common covariance tau for V_Q",
            "covariance PSD/scalar rules above",
            "exact"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_permutation <- function() {
    data.frame(
        design = c("independent", "blocked"),
        atomic_unit = c(
            "whole biological unit",
            "complete declared block; technical siblings move jointly"
        ),
        exchangeability_group = c(
            "acquisition stratum; further split by identical constrained-null design rows",
            "acquisition stratum; only canonical condition-position bundle swaps within block"
        ),
        variance_group = c(
            "acquisition stratum; no undeclared variance-group metadata invented",
            "acquisition stratum; exact nuisance-position compatibility required"
        ),
        invariant = c(
            "P Z0 equals Z0 within frozen tolerance; design labels stay fixed at target rows; source rows remain inside their constrained-null group",
            "P Z0 equals Z0 within frozen tolerance; source rows remain in the same block and nuisance signature; the tested categorical pivot may exchange"
        ),
        transformation_count = c(
            "product over compatible groups of factorial(group unit count)",
            "2^(complete swappable block count); incomplete blocks fixed"
        ),
        canonical_order = c(
            "radix groups/units; identity first; globally lexicographic canonical row-index maps",
            "radix blocks; identity first; ascending binary swap masks"
        ),
        incompatible = rep("association_incompatible_permutation", 2L),
        stringsAsFactors = FALSE
    )
}

m13a_permutation_counts <- function() {
    data.frame(
        case = c("insufficient", "exact", "monte_carlo"),
        total_transformations = c("T<20", "20<=T<=100000", "T>100000"),
        draw_rule = c(
            "none",
            "enumerate all distinct transforms canonically including identity",
            "9999 distinct nonidentity transforms from the canonical per-draw seeded streams; reject identity/duplicates and retry within that stream; 1000000 retries per draw then numerical failure"
        ),
        p_rule = c(
            "unavailable",
            "count(Wstar>=Wobs)/T",
            "(1+count(Wstar>=Wobs))/10000"
        ),
        diagnostic_count = c(
            "exact T",
            "exact T",
            "min(nearest finite-double T, double.xmax); branch uses log(T)"
        ),
        failure_code = c(
            "association_low_permutation_resolution",
            "any transformed fit failure => association_numerical_failure",
            "any transformed fit failure => association_numerical_failure"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_permutation_maps <- function() {
    data.frame(
        design = c("independent", "blocked"),
        candidate_map = c(
            "permute canonical biological-unit IDs within each constrained-null row group",
            "choose the radix-first categorical condition component as pivot; swap only its reference/target bundles independently in each complete block; other levels/incomplete blocks fixed; no pivot => incompatible"
        ),
        grouping = c(
            "Z0 rows share a group iff every pair has max_abs difference <= sqrt(double_epsilon)*max(1,max_abs(Z0)); nontransitive tolerance relation => incompatible",
            "bundle signature is canonical values of every nuisance column; each reference signature must have one target signature with equal sibling multiplicity"
        ),
        distinct = c(
            "distinct canonical row-index permutation vectors, even when transformed response/statistic values coincide",
            "distinct canonical row-index permutation vectors generated by binary complete-block swap masks"
        ),
        expansion = c(
            "one independent unit is one row",
            "within each matched signature, radix-ordered technical siblings map positionally and move jointly"
        ),
        validation = rep(
            "deduplicate row-index maps; identity first; require exact bijection and frozen P Z0 invariant; a byte-identical transformed residual refits original response bytes; tested-pivot equality is intentionally not required",
            2L
        ),
        stringsAsFactors = FALSE
    )
}

m13a_seed_protocol <- function() {
    data.frame(
        stream_id = "m13a_permutation_seed_stream_v4",
        protocol_id = .M13A_PROTOCOL_ID,
        input_hash = "association_mask_design_be_v1 SHA-256; never intensity values or R serialize bytes",
        key = "seed-v2 tag plus typed length-prefixed UTF-8/bytes vector(protocol_id,input_sha256,candidate_id,acquisition,hypothesis_id,decimal draw_id,decimal seed_nonce)",
        integer_map = "first seven SHA-256 hex digits plus one",
        collision = "canonical open addressing across the complete requested seed manifest; Monte Carlo map collisions use within-draw stream retries",
        order = "one complete seed call before maps; radix candidate/acquisition/hypothesis then integer draw ID; emitted result rows only for available Monte Carlo outcomes",
        rng = "Mersenne-Twister/Inversion/Rejection inside exact restored local scope; each draw seed starts one stream; identity/duplicate maps retry in that stream",
        stringsAsFactors = FALSE
    )
}

m13a_input_hash <- function() {
    data.frame(
        position = 1:8,
        field = c(
            "version_tag", "canonical_feature_names", "canonical_sample_names",
            "packed_missing_mask", "canonical_role_map", "canonical_interactions",
            "normalised_declared_values", "section_framing"
        ),
        encoding = c(
            "UTF-8 literal imputefinder:association-mask-design:be:v1 plus NUL",
            "signed big-endian i32 count; radix order; each name one-byte UTF8/bytes tag plus big-endian i32 byte length plus bytes",
            "same name encoding; order also defines design/mask sample alignment",
            "feature then sample traversal; base packBits LSB0; FALSE zero padding; big-endian i32 cell and byte counts",
            "fixed condition,nuisance,block,acquisition order; name-encoded role/column counts and values",
            "component names radix-sorted; interactions sorted by hex of their complete length-prefixed bytes",
            "declared columns in role order; categorical tagged/name-encoded; nonfactor integer/double tagged IEEE-754 binary64 big-endian with integer-to-double and negative-zero normalization",
            "one-byte section tags 1..7; every variable section carries explicit signed big-endian i32 counts"
        ),
        invariant = c(
            "version separation", "row permutation", "column permutation",
            "row/column permutation", "role selector order", "interaction order",
            "factor level order and integer/double equivalence", "unambiguous concatenation"
        ),
        stringsAsFactors = FALSE
    )
}

.M13A_HASH_ROLES <- c("condition", "nuisance", "block", "acquisition")

.m13a_be_i32 <- function(x) {
    if (!is.numeric(x) || anyNA(x) || any(x < 0) ||
        any(x > .Machine$integer.max) || any(x != floor(x))) {
        .m13a_fail("M13c hash length exceeds its signed-i32 protocol.")
    }
    writeBin(as.integer(x), raw(), size = 4L, endian = "big")
}

.m13a_canonical_text <- function(x) {
    output <- unname(as.character(x))
    text <- Encoding(output) != "bytes"
    output[text] <- enc2utf8(output[text])
    if (anyNA(output)) {
        .m13a_fail("M13c hash text must encode losslessly.")
    }
    output
}

.m13a_encode_text <- function(value) {
    value <- .m13a_canonical_text(value)
    if (length(value) != 1L) {
        .m13a_fail("M13c text records are scalar.")
    }
    tag <- if (identical(Encoding(value), "bytes")) as.raw(2L) else as.raw(1L)
    bytes <- charToRaw(value)
    c(tag, .m13a_be_i32(length(bytes)), bytes)
}

.m13a_encode_text_vector <- function(values) {
    values <- .m13a_canonical_text(values)
    c(
        .m13a_be_i32(length(values)),
        unlist(lapply(values, .m13a_encode_text), use.names = FALSE)
    )
}

.m13a_raw_hex <- function(bytes) {
    paste(sprintf("%02x", as.integer(bytes)), collapse = "")
}

.m13a_normalise_interactions <- function(interactions, declared) {
    if (is.null(interactions)) {
        interactions <- list()
    }
    if (!is.list(interactions)) {
        .m13a_fail("M13c hash interactions must be a list.")
    }
    output <- lapply(interactions, function(components) {
        components <- sort(.m13a_canonical_text(components), method = "radix")
        if (length(components) < 2L || anyDuplicated(components) ||
            !all(components %in% declared)) {
            .m13a_fail("M13c hash interaction components are malformed.")
        }
        components
    })
    encoded <- lapply(output, .m13a_encode_text_vector)
    if (length(encoded)) {
        keys <- vapply(encoded, .m13a_raw_hex, character(1L))
        if (anyDuplicated(keys)) {
            .m13a_fail("M13c hash interactions must be unique.")
        }
        output <- output[order(keys, method = "radix")]
    }
    output
}

.m13a_normalise_hash_roles <- function(roles, sample_data) {
    if (!is.list(roles) || !setequal(names(roles), .M13A_HASH_ROLES)) {
        .m13a_fail("M13c hash roles are malformed.")
    }
    roles <- roles[.M13A_HASH_ROLES]
    output <- lapply(seq_along(roles), function(index) {
        value <- .m13a_canonical_text(roles[[index]])
        if (identical(.M13A_HASH_ROLES[[index]], "nuisance")) {
            value <- sort(value, method = "radix")
        }
        value
    })
    names(output) <- .M13A_HASH_ROLES
    role_lengths <- lengths(output)
    if (role_lengths[[1L]] != 1L || any(role_lengths[c(3L, 4L)] > 1L)) {
        .m13a_fail("M13c hash role cardinality is malformed.")
    }
    declared <- unname(unlist(output, use.names = FALSE))
    if (anyDuplicated(declared) || !all(declared %in% names(sample_data))) {
        .m13a_fail("M13c hash roles must name unique sample-data columns.")
    }
    output
}

.m13a_encode_hash_roles <- function(roles) {
    c(
        .m13a_be_i32(length(roles)),
        unlist(lapply(names(roles), function(role) {
            c(.m13a_encode_text(role), .m13a_encode_text_vector(roles[[role]]))
        }), use.names = FALSE)
    )
}

.m13a_encode_hash_interactions <- function(interactions) {
    c(
        .m13a_be_i32(length(interactions)),
        unlist(lapply(interactions, .m13a_encode_text_vector), use.names = FALSE)
    )
}

.m13a_encode_design_column <- function(column, values) {
    numeric <- !is.factor(values) && (is.integer(values) || is.double(values))
    if (numeric) {
        values <- as.double(values)
        values[values == 0] <- 0
        if (anyNA(values) || any(!is.finite(values))) {
            .m13a_fail("M13c numeric hash values must be finite.")
        }
        payload <- c(
            as.raw(1L), .m13a_be_i32(length(values)),
            writeBin(values, raw(), size = 8L, endian = "big")
        )
    } else {
        values <- .m13a_canonical_text(values)
        if (!length(values) || any(!nzchar(values))) {
            .m13a_fail("M13c categorical hash values must be nonempty.")
        }
        payload <- c(as.raw(2L), .m13a_encode_text_vector(values))
    }
    c(.m13a_encode_text(column), payload)
}

m13a_mask_design_bytes <- function(
    missing_mask,
    sample_data,
    roles,
    interactions = list()
) {
    valid_mask <- is.matrix(missing_mask) && is.logical(missing_mask) &&
        !anyNA(missing_mask) && nrow(missing_mask) > 0L && ncol(missing_mask) > 0L &&
        !is.null(rownames(missing_mask)) && !is.null(colnames(missing_mask)) &&
        !anyNA(rownames(missing_mask)) && !anyNA(colnames(missing_mask)) &&
        all(nzchar(rownames(missing_mask))) && all(nzchar(colnames(missing_mask))) &&
        !anyDuplicated(rownames(missing_mask)) && !anyDuplicated(colnames(missing_mask))
    valid_data <- is.data.frame(sample_data) && nrow(sample_data) == ncol(missing_mask) &&
        !is.null(rownames(sample_data)) && !anyNA(rownames(sample_data)) &&
        !anyDuplicated(rownames(sample_data)) &&
        setequal(rownames(sample_data), colnames(missing_mask))
    if (!valid_mask || !valid_data) {
        .m13a_fail("M13c mask/design hash input is malformed.")
    }
    feature_names <- .m13a_canonical_text(rownames(missing_mask))
    sample_names <- .m13a_canonical_text(colnames(missing_mask))
    data_names <- .m13a_canonical_text(rownames(sample_data))
    column_names <- .m13a_canonical_text(names(sample_data))
    if (anyDuplicated(feature_names) || anyDuplicated(sample_names) ||
        anyDuplicated(data_names) || anyDuplicated(column_names) ||
        !setequal(sample_names, data_names)) {
        .m13a_fail("M13c canonical hash identifiers must remain unique/aligned.")
    }
    rownames(sample_data) <- data_names
    names(sample_data) <- column_names
    feature_order <- order(feature_names, method = "radix")
    sample_order <- order(sample_names, method = "radix")
    feature_names <- feature_names[feature_order]
    sample_names <- sample_names[sample_order]
    missing_mask <- missing_mask[feature_order, sample_order, drop = FALSE]
    rownames(missing_mask) <- feature_names
    colnames(missing_mask) <- sample_names
    sample_data <- sample_data[match(sample_names, data_names), , drop = FALSE]
    rownames(sample_data) <- sample_names
    roles <- .m13a_normalise_hash_roles(roles, sample_data)
    declared <- unname(unlist(roles, use.names = FALSE))
    interactions <- .m13a_normalise_interactions(interactions, declared)

    mask_bits <- as.vector(t(missing_mask))
    padding <- as.integer((-length(mask_bits)) %% 8L)
    mask_bytes <- packBits(c(mask_bits, rep(FALSE, padding)), type = "raw")
    design_bytes <- unlist(lapply(declared, function(column) {
        .m13a_encode_design_column(column, sample_data[[column]])
    }), use.names = FALSE)
    c(
        charToRaw("imputefinder:association-mask-design:be:v1"), as.raw(0L),
        as.raw(1L), .m13a_be_i32(dim(missing_mask)),
        as.raw(2L), .m13a_encode_text_vector(feature_names),
        as.raw(3L), .m13a_encode_text_vector(sample_names),
        as.raw(4L), .m13a_be_i32(length(mask_bits)),
        .m13a_be_i32(length(mask_bytes)), mask_bytes,
        as.raw(5L), .m13a_encode_hash_roles(roles),
        as.raw(6L), .m13a_encode_hash_interactions(interactions),
        as.raw(7L), .m13a_be_i32(length(declared)), design_bytes
    )
}

m13a_mask_design_sha256 <- function(
    missing_mask,
    sample_data,
    roles,
    interactions = list()
) {
    unname(tools::sha256sum(bytes = m13a_mask_design_bytes(
        missing_mask, sample_data, roles, interactions
    )))
}

m13a_permutation_map_sha256 <- function(map) {
    if (!is.integer(map) || length(map) == 0L || anyNA(map) ||
        !identical(sort(map, method = "radix"), seq_along(map))) {
        .m13a_fail("M13c permutation maps must be exact integer bijections.")
    }
    bytes <- c(
        charToRaw("imputefinder:association-permutation-map:v1"), as.raw(0L),
        .m13a_be_i32(length(map)), .m13a_be_i32(map)
    )
    unname(tools::sha256sum(bytes = bytes))
}

.m13a_canonical_groups <- function(groups, n) {
    if (!is.list(groups) || !length(groups) ||
        !all(vapply(groups, is.integer, logical(1L))) ||
        any(lengths(groups) == 0L)) {
        .m13a_fail("M13c independent Monte Carlo groups are malformed.")
    }
    groups <- lapply(groups, sort, method = "radix")
    first <- vapply(groups, `[[`, integer(1L), 1L)
    groups <- groups[order(first, method = "radix")]
    if (!identical(sort(unlist(groups), method = "radix"), seq_len(n))) {
        .m13a_fail("M13c independent Monte Carlo groups must partition rows.")
    }
    groups
}

.m13a_canonical_swaps <- function(swaps, n) {
    identity <- seq_len(n)
    valid <- is.list(swaps) && length(swaps) > 0L &&
        all(vapply(swaps, function(map) {
            is.integer(map) && length(map) == n && !anyNA(map) &&
                identical(sort(map, method = "radix"), identity) &&
                identical(map[map], identity) && !identical(map, identity)
        }, logical(1L)))
    if (!valid) {
        .m13a_fail("M13c blocked Monte Carlo swaps are malformed.")
    }
    moved <- lapply(swaps, function(map) which(map != identity))
    if (anyDuplicated(unlist(moved))) {
        .m13a_fail("M13c blocked Monte Carlo swaps must move disjoint rows.")
    }
    swaps
}

.m13a_restore_rng <- function(kind, had_seed, seed) {
    do.call(RNGkind, as.list(kind))
    if (had_seed) {
        assign(".Random.seed", seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
    }
}

m13a_monte_carlo_maps <- function(
    seed_manifest,
    design,
    n,
    groups = NULL,
    swaps = NULL
) {
    seed_fields <- c(
        "protocol_id", "candidate_id", "acquisition", "hypothesis_id",
        "draw_id", "seed", "seed_nonce"
    )
    valid <- is.data.frame(seed_manifest) &&
        identical(names(seed_manifest), seed_fields) && nrow(seed_manifest) > 0L &&
        all(vapply(seed_manifest[c(
            "protocol_id", "candidate_id", "acquisition", "hypothesis_id"
        )], is.character, logical(1L))) &&
        !anyNA(seed_manifest) &&
        all(seed_manifest$protocol_id == .M13A_PROTOCOL_ID) &&
        is.integer(seed_manifest$draw_id) &&
        identical(seed_manifest$draw_id, sort(seed_manifest$draw_id, method = "radix")) &&
        all(seed_manifest$draw_id > 0L) && !anyDuplicated(seed_manifest$draw_id) &&
        length(unique(seed_manifest$candidate_id)) == 1L &&
        identical(unique(seed_manifest$candidate_id), "a_fraction_freedman_lane") &&
        length(unique(seed_manifest$acquisition)) == 1L &&
        nzchar(unique(seed_manifest$acquisition)) &&
        length(unique(seed_manifest$hypothesis_id)) == 1L &&
        grepl("^a_[0-9a-f]{64}$", unique(seed_manifest$hypothesis_id)) &&
        is.integer(seed_manifest$seed) && all(seed_manifest$seed > 0L) &&
        !anyDuplicated(seed_manifest$seed) && is.integer(seed_manifest$seed_nonce) &&
        all(seed_manifest$seed_nonce >= 0L) && is.integer(n) && length(n) == 1L &&
        !is.na(n) && n > 0L && .m13a_scalar_character(design) &&
        design %in% c("independent", "blocked")
    if (!valid) {
        .m13a_fail("M13c Monte Carlo seed rows are malformed.")
    }
    if (identical(design, "independent")) {
        groups <- .m13a_canonical_groups(groups, n)
        if (sum(lgamma(lengths(groups) + 1)) <= log(100000)) {
            .m13a_fail("M13c Monte Carlo requires more than 100000 maps.")
        }
    } else {
        swaps <- .m13a_canonical_swaps(swaps, n)
        if (length(swaps) * log(2) <= log(100000)) {
            .m13a_fail("M13c Monte Carlo requires more than 100000 maps.")
        }
    }

    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(.m13a_restore_rng(old_kind, had_seed, old_seed), add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")

    identity <- seq_len(n)
    maps <- vector("list", nrow(seed_manifest))
    map_retry <- integer(nrow(seed_manifest))
    map_sha256 <- character(nrow(seed_manifest))
    accepted <- new.env(hash = TRUE, parent = emptyenv())
    for (index in seq_len(nrow(seed_manifest))) {
        set.seed(seed_manifest$seed[[index]])
        retry <- 0L
        repeat {
            map <- identity
            if (identical(design, "independent")) {
                for (group in groups) {
                    map[group] <- group[sample.int(length(group))]
                }
            } else {
                selected <- sample.int(2L, length(swaps), replace = TRUE) - 1L
                for (swap_index in which(selected == 1L)) {
                    map <- swaps[[swap_index]][map]
                }
            }
            hash <- m13a_permutation_map_sha256(map)
            duplicate <- exists(hash, envir = accepted, inherits = FALSE)
            if (duplicate && !identical(get(hash, envir = accepted), map)) {
                .m13a_fail("M13c permutation-map SHA-256 collision.")
            }
            if (!identical(map, identity) && !duplicate) {
                assign(hash, map, envir = accepted)
                break
            }
            retry <- retry + 1L
            if (retry > 1000000L) {
                .m13a_fail("M13c Monte Carlo map retry limit exceeded.")
            }
        }
        maps[[index]] <- map
        map_retry[[index]] <- retry
        map_sha256[[index]] <- hash
    }
    names(maps) <- as.character(seed_manifest$draw_id)
    manifest <- cbind(
        seed_manifest,
        map_retry = map_retry,
        map_sha256 = map_sha256,
        stringsAsFactors = FALSE
    )
    list(manifest = manifest, maps = maps)
}

m13a_multiplicity <- function() {
    data.frame(
        family_id = "a_encoded_coefficients_holm_v3",
        family = "one acquisition stratum x one candidate x all available eligible 1-df coefficient hypotheses",
        method = "Holm",
        level = 0.05,
        unavailable = "retained as structured output but absent from adjustment denominator",
        unsupported = "retained as structured output but absent from adjustment denominator",
        flag = "adjusted_p<=0.05",
        arithmetic = "stats::p.adjust(raw_p,method=holm) in canonical hypothesis order",
        stringsAsFactors = FALSE
    )
}

m13a_result_schema <- function() {
    rows <- list(
        association = c(
            schema = "character[1],exact:sentinel_association_v1",
            protocol = "list[9],exact:protocol",
            candidate = "character[1],candidate-id",
            methods = "named-character[5],exact:methods",
            response = "data.frame,response-schema,canonical-stratum/sample-order",
            hypotheses = "data.frame,hypotheses-schema,canonical-stratum/coefficient-order",
            support = "data.frame,support-schema,one-row-per-hypothesis",
            outcomes = "named-list,one-entry-per-hypothesis,same-order-and-names",
            multiplicity = "data.frame,multiplicity-schema,one-row-per-stratum",
            diagnostics = "list[4],exact:diagnostics"
        ),
        protocol = c(
            id = "character[1],exact:m13_a_association_protocol_v4",
            hash = "character[1],sha256:m13a-protocol-hash",
            contract_hash = "character[1],sha256:m13a-contract-hash",
            gate_registry_hash = "character[1],sha256:m13a-gate-registry-hash",
            implementation_hash = "character[1],sha256:locked-implementation-manifest",
            candidate_evidence_hash = "character[1],sha256:locked-candidate-evidence",
            winner_state = "character[1]:winner_locked|development_pass|confirmation_pass",
            m12_contract_hash = "character[1],sha256:bound-upstream",
            m12_gate_registry_hash = "character[1],sha256:bound-upstream"
        ),
        association_unavailable = c(
            status = "character[1],exact:unavailable",
            quantity = "character[1],exact:association",
            code = "character[1]:association_no_observable_features|association_no_testable_hypotheses",
            message = "character[1],nonempty",
            requires = "character,nonempty-unique"
        ),
        methods = c(
            response = "character[1],exact:a_sample_detection_fraction_v3",
            encoding = "character[1],exact:canonical_treatment_contrasts_v1",
            rank = "character[1],exact:svd_relative_v1",
            candidate = "character[1],same-as-association-candidate",
            multiplicity = "character[1],exact:Holm"
        ),
        response = c(
            stratum = "character,unique-level:all|s_[0-9a-f]{64}",
            acquisition = "character,nonempty;literal:undeclared-when-absent",
            sample = "character,unique-within-analysis",
            globally_observable_count = "integer,constant-positive",
            detected_count = "integer,0..globally_observable_count",
            detection_fraction = "double,exact-detected/count"
        ),
        hypotheses = c(
            hypothesis = "character,unique,pattern:a_[0-9a-f]{64}",
            stratum = "character,response-foreign-key",
            coefficient = "character,stratum-model-coefficient-foreign-key",
            label = "character,canonical-model-label",
            term_id = "character,stratum-model-term-foreign-key",
            term = "character,canonical-model-term",
            kind = "character:main|interaction",
            components = "AsIs-list-of-character,canonical-component-order",
            component_encodings = "AsIs-list-of-named-character:numeric|treatment,exact-components-and-order",
            role = "character:condition|nuisance|interaction",
            eligible = "logical,role-rule;always-TRUE-in-this-table",
            estimable = "logical,row-space-rule",
            family_member = "logical,TRUE-iff-available-finite-raw-p"
        ),
        support = c(
            hypothesis = "character,hypotheses-foreign-key,same-order",
            design = "character:independent|blocked",
            side_definition = "character,nonempty-canonical-description",
            side_count_low = "integer,nonnegative",
            side_count_high = "integer,nonnegative",
            complete_block_count = "integer,nonnegative-or-NA-independent",
            cell_count_min = "integer,nonnegative-or-NA-main",
            numeric_components = "AsIs-list-of-data.frame(component,median,minimum,maximum,extrapolates_0_1),one-row-per-numeric-component-canonical",
            eligible = "logical,common-support-rule",
            code = "character,NA-iff-eligible-else-support-code"
        ),
        available_outcome = c(
            status = "character[1],exact:available",
            quantity = "character[1],same-as-outcome-list-key",
            candidate = "character[1],same-as-association-candidate",
            stratum = "character[1],hypothesis-row-value",
            hypothesis = "character[1],same-as-outcome-list-key",
            effect = "double[1],finite-probability-units",
            standard_error = "double[1],finite-positive",
            conf_low = "double[1],finite<=effect",
            conf_high = "double[1],finite>=effect",
            statistic = "double[1],finite",
            reference_df = "double[1],finite>=3",
            raw_p = "double[1],finite-[0,1]",
            adjusted_p = "double[1],finite-[0,1]",
            flag = "logical[1],exact-adjusted_p<=0.05",
            log_odds = "double[1],quasi-finite-else-NA",
            log_odds_standard_error = "double[1],quasi-finite-positive-else-NA",
            log_odds_conf_low = "double[1],quasi-finite-else-NA",
            log_odds_conf_high = "double[1],quasi-finite-else-NA",
            permutation_count = "integer[1],FL>=20-else-NA",
            diagnostics = "list[12],exact:outcome_diagnostics"
        ),
        outcome_diagnostics = c(
            variance_method = "character[1]:HC3|CR2|quasibinomial",
            rank = "integer[1],positive",
            residual_df = "integer[1],>=3",
            leverage_max = "double[1],HC3/CR2-[0,1]-else-NA",
            satterthwaite_df = "double[1],CR2>=5-else-NA",
            dispersion = "double[1],quasi-positive-else-NA",
            converged = "logical[1],quasi-TRUE-else-NA",
            boundary = "logical[1],quasi-FALSE-else-NA",
            allowable_transformations = "double[1],FL-integer>=20-else-NA",
            evaluated_transformations = "integer[1],FL-positive-else-NA",
            exceedance_count = "integer[1],FL-0..evaluated-else-NA",
            permutation_mode = "character[1]:none|exact|monte_carlo"
        ),
        multiplicity = c(
            stratum = "character,unique,response-foreign-key",
            method = "character,exact:Holm",
            level = "double,exact:0.05",
            family_size = "integer,exact-available-count",
            available_count = "integer,nonnegative",
            unavailable_count = "integer,exact-hypothesis-minus-available"
        ),
        diagnostics = c(
            input_sha256 = "character[1],sha256:association_mask_design_be_v1",
            strata = "data.frame,stratum_diagnostics-schema",
            seed_manifest = "data.frame,seed_manifest-schema;typed-empty-unless-MC-FL",
            warnings = "character,unique-canonical-codes"
        ),
        stratum_diagnostics = c(
            stratum = "character,unique",
            acquisition = "character,nonempty",
            sample_count = "integer,positive",
            rank = "integer,positive",
            coefficient_count = "integer,positive",
            residual_df = "integer,may-be<3-before-abstention",
            testable_count = "integer,nonnegative",
            unavailable_count = "integer,nonnegative"
        ),
        seed_manifest = c(
            protocol_id = "character,exact:m13_a_association_protocol_v4",
            candidate_id = "character,exact:a_fraction_freedman_lane",
            acquisition = "character,nonempty",
            hypothesis_id = "character,hypotheses-foreign-key",
            draw_id = "integer,positive",
            seed = "integer,positive-unique-over-complete-manifest",
            seed_nonce = "integer,nonnegative-global-seed-collision-count",
            map_retry = "integer,0-to-1000000-within-draw-stream-retry-count",
            map_sha256 = "character,sha256:canonical-row-index-map"
        )
    )
    do.call(rbind, lapply(names(rows), function(record) {
        fields <- rows[[record]]
        data.frame(
            record = record,
            position = seq_along(fields),
            field = names(fields),
            rule = unname(fields),
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    }))
}

m13a_result_records <- function() {
    data.frame(
        record = c(
            "association", "association_unavailable", "protocol", "methods", "response", "hypotheses",
            "support", "outcomes", "available_outcome", "unavailable_outcome",
            "outcome_diagnostics", "multiplicity", "diagnostics",
            "stratum_diagnostics", "seed_manifest"
        ),
        container = c(
            "list", "list", "list", "named_character", rep("data.frame", 3L), "named_list",
            "list", "list", "list", "data.frame", "list", "data.frame", "data.frame"
        ),
        class = c(
            "imputefinder_association_panel", "imputefinder_unavailable", "none", "none", rep("data.frame", 3L),
            "none", "imputefinder_association", "imputefinder_unavailable",
            "none", "data.frame", "none", "data.frame", "data.frame"
        ),
        key = c(
            "one panel", "quantity=association", "one panel", "one panel", "stratum+sample",
            "hypothesis", "hypothesis", "hypothesis-list-name",
            "quantity=hypothesis-list-name", "quantity=hypothesis-list-name",
            "one per available outcome", "stratum", "one panel", "stratum",
            "candidate+acquisition+hypothesis+draw"
        ),
        order = c(
            "exact fields", "existing exact unavailable fields", "exact fields", "exact names", "canonical",
            "canonical", "hypothesis order", "hypothesis order", "exact fields",
            "existing exact unavailable fields", "exact fields", "stratum order",
            "exact fields", "stratum order", "canonical manifest order"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_lifecycle <- function() {
    data.frame(
        sentinel_fields = c(
            "pre_rescue", "pre_rescue;coverage",
            "pre_rescue;coverage;association"
        ),
        readability = c(
            "legacy_read", "legacy_read", "current_read_write"
        ),
        association_rule = c(
            "absent", "absent",
            "exact sentinel_association_v1 or exact top-level association unavailable; reject unknown schema; candidate comparison tables forbidden"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_hypothesis_id <- function(stratum, coefficient, label, term_id) {
    values <- c(stratum, coefficient, label, term_id)
    if (!is.character(values) || length(values) != 4L || anyNA(values) ||
        any(!nzchar(values))) {
        .m13a_fail("M13c hypothesis identity fields must be nonempty scalars.")
    }
    bytes <- c(
        charToRaw("imputefinder:association-hypothesis:v1"), as.raw(0L),
        .m13a_encode_text_vector(values)
    )
    paste0("a_", unname(tools::sha256sum(bytes = bytes)))
}

.m13a_schema_fields <- function(record) {
    schema <- m13a_result_schema()
    schema$field[schema$record == record]
}

.m13a_scalar_character <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.m13a_sha256 <- function(x) {
    .m13a_scalar_character(x) && grepl("^[0-9a-f]{64}$", x)
}

.m13a_double_scalar <- function(x, finite = TRUE) {
    is.double(x) && length(x) == 1L &&
        if (finite) is.finite(x) else is.na(x)
}

.m13a_integer_scalar <- function(x, minimum = NULL, na = FALSE) {
    valid <- is.integer(x) && length(x) == 1L
    if (!valid) return(FALSE)
    if (na) return(is.na(x))
    !is.na(x) && (is.null(minimum) || x >= minimum)
}

.m13a_numeric_close <- function(x, y) {
    is.finite(x) && is.finite(y) &&
        abs(x - y) <= sqrt(.Machine$double.eps) * max(1, abs(x), abs(y))
}

.m13a_validate_outcome_diagnostics <- function(outcome, candidate) {
    diagnostics <- outcome$diagnostics
    valid <- is.list(diagnostics) &&
        identical(names(diagnostics), .m13a_schema_fields("outcome_diagnostics")) &&
        .m13a_scalar_character(diagnostics$variance_method) &&
        diagnostics$variance_method %in% c("HC3", "CR2", "quasibinomial") &&
        .m13a_integer_scalar(diagnostics$rank, 1L) &&
        .m13a_integer_scalar(diagnostics$residual_df, 3L) &&
        .m13a_double_scalar(diagnostics$leverage_max) ==
            diagnostics$variance_method %in% c("HC3", "CR2") &&
        (!is.finite(diagnostics$leverage_max) ||
            diagnostics$leverage_max >= 0 && diagnostics$leverage_max <= 1) &&
        (.m13a_double_scalar(diagnostics$satterthwaite_df) &&
            diagnostics$satterthwaite_df >= 5) ==
            identical(diagnostics$variance_method, "CR2") &&
        (.m13a_double_scalar(diagnostics$dispersion) &&
            diagnostics$dispersion > 0) ==
            identical(diagnostics$variance_method, "quasibinomial") &&
        is.logical(diagnostics$converged) && length(diagnostics$converged) == 1L &&
        is.logical(diagnostics$boundary) && length(diagnostics$boundary) == 1L &&
        .m13a_scalar_character(diagnostics$permutation_mode) &&
        diagnostics$permutation_mode %in% c("none", "exact", "monte_carlo")
    if (!valid) return(FALSE)

    quasi <- identical(candidate, "a_fraction_quasibinomial")
    freedman_lane <- identical(candidate, "a_fraction_freedman_lane")
    if (quasi) {
        if (!identical(diagnostics$variance_method, "quasibinomial") ||
            !identical(diagnostics$converged, TRUE) ||
            !identical(diagnostics$boundary, FALSE) ||
            !.m13a_double_scalar(diagnostics$leverage_max, FALSE) ||
            !.m13a_double_scalar(diagnostics$satterthwaite_df, FALSE)) {
            return(FALSE)
        }
    } else if (!diagnostics$variance_method %in% c("HC3", "CR2") ||
        !is.na(diagnostics$converged) || !is.na(diagnostics$boundary) ||
        !.m13a_double_scalar(diagnostics$dispersion, FALSE)) {
        return(FALSE)
    }

    expected_df <- if (identical(diagnostics$variance_method, "CR2")) {
        diagnostics$satterthwaite_df
    } else as.double(diagnostics$residual_df)
    if (!.m13a_numeric_close(outcome$reference_df, expected_df)) return(FALSE)

    if (!freedman_lane) {
        return(
            .m13a_integer_scalar(outcome$permutation_count, na = TRUE) &&
            .m13a_double_scalar(diagnostics$allowable_transformations, FALSE) &&
            .m13a_integer_scalar(diagnostics$evaluated_transformations, na = TRUE) &&
            .m13a_integer_scalar(diagnostics$exceedance_count, na = TRUE) &&
            identical(diagnostics$permutation_mode, "none")
        )
    }

    allowable <- diagnostics$allowable_transformations
    if (!.m13a_double_scalar(allowable) || allowable < 20 ||
        allowable != floor(allowable) ||
        !.m13a_integer_scalar(diagnostics$evaluated_transformations, 1L) ||
        !.m13a_integer_scalar(diagnostics$exceedance_count, 0L) ||
        !.m13a_integer_scalar(outcome$permutation_count, 20L)) {
        return(FALSE)
    }
    if (allowable <= 100000) {
        valid <- identical(diagnostics$permutation_mode, "exact") &&
            diagnostics$evaluated_transformations == allowable &&
            outcome$permutation_count == diagnostics$evaluated_transformations &&
            diagnostics$exceedance_count >= 1L &&
            diagnostics$exceedance_count <= diagnostics$evaluated_transformations
        expected_p <- diagnostics$exceedance_count /
            diagnostics$evaluated_transformations
    } else {
        valid <- identical(diagnostics$permutation_mode, "monte_carlo") &&
            diagnostics$evaluated_transformations == 9999L &&
            outcome$permutation_count == 10000L &&
            diagnostics$exceedance_count <= 9999L
        expected_p <- (1 + diagnostics$exceedance_count) / 10000
    }
    valid && identical(outcome$raw_p, as.double(expected_p))
}

.m13a_validate_available_outcome <- function(
    outcome,
    key,
    candidate,
    hypothesis_row
) {
    valid <- is.list(outcome) &&
        identical(class(outcome), "imputefinder_association") &&
        identical(names(outcome), .m13a_schema_fields("available_outcome")) &&
        identical(outcome$status, "available") &&
        identical(outcome$quantity, key) && identical(outcome$hypothesis, key) &&
        identical(outcome$candidate, candidate) &&
        identical(outcome$stratum, hypothesis_row$stratum) &&
        all(vapply(outcome[c(
            "effect", "standard_error", "conf_low", "conf_high", "statistic",
            "reference_df", "raw_p", "adjusted_p"
        )], function(x) is.double(x) && length(x) == 1L && is.finite(x), logical(1L))) &&
        outcome$standard_error > 0 && outcome$conf_low <= outcome$effect &&
        outcome$conf_high >= outcome$effect && outcome$reference_df >= 3 &&
        outcome$raw_p >= 0 && outcome$raw_p <= 1 &&
        outcome$adjusted_p >= 0 && outcome$adjusted_p <= 1 &&
        is.logical(outcome$flag) && length(outcome$flag) == 1L &&
        !is.na(outcome$flag) && identical(outcome$flag, outcome$adjusted_p <= 0.05) &&
        .m13a_validate_outcome_diagnostics(outcome, candidate)
    secondary <- unlist(outcome[c(
        "log_odds", "log_odds_standard_error", "log_odds_conf_low",
        "log_odds_conf_high"
    )], use.names = FALSE)
    quasi <- identical(candidate, "a_fraction_quasibinomial")
    secondary_valid <- is.double(secondary) &&
        if (quasi) {
            all(is.finite(secondary)) && outcome$log_odds_standard_error > 0 &&
                outcome$log_odds_conf_low <= outcome$log_odds &&
                outcome$log_odds_conf_high >= outcome$log_odds &&
                .m13a_numeric_close(
                    outcome$statistic,
                    outcome$log_odds / outcome$log_odds_standard_error
                )
        } else all(is.na(secondary))
    statistic_valid <- if (identical(candidate, "a_fraction_freedman_lane")) {
        outcome$statistic >= 0 && .m13a_numeric_close(
            outcome$statistic,
            (outcome$effect / outcome$standard_error)^2
        )
    } else if (!quasi) {
        .m13a_numeric_close(
            outcome$statistic,
            outcome$effect / outcome$standard_error
        )
    } else TRUE
    parametric_p_valid <- if (identical(candidate, "a_fraction_freedman_lane")) {
        TRUE
    } else {
        .m13a_numeric_close(
            outcome$raw_p,
            2 * stats::pt(-abs(outcome$statistic), outcome$reference_df)
        )
    }
    valid && secondary_valid && statistic_valid && parametric_p_valid
}

.m13a_validate_unavailable_outcome <- function(outcome, key, codes) {
    is.list(outcome) &&
        identical(class(outcome), "imputefinder_unavailable") &&
        identical(
            names(outcome),
            c("status", "quantity", "code", "message", "requires")
        ) &&
        identical(outcome$status, "unavailable") &&
        identical(outcome$quantity, key) && outcome$code %in% codes &&
        .m13a_scalar_character(outcome$message) &&
        is.character(outcome$requires) && length(outcome$requires) > 0L &&
        !anyNA(outcome$requires) && all(nzchar(outcome$requires)) &&
        !anyDuplicated(outcome$requires)
}

m13a_validate_result_shape <- function(result) {
    panel_codes <- c(
        "association_no_observable_features",
        "association_no_testable_hypotheses"
    )
    if (inherits(result, "imputefinder_unavailable")) {
        if (!.m13a_validate_unavailable_outcome(
            result,
            "association",
            panel_codes
        )) {
            .m13a_fail("M13c top-level association abstention is malformed.")
        }
        return(invisible(result))
    }
    candidates <- m13a_candidate_engines()$candidate_id
    expected_methods <- c(
        response = "a_sample_detection_fraction_v3",
        encoding = "canonical_treatment_contrasts_v1",
        rank = "svd_relative_v1",
        candidate = result$candidate,
        multiplicity = "Holm"
    )
    valid <- is.list(result) &&
        identical(class(result), "imputefinder_association_panel") &&
        identical(names(result), .m13a_schema_fields("association")) &&
        identical(result$schema, "sentinel_association_v1") &&
        .m13a_scalar_character(result$candidate) && result$candidate %in% candidates &&
        is.list(result$protocol) &&
        identical(names(result$protocol), .m13a_schema_fields("protocol")) &&
        identical(result$protocol$id, .M13A_PROTOCOL_ID) &&
        identical(result$protocol$hash, m13a_protocol_hash()) &&
        identical(result$protocol$contract_hash,
            unname(.M13A_EXPECTED_HASHES[["contract"]])) &&
        identical(result$protocol$gate_registry_hash,
            unname(.M13A_EXPECTED_HASHES[["gate_registry"]])) &&
        .m13a_sha256(result$protocol$implementation_hash) &&
        .m13a_sha256(result$protocol$candidate_evidence_hash) &&
        .m13a_scalar_character(result$protocol$winner_state) &&
        result$protocol$winner_state %in% c(
            "winner_locked", "development_pass", "confirmation_pass"
        ) &&
        identical(result$protocol$m12_contract_hash, .M13A_EXPECTED_UPSTREAM[["m12_contract"]]) &&
        identical(result$protocol$m12_gate_registry_hash, .M13A_EXPECTED_UPSTREAM[["m12_gate_registry"]]) &&
        is.character(result$methods) &&
        identical(result$methods, expected_methods)
    if (!valid) {
        .m13a_fail("M13c association panel header is malformed.")
    }

    response <- result$response
    response_valid <- is.data.frame(response) && nrow(response) > 0L &&
        identical(names(response), .m13a_schema_fields("response")) &&
        all(vapply(response[c("stratum", "acquisition", "sample")], is.character, logical(1L))) &&
        is.integer(response$globally_observable_count) &&
        is.integer(response$detected_count) && is.double(response$detection_fraction) &&
        !anyNA(response) &&
        all(nzchar(response$stratum)) && all(nzchar(response$acquisition)) &&
        all(nzchar(response$sample)) && !anyDuplicated(response$sample) &&
        all(response$stratum == "all" |
            grepl("^s_[0-9a-f]{64}$", response$stratum)) &&
        all(response$globally_observable_count > 0L) &&
        identical(
            response$globally_observable_count,
            rep(response$globally_observable_count[[1L]], nrow(response))
        ) &&
        all(response$detected_count >= 0L) &&
        all(response$detected_count <= response$globally_observable_count) &&
        identical(
            response$detection_fraction,
            response$detected_count / response$globally_observable_count
        )
    if (response_valid) {
        stratum_rows <- !duplicated(response$stratum)
        stratum_map <- response[stratum_rows, c("stratum", "acquisition"), drop = FALSE]
        one_acquisition <- vapply(split(
            response$acquisition,
            response$stratum
        ), function(x) length(unique(x)) == 1L, logical(1L))
        undeclared <- identical(stratum_map$stratum, "all") &&
            identical(stratum_map$acquisition, "undeclared")
        declared <- !any(stratum_map$stratum == "all") &&
            identical(
                stratum_map$stratum,
                unname(vapply(
                    stratum_map$acquisition,
                    m13a_stratum_id,
                    character(1L)
                ))
            )
        canonical <- do.call(
            order,
            list(response$acquisition, response$sample, method = "radix")
        )
        response_valid <- all(one_acquisition) && (undeclared || declared) &&
            identical(canonical, seq_len(nrow(response)))
    }
    if (!response_valid) {
        .m13a_fail("M13c association response is malformed.")
    }

    hypotheses <- result$hypotheses
    hypothesis_valid <- is.data.frame(hypotheses) && nrow(hypotheses) > 0L &&
        identical(names(hypotheses), .m13a_schema_fields("hypotheses")) &&
        is.character(hypotheses$hypothesis) &&
        all(grepl("^a_[0-9a-f]{64}$", hypotheses$hypothesis)) &&
        !anyDuplicated(hypotheses$hypothesis) &&
        all(vapply(hypotheses[c(
            "stratum", "coefficient", "label", "term_id", "term", "kind", "role"
        )], is.character, logical(1L))) &&
        !anyNA(hypotheses[c(
            "stratum", "coefficient", "label", "term_id", "term", "kind", "role"
        )]) &&
        all(vapply(hypotheses[c(
            "stratum", "coefficient", "label", "term_id", "term", "kind", "role"
        )], function(x) all(nzchar(x)), logical(1L))) &&
        is.list(hypotheses$components) && inherits(hypotheses$components, "AsIs") &&
        all(vapply(hypotheses$components, function(x) {
            is.character(x) && length(x) > 0L && !anyNA(x) &&
                all(nzchar(x)) && !anyDuplicated(x) &&
                identical(x, sort(x, method = "radix"))
        }, logical(1L))) &&
        is.list(hypotheses$component_encodings) &&
        inherits(hypotheses$component_encodings, "AsIs") &&
        all(vapply(seq_len(nrow(hypotheses)), function(index) {
            encoding <- hypotheses$component_encodings[[index]]
            is.character(encoding) &&
                identical(names(encoding), hypotheses$components[[index]]) &&
                all(encoding %in% c("numeric", "treatment"))
        }, logical(1L))) &&
        all(hypotheses$kind %in% c("main", "interaction")) &&
        all(hypotheses$role %in% c("condition", "nuisance", "interaction")) &&
        all((hypotheses$kind == "main") ==
            (hypotheses$role %in% c("condition", "nuisance"))) &&
        all(lengths(hypotheses$components)[hypotheses$kind == "main"] == 1L) &&
        all(lengths(hypotheses$components)[hypotheses$kind == "interaction"] >= 2L) &&
        all(grepl("^coef_[0-9]{4,}$", hypotheses$coefficient)) &&
        all(grepl("^term_[0-9]{4,}$", hypotheses$term_id)) &&
        is.logical(hypotheses$eligible) && all(hypotheses$eligible) &&
        is.logical(hypotheses$estimable) && is.logical(hypotheses$family_member) &&
        !anyNA(hypotheses$estimable) && !anyNA(hypotheses$family_member) &&
        all(hypotheses$stratum %in% response$stratum) &&
        !anyDuplicated(hypotheses[c("stratum", "coefficient")])
    if (hypothesis_valid && nrow(hypotheses)) {
        expected_ids <- vapply(seq_len(nrow(hypotheses)), function(index) {
            m13a_hypothesis_id(
                hypotheses$stratum[[index]], hypotheses$coefficient[[index]],
                hypotheses$label[[index]], hypotheses$term_id[[index]]
            )
        }, character(1L))
        strata <- unique(response$stratum)
        canonical <- do.call(order, list(
            match(hypotheses$stratum, strata),
            hypotheses$coefficient,
            method = "radix"
        ))
        hypothesis_valid <- identical(hypotheses$hypothesis, expected_ids) &&
            identical(canonical, seq_len(nrow(hypotheses)))
    }
    support <- result$support
    code_catalog <- m13a_unavailable_codes()
    support_codes <- code_catalog$code[code_catalog$stage == "support"]
    support_valid <- is.data.frame(support) &&
        identical(names(support), .m13a_schema_fields("support")) &&
        identical(support$hypothesis, hypotheses$hypothesis) &&
        all(vapply(support[c("hypothesis", "design", "side_definition")],
            is.character, logical(1L))) &&
        !anyNA(support[c("hypothesis", "design", "side_definition")]) &&
        all(nzchar(support$side_definition)) &&
        all(support$design %in% c("independent", "blocked")) &&
        all(vapply(support[c(
            "side_count_low", "side_count_high", "complete_block_count",
            "cell_count_min"
        )], is.integer, logical(1L))) &&
        all(support$side_count_low >= 0L) && all(support$side_count_high >= 0L) &&
        all(is.na(support$complete_block_count) |
            support$complete_block_count >= 0L) &&
        all(is.na(support$cell_count_min) | support$cell_count_min >= 0L) &&
        all(is.na(support$complete_block_count[support$design == "independent"])) &&
        all(!is.na(support$complete_block_count[support$design == "blocked"])) &&
        all(is.na(support$cell_count_min[hypotheses$kind == "main"])) &&
        all(!is.na(support$cell_count_min[hypotheses$kind == "interaction"])) &&
        is.list(support$numeric_components) &&
        inherits(support$numeric_components, "AsIs") &&
        all(vapply(seq_len(nrow(support)), function(index) {
            numeric <- support$numeric_components[[index]]
            is.data.frame(numeric) && identical(
                names(numeric),
                c(
                    "component", "median", "minimum", "maximum",
                    "extrapolates_0_1"
                )
            ) && is.character(numeric$component) &&
                is.double(numeric$median) && is.double(numeric$minimum) &&
                is.double(numeric$maximum) &&
                is.logical(numeric$extrapolates_0_1) && !anyNA(numeric) &&
                all(nzchar(numeric$component)) &&
                !anyDuplicated(numeric$component) &&
                identical(
                    numeric$component,
                    names(hypotheses$component_encodings[[index]])[
                        hypotheses$component_encodings[[index]] == "numeric"
                    ]
                ) &&
                all(is.finite(numeric$median)) &&
                all(is.finite(numeric$minimum)) &&
                all(is.finite(numeric$maximum)) &&
                all(numeric$minimum <= numeric$median) &&
                all(numeric$median <= numeric$maximum) &&
                identical(
                    numeric$extrapolates_0_1,
                    numeric$minimum > 0 | numeric$maximum < 1
                )
        }, logical(1L))) &&
        is.logical(support$eligible) && !anyNA(support$eligible) &&
        is.character(support$code) &&
        identical(is.na(support$code), support$eligible) &&
        all(is.na(support$code) | support$code %in% support_codes)
    outcomes <- result$outcomes
    outcome_valid <- is.list(outcomes) &&
        identical(names(outcomes), hypotheses$hypothesis)
    if (!hypothesis_valid || !support_valid || !outcome_valid) {
        .m13a_fail("M13c hypothesis/support/outcome keys are malformed.")
    }

    available <- logical(length(outcomes))
    codes <- m13a_unavailable_codes()$code
    for (index in seq_along(outcomes)) {
        outcome <- outcomes[[index]]
        key <- names(outcomes)[[index]]
        available[[index]] <- inherits(outcome, "imputefinder_association")
        valid_outcome <- if (available[[index]]) {
            .m13a_validate_available_outcome(
                outcome, key, result$candidate,
                hypotheses[index, , drop = FALSE]
            )
        } else {
            .m13a_validate_unavailable_outcome(outcome, key, codes)
        }
        if (!valid_outcome) {
            .m13a_fail("M13c association outcome is malformed: ", key)
        }
    }
    if (!identical(hypotheses$family_member, available) ||
        any(available & (!hypotheses$estimable | !support$eligible))) {
        .m13a_fail("M13c Holm membership must equal finite available outcomes.")
    }
    for (index in which(!available)) {
        code <- outcomes[[index]]$code
        expected <- if (!hypotheses$estimable[[index]]) {
            "association_nonestimable"
        } else if (!support$eligible[[index]]) {
            support$code[[index]]
        } else NULL
        if (!is.null(expected) && !identical(code, expected)) {
            .m13a_fail("M13c unavailable precedence is malformed.")
        }
    }
    for (stratum in unique(hypotheses$stratum)) {
        selected <- hypotheses$stratum == stratum & available
        if (any(selected)) {
            raw <- vapply(outcomes[selected], `[[`, numeric(1L), "raw_p")
            adjusted <- vapply(outcomes[selected], `[[`, numeric(1L), "adjusted_p")
            if (!identical(adjusted, stats::p.adjust(raw, method = "holm"))) {
                .m13a_fail("M13c stored Holm arithmetic is malformed.")
            }
        }
    }

    multiplicity <- result$multiplicity
    strata <- unique(response$stratum)
    multiplicity_valid <- is.data.frame(multiplicity) &&
        identical(names(multiplicity), .m13a_schema_fields("multiplicity")) &&
        identical(multiplicity$stratum, strata) &&
        !anyNA(multiplicity) &&
        all(multiplicity$method == "Holm") && all(multiplicity$level == 0.05) &&
        is.integer(multiplicity$family_size) &&
        is.integer(multiplicity$available_count) &&
        is.integer(multiplicity$unavailable_count) &&
        all(multiplicity$family_size >= 0L) &&
        all(multiplicity$available_count >= 0L) &&
        all(multiplicity$unavailable_count >= 0L)
    if (multiplicity_valid) {
        for (index in seq_along(strata)) {
            selected <- hypotheses$stratum == strata[[index]]
            count <- sum(available[selected])
            multiplicity_valid <- multiplicity_valid &&
                multiplicity$family_size[[index]] == count &&
                multiplicity$available_count[[index]] == count &&
                multiplicity$unavailable_count[[index]] == sum(!available[selected])
        }
    }
    diagnostics <- result$diagnostics
    seed <- diagnostics$seed_manifest
    stratum_diagnostics <- diagnostics$strata
    diagnostic_valid <- is.list(diagnostics) &&
        identical(names(diagnostics), .m13a_schema_fields("diagnostics")) &&
        .m13a_sha256(diagnostics$input_sha256) &&
        is.data.frame(stratum_diagnostics) &&
        identical(names(stratum_diagnostics), .m13a_schema_fields("stratum_diagnostics")) &&
        identical(stratum_diagnostics$stratum, strata) &&
        all(vapply(stratum_diagnostics[c("stratum", "acquisition")],
            is.character, logical(1L))) &&
        all(vapply(stratum_diagnostics[c(
            "sample_count", "rank", "coefficient_count", "residual_df",
            "testable_count", "unavailable_count"
        )], is.integer, logical(1L))) &&
        !anyNA(stratum_diagnostics) &&
        all(stratum_diagnostics$sample_count > 0L) &&
        all(stratum_diagnostics$rank > 0L) &&
        all(stratum_diagnostics$coefficient_count > 0L) &&
        all(stratum_diagnostics$residual_df ==
            stratum_diagnostics$sample_count - stratum_diagnostics$rank) &&
        all(stratum_diagnostics$testable_count >= 0L) &&
        all(stratum_diagnostics$unavailable_count >= 0L) &&
        is.data.frame(seed) &&
        identical(names(seed), .m13a_schema_fields("seed_manifest")) &&
        is.character(diagnostics$warnings) && !anyNA(diagnostics$warnings) &&
        all(nzchar(diagnostics$warnings)) && !anyDuplicated(diagnostics$warnings) &&
        identical(diagnostics$warnings, sort(diagnostics$warnings, method = "radix"))
    if (diagnostic_valid) {
        for (index in seq_along(strata)) {
            stratum <- strata[[index]]
            response_selected <- response$stratum == stratum
            hypothesis_selected <- hypotheses$stratum == stratum
            diagnostic_valid <- diagnostic_valid &&
                identical(
                    stratum_diagnostics$acquisition[[index]],
                    unique(response$acquisition[response_selected])
                ) &&
                stratum_diagnostics$sample_count[[index]] == sum(response_selected) &&
                stratum_diagnostics$coefficient_count[[index]] >=
                    sum(hypothesis_selected) &&
                stratum_diagnostics$testable_count[[index]] ==
                    sum(available[hypothesis_selected]) &&
                stratum_diagnostics$unavailable_count[[index]] ==
                    sum(!available[hypothesis_selected])
        }
    }
    if (diagnostic_valid && nrow(seed)) {
        seed_key <- seed[c(
            "candidate_id", "acquisition", "hypothesis_id", "draw_id"
        )]
        canonical <- do.call(order, c(unname(seed_key), list(method = "radix")))
        monte_carlo <- vapply(outcomes, function(outcome) {
            inherits(outcome, "imputefinder_association") &&
                identical(outcome$diagnostics$permutation_mode, "monte_carlo")
        }, logical(1L))
        mc_ids <- hypotheses$hypothesis[monte_carlo]
        counts <- table(factor(seed$hypothesis_id, levels = mc_ids))
        response_acquisition <- stats::setNames(
            stratum_diagnostics$acquisition,
            stratum_diagnostics$stratum
        )
        expected_acquisition <- unname(response_acquisition[
            hypotheses$stratum[match(seed$hypothesis_id, hypotheses$hypothesis)]
        ])
        diagnostic_valid <- identical(result$candidate, "a_fraction_freedman_lane") &&
            !anyNA(seed) &&
            all(vapply(seed[c(
                "protocol_id", "candidate_id", "acquisition", "hypothesis_id",
                "map_sha256"
            )], is.character, logical(1L))) &&
            all(seed$protocol_id == .M13A_PROTOCOL_ID) &&
            all(seed$candidate_id == result$candidate) &&
            all(seed$hypothesis_id %in% mc_ids) &&
            identical(seed$acquisition, expected_acquisition) &&
            is.integer(seed$draw_id) && all(seed$draw_id > 0L) &&
            is.integer(seed$seed) && all(seed$seed > 0L) &&
            !anyDuplicated(seed$seed) && is.integer(seed$seed_nonce) &&
            all(seed$seed_nonce >= 0L) && is.integer(seed$map_retry) &&
            all(seed$map_retry >= 0L) && all(seed$map_retry <= 1000000L) &&
            all(vapply(seed$map_sha256, .m13a_sha256, logical(1L))) &&
            all(vapply(split(
                seed$map_sha256,
                seed$hypothesis_id
            ), function(hashes) {
                !anyDuplicated(hashes)
            }, logical(1L))) &&
            !anyDuplicated(seed_key) &&
            identical(canonical, seq_len(nrow(seed))) &&
            length(mc_ids) > 0L && all(counts == 9999L) &&
            all(vapply(split(seed$draw_id, seed$hypothesis_id), function(x) {
                identical(x, 1:9999)
            }, logical(1L)))
    } else if (diagnostic_valid) {
        monte_carlo <- vapply(outcomes, function(outcome) {
            inherits(outcome, "imputefinder_association") &&
                identical(outcome$diagnostics$permutation_mode, "monte_carlo")
        }, logical(1L))
        diagnostic_valid <- !any(monte_carlo)
    }
    if (!multiplicity_valid || !diagnostic_valid) {
        .m13a_fail("M13c multiplicity/diagnostic records are malformed.")
    }
    invisible(result)
}

.m13a_result_fixture <- function(available = TRUE) {
    hypothesis <- m13a_hypothesis_id(
        "all", "coef_0002", "condition[B]", "term_0002"
    )
    response <- data.frame(
        stratum = rep("all", 8L),
        acquisition = rep("undeclared", 8L),
        sample = sprintf("s%02d", 1:8),
        globally_observable_count = rep(2L, 8L),
        detected_count = as.integer(c(1, 1, 1, 1, 2, 2, 2, 2)),
        detection_fraction = c(0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1),
        stringsAsFactors = FALSE
    )
    hypotheses <- data.frame(
        hypothesis = hypothesis,
        stratum = "all",
        coefficient = "coef_0002",
        label = "condition[B]",
        term_id = "term_0002",
        term = "condition",
        kind = "main",
        role = "condition",
        eligible = TRUE,
        estimable = TRUE,
        family_member = available,
        stringsAsFactors = FALSE
    )
    hypotheses$components <- I(list("condition"))
    hypotheses$component_encodings <- I(list(c(condition = "treatment")))
    hypotheses <- hypotheses[.m13a_schema_fields("hypotheses")]
    support <- data.frame(
        hypothesis = hypothesis,
        design = "independent",
        side_definition = "condition[A] versus condition[B]",
        side_count_low = 4L,
        side_count_high = 4L,
        complete_block_count = NA_integer_,
        cell_count_min = NA_integer_,
        eligible = TRUE,
        code = NA_character_,
        stringsAsFactors = FALSE
    )
    support$numeric_components <- I(list(data.frame(
        component = character(), median = double(), minimum = double(),
        maximum = double(), extrapolates_0_1 = logical(),
        stringsAsFactors = FALSE
    )))
    support <- support[.m13a_schema_fields("support")]
    outcome <- if (available) {
        structure(
            list(
                status = "available",
                quantity = hypothesis,
                candidate = "a_fraction_ols_hc3_cr2",
                stratum = "all",
                hypothesis = hypothesis,
                effect = 0.5,
                standard_error = 0.1,
                conf_low = 0.255,
                conf_high = 0.745,
                statistic = 5.0,
                reference_df = 6.0,
                raw_p = 2 * stats::pt(-5, 6),
                adjusted_p = 2 * stats::pt(-5, 6),
                flag = TRUE,
                log_odds = NA_real_,
                log_odds_standard_error = NA_real_,
                log_odds_conf_low = NA_real_,
                log_odds_conf_high = NA_real_,
                permutation_count = NA_integer_,
                diagnostics = list(
                    variance_method = "HC3", rank = 2L, residual_df = 6L,
                    leverage_max = 0.25, satterthwaite_df = NA_real_,
                    dispersion = NA_real_, converged = NA, boundary = NA,
                    allowable_transformations = NA_real_,
                    evaluated_transformations = NA_integer_,
                    exceedance_count = NA_integer_, permutation_mode = "none"
                )
            ),
            class = "imputefinder_association"
        )
    } else {
        structure(
            list(
                status = "unavailable", quantity = hypothesis,
                code = "association_singular_covariance",
                message = "The scalar robust covariance is singular.",
                requires = "positive scalar covariance"
            ),
            class = "imputefinder_unavailable"
        )
    }
    multiplicity <- data.frame(
        stratum = "all", method = "Holm", level = 0.05,
        family_size = as.integer(available),
        available_count = as.integer(available),
        unavailable_count = as.integer(!available),
        stringsAsFactors = FALSE
    )
    seed_manifest <- data.frame(
        protocol_id = character(), candidate_id = character(),
        acquisition = character(), hypothesis_id = character(),
        draw_id = integer(), seed = integer(), seed_nonce = integer(),
        map_retry = integer(), map_sha256 = character(),
        stringsAsFactors = FALSE
    )
    diagnostics <- list(
        input_sha256 = paste(rep("0", 64L), collapse = ""),
        strata = data.frame(
            stratum = "all", acquisition = "undeclared", sample_count = 8L,
            rank = 2L, coefficient_count = 2L, residual_df = 6L,
            testable_count = as.integer(available),
            unavailable_count = as.integer(!available),
            stringsAsFactors = FALSE
        ),
        seed_manifest = seed_manifest,
        warnings = character()
    )
    structure(
        list(
            schema = "sentinel_association_v1",
            protocol = list(
                id = .M13A_PROTOCOL_ID,
                hash = m13a_protocol_hash(),
                contract_hash = unname(.M13A_EXPECTED_HASHES[["contract"]]),
                gate_registry_hash = unname(.M13A_EXPECTED_HASHES[["gate_registry"]]),
                implementation_hash = paste(rep("1", 64L), collapse = ""),
                candidate_evidence_hash = paste(rep("2", 64L), collapse = ""),
                winner_state = "winner_locked",
                m12_contract_hash = .M13A_EXPECTED_UPSTREAM[["m12_contract"]],
                m12_gate_registry_hash = .M13A_EXPECTED_UPSTREAM[["m12_gate_registry"]]
            ),
            candidate = "a_fraction_ols_hc3_cr2",
            methods = c(
                response = "a_sample_detection_fraction_v3",
                encoding = "canonical_treatment_contrasts_v1",
                rank = "svd_relative_v1",
                candidate = "a_fraction_ols_hc3_cr2",
                multiplicity = "Holm"
            ),
            response = response,
            hypotheses = hypotheses,
            support = support,
            outcomes = stats::setNames(list(outcome), hypothesis),
            multiplicity = multiplicity,
            diagnostics = diagnostics
        ),
        class = "imputefinder_association_panel"
    )
}

.m13a_monte_carlo_result_fixture <- function() {
    result <- .m13a_result_fixture(TRUE)
    result$candidate <- "a_fraction_freedman_lane"
    result$methods[["candidate"]] <- "a_fraction_freedman_lane"
    outcome <- result$outcomes[[1L]]
    outcome$candidate <- "a_fraction_freedman_lane"
    outcome$statistic <- (outcome$effect / outcome$standard_error)^2
    outcome$raw_p <- 0.0001
    outcome$adjusted_p <- 0.0001
    outcome$permutation_count <- 10000L
    outcome$diagnostics$allowable_transformations <- 362880.0
    outcome$diagnostics$evaluated_transformations <- 9999L
    outcome$diagnostics$exceedance_count <- 0L
    outcome$diagnostics$permutation_mode <- "monte_carlo"
    result$outcomes[[1L]] <- outcome

    draw_id <- seq_len(9999L)
    result$diagnostics$seed_manifest <- data.frame(
        protocol_id = rep(.M13A_PROTOCOL_ID, 9999L),
        candidate_id = rep("a_fraction_freedman_lane", 9999L),
        acquisition = rep("undeclared", 9999L),
        hypothesis_id = rep(result$hypotheses$hypothesis, 9999L),
        draw_id = draw_id,
        seed = draw_id,
        seed_nonce = integer(9999L),
        map_retry = integer(9999L),
        map_sha256 = paste0(
            strrep("0", 56L),
            sprintf("%08x", draw_id)
        ),
        stringsAsFactors = FALSE
    )
    result
}

.m13a_panel_unavailable_fixture <- function() {
    structure(
        list(
            status = "unavailable",
            quantity = "association",
            code = "association_no_observable_features",
            message = "No feature is observed in the supplied input.",
            requires = "one globally observable feature"
        ),
        class = "imputefinder_unavailable"
    )
}

m13a_unavailable_codes <- function() {
    data.frame(
        code = c(
            "association_no_observable_features",
            "association_nonestimable",
            "association_low_independent_support",
            "association_low_block_support",
            "association_low_interaction_support",
            "association_incompatible_unit_positions",
            "association_low_residual_df",
            "association_low_reference_df",
            "association_low_permutation_resolution",
            "association_incompatible_permutation",
            "association_degenerate_response",
            "association_singular_covariance",
            "association_quasibinomial_boundary",
            "association_quasibinomial_nonconvergence",
            "association_quasibinomial_paired_scope",
            "association_numerical_failure",
            "association_no_testable_hypotheses"
        ),
        stage = c(
            "response", "estimability", "support", "support", "support",
            "support", "fit", "reference", "permutation", "permutation", "fit",
            "covariance", "quasibinomial", "quasibinomial", "scope", "numerics",
            "family"
        ),
        denominator = c(
            "disposition_only", "disposition_only", "disposition_only",
            "disposition_only", "disposition_only", "disposition_only",
            "disposition_only", "opportunity_miss", "opportunity_miss",
            "opportunity_miss", "opportunity_miss", "opportunity_miss",
            "opportunity_miss", "opportunity_miss", "opportunity_miss",
            "opportunity_miss", "diagnostic_only"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_truth <- function() {
    data.frame(
        truth_id = c(
            "null_families", "dda_batch_crossed_target",
            "dda_monotone_unequal_target", "dia_batch_partial_target",
            "dia_monotone_paired_target"
        ),
        scenarios = c(
            .m13a_tokens("dda_null_balanced", "dia_null_unequal"),
            "dda_batch_crossed", "dda_monotone_unequal",
            "dia_batch_partial", "dia_monotone_paired"
        ),
        target_label = c(
            "whole_holm_family", "batch[batch_2]", "condition[D]",
            "batch[batch_3]", "condition[B]"
        ),
        reference = c(
            "all eligible coefficients zero", "batch_1", "A", "batch_1", "A"
        ),
        response = c(
            "all eligible coefficients have zero generator plug-in detection-fraction association",
            rep("per sample mean(1-missing_probability) over the realized globally-observable feature set; not a conditional expectation", 4L)
        ),
        unavailable = c(
            "family false flag if any adjusted p<=0.05; zero available p-values rejects screening/selected-winner gate",
            rep("in licensed screening scope counts as calibration miss; selected-winner unavailability fails full gate; screening-out-of-scope excluded only from performance", 4L)
        ),
        sample_size_candidate = c(64L, 32L, 32L, 32L, 32L),
        stringsAsFactors = FALSE
    )
}

m13a_truth_methods <- function() {
    data.frame(
        candidate_id = c(
            "a_fraction_ols_hc3_cr2",
            "a_fraction_freedman_lane",
            "a_fraction_quasibinomial"
        ),
        oracle_projection = c(
            "OLS rank-coordinate projection of plug-in sample fractions",
            "same OLS rank-coordinate projection as its reported coefficient",
            "quasibinomial rank-coordinate projection of fractional plug-in counts"
        ),
        gated_effect = c(
            "canonical linear-projection coefficient",
            "canonical linear-projection coefficient",
            "standardized encoded-column 0-to-1 probability difference from oracle projection"
        ),
        interpretation = c(
            rep("candidate-specific projection calibration on one common plug-in response; effect magnitudes are not assumed identical across model classes", 3L)
        ),
        stringsAsFactors = FALSE
    )
}

m13a_opportunities <- function() {
    data.frame(
        manifest_id = c(
            "a_role_eligible_manifest_v3", "a_candidate_opportunities_v3"
        ),
        construction = c(
            "candidate-agnostic after acquisition rebuild and role eligibility, before estimability/support",
            "role-eligible rows with observable response passing algebraic estimability common support and residual-df floor"
        ),
        row = rep("scenario x replicate x stratum x eligible canonical coefficient", 2L),
        denominator = c(
            "disposition audit: every row must retain available or exact unavailable output",
            "scope ranking: every opportunity; no successful-fit denominator"
        ),
        candidate_failure = c(
            "nonestimable unsupported and out-of-scope remain visible",
            "in-scope fit failure and out-of-scope both remain denominator misses"
        ),
        numerator = c(
            "available or exact structured unavailable output",
            "complete finite available outcome"
        ),
        null_denominator = c(
            "not_applicable",
            "32 scenario-replicates per acquisition; zero available p-values makes the selected-winner gate fail"
        ),
        target_denominator = c(
            "not_applicable",
            "four fixed targets x 32 replicates for the selected-winner gates"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_evaluation_allocation <- function() {
    manifest <- .m13a_m12$.m12_generator_manifest()
    effective_tracks <- manifest$required_tracks
    effective_tracks[manifest$scenario_id == "dda_monotone_unequal"] <- "A;B;C"
    data.frame(
        scenario_id = manifest$scenario_id,
        acquisition = manifest$acquisition,
        candidate_replicates = "1-32",
        candidate_count = 32L,
        development_replicates = "33-64",
        development_count = 32L,
        candidate_state = "authorized_after_v4_and_implementation_hash_freeze",
        development_state = "sealed_until_winner_lock",
        m12_required_tracks = manifest$required_tracks,
        effective_required_tracks = effective_tracks,
        stringsAsFactors = FALSE
    )
}

m13a_selection <- function() {
    data.frame(
        selection_id = "a_association_selection_v4",
        candidates = .m13a_tokens(
            "a_fraction_ols_hc3_cr2", "a_fraction_freedman_lane",
            "a_fraction_quasibinomial"
        ),
        computation_scope = "all 13 scenarios x replicates 1-32 only",
        hard_rejection = "candidate screening: within-licensed-scope null/calibration miss malformed disposition causal wording or any sealed-engine execution failure; out-of-scope is not a screening rejection; every screening-qualified candidate then faces the full v4 synthetic gate",
        scope_denominator = "candidate-agnostic opportunity manifest after common estimability/support",
        scope_numerator = "denominator hypotheses returning complete finite available outcomes",
        out_of_scope = "named unavailable retained in denominator and therefore lowers scope coverage",
        screening_metrics = "same thresholds within licensed target/null subsets; in-scope unavailable stays a miss; exact denominators reported",
        ranking = "rank only full-v4-gate passers by the deterministic anchored-cluster algorithm in m13a_rank_candidates",
        ranking_near_tie = "absolute metric difference <= 0.02 + 8*.Machine$double.eps*pmax(1,abs(value),abs(anchor),0.02); mathematically exact decimal boundary included and materially outside excluded",
        selected_winner_gate = "apply the unchanged-threshold full 64-null/128-target v4 registry to every screening-qualified candidate before ranking; unavailable target or zero-p family eliminates that candidate",
        no_winner = "no screening-qualified full-gate passer kills/parks the panel; one candidate failure never suppresses another passer; design core remains mandatory",
        development_unlock = "lock one winner and candidate evidence hash before opening replicates 33-64 or HarmonizR bytes",
        confirmation_unlock = "M15P synchronized opening only; MultiPro remains sealed through M13",
        stringsAsFactors = FALSE
    )
}

.m13a_ranking_near <- function(value, anchor) {
    slack <- 8 * .Machine$double.eps * pmax(
        1,
        abs(value),
        abs(anchor),
        0.02
    )
    abs(value - anchor) <= 0.02 + slack
}

m13a_rank_candidates <- function(metrics) {
    fields <- c(
        "candidate_id", "scope_coverage", "power", "absolute_bias",
        "median_runtime"
    )
    candidates <- m13a_candidate_engines()
    valid <- is.data.frame(metrics) && identical(names(metrics), fields) &&
        nrow(metrics) > 0L && !anyDuplicated(metrics$candidate_id) &&
        all(metrics$candidate_id %in% candidates$candidate_id) &&
        all(vapply(metrics[-1L], is.double, logical(1L))) &&
        !anyNA(metrics) && all(is.finite(as.matrix(metrics[-1L]))) &&
        all(metrics$scope_coverage >= 0 & metrics$scope_coverage <= 1) &&
        all(metrics$power >= 0 & metrics$power <= 1) &&
        all(metrics$absolute_bias >= 0) && all(metrics$median_runtime >= 0)
    if (!valid) {
        .m13a_fail("M13c candidate-ranking metrics are malformed.")
    }
    metrics$simplicity_rank <- candidates$simplicity_rank[
        match(metrics$candidate_id, candidates$candidate_id)
    ]
    output <- character()
    remaining <- metrics
    while (nrow(remaining)) {
        anchor_order <- order(
            -remaining$scope_coverage,
            remaining$candidate_id,
            method = "radix"
        )
        anchor <- remaining[anchor_order[[1L]], , drop = FALSE]
        near <- .m13a_ranking_near(
            remaining$scope_coverage,
            anchor$scope_coverage
        ) & .m13a_ranking_near(
            remaining$power,
            anchor$power
        ) & .m13a_ranking_near(
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

m13a_screening_allocation <- function() {
    data.frame(
        candidate_id = c(
            "a_fraction_ols_hc3_cr2",
            "a_fraction_freedman_lane",
            "a_fraction_quasibinomial"
        ),
        null_scenarios = rep(
            .m13a_tokens("dda_null_balanced", "dia_null_unequal"),
            3L
        ),
        null_family_count = c(64L, 64L, 64L),
        alternative_scenarios = c(
            rep(.m13a_tokens(
                "dda_batch_crossed", "dda_monotone_unequal",
                "dia_batch_partial", "dia_monotone_paired"
            ), 2L),
            .m13a_tokens(
                "dda_batch_crossed", "dda_monotone_unequal",
                "dia_batch_partial"
            )
        ),
        alternative_target_count = c(128L, 128L, 96L),
        in_scope_failure = rep(
            "retained as calibration/null miss; never removed from denominator",
            3L
        ),
        out_of_scope = c(
            "none", "none",
            "dia_monotone_paired excluded from screening performance only; retained as common scope-opportunity miss"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_screening_thresholds <- function() {
    data.frame(
        metric = c(
            "family_false_flag_clopper_pearson_upper95",
            "maximum_acquisition_false_flag_fraction",
            "candidate_specific_projection_95_interval_coverage",
            "holm_target_detection_fraction",
            "median_absolute_candidate_specific_projection_bias"
        ),
        operator = c("<=", "<=", ">=", ">=", "<="),
        threshold = c(0.15, 0.125, 0.90, 0.70, 0.03),
        scale = c("probability", "fraction", "fraction", "fraction", "probability points"),
        denominator = c(
            "licensed null families; exact candidate allocation; zero-p family rejects screening",
            "licensed null families by acquisition; exact candidate allocation",
            "licensed fixed alternative targets; in-scope unavailable is miss",
            "licensed fixed alternative targets; in-scope unavailable is miss",
            "licensed fixed alternative targets; in-scope unavailable makes metric failed"
        ),
        role = "candidate screening only; every screening-qualified candidate then faces the full v4 gate registry before ranking",
        stringsAsFactors = FALSE
    )
}

m13a_gate_statistics <- function() {
    data.frame(
        m12_gate_id = c(
            "a_assoc_null_upper", "a_assoc_null_stratum",
            "a_assoc_interval_coverage", "a_assoc_alternative_power",
            "a_assoc_effect_bias"
        ),
        estimate = c(
            "one-sided 95% Clopper-Pearson upper bound for the pooled false-family fraction",
            "maximum DDA/DIA false-family fraction; lexical acquisition breaks an exact tie",
            "mean of 128 fixed-target binary indicators that conf_low<=candidate-specific nonzero oracle truth<=conf_high, with both endpoints closed",
            "mean of 128 fixed-target binary indicators adjusted_p<=0.05",
            "stats::median of 128 fixed-target absolute candidate-specific projection errors; arithmetic mean of ordered positions 64 and 65"
        ),
        numerator = c(
            "false-family count", "false-family count in selected acquisition",
            "covered-target count", "detected-target count",
            "sum absolute projection error; audit field only"
        ),
        denominator = c(
            "64 fixed null scenario-replicates",
            "32 fixed null scenario-replicates in selected acquisition; gate-registry sample_size=64 denotes total pooled null allocation",
            rep("128 fixed targets: four scenarios x replicates 1-32", 3L)
        ),
        uncertainty = c(
            "one-sided exact binomial lower=0 upper=qbeta(.95,x+1,n-x); x=n gives 1; observed fraction retained as numerator/denominator",
            "two-sided exact binomial qbeta(.025/.975) for the selected acquisition; boundary endpoints 0/1",
            rep("9999-draw scenario-stratified whole-generator-replicate bootstrap; canonical target-table scenario order and replicate order 1:32; sample.int(32L,32L,replace=TRUE) once per scenario per draw; type-8 .025/.975 quantiles expanded only to contain the point estimate", 3L)
        ),
        seed = c(
            rep("not_applicable", 2L),
            rep("SHA-256 UTF-8 paste(protocol_id,candidate_id,versioned_gate_id,scenario_stratified_whole_replicate_bootstrap_v1,sep=TAB); integer first seven hex digits plus one", 3L)
        ),
        rng = c(
            rep("not_applicable", 2L),
            rep("R RNGkind Mersenne-Twister/Inversion/Rejection; caller kind seed and seed absence restored exactly", 3L)
        ),
        dependence = c(
            rep("not_applicable", 2L),
            rep("the resampling unit is one complete generated matrix/result; all within-replicate feature-module dependence stays intact and features/modules are never resampled separately", 3L)
        ),
        passage = c(
            "threshold compares the one-sided Clopper-Pearson upper endpoint stored as estimate",
            rep(
                "threshold compares the point statistic; interval is audit evidence",
                4L
            )
        ),
        bootstrap_statistic = c(
            rep("not_applicable", 2L),
            "mean", "mean",
            "stats::median; even-length arithmetic mean of the two central ordered values in every draw"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_effective_manifest <- function() {
    manifest <- .m13a_m12$.m12_generator_manifest()
    manifest$required_tracks[manifest$scenario_id == "dda_monotone_unequal"] <-
        "A;B;C"
    manifest
}

m13a_normative_digest <- function() {
    path <- "dev/m13-association-contract.md"
    if (!file.exists(path)) {
        .m13a_fail("M13c normative detail file is absent.")
    }
    data.frame(
        path = path,
        sha256 = unname(tools::sha256sum(path)),
        role = "normative_protocol_detail",
        stringsAsFactors = FALSE
    )
}

.m13a_protocol_frames <- function() {
    list(
        normative_digest = m13a_normative_digest(),
        normative_source_digest = m13a_normative_source_digest(),
        implementation_bindings = m13a_implementation_bindings(),
        response = m13a_response(),
        strata = m13a_strata(),
        hypotheses = m13a_hypotheses(),
        support = m13a_support(),
        unit_positions = m13a_unit_positions(),
        candidate_engines = m13a_candidate_engines(),
        numerics = m13a_numerics(),
        permutation = m13a_permutation(),
        permutation_counts = m13a_permutation_counts(),
        permutation_maps = m13a_permutation_maps(),
        seed_protocol = m13a_seed_protocol(),
        input_hash = m13a_input_hash(),
        multiplicity = m13a_multiplicity(),
        result_schema = m13a_result_schema(),
        result_records = m13a_result_records(),
        lifecycle = m13a_lifecycle(),
        unavailable_codes = m13a_unavailable_codes(),
        truth = m13a_truth(),
        truth_methods = m13a_truth_methods(),
        opportunities = m13a_opportunities(),
        selection = m13a_selection(),
        screening_allocation = m13a_screening_allocation(),
        screening_thresholds = m13a_screening_thresholds(),
        gate_statistics = m13a_gate_statistics(),
        structural_gates = m13a_structural_gates(),
        execution_seal = m13a_execution_seal(),
        evidence_hash_protocol = m13a_evidence_hash_protocol(),
        implementation_manifest_schema = m13a_implementation_manifest_schema(),
        candidate_evidence_schema = m13a_candidate_evidence_schema()
    )
}

m13a_protocol_hash <- function() {
    hashes <- .m13a_hash_frames(.m13a_protocol_frames())
    aggregate <- paste(names(hashes), hashes, sep = "\t", collapse = "\n")
    unname(tools::sha256sum(bytes = charToRaw(enc2utf8(aggregate))))
}

.m13a_effective_evidence_hash <- function(data_ids) {
    contract <- .m13a_m12$m12_validation_contract()
    .m13a_m12$.m12_evidence_hash(
        data_ids,
        generator = m13a_effective_manifest(),
        generator_protocol = contract$generator_protocol,
        internal = contract$internal_evidence,
        cards = contract$data_cards,
        roles = contract$data_roles,
        artifacts = contract$artifact_manifest
    )
}

.m13a_association_stage <- function(gate_id) {
    if (grepl("^a_public_confirmation_", gate_id)) {
        "m15p_confirmation"
    } else if (grepl("^a_public_development_", gate_id)) {
        "post_selection_development"
    } else {
        "selected_winner_synthetic"
    }
}

.m13a_association_replicates <- function(stage) {
    ifelse(
        stage == "selected_winner_synthetic",
        "1-32",
        ifelse(stage == "post_selection_development", "external_development", "not_applicable")
    )
}

.m13a_v4_comparator <- function(gate_id, comparator) {
    replacements <- c(
        a_public_confirmation_coverage = paste0(
            "all estimable per-acquisition condition/instrument/nuisance ",
            "coefficients return complete evidence; acquisition descriptive"
        ),
        a_assoc_interval_coverage = paste0(
            "candidate-specific projection of generator plug-in sample-",
            "detection fractions"
        ),
        a_assoc_alternative_power = paste0(
            "nonzero fixed candidate-specific generator plug-in projection target"
        ),
        a_assoc_effect_bias = paste0(
            "candidate-specific projection of generator plug-in sample-",
            "detection fractions"
        )
    )
    if (gate_id %in% names(replacements)) {
        unname(replacements[[gate_id]])
    } else comparator
}

.m13a_v4_metric <- function(gate_id, metric) {
    if (identical(gate_id, "a_assoc_interval_coverage")) {
        "fixed_target_95_interval_coverage_fraction"
    } else {
        metric
    }
}

.m13a_v4_uncertainty <- function(gate_id, uncertainty_method) {
    if (gate_id %in% c(
        "a_assoc_interval_coverage",
        "a_assoc_alternative_power",
        "a_assoc_effect_bias"
    )) {
        paste0(
            "9999-draw scenario-stratified whole-generator-replicate ",
            "bootstrap; type-8 .025/.975; point-contained"
        )
    } else {
        uncertainty_method
    }
}

m13a_gate_registry <- function() {
    base <- .m13a_m12$m12_validation_contract()$gate_registry
    selected <- base$track == "A" &
        (grepl("^a_assoc_", base$gate_id) | grepl("^a_public_", base$gate_id))
    base <- base[selected, , drop = FALSE]
    base_bindings <- .m13a_m12$m12_gate_binding_hashes()
    stage <- vapply(base$gate_id, .m13a_association_stage, character(1L))
    output <- data.frame(
        gate_id = paste0(base$gate_id, "_v4"),
        registry_version = "m13a_gate_registry_v4",
        claim_id = base$claim_id,
        track = base$track,
        metric = unname(mapply(
            .m13a_v4_metric,
            base$gate_id,
            base$metric,
            USE.NAMES = FALSE
        )),
        unit = base$unit,
        comparator = vapply(
            seq_len(nrow(base)),
            function(index) .m13a_v4_comparator(
                base$gate_id[[index]], base$comparator[[index]]
            ),
            character(1L)
        ),
        strata = base$strata,
        sample_size = base$sample_size,
        uncertainty_method = unname(mapply(
            .m13a_v4_uncertainty,
            base$gate_id,
            base$uncertainty_method,
            USE.NAMES = FALSE
        )),
        threshold_operator = base$threshold_operator,
        threshold_value = base$threshold_value,
        threshold_scale = base$threshold_scale,
        failure_treatment = base$failure_treatment,
        data_ids = base$data_ids,
        protocol_id = .M13A_PROTOCOL_ID,
        data_hash = vapply(
            base$data_ids,
            .m13a_effective_evidence_hash,
            character(1L)
        ),
        protocol_hash = m13a_protocol_hash(),
        stage = stage,
        replicate_range = .m13a_association_replicates(stage),
        m12_gate_id = base$gate_id,
        m12_gate_binding_hash = unname(base_bindings[base$gate_id]),
        result_state = ifelse(stage == "m15p_confirmation", "sealed", "unrun"),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    output$gate_binding_hash <- vapply(seq_len(nrow(output)), function(index) {
        .m13a_m12$.m12_hash_frame(output[index, , drop = FALSE])
    }, character(1L))
    output
}

m13a_completed_gate_lineage <- function() {
    base <- .m13a_m12$m12_validation_contract()$gate_registry
    selected <- base$track == "A" &
        !grepl("^a_assoc_", base$gate_id) & !grepl("^a_public_", base$gate_id)
    base <- base[selected, , drop = FALSE]
    bindings <- .m13a_m12$m12_gate_binding_hashes()
    data.frame(
        gate_id = base$gate_id,
        registry_version = base$registry_version,
        gate_binding_hash = unname(bindings[base$gate_id]),
        evidence_stage = "m13a_completed_design_core",
        result_state = "completed",
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

m13a_gate_binding_hashes <- function(registry = m13a_gate_registry()) {
    canonical <- m13a_gate_registry()
    if (!identical(registry, canonical)) {
        .m13a_fail("M13c gate bindings require the canonical v4 registry.")
    }
    stats::setNames(registry$gate_binding_hash, registry$gate_id)
}

m13a_gate_registry_hash <- function(registry = m13a_gate_registry()) {
    canonical <- m13a_gate_registry()
    if (!identical(registry, canonical)) {
        .m13a_fail("M13c gate-registry hash requires the canonical v4 registry.")
    }
    .m13a_m12$.m12_hash_frame(registry)
}

m13a_evaluate_gate_results <- function(
    results,
    stage
) {
    registry <- m13a_gate_registry()
    stages <- unique(registry$stage)
    if (!.m13a_scalar_character(stage) || !stage %in% stages) {
        .m13a_fail("M13c gate evaluation requires one exact evidence stage.")
    }
    registry <- registry[registry$stage == stage, , drop = FALSE]
    gate_ids <- registry$gate_id
    expected_fields <- c(
        "gate_id", "registry_version", "metric", "estimate", "numerator",
        "denominator", "lower", "upper", "status", "gate_binding_hash",
        "evidence_hash", "note"
    )
    if (!is.data.frame(results) ||
        !identical(names(results), expected_fields) ||
        nrow(results) != length(gate_ids) || anyDuplicated(results$gate_id) ||
        !setequal(results$gate_id, gate_ids)) {
        .m13a_fail("M13c gate results have malformed rows or keys.")
    }
    results <- results[match(gate_ids, results$gate_id), , drop = FALSE]
    row.names(results) <- NULL
    row.names(registry) <- NULL
    measured <- results$status == "measured"
    numeric_fields <- c("estimate", "numerator", "denominator", "lower", "upper")
    valid <- all(vapply(results[numeric_fields], is.double, logical(1L))) &&
        all(results$status %in% c("measured", "failed", "unavailable")) &&
        identical(results$registry_version, registry$registry_version) &&
        identical(results$metric, registry$metric) &&
        identical(results$gate_binding_hash, registry$gate_binding_hash) &&
        all(grepl("^[0-9a-f]{64}$", results$evidence_hash)) &&
        is.character(results$note) && !anyNA(results$note) &&
        all(nzchar(results$note)) &&
        all(is.finite(results$estimate[measured])) &&
        all(is.finite(results$numerator[measured])) &&
        all(is.finite(results$denominator[measured])) &&
        all(results$denominator[measured] > 0) &&
        all(is.finite(results$lower[measured])) &&
        all(is.finite(results$upper[measured])) &&
        all(results$lower[measured] <= results$estimate[measured]) &&
        all(results$estimate[measured] <= results$upper[measured]) &&
        all(is.na(results$estimate[!measured])) &&
        all(is.na(results$numerator[!measured])) &&
        all(is.na(results$denominator[!measured])) &&
        all(is.na(results$lower[!measured])) &&
        all(is.na(results$upper[!measured]))
    if (!valid) {
        .m13a_fail("M13c gate result values or bindings are malformed.")
    }
    compare <- function(value, operator, threshold) {
        switch(
            operator,
            "<=" = value <= threshold, ">=" = value >= threshold,
            "==" = value == threshold, "<" = value < threshold,
            ">" = value > threshold,
            FALSE
        )
    }
    passed <- measured & mapply(
        compare,
        results$estimate,
        registry$threshold_operator,
        registry$threshold_value,
        USE.NAMES = FALSE
    )
    data.frame(
        gate_id = registry$gate_id,
        m12_gate_id = registry$m12_gate_id,
        stage = registry$stage,
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

m13a_structural_gates <- function() {
    data.frame(
        rail_id = c(
            "a_candidate_disposition_exact", "a_candidate_holm_exact",
            "a_candidate_engine_execution_complete",
            "a_screening_null_available",
            "a_full_gate_candidate_null_available",
            "a_full_gate_candidate_target_finite"
        ),
        rule = c(
            "finite_available_or_exact_structured_unavailable_fraction",
            "canonical_family_holm_recalculation_match_fraction",
            "every sealed candidate-engine run completes without a top-level execution failure; a caught error emits association_numerical_failure for every prepared hypothesis but remains a global structural rejection",
            "each licensed-scope null family has at least one available raw p-value",
            "each screening-qualified candidate full-denominator null family has at least one available raw p-value",
            "each screening-qualified candidate has complete finite outcomes for every full-denominator alternative target"
        ),
        denominator = c(
            "all eligible coefficient outputs including nonestimable unsupported and out-of-scope",
            "all nonempty available families",
            "all 416 scenario-replicate executions per candidate",
            "64 null scenario-replicates per candidate",
            "64 null scenario-replicates per screening-qualified candidate",
            "128 fixed alternative targets per screening-qualified candidate"
        ),
        required = rep(TRUE, 6L),
        failure = c(
            "candidate screening rejection", "candidate screening rejection",
            "candidate screening rejection",
            "candidate screening rejection",
            "candidate eliminated before ranking", "candidate eliminated before ranking"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_execution_seal <- function() {
    data.frame(
        seal_id = "m13a_candidate_implementation_v1",
        timing = "after all implementation tests pass and before any replicate 1-32 candidate computation",
        source_manifest = "exact source_file schema; safe canonical repository paths plus current-byte SHA-256; roles contract,engine,evaluator,schema,test,harness all required",
        environment_manifest = "exact sorted name/value/required schema including R.version.string R.platform imputefinder.version and every loaded non-base dependency",
        test_manifest = "exact sorted command/positive-expectation-count/passed-status/evidence-SHA-256 schema",
        binding = "candidate evidence stores v4 contract hash protocol hash implementation-manifest hash effective-manifest hash and exact scenario/replicate IDs",
        replay = "production resolver stores no caller callbacks; package-owned readers regenerate inputs from the manifest-bound M12 generator and re-execute each engine from candidate-visible preparation only",
        authorization = "fixture replay has a distinct schema/class; production resolver requires the exact on-disk pre-allocation manifest; authorization requires the exact canonical manifest+416-input+1248-result inventory with no extras and matching byte hashes before/after complete replay; candidate evidence persists that allocation-inventory aggregate and verification additionally binds the stored evidence object",
        mutation_rule = "any source environment allocation input result or replay-binding drift rejects execution and requires a new pre-result implementation seal",
        current_state = "required_uncreated",
        stringsAsFactors = FALSE
    )
}

m13a_evidence_hash_protocol <- function() {
    data.frame(
        position = 1:5,
        component = c(
            "scalar_fields", "ordered_vectors", "nested_frames",
            "component_hashes", "self_hash"
        ),
        rule = c(
            "one position/field/value/type data frame in exact enclosing-list order; NA uses the bound frame canonicalizer token",
            "one position/value/type data frame per vector; exact order retained",
            "exact schema field order/types and schema-key row order required before hashing",
            "M12 percent-escaped UTF-8 tab/newline frame SHA-256; aggregate component-name TAB hash lines in schema order",
            "manifest_hash/evidence_hash omitted from its own components then set to aggregate SHA-256"
        ),
        implementation = "m13a_evidence_bundle_hash_v1; validator recomputes every component and rejects mismatch",
        stringsAsFactors = FALSE
    )
}

m13a_evidence_bundle_hash <- function(components) {
    valid <- is.list(components) && length(components) > 0L &&
        !is.null(names(components)) && !anyNA(names(components)) &&
        all(nzchar(names(components))) && !anyDuplicated(names(components)) &&
        all(vapply(components, is.data.frame, logical(1L)))
    if (!valid) {
        .m13a_fail("M13c evidence-hash components are malformed.")
    }
    hashes <- vapply(
        components,
        .m13a_m12$.m12_hash_frame,
        character(1L)
    )
    aggregate <- paste(names(hashes), hashes, sep = "\t", collapse = "\n")
    unname(tools::sha256sum(bytes = charToRaw(enc2utf8(aggregate))))
}

m13a_implementation_manifest_schema <- function() {
    rows <- list(
        implementation_manifest = c(
            "schema=character[1]:m13a_implementation_manifest_v1",
            "contract_hash=sha256:current-v4", "protocol_hash=sha256:current-v4",
            "gate_registry_hash=sha256:current-v4",
            "effective_manifest_hash=sha256:current-v4",
            "source_files=data.frame:source_file",
            "environment=data.frame:environment_value",
            "tests=data.frame:test_evidence", "state=character[1]:locked",
            "manifest_hash=sha256:self-excluding:evidence-hash-protocol"
        ),
        source_file = c(
            "path=character:unique-safe-relative-radix-sorted",
            "role=character:contract|engine|evaluator|schema|test|harness",
            "sha256=sha256:current-file-bytes"
        ),
        environment_value = c(
            "name=character:unique-radix-sorted", "value=character:nonempty",
            "required=logical;TRUE-for-R.version.string,R.platform,imputefinder.version-and-every-loaded-nonbase-dependency"
        ),
        test_evidence = c(
            "command=character:unique-radix-sorted",
            "expectation_count=integer:positive", "status=character:passed",
            "evidence_sha256=sha256:exact-command-output"
        ),
        implementation_validator = c(
            "timing=before-any-candidate-allocation",
            "checks=exact-fields-types-cardinality-order-live-file-hashes-required-environment-tests-and-self-hash",
            "implementation=required-red-to-green-rail-in-locked-implementation-source"
        )
    )
    do.call(rbind, lapply(names(rows), function(record) {
        data.frame(
            record = record,
            position = seq_along(rows[[record]]),
            rule = rows[[record]],
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    }))
}

m13a_candidate_evidence_schema <- function() {
    rows <- list(
        candidate_evidence = c(
            "schema=character[1]:m13a_candidate_evidence_v1",
            "contract_hash=sha256", "protocol_hash=sha256",
            "gate_registry_hash=sha256", "implementation_hash=sha256",
            "effective_manifest_hash=sha256",
            "generator_protocol_hash=sha256",
            "artifact_inventory_hash=sha256:canonical-manifest-plus-416-input-plus-1248-result-byte-inventory; production-only; fixture-NA",
            "candidate_ids=character[3]",
            "scenario_ids=character[13]", "replicate_ids=integer[32]:1..32",
            "run_bindings=data.frame:run_binding",
            "outcome_bindings=data.frame:outcome_binding",
            "screening_metrics=data.frame:screening_metric",
            "full_gate_results=data.frame:candidate-plus-exact-stage-gate-result",
            "ranking_metrics=data.frame:ranking_metric-full-gate-passers-only",
            "candidate_order=character:exact-m13a_rank_candidates-output",
            "selected_candidate=character[1]-or-NA:first-candidate-order",
            "state=character[1]:winner_locked|no_winner",
            "evidence_hash=sha256:self-excluding:evidence-hash-protocol"
        ),
        run_binding = c(
            "candidate_id=character:engine-order",
            "scenario_id=character:effective-manifest-order",
            "replicate=integer:1..32", "input_sha256=sha256",
            "result_sha256=sha256-internal-candidate-artifact-without-public-winner-wrapper-measured-else-NA",
            "result_kind=character:panel|unavailable|execution_failure",
            "elapsed_seconds=double:finite-nonnegative",
            "status=character:measured|failed",
            "failure_code=character:NA-measured-else-nonempty",
            "key=complete-Cartesian-3x13x32-in-canonical-order"
        ),
        outcome_binding = c(
            "candidate_id=character", "scenario_id=character",
            "replicate=integer", "stratum=character-or-NA-panel-abstention",
            "hypothesis=character:a_sha-or-literal-association-panel-abstention",
            "status=character:available|unavailable",
            "code=character:NA-available-else-frozen-unavailable-code",
            "key=candidate+scenario+replicate+stratum+hypothesis-canonical"
        ),
        screening_metric = c(
            "candidate_id=character:three-candidate-order",
            "metric=character:common_opportunity_coverage-plus-five-screening-thresholds",
            "numerator=double-or-NA-failed",
            "denominator=double-positive-or-NA-failed",
            "estimate=double-finite-or-NA-failed",
            "licensed_scope=character:exact-screening-allocation",
            "status=character:measured|failed", "passed=logical:exact-comparator"
        ),
        full_gate_result = c(
            "candidate_id=character:screening-qualified-candidate-order",
            "gate_id=character:exact-five-v4-selected-winner-gates",
            "m12_gate_id=character:bound-upstream-gate",
            "stage=character:selected_winner_synthetic",
            "metric=character:exact-gate-registry-metric",
            "estimate=double-finite-measured-else-NA",
            "numerator=double-nonnegative-measured-else-NA",
            "denominator=double-positive-measured-else-NA",
            "lower=double<=estimate-measured-else-NA",
            "upper=double>=estimate-measured-else-NA",
            "operator=character:exact-gate-registry-operator",
            "threshold=double:exact-gate-registry-threshold",
            "status=character:measured|failed|unavailable",
            "passed=logical:measured-and-exact-comparator",
            "failure_treatment=character:exact-gate-registry-treatment",
            "evidence_hash=sha256:exact-gate-audit-frame",
            "note=character:nonempty-statistic-or-failure-description"
        ),
        ranking_metric = c(
            "candidate_id=character:full-v4-gate-passer",
            "scope_coverage=double:[0,1]", "power=double:[0,1]",
            "absolute_bias=double:nonnegative",
            "median_runtime=double:finite-nonnegative"
        ),
        candidate_evidence_validator = c(
            "checks=exact-fields-types-cardinality-canonical-orders-self-hash-and-3x13x32-run-grid",
            "screening=derive-five-threshold-passers; require-full-selected-winner-stage-gate-set-for-each",
            "selection=rank-only-full-gate-passers-with-m13a_rank_candidates; selected-is-first-or-NA",
            "implementation=required-red-to-green-rail-before-result-allocation"
        )
    )
    do.call(rbind, lapply(names(rows), function(record) {
        data.frame(
            record = record,
            position = seq_along(rows[[record]]),
            rule = rows[[record]],
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    }))
}

m13a_release_boundary <- function() {
    data.frame(
        state = c(
            "contract_frozen", "candidate_study", "winner_locked",
            "development_pass", "confirmation_pass"
        ),
        production_behavior = c(
            "implement/test engines and lock implementation manifest; candidate execution still forbidden",
            "candidate engines remain internal dev code; sentinel association absent",
            "only selected candidate may populate sentinel association",
            "association panel may become confirmation_candidate; no release claim",
            "association panel may ship if synchronized M15P gates pass"
        ),
        evidence_open = c(
            "no association result allocation",
            "synthetic replicates 1-32 only",
            "synthetic replicates 33-64 plus HarmonizR development",
            "no confirmation evidence",
            "predeclared synchronized confirmation artifacts"
        ),
        stringsAsFactors = FALSE
    )
}

m13a_contract <- function() {
    list(
        descriptor = m13a_descriptor(),
        upstream_bindings = m13a_upstream_bindings(),
        normative_digest = m13a_normative_digest(),
        normative_source_digest = m13a_normative_source_digest(),
        implementation_bindings = m13a_implementation_bindings(),
        effective_manifest = m13a_effective_manifest(),
        response = m13a_response(),
        strata = m13a_strata(),
        hypotheses = m13a_hypotheses(),
        support = m13a_support(),
        unit_positions = m13a_unit_positions(),
        candidate_engines = m13a_candidate_engines(),
        numerics = m13a_numerics(),
        permutation = m13a_permutation(),
        permutation_counts = m13a_permutation_counts(),
        permutation_maps = m13a_permutation_maps(),
        seed_protocol = m13a_seed_protocol(),
        input_hash = m13a_input_hash(),
        multiplicity = m13a_multiplicity(),
        result_schema = m13a_result_schema(),
        result_records = m13a_result_records(),
        lifecycle = m13a_lifecycle(),
        unavailable_codes = m13a_unavailable_codes(),
        truth = m13a_truth(),
        truth_methods = m13a_truth_methods(),
        opportunities = m13a_opportunities(),
        evaluation_allocation = m13a_evaluation_allocation(),
        selection = m13a_selection(),
        screening_allocation = m13a_screening_allocation(),
        screening_thresholds = m13a_screening_thresholds(),
        gate_statistics = m13a_gate_statistics(),
        gate_registry = m13a_gate_registry(),
        completed_gate_lineage = m13a_completed_gate_lineage(),
        structural_gates = m13a_structural_gates(),
        execution_seal = m13a_execution_seal(),
        evidence_hash_protocol = m13a_evidence_hash_protocol(),
        implementation_manifest_schema = m13a_implementation_manifest_schema(),
        candidate_evidence_schema = m13a_candidate_evidence_schema(),
        release_boundary = m13a_release_boundary()
    )
}

m13a_seed_manifest <- function(
    input_sha256,
    instances,
    draw_ids
) {
    valid <- is.character(input_sha256) && length(input_sha256) == 1L &&
        grepl("^[0-9a-f]{64}$", input_sha256) &&
        is.data.frame(instances) &&
        identical(
            names(instances),
            c("candidate_id", "acquisition", "hypothesis_id")
        ) &&
        nrow(instances) > 0L &&
        all(vapply(instances, is.character, logical(1L))) &&
        !anyNA(instances) &&
        all(instances$candidate_id == "a_fraction_freedman_lane") &&
        all(nzchar(instances$acquisition)) &&
        all(grepl("^a_[0-9a-f]{64}$", instances$hypothesis_id)) &&
        !anyDuplicated(instances) &&
        is.integer(draw_ids) && length(draw_ids) > 0L && !anyNA(draw_ids) &&
        all(draw_ids > 0L) && !anyDuplicated(draw_ids)
    if (!valid) {
        .m13a_fail("M13c seed-manifest input is malformed.")
    }
    order_index <- do.call(
        order,
        c(unname(instances), list(method = "radix"))
    )
    instances <- instances[order_index, , drop = FALSE]
    row.names(instances) <- NULL
    grid <- instances[rep(seq_len(nrow(instances)), each = length(draw_ids)), , drop = FALSE]
    grid$draw_id <- rep(sort(draw_ids, method = "radix"), nrow(instances))
    row.names(grid) <- NULL
    used <- new.env(hash = TRUE, parent = emptyenv())
    seed <- seed_nonce <- integer(nrow(grid))
    for (index in seq_len(nrow(grid))) {
        candidate_nonce <- 0L
        repeat {
            key <- c(
                .M13A_PROTOCOL_ID, input_sha256,
                grid$candidate_id[[index]], grid$acquisition[[index]],
                grid$hypothesis_id[[index]], grid$draw_id[[index]],
                candidate_nonce
            )
            digest <- unname(tools::sha256sum(
                bytes = c(
                    charToRaw(
                        "imputefinder:association-permutation-seed:v2"
                    ),
                    as.raw(0L),
                    .m13a_encode_text_vector(as.character(key))
                )
            ))
            candidate_seed <- as.integer(
                strtoi(substr(digest, 1L, 7L), 16L) + 1L
            )
            seed_key <- as.character(candidate_seed)
            if (!exists(seed_key, envir = used, inherits = FALSE)) {
                assign(seed_key, TRUE, envir = used)
                break
            }
            candidate_nonce <- candidate_nonce + 1L
        }
        seed[[index]] <- candidate_seed
        seed_nonce[[index]] <- candidate_nonce
    }
    data.frame(
        protocol_id = .M13A_PROTOCOL_ID,
        candidate_id = grid$candidate_id,
        acquisition = grid$acquisition,
        hypothesis_id = grid$hypothesis_id,
        draw_id = as.integer(grid$draw_id),
        seed = seed,
        seed_nonce = seed_nonce,
        stringsAsFactors = FALSE
    )
}

# END M13A NORMATIVE SOURCE

.m13a_quasibinomial_result_fixture <- function() {
    result <- .m13a_result_fixture(TRUE)
    result$candidate <- "a_fraction_quasibinomial"
    result$methods[["candidate"]] <- "a_fraction_quasibinomial"
    outcome <- result$outcomes[[1L]]
    outcome$candidate <- "a_fraction_quasibinomial"
    outcome$log_odds <- 1.0
    outcome$log_odds_standard_error <- 0.2
    outcome$log_odds_conf_low <- 0.51
    outcome$log_odds_conf_high <- 1.49
    outcome$statistic <- 5.0
    outcome$raw_p <- 2 * stats::pt(-5, 6)
    outcome$adjusted_p <- outcome$raw_p
    outcome$diagnostics$variance_method <- "quasibinomial"
    outcome$diagnostics$leverage_max <- NA_real_
    outcome$diagnostics$dispersion <- 1.25
    outcome$diagnostics$converged <- TRUE
    outcome$diagnostics$boundary <- FALSE
    result$outcomes[[1L]] <- outcome
    result
}

m13a_normative_source_digest <- function() {
    path <- "dev/m13-association-contract.R"
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    begin <- which(lines == "# BEGIN M13A NORMATIVE SOURCE")
    end <- which(lines == "# END M13A NORMATIVE SOURCE")
    if (length(begin) != 1L || length(end) != 1L || end <= begin + 1L) {
        .m13a_fail("M13c normative source markers are malformed.")
    }
    bytes <- charToRaw(enc2utf8(paste(
        lines[seq.int(begin + 1L, end - 1L)],
        collapse = "\n"
    )))
    data.frame(
        path = path,
        sha256 = unname(tools::sha256sum(bytes = bytes)),
        role = "normative_contract_source",
        stringsAsFactors = FALSE
    )
}

m13a_validate_contract <- function(contract = m13a_contract()) {
    live <- .m13a_live_upstream()
    if (!identical(live, .M13A_EXPECTED_UPSTREAM)) {
        changed <- names(live)[live != .M13A_EXPECTED_UPSTREAM]
        .m13a_fail("M13c upstream binding drift: ", paste(changed, collapse = ", "))
    }
    expected <- m13a_contract()
    if (!is.list(contract) || !identical(names(contract), names(expected))) {
        .m13a_fail("M13c contract components/order differ from v4.")
    }
    valid <- vapply(names(expected), function(component) {
        identical(contract[[component]], expected[[component]])
    }, logical(1L))
    if (!all(valid)) {
        .m13a_fail(
            "M13c contract component drift: ",
            paste(names(valid)[!valid], collapse = ", ")
        )
    }
    if (!identical(
        contract$evaluation_allocation$candidate_replicates,
        rep("1-32", 13L)
    ) || !identical(
        contract$evaluation_allocation$development_replicates,
        rep("33-64", 13L)
    ) || !identical(
        contract$evaluation_allocation$effective_required_tracks[
            contract$evaluation_allocation$scenario_id == "dda_monotone_unequal"
        ],
        "A;B;C"
    )) {
        .m13a_fail("M13c allocation/bookkeeping correction is not sealed.")
    }
    registry <- contract$gate_registry
    binding_fields <- setdiff(names(registry), "gate_binding_hash")
    recomputed <- vapply(seq_len(nrow(registry)), function(index) {
        .m13a_m12$.m12_hash_frame(
            registry[index, binding_fields, drop = FALSE]
        )
    }, character(1L))
    if (nrow(registry) != 9L ||
        any(registry$stage == "m15p_confirmation" &
            registry$result_state != "sealed") ||
        any(registry$stage == "selected_winner_synthetic" &
            registry$replicate_range != "1-32") ||
        any(registry$protocol_hash != m13a_protocol_hash()) ||
        !identical(registry$gate_binding_hash, recomputed) ||
        !all(grepl("^[0-9a-f]{64}$", registry$data_hash)) ||
        !all(grepl("^[0-9a-f]{64}$", registry$m12_gate_binding_hash))) {
        .m13a_fail("M13c gate stages or replicate bindings are malformed.")
    }
    corrected <- contract$effective_manifest$scenario_id ==
        "dda_monotone_unequal"
    if (sum(corrected) != 1L ||
        contract$effective_manifest$required_tracks[corrected] != "A;B;C") {
        .m13a_fail("M13c effective generator manifest is not corrected.")
    }
    invisible(contract)
}

.m13a_hash_frames <- function(frames) {
    vapply(frames, .m13a_m12$.m12_hash_frame, character(1L))
}

m13a_contract_hashes <- function(contract = m13a_contract()) {
    m13a_validate_contract(contract)
    hashes <- .m13a_hash_frames(contract)
    aggregate <- paste(names(hashes), hashes, sep = "\t", collapse = "\n")
    c(
        hashes,
        contract = unname(tools::sha256sum(
            bytes = charToRaw(enc2utf8(aggregate))
        ))
    )
}

m13a_verify_hashes <- function(contract = m13a_contract()) {
    hashes <- m13a_contract_hashes(contract)
    if (!identical(hashes, .M13A_EXPECTED_HASHES)) {
        changed <- union(
            names(hashes)[hashes != .M13A_EXPECTED_HASHES[names(hashes)]],
            names(.M13A_EXPECTED_HASHES)[
                .M13A_EXPECTED_HASHES != hashes[names(.M13A_EXPECTED_HASHES)]
            ]
        )
        .m13a_fail("M13c frozen hash drift: ", paste(changed, collapse = ", "))
    }
    hashes
}

.m13a_expect_error <- function(code) {
    inherits(try(force(code), silent = TRUE), "try-error")
}

.m13a_gate_self_tests <- function() {
    registry <- m13a_gate_registry()
    results <- data.frame(
        gate_id = registry$gate_id,
        registry_version = registry$registry_version,
        metric = registry$metric,
        estimate = as.double(registry$threshold_value),
        numerator = as.double(registry$threshold_value),
        denominator = rep(1.0, nrow(registry)),
        lower = as.double(registry$threshold_value),
        upper = as.double(registry$threshold_value),
        status = rep("measured", nrow(registry)),
        gate_binding_hash = registry$gate_binding_hash,
        evidence_hash = rep(paste(rep("0", 64L), collapse = ""), nrow(registry)),
        note = rep("boundary self-test", nrow(registry)),
        stringsAsFactors = FALSE
    )
    reports <- lapply(unique(registry$stage), function(stage) {
        selected <- registry$stage == stage
        m13a_evaluate_gate_results(results[selected, , drop = FALSE], stage)
    })
    report <- do.call(rbind, reports)
    detached <- results
    first <- substr(detached$gate_binding_hash[[1L]], 1L, 1L)
    detached$gate_binding_hash[[1L]] <- paste0(
        if (identical(first, "0")) "1" else "0",
        substr(detached$gate_binding_hash[[1L]], 2L, 64L)
    )
    selected_stage <- registry$stage[[1L]]
    selected <- registry$stage == selected_stage
    stage_results <- detached[selected, , drop = FALSE]
    c(
        boundary_pass = all(report$passed),
        binding_rejected = .m13a_expect_error(
            m13a_evaluate_gate_results(
                stage_results,
                selected_stage
            )
        ),
        omission_rejected = .m13a_expect_error(
            m13a_evaluate_gate_results(
                stage_results[-1L, , drop = FALSE],
                selected_stage
            )
        )
    )
}

.m13a_hash_self_tests <- function() {
    mask <- matrix(
        c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
        nrow = 2L,
        byrow = TRUE,
        dimnames = list(c("feature_b", "feature_a"), c("s4", "s2", "s1", "s3"))
    )
    design <- data.frame(
        condition = factor(c("B", "A", "A", "B"), levels = c("B", "A")),
        batch = c("two", "one", "two", "one"),
        run = c(4L, 2L, 1L, 3L),
        acquisition = c("DDA", "DDA", "DDA", "DDA"),
        row.names = colnames(mask),
        stringsAsFactors = FALSE
    )
    roles <- list(
        condition = "condition", nuisance = c("run", "batch"),
        block = character(), acquisition = "acquisition"
    )
    interactions <- list(c("condition", "run"), c("batch", "condition"))
    baseline <- m13a_mask_design_sha256(mask, design, roles, interactions)

    permuted <- mask[c(2L, 1L), c(3L, 1L, 4L, 2L), drop = FALSE]
    recoded <- design[colnames(permuted), , drop = FALSE]
    recoded$condition <- as.character(recoded$condition)
    recoded$run <- as.double(recoded$run)
    reordered_roles <- roles[c("acquisition", "block", "nuisance", "condition")]
    reordered_roles$nuisance <- rev(reordered_roles$nuisance)
    invariant <- m13a_mask_design_sha256(
        permuted, recoded, reordered_roles,
        lapply(rev(interactions), rev)
    )

    changed_mask <- mask
    changed_mask[[1L]] <- !changed_mask[[1L]]
    changed_design <- design
    changed_design$run[[1L]] <- changed_design$run[[1L]] + 1L
    changed_roles <- roles
    changed_roles$nuisance <- c("run", "acquisition")
    changed_roles$acquisition <- "batch"
    c(
        order_and_reencoding = identical(baseline, invariant),
        mask_sensitive = !identical(
            baseline,
            m13a_mask_design_sha256(changed_mask, design, roles, interactions)
        ),
        design_sensitive = !identical(
            baseline,
            m13a_mask_design_sha256(mask, changed_design, roles, interactions)
        ),
        role_sensitive = !identical(
            baseline,
            m13a_mask_design_sha256(mask, design, changed_roles, interactions)
        ),
        interaction_sensitive = !identical(
            baseline,
            m13a_mask_design_sha256(mask, design, roles, interactions[-1L])
        )
    )
}

m13a_self_tests <- function(contract = m13a_contract()) {
    m13a_validate_contract(contract)
    changed <- contract
    changed$multiplicity$level <- 0.10
    seeds_a <- m13a_seed_manifest(
        paste(rep("0", 64L), collapse = ""),
        data.frame(
            candidate_id = rep("a_fraction_freedman_lane", 2L),
            acquisition = c("DIA", "DDA"),
            hypothesis_id = paste0(
                "a_", c(paste(rep("1", 64L), collapse = ""), paste(rep("2", 64L), collapse = ""))
            ),
            stringsAsFactors = FALSE
        ),
        1:3
    )
    seeds_b <- m13a_seed_manifest(
        paste(rep("0", 64L), collapse = ""),
        data.frame(
            candidate_id = rep("a_fraction_freedman_lane", 2L),
            acquisition = c("DDA", "DIA"),
            hypothesis_id = paste0(
                "a_", c(paste(rep("2", 64L), collapse = ""), paste(rep("1", 64L), collapse = ""))
            ),
            stringsAsFactors = FALSE
        ),
        as.integer(c(3, 1, 2))
    )
    result_fixture <- .m13a_result_fixture(TRUE)
    monte_carlo_fixture <- .m13a_monte_carlo_result_fixture()
    quasibinomial_fixture <- .m13a_quasibinomial_result_fixture()
    duplicate_map_fixture <- monte_carlo_fixture
    duplicate_map_fixture$diagnostics$seed_manifest$map_sha256[[2L]] <-
        duplicate_map_fixture$diagnostics$seed_manifest$map_sha256[[1L]]
    unavailable_fixture <- .m13a_result_fixture(FALSE)
    panel_unavailable_fixture <- .m13a_panel_unavailable_fixture()
    declared_fixture <- result_fixture
    declared_stratum <- m13a_stratum_id("DDA")
    declared_fixture$response$stratum <- declared_stratum
    declared_fixture$response$acquisition <- "DDA"
    declared_fixture$hypotheses$stratum <- declared_stratum
    declared_hypothesis <- m13a_hypothesis_id(
        declared_stratum,
        declared_fixture$hypotheses$coefficient,
        declared_fixture$hypotheses$label,
        declared_fixture$hypotheses$term_id
    )
    declared_fixture$hypotheses$hypothesis <- declared_hypothesis
    declared_fixture$support$hypothesis <- declared_hypothesis
    declared_fixture$outcomes[[1L]]$quantity <- declared_hypothesis
    declared_fixture$outcomes[[1L]]$stratum <- declared_stratum
    declared_fixture$outcomes[[1L]]$hypothesis <- declared_hypothesis
    names(declared_fixture$outcomes) <- declared_hypothesis
    declared_fixture$multiplicity$stratum <- declared_stratum
    declared_fixture$diagnostics$strata$stratum <- declared_stratum
    declared_fixture$diagnostics$strata$acquisition <- "DDA"
    corrupted_result <- result_fixture
    corrupted_result$outcomes[[1L]]$quantity <- "wrong_key"
    missing_numeric_support <- result_fixture
    missing_numeric_support$hypotheses$component_encodings[[1L]] <- c(
        condition = "numeric"
    )
    wrong_extrapolation <- missing_numeric_support
    wrong_extrapolation$support$numeric_components[[1L]] <- data.frame(
        component = "condition", median = 0.5, minimum = 0.2,
        maximum = 0.8, extrapolates_0_1 = FALSE,
        stringsAsFactors = FALSE
    )
    map_seeds <- seeds_a[
        seeds_a$acquisition == seeds_a$acquisition[[1L]] &
            seeds_a$hypothesis_id == seeds_a$hypothesis_id[[1L]],
        ,
        drop = FALSE
    ]
    wrong_protocol_seeds <- map_seeds
    wrong_protocol_seeds$protocol_id <- "wrong_protocol"
    negative_nonce_seeds <- map_seeds
    negative_nonce_seeds$seed_nonce[[1L]] <- -1L
    double_draw_seeds <- map_seeds
    double_draw_seeds$draw_id <- as.double(double_draw_seeds$draw_id)
    rng_kind_before <- RNGkind()
    rng_exists_before <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    rng_seed_before <- if (rng_exists_before) {
        get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    maps_a <- m13a_monte_carlo_maps(
        map_seeds,
        "independent",
        9L,
        groups = list(1:9)
    )
    maps_b <- m13a_monte_carlo_maps(
        map_seeds,
        "independent",
        9L,
        groups = list(1:9)
    )
    ranking <- data.frame(
        candidate_id = m13a_candidate_engines()$candidate_id,
        scope_coverage = c(0.90, 0.89, 0.88),
        power = c(0.80, 0.79, 0.78),
        absolute_bias = c(0.01, 0.02, 0.03),
        median_runtime = c(2, 1, 0.5),
        stringsAsFactors = FALSE
    )
    evidence_components <- list(
        scalar_fields = data.frame(
            position = 1:2,
            field = c("schema", "state"),
            value = c("fixture", "locked"),
            type = "character",
            stringsAsFactors = FALSE
        ),
        ordered_vector = data.frame(
            position = 1:3,
            value = c("a", "b", "c"),
            type = "character",
            stringsAsFactors = FALSE
        )
    )
    evidence_hash <- m13a_evidence_bundle_hash(evidence_components)
    changed_evidence <- evidence_components
    changed_evidence$ordered_vector$value[[2L]] <- "x"
    c(
        mutation_rejected = .m13a_expect_error(m13a_validate_contract(changed)),
        hashes_sealed = identical(
            m13a_contract_hashes(contract),
            .M13A_EXPECTED_HASHES
        ),
        seed_order_invariant = identical(seeds_a, seeds_b),
        seed_unique = !anyDuplicated(seeds_a$seed),
        candidate_range_sealed = all(
            contract$gate_registry$replicate_range[
                contract$gate_registry$stage == "selected_winner_synthetic"
            ] == "1-32"
        ),
        confirmation_sealed = all(
            contract$gate_registry$result_state[
                contract$gate_registry$stage == "m15p_confirmation"
            ] == "sealed"
        ),
        historical_mismatch_corrected = identical(
            contract$evaluation_allocation$effective_required_tracks[
                contract$evaluation_allocation$scenario_id == "dda_monotone_unequal"
            ],
            "A;B;C"
        ),
        result_schema_available = !.m13a_expect_error(
            m13a_validate_result_shape(result_fixture)
        ),
        result_schema_monte_carlo = !.m13a_expect_error(
            m13a_validate_result_shape(monte_carlo_fixture)
        ),
        result_schema_quasibinomial = !.m13a_expect_error(
            m13a_validate_result_shape(quasibinomial_fixture)
        ),
        result_seed_map_duplicate_rejected = .m13a_expect_error(
            m13a_validate_result_shape(duplicate_map_fixture)
        ),
        result_schema_unavailable = !.m13a_expect_error(
            m13a_validate_result_shape(unavailable_fixture)
        ),
        result_schema_declared = !.m13a_expect_error(
            m13a_validate_result_shape(declared_fixture)
        ),
        panel_abstention_available = !.m13a_expect_error(
            m13a_validate_result_shape(panel_unavailable_fixture)
        ),
        result_join_rejected = .m13a_expect_error(
            m13a_validate_result_shape(corrupted_result)
        ),
        numeric_support_subset_rejected = .m13a_expect_error(
            m13a_validate_result_shape(missing_numeric_support)
        ),
        numeric_extrapolation_rejected = .m13a_expect_error(
            m13a_validate_result_shape(wrong_extrapolation)
        ),
        monte_carlo_maps_reproducible = identical(maps_a, maps_b),
        monte_carlo_seed_rows_rejected = .m13a_expect_error(
            m13a_monte_carlo_maps(
                wrong_protocol_seeds,
                "independent",
                9L,
                groups = list(1:9)
            )
        ) && .m13a_expect_error(
            m13a_monte_carlo_maps(
                negative_nonce_seeds,
                "independent",
                9L,
                groups = list(1:9)
            )
        ) && .m13a_expect_error(
            m13a_monte_carlo_maps(
                double_draw_seeds,
                "independent",
                9L,
                groups = list(1:9)
            )
        ),
        monte_carlo_maps_distinct = !anyDuplicated(
            maps_a$manifest$map_sha256
        ) && all(vapply(maps_a$maps, function(map) {
            !identical(map, seq_along(map))
        }, logical(1L))),
        monte_carlo_rng_restored = identical(RNGkind(), rng_kind_before) &&
            identical(
                exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
                rng_exists_before
            ) && (
                !rng_exists_before || identical(
                    get(".Random.seed", envir = .GlobalEnv),
                    rng_seed_before
                )
            ),
        ranking_order_invariant = identical(
            m13a_rank_candidates(ranking),
            m13a_rank_candidates(ranking[c(3L, 1L, 2L), , drop = FALSE])
        ),
        evidence_hash_sensitive = !identical(
            evidence_hash,
            m13a_evidence_bundle_hash(changed_evidence)
        ),
        .m13a_gate_self_tests(),
        .m13a_hash_self_tests()
    )
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if (!identical(args, "--verify")) {
        .m13a_fail("usage: Rscript --vanilla dev/m13-association-contract.R --verify")
    }
    contract <- m13a_contract()
    tests <- m13a_self_tests(contract)
    if (!all(tests)) {
        .m13a_fail(
            "M13c self-test failure: ",
            paste(names(tests)[!tests], collapse = ", ")
        )
    }
    hashes <- m13a_verify_hashes(contract)
    cat("contract_version: ", .M13A_VERSION, "\n", sep = "")
    cat("protocol_hash: ", m13a_protocol_hash(), "\n", sep = "")
    cat("state: frozen_unrun; candidate=1-32; development=sealed; confirmation=sealed\n")
    cat(
        "self_tests: ",
        paste(names(tests), tests, sep = "=", collapse = "; "),
        "\n",
        sep = ""
    )
    for (name in names(hashes)) {
        cat(name, ": ", hashes[[name]], "\n", sep = "")
    }
}
