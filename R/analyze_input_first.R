#' Build an Experimental Missingness Analysis Sidecar
#'
#' Construct an input-first companion analysis without changing the stable
#' [classify_missingness()] result contract. The sidecar records the original
#' missingness mask, a deterministic input fingerprint, aligned typed design,
#' algebraic design/estimability evidence, pre-rescue evidence, classic result
#' or structured failure, and provenance.
#'
#' @param x An ordinary numeric matrix or a `SummarizedExperiment` containing
#'   unnormalised log2 protein intensities with features in rows and samples in
#'   columns.
#' @param design A [missingness_design()] with sample names exactly matching
#'   `x`. For `SummarizedExperiment` input, its condition role must also name an
#'   equal `colData` column.
#' @param assay Assay name for `SummarizedExperiment` input. It may be omitted
#'   only when exactly one assay exists. It must be `NULL` for matrix input.
#' @param base_fit Optional complete `imputefinder_result`. Its cutoff-source
#'   policy is recomputed against `x`, `design`, and `seed`; reuse occurs only
#'   when exact, canonical, and tolerance comparisons all pass.
#' @param cutoffs Optional named numeric vector of manual condition cutoffs,
#'   with automatic selection filling omitted conditions. `base_fit` and
#'   `cutoffs` are conflicting specifications.
#' @param seed Integer seed for the classic condition-rescue path.
#' @param modules Unique subset of `"sentinel"` and `"stability"`. The default
#'   constructs pre-rescue sentinel evidence and requests stability. Stability
#'   currently returns a structured `imputefinder_unavailable` record pending
#'   its frozen validation milestone. Use `character()` to skip both modules.
#'
#' @return An experimental `imputefinder_analysis` with seven fields:
#'   `classic`, `spec`, `design`, `input`, `sentinel`, `stability`, and
#'   `provenance`. `design$estimability` stores the canonical declared model,
#'   SVD rank/null-space aliases, and block-aware independent-unit accounting.
#'   The sidecar retains computed evidence and a packed missingness mask, but
#'   does not retain a second copy of the original numeric matrix.
#'
#' @details
#' The sidecar reports observed evidence, compatibility, and unavailable
#' quantities. It does not identify a causal missingness mechanism. Serialized
#' objects carry an experimental versioned schema and unknown schemas fail
#' explicitly instead of being silently upgraded.
#'
#' The mandatory design core runs independently of `modules`. It reports
#' algebraic separability and declared units; it performs no association test,
#' automatic correction, sample exclusion, or causal attribution.
#'
#' @examples
#' x <- matrix(
#'     seq_len(16),
#'     nrow = 4,
#'     dimnames = list(paste0("protein", 1:4), paste0("sample", 1:4))
#' )
#' design <- missingness_design(
#'     data.frame(
#'         condition = c("control", "control", "treated", "treated"),
#'         row.names = colnames(x)
#'     ),
#'     condition = "condition"
#' )
#'
#' analysis <- analyze_missingness(x, design, modules = "sentinel")
#' analysis
#'
#' @seealso [missingness_design()] and [classify_missingness()].
#' @export
analyze_missingness <- function(
    x,
    design,
    assay = NULL,
    base_fit = NULL,
    cutoffs = NULL,
    seed = 1L,
    modules = c("sentinel", "stability")
) {
    .analyze_missingness_input_first(
        x = x,
        design = design,
        assay = assay,
        base_fit = base_fit,
        cutoffs = cutoffs,
        seed = seed,
        modules = modules,
        call = match.call()
    )
}

.abort_design_condition <- function(
    message,
    condition_role,
    mismatched_samples = character()
) {
    .abort_missingness_design(
        message,
        "imputefinder_design_condition_error",
        condition_role = condition_role,
        mismatched_samples = mismatched_samples
    )
}

.require_se_design_condition <- function(x, condition_role) {
    matches <- which(names(SummarizedExperiment::colData(x)) == condition_role)
    if (length(matches) != 1L) {
        .abort_design_condition(
            paste0(
                "The design condition role must name exactly one ",
                "SummarizedExperiment colData column."
            ),
            condition_role
        )
    }

    invisible(condition_role)
}

.aligned_design_groups <- function(design) {
    condition_role <- design$roles$condition
    stats::setNames(
        design$sample_data[[condition_role]],
        rownames(design$sample_data)
    )
}

.validate_se_design_condition <- function(prepared, design) {
    expected <- .aligned_design_groups(design)
    actual <- prepared$groups_by_sample
    mismatched <- names(expected)[expected != actual]
    if (length(mismatched) > 0L) {
        .abort_design_condition(
            paste0(
                "The aligned design condition must equal the selected ",
                "SummarizedExperiment colData condition for every sample."
            ),
            design$roles$condition,
            mismatched_samples = mismatched
        )
    }

    invisible(expected)
}

.prepare_input_first_matrix <- function(x, design, assay) {
    .validate_matrix_data(x)
    aligned <- .align_missingness_design(design, colnames(x))
    prepared <- .prepare_input(
        x,
        group = .aligned_design_groups(aligned),
        assay = assay
    )

    list(prepared = prepared, design = aligned)
}

.prepare_input_first_se <- function(x, design, assay) {
    condition_role <- design$roles$condition
    .require_se_design_condition(x, condition_role)
    prepared <- .prepare_input(
        x,
        group_col = condition_role,
        assay = assay
    )
    aligned <- .align_missingness_design(design, colnames(prepared$data))
    .validate_se_design_condition(prepared, aligned)

    list(prepared = prepared, design = aligned)
}

.prepare_input_first_analysis <- function(x, design, assay) {
    .validate_missingness_design(design)
    if (is.matrix(x)) {
        return(.prepare_input_first_matrix(x, design, assay))
    }
    if (methods::is(x, "SummarizedExperiment")) {
        return(.prepare_input_first_se(x, design, assay))
    }

    stop(
        "`x` must be an ordinary numeric matrix or a ",
        "SummarizedExperiment object.",
        call. = FALSE
    )
}

.prepared_sidecar_assay <- function(prepared) {
    if (identical(prepared$representation, "matrix")) {
        NA_character_
    } else {
        prepared$assay_name
    }
}

.attempt_input_first_classic <- function(
    prepared,
    original,
    call,
    cutoffs,
    seed
) {
    stage <- new.env(parent = emptyenv())
    stage$current <- "classic"
    observe_stage <- function(value) {
        stage$current <- value
        invisible(value)
    }

    tryCatch(
        .classify_prepared_missingness(
            prepared = prepared,
            original = original,
            matched_call = call,
            cutoffs = cutoffs,
            seed = seed,
            stage_observer = observe_stage
        ),
        error = function(condition) {
            .new_classic_failure(
                condition,
                stage = stage$current,
                call = call
            )
        }
    )
}

.validate_input_first_call <- function(call) {
    if (!is.language(call)) {
        .abort_sidecar(
            "Analysis provenance call must be language.",
            "imputefinder_analysis_schema_error",
            field = "call"
        )
    }
    invisible(call)
}

.validate_base_fit_cutoff_conflict <- function(base_fit, cutoffs) {
    if (!is.null(base_fit) && !is.null(cutoffs)) {
        .abort_sidecar(
            "`base_fit` and `cutoffs` are conflicting specifications.",
            "imputefinder_base_fit_conflict_error",
            arguments = c("base_fit", "cutoffs")
        )
    }
    invisible(base_fit)
}

.resolve_input_first_classic <- function(
    prepared,
    original,
    call,
    base_fit,
    cutoffs,
    seed
) {
    if (is.null(base_fit)) {
        return(.attempt_input_first_classic(
            prepared,
            original,
            call,
            cutoffs,
            seed
        ))
    }
    .accept_compatible_base_fit(
        base_fit,
        prepared,
        original,
        seed
    )
}

.analyze_missingness_input_first <- function(
    x,
    design,
    assay = NULL,
    base_fit = NULL,
    cutoffs = NULL,
    seed = 1L,
    modules = c("sentinel", "stability"),
    call = NULL
) {
    analysis_call <- if (is.null(call)) match.call() else call
    .validate_input_first_call(analysis_call)
    .validate_base_fit_cutoff_conflict(base_fit, cutoffs)
    modules <- .normalise_sidecar_modules(modules)
    seed <- .normalise_rescue_seed(seed)
    resolved <- .prepare_input_first_analysis(x, design, assay)
    prepared <- resolved$prepared
    acquisition <- .expected_sidecar_acquisition(resolved$design)
    input <- .new_sidecar_input(
        prepared$data,
        representation = prepared$representation,
        assay = .prepared_sidecar_assay(prepared),
        acquisition = acquisition
    )
    sentinel <- .run_pre_rescue_sentinel(
        modules,
        prepared$data,
        prepared$groups_by_sample
    )
    classic <- .resolve_input_first_classic(
        prepared,
        x,
        analysis_call,
        base_fit,
        cutoffs,
        seed
    )
    stability <- .run_pending_stability(modules)

    .new_imputefinder_analysis(
        classic = classic,
        design = resolved$design,
        input = input,
        call = analysis_call,
        sentinel = sentinel,
        stability = stability,
        seeds = list(classic_rescue = seed)
    )
}
