#' ImputeFinder: Condition-Aware Missingness Classification
#'
#' ImputeFinder classifies protein-intensity missingness independently within
#' each experimental condition. It preserves evidence-backed,
#' condition-specific on/off features, then returns retained features in
#' condition-specific groups for downstream mixed imputation.
#'
#' The `MNAR` and `MAR` labels describe likely mechanisms supported by the
#' observed intensity and missingness pattern; they do not prove the true
#' missingness mechanism. Results retain the classifications, cutoff profiles
#' and diagnostics, rescue-seed provenance, and feature-level drop reasons
#' needed to audit those decisions.
#'
#' @section Expected input:
#' [classify_missingness()] accepts unnormalised log2 protein intensities with
#' features in rows, named samples in columns, and one condition label per
#' sample. Missing observations must be represented by `NA`. Ordinary numeric
#' matrices and [SummarizedExperiment::SummarizedExperiment] objects share the
#' same matrix core.
#'
#' @section Workflow boundary:
#' ImputeFinder filters and seed-modifies the input; it does not normalise or
#' impute it. Use the returned data for normalisation, split it by condition,
#' and apply the chosen MNAR and MAR methods using the corresponding retained
#' feature groups.
#'
#' The experimental [analyze_missingness()] sidecar starts from the original
#' input and a typed [missingness_design()]. It records identity, canonical
#' design rank/aliasing and declared units, and pre-rescue sample, condition,
#' role-coverage, and feature-overlap evidence without changing the stable
#' classifier result.
#'
#' @seealso [classify_missingness()], [analyze_missingness()],
#'   [plot_missingness()], and `vignette("imputefinder")`.
#' @keywords internal
"_PACKAGE"
