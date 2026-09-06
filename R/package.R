#' ImputeFinder: Condition-Aware Missingness Classification
#'
#' ImputeFinder classifies protein-intensity missingness by condition, preserving
#' on/off patterns and preparing retained features for downstream imputation.
#'
#' `MNAR` and `MAR` are intensity-based classifications, not proven mechanisms.
#' Results record states, cutoffs, rescue values, and retention or drop reasons.
#'
#' @section Expected input:
#' [classify_missingness()] accepts unnormalised log2 protein intensities with
#' features in rows, named samples in columns, and one condition per sample.
#' Use `NA` for missingness. Numeric matrices and
#' [SummarizedExperiment::SummarizedExperiment] objects share the same matrix core.
#'
#' @section Workflow boundary:
#' The default minimum-rescue policy adds logged anchors to filtered data;
#' `rescue = "none"` keeps unobserved blocks `absent`. ImputeFinder neither
#' normalises nor imputes. Normalise returned data, split by condition, and
#' apply chosen methods to the supplied groups. Handle absent blocks separately
#' from both imputation branches.
#'
#' @seealso [classify_missingness()], [plot_missingness()],
#'   `vignette("imputefinder")`, and `citation("imputefinder")`.
#' @keywords internal
"_PACKAGE"
