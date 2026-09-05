#' ImputeFinder: Condition-Aware Missingness Classification
#'
#' ImputeFinder classifies protein-intensity missingness independently within
#' each experimental condition. This keeps condition-specific on/off patterns
#' visible and prepares retained features for downstream imputation.
#'
#' The `MNAR` and `MAR` labels are working classifications based on observed
#' intensity and missingness. They do not prove the underlying mechanism. The
#' result includes classifications, cutoff information, a rescue log, and the
#' reason each feature was retained or dropped.
#'
#' @section Expected input:
#' [classify_missingness()] accepts unnormalised log2 protein intensities with
#' features in rows, named samples in columns, and one condition label per
#' sample. Missing observations must be represented by `NA`. Ordinary numeric
#' matrices and [SummarizedExperiment::SummarizedExperiment] objects share the
#' same matrix core.
#'
#' @section Workflow boundary:
#' ImputeFinder filters the input. The default minimum-rescue policy adds
#' recorded anchors; `rescue = "none"` keeps fully unobserved blocks `absent`.
#' It does not normalise or impute the data. Normalise the returned data, split
#' it by condition, and apply the chosen methods to the supplied feature groups.
#' Absent blocks require a separate downstream decision; they are not assigned
#' to either imputation branch.
#'
#' @seealso [classify_missingness()], [plot_missingness()],
#'   `vignette("imputefinder")`, and `citation("imputefinder")`.
#' @keywords internal
"_PACKAGE"
