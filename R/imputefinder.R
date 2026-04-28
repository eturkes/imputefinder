#    This file is part of ImputeFinder.
#    Copyright (C) 2025  Emir Turkes, UK DRI at UCL
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

# Names used inside ggplot2 aes() and the running-accumulator inside
# Map() are flagged by R CMD check as "undefined globals"; declare them
# here so the static check stays clean.
utils::globalVariables(c("x", "y", "final_MAR"))

# Internal helper that builds the cumulative intensity plot consumed by
# `classify_missingness()` for elbow detection. Not exported.
plot_detect_custom <- function(se, elbow = TRUE, threshold = 0.35) {
  m <- SummarizedExperiment::assay(se)
  means <- rowMeans(m, na.rm = TRUE)
  means <- means[is.finite(means)]
  if (length(means) < 3L) {
    stop("Not enough finite mean intensities for elbow detection.")
  }
  sorted_means <- sort(means)
  df <- data.frame(
    x = sorted_means,
    y = seq_along(sorted_means) / length(sorted_means)
  )
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line()
}

#' Classify Missingness as MAR or MNAR Based on Intensity Distribution
#'
#' Identifies Missing Not At Random (MNAR) and Missing At Random (MAR)
#' features in a [SummarizedExperiment::SummarizedExperiment] of
#' quantitative proteomics intensities. For each experimental condition
#' the empirical distribution of feature mean intensities is inspected
#' and an elbow on the cumulative curve is used to derive an
#' intensity cutoff: features with mean intensity below the cutoff are
#' classified as MNAR, while the remaining features with at least
#' `min_non_na` non-missing values are classified as MAR.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] object whose
#'   first assay contains the quantitative intensities and whose
#'   `colData` has a `condition` column describing the experimental
#'   groups.
#' @param threshold A numeric value for the second derivative threshold
#'   used to locate the elbow on the cumulative intensity curve
#'   (default `0.35`).
#' @param min_non_na Integer. Minimum number of non-missing values
#'   required for a feature to be retained as MAR (default `5`).
#' @param return_plot Logical. If `TRUE`, a list of `ggplot` objects
#'   (one per condition) is returned alongside the classification
#'   results.
#'
#' @return A named `list` with the following elements:
#' \describe{
#'   \item{`data`}{The input `SummarizedExperiment` filtered to the
#'     union of MAR and MNAR features.}
#'   \item{`MAR`}{Character vector of feature names classified as MAR
#'     in at least one condition.}
#'   \item{`MNAR`}{Character vector of feature names classified as MNAR
#'     in at least one condition.}
#'   \item{`cutoffs`}{Named numeric vector giving the intensity cutoff
#'     used for each condition.}
#'   \item{`plot`}{Optional list of `ggplot` objects, one per
#'     condition, returned only when `return_plot = TRUE`.}
#' }
#'
#' @examples
#' library(SummarizedExperiment)
#'
#' # `classify_missingness()` requires a SummarizedExperiment as input;
#' # supplying any other object raises an informative error.
#' res <- try(classify_missingness(matrix(1:6, nrow = 2)), silent = TRUE)
#' inherits(res, "try-error")
#'
#' # It also requires that the assay actually contains missing values,
#' # otherwise classification is not meaningful.
#' counts <- matrix(
#'   seq_len(6),
#'   nrow = 2,
#'   dimnames = list(c("p1", "p2"), c("s1", "s2", "s3"))
#' )
#' se <- SummarizedExperiment(
#'   assays = list(intensity = counts),
#'   colData = DataFrame(condition = c("A", "A", "B"))
#' )
#' res <- try(classify_missingness(se), silent = TRUE)
#' inherits(res, "try-error")
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ggplot2 ggplot_build
#' @export
classify_missingness <- function(se, threshold = 0.35, min_non_na = 5, return_plot = FALSE) {
  stopifnot(inherits(se, "SummarizedExperiment"))
  if (!any(is.na(assay(se)))) {
    stop("No missing values detected in assay.")
  }

  data <- se
  conditions <- unique(colData(data)$condition)
  MNAR <- list()
  MAR <- list()
  cutoffs <- numeric()
  plots <- list()

  for (cond in conditions) {
    idx <- which(colData(data)$condition == cond)
    subset_data <- data[, idx]

    plot_obj <- plot_detect_custom(subset_data, elbow = TRUE, threshold = threshold)
    plot_data <- ggplot_build(plot_obj)
    x_vals <- plot_data$data[[1]]$x
    y_vals <- plot_data$data[[1]]$y

    # Find elbow
    d1 <- diff(y_vals) / diff(x_vals)
    d2 <- diff(d1) / diff(x_vals[-1])
    elbow_x <- min(x_vals[which(abs(d2) > threshold)])
    cutoffs[cond] <- elbow_x

    # Determine MNAR
    means <- rowMeans(assay(data)[, idx], na.rm = TRUE)
    order_idx <- order(means)
    data <- data[order_idx, ]
    mnar_features <- names(which(rowMeans(assay(data)[, idx], na.rm = TRUE) < elbow_x))
    MNAR[[cond]] <- mnar_features

    # Remove MAR proteins with mostly missing values
    mar_data <- data[!(rownames(data) %in% mnar_features), idx]
    keep <- rowSums(!is.na(assay(mar_data))) >= min_non_na
    MAR[[cond]] <- rownames(mar_data)[keep]

    if (return_plot) {
      plots[[cond]] <- plot_obj + ggplot2::geom_vline(xintercept = elbow_x) +
        ggplot2::ggtitle(cond)
    }
  }

  # Create unified lists
  final_MAR <- Reduce(union, Map(function(x, y) intersect(x, union(y, final_MAR %||% character(0))), MAR, MNAR))
  final_MNAR <- Reduce(union, MNAR)

  keep <- unique(c(final_MAR, final_MNAR))
  data <- data[rownames(data) %in% keep, ]

  output <- list(
    data = data,
    MAR = final_MAR,
    MNAR = final_MNAR,
    cutoffs = cutoffs
  )
  if (return_plot) {
    output$plot <- plots
  }
  return(output)
}
