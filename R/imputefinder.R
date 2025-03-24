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

#' Classify Missingness as MAR or MNAR Based on Intensity Distribution
#'
#' This function identifies MNAR and MAR features using intensity distributions and missingness patterns.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param threshold A numeric value for the second derivative threshold (default = 0.35).
#' @param min_non_na Integer. Minimum number of non-missing values required for MAR retention (default = 5).
#' @param return_plot Logical. If TRUE, a ggplot object is returned alongside results.
#' @return A list with:
#'   \item{data}{Filtered \code{SummarizedExperiment} object.}
#'   \item{MAR}{Character vector of MAR feature names.}
#'   \item{MNAR}{Character vector of MNAR feature names.}
#'   \item{cutoffs}{Named numeric vector of intensity cutoffs per condition.}
#'   \item{plot}{List of ggplot objects, if \code{return_plot = TRUE}.}
#'
#' @importFrom SummarizedExperiment assay colData
#' @export
classify_missingness <- function(se, threshold = 0.35, min_non_na = 5, return_plot = FALSE) {
  stopifnot(inherits(se, "SummarizedExperiment"))
  if (!any(is.na(assay(se)))) stop("No missing values detected in assay.")

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
