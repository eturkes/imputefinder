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

#' Classify Missingness by Condition
#'
#' Prepare condition-aligned quantitative proteomics data for condition-specific
#' missingness classification. The scientific classification core is introduced
#' in the subsequent implementation milestones.
#'
#' @param x Numeric matrix or a \code{SummarizedExperiment} object.
#' @param group Atomic condition vector or aligned design data frame for matrix
#'   input; omitted for \code{SummarizedExperiment} input.
#' @param group_col Condition column for design or
#'   \code{SummarizedExperiment} input.
#' @param assay Assay name for \code{SummarizedExperiment} input.
#' @param cutoffs Optional named numeric condition cutoffs.
#' @param seed Integer rescue seed.
#' @return An \code{imputefinder_result} object once classification is complete.
#'
#' @export
classify_missingness <- function(
    x,
    group = NULL,
    group_col = NULL,
    assay = NULL,
    cutoffs = NULL,
    seed = 1L
) {
    .prepare_input(x, group, group_col, assay)
    stop(
        "The missingness classification core is not implemented yet.",
        call. = FALSE
    )
}
