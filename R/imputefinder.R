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
#' Rescue condition-specific absences, classify missingness independently within
#' each condition, and retain features suitable for condition-aware downstream
#' imputation. Input values must be unnormalised log2 protein intensities.
#'
#' @param x Numeric matrix or a \code{SummarizedExperiment} object.
#' @param group Atomic condition vector or aligned design data frame for matrix
#'   input; omitted for \code{SummarizedExperiment} input.
#' @param group_col Condition column for design or
#'   \code{SummarizedExperiment} input.
#' @param assay Assay name for \code{SummarizedExperiment} input.
#' @param cutoffs Optional named numeric vector. A cutoff is required for each
#'   condition containing missing values; conditions without missing values need
#'   no cutoff.
#' @param seed Integer rescue seed.
#' @return An \code{imputefinder_result} containing filtered seed-modified data,
#'   per-condition classifications and groups, feature audit status, cutoffs,
#'   diagnostics, and seed provenance.
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
    matched_call <- match.call()
    prepared <- .prepare_input(x, group, group_col, assay)
    rescued <- .seed_missing_conditions(
        prepared$data,
        prepared$groups_by_sample,
        seed
    )
    statistics <- .feature_condition_statistics(
        rescued$data,
        prepared$groups_by_sample,
        rescued$seed_log
    )
    conditions <- sort(
        unique(unname(prepared$groups_by_sample)),
        method = "radix"
    )
    resolved <- .resolve_manual_cutoffs(
        statistics,
        cutoffs,
        conditions
    )
    classified <- .assign_condition_states(statistics, resolved$cutoffs)
    reconciled <- .reconcile_condition_states(
        classified,
        rescued$feature_status
    )
    groups <- .retained_condition_groups(
        reconciled$classifications,
        reconciled$feature_status$feature
    )
    retained_features <- reconciled$feature_status$feature[
        reconciled$feature_status$retained
    ]
    filtered_data <- rescued$data[
        retained_features,
        ,
        drop = FALSE
    ]
    output_data <- .restore_output_data(filtered_data, prepared, x)

    .new_imputefinder_result(
        data = output_data,
        classifications = reconciled$classifications,
        groups = groups,
        feature_status = reconciled$feature_status,
        cutoffs = resolved$cutoffs,
        cutoff_diagnostics = resolved$diagnostics,
        seed_log = rescued$seed_log,
        groups_by_sample = prepared$groups_by_sample,
        call = matched_call
    )
}
