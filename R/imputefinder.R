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

#' Classify missing values by condition
#'
#' Classify each feature separately within each experimental condition. Fully
#' missing condition blocks receive one recorded low-intensity rescue value so
#' that condition-specific on/off patterns can be retained. Input values must be
#' unnormalised log2 protein intensities.
#'
#' @param x An ordinary numeric matrix or a \code{SummarizedExperiment}. Put
#'   features in rows and samples in columns.
#' @param group A condition vector or an aligned design data frame for matrix
#'   input. Do not supply it for \code{SummarizedExperiment} input.
#' @param group_col The condition column in a design data frame or
#'   \code{SummarizedExperiment} object.
#' @param assay The assay name for \code{SummarizedExperiment} input.
#' @param cutoffs Optional named numeric vector of manual condition cutoffs.
#'   Automatic estimation is used for omitted conditions that contain missing
#'   values. Conditions without missing values need no cutoff.
#' @param seed An integer that controls which sample receives each rescue value.
#' @return An \code{imputefinder_result} containing filtered data,
#'   feature-condition classifications, downstream groups, feature-level
#'   decisions, cutoffs, diagnostic profiles, and a log of rescue values.
#' @seealso \code{\link{plot_missingness}},
#'   \code{vignette("imputefinder")}, and
#'   \code{citation("imputefinder")}.
#' @examples
#' x <- rbind(
#'     on_off = c(NA, NA, NA, NA, 15, 16, 15, 17),
#'     mar = c(14, 15, NA, 16, 14, NA, 15, 16),
#'     low_anchor = c(8, NA, NA, NA, 14, 15, 16, 17)
#' )
#' colnames(x) <- paste0("sample", seq_len(ncol(x)))
#' group <- rep(c("A", "B"), each = 4)
#'
#' fit <- classify_missingness(
#'     x,
#'     group = group,
#'     cutoffs = c(A = 12, B = 12)
#' )
#' fit
#' fit$feature_status
#' fit$groups$A
#' fit$seed_log
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
    profiles <- .build_missingness_profiles(statistics)
    conditions <- sort(
        unique(unname(prepared$groups_by_sample)),
        method = "radix"
    )
    resolved <- .resolve_cutoffs(
        statistics,
        profiles,
        cutoffs,
        conditions
    )
    for (condition in conditions) {
        resolved$diagnostics[[condition]]$profile <-
            profiles[[condition]]$metadata
    }
    classified <- .assign_condition_states(statistics, resolved$cutoffs)
    reconciled <- .reconcile_condition_states(
        classified,
        rescued$feature_status
    )
    .finalise_missingness_result(
        prepared,
        x,
        matched_call,
        rescued,
        resolved,
        profiles,
        reconciled
    )
}

.finalise_missingness_result <- function(
    prepared,
    original,
    matched_call,
    rescued,
    resolved,
    profiles,
    reconciled
) {
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
    output_data <- .restore_output_data(filtered_data, prepared, original)

    .new_imputefinder_result(
        data = output_data,
        classifications = reconciled$classifications,
        groups = groups,
        feature_status = reconciled$feature_status,
        cutoffs = resolved$cutoffs,
        cutoff_diagnostics = resolved$diagnostics,
        profiles = profiles,
        seed_log = rescued$seed_log,
        groups_by_sample = prepared$groups_by_sample,
        call = matched_call
    )
}
