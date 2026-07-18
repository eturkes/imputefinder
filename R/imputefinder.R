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
#' @param cutoffs Optional named numeric vector of manual condition cutoffs.
#'   Automatic selection fills omitted conditions containing missing values;
#'   conditions without missing values need no cutoff.
#' @param seed Integer rescue seed.
#' @return An \code{imputefinder_result} containing filtered seed-modified data,
#'   per-condition classifications and groups, feature audit status, cutoffs,
#'   explicit missingness profiles, diagnostics, and seed provenance.
#' @seealso \code{\link{plot_missingness}} and
#'   \code{vignette("imputefinder")}.
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
    .classify_prepared_missingness(
        prepared,
        x,
        matched_call,
        cutoffs,
        seed
    )
}

.observe_classic_stage <- function(stage_observer, stage) {
    if (!is.null(stage_observer)) {
        stage_observer(stage)
    }
    invisible(stage)
}

.resolve_prepared_cutoffs <- function(
    prepared,
    statistics,
    profiles,
    cutoffs
) {
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

    resolved
}

.classify_prepared_missingness <- function(
    prepared,
    original,
    matched_call,
    cutoffs,
    seed,
    stage_observer = NULL
) {
    .observe_classic_stage(stage_observer, "rescue")
    rescued <- .seed_missing_conditions(
        prepared$data,
        prepared$groups_by_sample,
        seed
    )
    .observe_classic_stage(stage_observer, "block_statistics")
    statistics <- .feature_condition_statistics(
        rescued$data,
        prepared$groups_by_sample,
        rescued$seed_log
    )
    .observe_classic_stage(stage_observer, "profile")
    profiles <- .build_missingness_profiles(statistics)
    .observe_classic_stage(stage_observer, "cutoff")
    resolved <- .resolve_prepared_cutoffs(
        prepared,
        statistics,
        profiles,
        cutoffs
    )
    .observe_classic_stage(stage_observer, "state")
    classified <- .assign_condition_states(statistics, resolved$cutoffs)
    .observe_classic_stage(stage_observer, "reconciliation")
    reconciled <- .reconcile_condition_states(
        classified,
        rescued$feature_status
    )
    .observe_classic_stage(stage_observer, "result")
    .finalise_missingness_result(
        prepared,
        original,
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
