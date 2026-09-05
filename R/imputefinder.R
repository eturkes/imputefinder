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
#' missing condition blocks can be retained without inserting values. The
#' legacy minimum-rescue mode remains the default. Input values must be
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
#'   Automatic estimation is used for omitted conditions with partially observed
#'   blocks. Complete and fully absent blocks need no cutoff.
#' @param seed An integer that controls which sample receives each rescue value.
#' @param rescue Either \code{"minimum"}, the legacy default, or \code{"none"}.
#'   Minimum rescue inserts one condition-minimum value into each fully missing
#'   block. With no rescue, these blocks remain missing and receive the state
#'   \code{absent}. This state does not identify a missingness mechanism.
#' @param ambiguity Policy for multiple credible automatic transitions:
#'   \code{"warn"} selects the lowest transition and signals a warning;
#'   \code{"error"} requires a manual cutoff; \code{"lowest"} explicitly accepts
#'   the lowest transition with stored diagnostics only. The default is
#'   \code{"warn"}.
#' @details
#' When several credible downward transitions exist, automatic estimation
#' follows \code{ambiguity}. Candidate transitions remain in the diagnostics.
#' Inspect the stored profile and supply a manual cutoff when appropriate.
#'
#' With \code{rescue = "none"}, absent blocks do not contribute to intensity
#' profiles or cutoff estimation. They need no cutoff. A condition can contain
#' no observations. Features are dropped when all values are missing, any
#' condition is insufficient, or every condition is MNAR or absent. Other
#' on/off features remain available through each condition's \code{absent}
#' group. No missing value is replaced in this mode.
#'
#' With one sample in a condition, minimum rescue fills the missing block, which
#' then becomes \code{complete}. Inspect \code{seeded} and \code{seed_log} to
#' distinguish these blocks from originally complete measurements. Replicated
#' conditions are needed to assess within-condition missingness.
#' @return An \code{imputefinder_result} containing filtered data,
#'   feature-condition classifications, downstream groups, feature-level
#'   decisions, cutoffs, diagnostic profiles, and a log of rescue values.
#'   Classifications include \code{original_observed_count},
#'   \code{original_missing_count}, and \code{original_mean_intensity} before
#'   rescue. The resolved \code{rescue} and \code{ambiguity} policies are stored
#'   in \code{policies}.
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
    seed = 1L,
    rescue = "minimum",
    ambiguity = "warn"
) {
    matched_call <- match.call()
    rescue <- .missingness_policy(rescue, c("minimum", "none"), "rescue")
    ambiguity <- .missingness_policy(
        ambiguity, c("warn", "error", "lowest"), "ambiguity"
    )
    prepared <- .prepare_input(x, group, group_col, assay)
    rescued <- .seed_missing_conditions(
        prepared$data,
        prepared$groups_by_sample,
        seed,
        rescue
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
    profiles <- .build_missingness_profiles(statistics, conditions)
    resolved <- .resolve_cutoffs(
        statistics,
        profiles,
        cutoffs,
        conditions,
        ambiguity
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
        reconciled,
        policies = list(rescue = rescue, ambiguity = ambiguity)
    )
}

.missingness_policy <- function(value, choices, argument) {
    valid <- is.character(value) && !is.object(value) &&
        is.null(dim(value)) && length(value) == 1L &&
        !is.na(value) && value %in% choices
    if (!valid) {
        stop(sprintf("`%s` must be one of: %s.", argument,
            paste(choices, collapse = ", ")), call. = FALSE)
    }
    unname(value)
}

.finalise_missingness_result <- function(
    prepared,
    original,
    matched_call,
    rescued,
    resolved,
    profiles,
    reconciled,
    policies
) {
    groups <- .retained_condition_groups(
        reconciled$classifications,
        reconciled$feature_status$feature,
        conditions = names(resolved$cutoffs)
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

    classifications <- reconciled$classifications
    classifications$original_observed_count <-
        classifications$observed_count - as.integer(classifications$seeded)
    classifications$original_missing_count <-
        classifications$missing_count + as.integer(classifications$seeded)
    classifications$original_mean_intensity <- classifications$mean_intensity
    classifications$original_mean_intensity[classifications$seeded] <- NA_real_

    .new_imputefinder_result(
        data = output_data,
        classifications = classifications,
        groups = groups,
        feature_status = reconciled$feature_status,
        cutoffs = resolved$cutoffs,
        cutoff_diagnostics = resolved$diagnostics,
        profiles = profiles,
        seed_log = rescued$seed_log,
        groups_by_sample = prepared$groups_by_sample,
        call = matched_call,
        policies = policies
    )
}
