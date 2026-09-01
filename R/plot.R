utils::globalVariables(c("intensity", "mean_intensity", "missing_proportion"))

#' Plot a Stored Missingness Profile
#'
#' Plot the count-weighted proportion of features with missing observations for
#' one condition. The curve is read directly from the profile stored in
#' \code{result}; the input data and classifications are not recomputed.
#'
#' @param result An \code{imputefinder_result} returned by
#'   \code{\link{classify_missingness}}.
#' @param condition One condition label in \code{result}.
#' @return A \code{ggplot} object. A recorded cutoff is shown as a vertical
#'   line; bottom ticks identify seeded feature-condition blocks. Stored profile
#'   and cutoff warnings appear in the caption.
#' @seealso \code{\link{classify_missingness}} and
#'   \code{vignette("imputefinder")}.
#' @examples
#' x <- rbind(
#'     on_off = c(NA, NA, NA, NA, 15, 16, 15, 17),
#'     mar = c(14, 15, NA, 16, 14, NA, 15, 16),
#'     low_anchor = c(8, NA, NA, NA, 14, 15, 16, 17)
#' )
#' colnames(x) <- paste0("sample", seq_len(ncol(x)))
#' group <- rep(c("A", "B"), each = 4)
#' fit <- classify_missingness(x, group, cutoffs = c(A = 12, B = 12))
#'
#' plot_missingness(fit, "A")
#'
#' @export
plot_missingness <- function(result, condition) {
    profile <- .stored_plot_profile(result, condition)
    cutoff <- .stored_plot_cutoff(result, condition)
    diagnostic <- .stored_plot_diagnostic(result, condition)
    plot <- .base_missingness_plot(profile)
    plot <- .add_missingness_markers(plot, profile, cutoff)
    .label_missingness_plot(plot, condition, cutoff, diagnostic, profile)
}

.base_missingness_plot <- function(profile) {
    ggplot2::ggplot(
        profile$grid,
        ggplot2::aes(x = intensity, y = missing_proportion)
    ) +
        ggplot2::geom_area(
            fill = "#0F766E",
            alpha = 0.18,
            na.rm = TRUE
        ) +
        ggplot2::geom_line(
            colour = "#0F5E62",
            linewidth = 0.8,
            na.rm = TRUE
        )
}

.add_missingness_markers <- function(plot, profile, cutoff) {
    seeded <- profile$raw[profile$raw$seeded, , drop = FALSE]
    if (nrow(seeded) > 0L) {
        plot <- plot + ggplot2::geom_rug(
            data = seeded,
            mapping = ggplot2::aes(x = mean_intensity),
            inherit.aes = FALSE,
            sides = "b",
            colour = "#C2410C",
            linewidth = 0.55
        )
    }
    if (is.finite(cutoff)) {
        plot <- plot + ggplot2::geom_vline(
            xintercept = cutoff,
            colour = "#B91C1C",
            linetype = "22",
            linewidth = 0.7
        )
    }

    plot
}

.label_missingness_plot <- function(
    plot,
    condition,
    cutoff,
    diagnostic,
    profile
) {
    plot +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.01, 0.03))
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.25),
            labels = .percentage_labels,
            expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
            title = sprintf("Missingness profile - condition %s", condition),
            subtitle = .cutoff_subtitle(cutoff, diagnostic),
            x = "Mean log2 intensity",
            y = "Features with missing values",
            caption = .plot_caption(profile, diagnostic)
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.title.position = "plot",
            plot.caption.position = "plot",
            plot.caption = ggplot2::element_text(
                colour = "#4B5563",
                hjust = 0
            )
        )
}

.stored_plot_profile <- function(result, condition) {
    if (!inherits(result, "imputefinder_result")) {
        stop(
            "`result` must inherit from `imputefinder_result`.",
            call. = FALSE
        )
    }
    if (!is.character(condition) || length(condition) != 1L ||
        is.na(condition) || !nzchar(condition)) {
        stop("`condition` must be one non-empty string.", call. = FALSE)
    }

    conditions <- names(result$profiles)
    valid_conditions <- is.character(conditions) &&
        length(conditions) > 0L &&
        !anyNA(conditions) &&
        all(nzchar(conditions)) &&
        !anyDuplicated(conditions)
    if (!valid_conditions) {
        stop("Result has no valid condition profiles.", call. = FALSE)
    }
    if (!condition %in% conditions) {
        stop(
            sprintf(
                "Unknown condition `%s`; available conditions: %s.",
                condition,
                paste(sprintf("`%s`", conditions), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    profile <- result$profiles[[condition]]
    if (is.null(profile)) {
        stop(
            sprintf(
                "Result has no stored missingness profile for condition `%s`.",
                condition
            ),
            call. = FALSE
        )
    }
    .validate_stored_plot_profile(profile, condition)
    profile
}

.validate_stored_plot_profile <- function(profile, condition) {
    valid_container <- is.list(profile) &&
        all(c("raw", "grid", "metadata") %in% names(profile)) &&
        is.data.frame(profile$raw) &&
        is.data.frame(profile$grid) &&
        is.list(profile$metadata)
    if (!valid_container) {
        stop(
            sprintf("Stored profile for condition `%s` is invalid.", condition),
            call. = FALSE
        )
    }

    grid <- profile$grid
    raw <- profile$raw
    required_grid <- c("intensity", "missing_proportion", "supported")
    valid_grid <- all(required_grid %in% names(grid)) &&
        nrow(grid) > 0L &&
        is.numeric(grid$intensity) &&
        is.numeric(grid$missing_proportion) &&
        is.logical(grid$supported) &&
        !anyNA(grid$intensity) &&
        !anyNA(grid$supported) &&
        all(is.finite(grid$intensity)) &&
        all(
            !grid$supported |
                (
                    is.finite(grid$missing_proportion) &
                        grid$missing_proportion >= 0 &
                        grid$missing_proportion <= 1
                )
        )
    valid_raw <- all(c("mean_intensity", "seeded") %in% names(raw)) &&
        is.numeric(raw$mean_intensity) &&
        is.logical(raw$seeded) &&
        !anyNA(raw$mean_intensity) &&
        !anyNA(raw$seeded) &&
        all(is.finite(raw$mean_intensity))
    if (!valid_grid || !valid_raw) {
        stop(
            sprintf("Stored profile for condition `%s` is invalid.", condition),
            call. = FALSE
        )
    }

    invisible(profile)
}

.stored_plot_cutoff <- function(result, condition) {
    cutoffs <- result$cutoffs
    valid <- is.numeric(cutoffs) &&
        is.null(dim(cutoffs)) &&
        !is.null(names(cutoffs)) &&
        condition %in% names(cutoffs)
    if (!valid) {
        stop(
            sprintf(
                "Result has no recorded cutoff for condition `%s`.",
                condition
            ),
            call. = FALSE
        )
    }

    cutoff <- unname(cutoffs[[condition]])
    if (length(cutoff) != 1L || is.nan(cutoff) ||
        (!is.na(cutoff) && !is.finite(cutoff))) {
        stop(
            sprintf(
                "Recorded cutoff for condition `%s` is invalid.",
                condition
            ),
            call. = FALSE
        )
    }

    cutoff
}

.stored_plot_diagnostic <- function(result, condition) {
    diagnostics <- result$cutoff_diagnostics
    if (!is.list(diagnostics) || is.null(names(diagnostics)) ||
        !condition %in% names(diagnostics)) {
        return(NULL)
    }

    diagnostic <- diagnostics[[condition]]
    if (!is.null(diagnostic) && !is.list(diagnostic)) {
        stop(
            sprintf(
                "Stored cutoff diagnostic for condition `%s` is invalid.",
                condition
            ),
            call. = FALSE
        )
    }
    diagnostic
}

.cutoff_subtitle <- function(cutoff, diagnostic) {
    if (!is.finite(cutoff)) {
        return(NULL)
    }

    source <- if (is.list(diagnostic) &&
        is.character(diagnostic$source) &&
        length(diagnostic$source) == 1L &&
        !is.na(diagnostic$source) &&
        nzchar(diagnostic$source)) {
        diagnostic$source
    } else {
        "recorded"
    }
    sprintf(
        "Cutoff: %s (%s)",
        format(cutoff, trim = TRUE, scientific = FALSE),
        source
    )
}

.percentage_labels <- function(values) {
    paste0(100 * values, "%")
}

.plot_caption <- function(profile, diagnostic) {
    notes <- if (any(profile$raw$seeded)) {
        "Bottom ticks mark seeded feature-condition blocks."
    } else {
        character()
    }
    profile_warnings <- .plot_warnings(profile$metadata$warnings, "profile")
    cutoff_warnings <- .plot_warnings(
        if (is.list(diagnostic)) diagnostic$warnings else NULL,
        "cutoff diagnostic"
    )
    warnings <- unique(c(profile_warnings, cutoff_warnings))
    if (length(warnings) > 0L) {
        notes <- c(notes, paste0("Warning: ", warnings))
    }

    if (length(notes) == 0L) NULL else paste(notes, collapse = "\n")
}

.plot_warnings <- function(warnings, source) {
    if (is.null(warnings)) {
        return(character())
    }
    if (!is.character(warnings) || anyNA(warnings)) {
        stop("Stored ", source, " notes are invalid.", call. = FALSE)
    }

    warnings[nzchar(warnings)]
}
