plot_geom_classes <- function(plot) {
    vapply(
        plot$layers,
        function(layer) class(layer$geom)[[1L]],
        character(1L)
    )
}

test_that("missingness plots use the stored profile and recorded cutoff", {
    result <- classify_normative_fixture()
    stored_grid <- result$profiles$A$grid

    # Conflicting live values make profile recomputation detectable.
    result$data[,] <- 999
    result$classifications$mean_intensity <- -999
    result$classifications$missing_fraction <- 0

    plot <- plot_missingness(result, "A")

    expect_s3_class(plot, "ggplot")
    expect_identical(plot$data, stored_grid)
    expect_identical(plot$labels$title, "Missingness profile - condition A")
    expect_identical(plot$labels$x, "Mean log2 intensity")
    expect_identical(plot$labels$y, "Features with missing values")

    y_scale <- plot$scales$get_scales("y")
    expect_identical(y_scale$limits, c(0, 1))
    expect_identical(y_scale$breaks, seq(0, 1, by = 0.25))
    expect_identical(
        y_scale$labels(c(0, 0.25, 0.5, 0.75, 1)),
        c("0%", "25%", "50%", "75%", "100%")
    )

    geom_classes <- plot_geom_classes(plot)
    expect_true("GeomVline" %in% geom_classes)
    built <- expect_silent(ggplot2::ggplot_build(plot))
    cutoff_layer <- built$data[[match("GeomVline", geom_classes)]]
    expect_identical(unique(cutoff_layer$xintercept), 12)
})

test_that("missingness plots mark seeded feature-condition blocks", {
    result <- classify_normative_fixture()
    plot <- plot_missingness(result, "A")
    geom_classes <- plot_geom_classes(plot)

    expect_true("GeomRug" %in% geom_classes)
    built <- expect_silent(ggplot2::ggplot_build(plot))
    seed_layer <- built$data[[match("GeomRug", geom_classes)]]
    expect_identical(unique(seed_layer$x), 8)
    expect_match(
        plot$labels$caption,
        "Bottom ticks mark seeded feature-condition blocks.",
        fixed = TRUE
    )
})

test_that("missingness plots surface stored warnings without a default warning", {
    result <- classify_normative_fixture()
    expect_null(plot_missingness(result, "B")$labels$caption)

    result$profiles$B$metadata$warnings <- "Profile support warning."
    result$cutoff_diagnostics$B$warnings <- "Cutoff quality warning."
    plot <- plot_missingness(result, "B")

    expect_match(
        plot$labels$caption,
        "Warning: Profile support warning.",
        fixed = TRUE
    )
    expect_match(
        plot$labels$caption,
        "Warning: Cutoff quality warning.",
        fixed = TRUE
    )
})

test_that("complete-only plots omit irrelevant cutoff and seed layers", {
    x <- rbind(
        low = c(1, 2, 3, 4),
        high = c(5, 6, 7, 8)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    result <- classify_missingness(x, rep(c("A", "B"), each = 2L))

    plot <- plot_missingness(result, "A")
    geom_classes <- plot_geom_classes(plot)

    expect_false("GeomVline" %in% geom_classes)
    expect_false("GeomRug" %in% geom_classes)
    expect_null(plot$labels$caption)
    expect_silent(ggplot2::ggplot_build(plot))
})

test_that("missingness plots reject invalid results and conditions", {
    result <- classify_normative_fixture()

    expect_error(
        plot_missingness(list(), "A"),
        "`result` must inherit from `imputefinder_result`.",
        fixed = TRUE
    )
    for (condition in list(character(), c("A", "B"), NA_character_, "")) {
        expect_error(
            plot_missingness(result, condition),
            "`condition` must be one non-empty string.",
            fixed = TRUE
        )
    }
    expect_error(
        plot_missingness(result, "unknown"),
        "Unknown condition `unknown`; available conditions: `A`, `B`.",
        fixed = TRUE
    )

    result$profiles["A"] <- list(NULL)
    expect_error(
        plot_missingness(result, "A"),
        "Result has no stored missingness profile for condition `A`.",
        fixed = TRUE
    )
})
