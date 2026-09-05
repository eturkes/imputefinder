observation_policy_fixture <- function() {
    x <- rbind(
        on_off = c(NA, NA, NA, 15, 16, 17),
        mar = c(15, NA, 16, 15, 16, 17),
        low = c(8, NA, NA, 15, 16, 17),
        complete = c(16, 16, 16, 16, 16, 16),
        all_missing = rep(NA_real_, 6)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    list(x = x, group = rep(c("A", "B"), each = 3))
}

render_policy_profile <- function(result, condition) {
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(grDevices::dev.off(device), add = TRUE)
    ggplot2::ggplotGrob(plot_missingness(result, condition))
}

test_that("no-rescue mode retains absence without inventing observations", {
    fixture <- observation_policy_fixture()
    result <- classify_missingness(fixture$x, fixture$group,
        cutoffs = c(A = 12), rescue = "none")
    expect_identical(result$data, fixture$x[1:4, , drop = FALSE])
    expect_identical(nrow(result$seed_log), 0L)
    expect_identical(result$policies, list(rescue = "none", ambiguity = "warn"))
    absent <- classification_row(result, "on_off", "A")
    expect_identical(absent$state, "absent")
    expect_identical(absent$original_observed_count, 0L)
    expect_identical(absent$original_missing_count, 3L)
    expect_true(is.na(absent$mean_intensity))
    expect_true(is.na(absent$original_mean_intensity))
    expect_false(absent$seeded)
    expect_true(absent$retained)
    expect_identical(result$groups$A$absent, "on_off")
    expect_identical(result$groups$A$MNAR, "low")
    expect_identical(result$groups$A$MAR, "mar")
    expect_identical(result$profiles$A$metadata$excluded_absent_feature_count, 1L)
    expect_identical(result$profiles$A$metadata$estimated_feature_count, 3L)
    expect_match(plot_missingness(result, "A")$labels$caption, "unobserved")
    device_before <- grDevices::dev.cur()
    expect_silent(render_policy_profile(result, "A"))
    expect_identical(grDevices::dev.cur(), device_before)
    expect_identical(summary(result)$states$absent, c(1L, 0L))
})

test_that("absence does not contribute to automatic cutoff evidence", {
    fixture <- automatic_cutoff_matrix_fixture()
    base <- classify_missingness(fixture$x, fixture$group,
        rescue = "none", ambiguity = "lowest")
    absent <- matrix(15, 10, ncol(fixture$x),
        dimnames = list(paste0("absent", 1:10), colnames(fixture$x)))
    absent[, fixture$group == "A"] <- NA_real_
    extended <- classify_missingness(rbind(fixture$x, absent), fixture$group,
        rescue = "none", ambiguity = "lowest")
    expect_identical(extended$cutoffs[["A"]], base$cutoffs[["A"]])
    expect_identical(extended$profiles$A$grid, base$profiles$A$grid)
    expect_identical(extended$cutoff_diagnostics$A$quality,
        base$cutoff_diagnostics$A$quality)
    expect_identical(extended$profiles$A$metadata$excluded_absent_feature_count, 10L)
})

test_that("unobserved conditions need no invented minimum or cutoff", {
    x <- rbind(on = c(NA, NA, 10, 11), off = c(NA, NA, 12, 13))
    colnames(x) <- paste0("s", 1:4)
    result <- classify_missingness(x, rep(c("A", "B"), each = 2),
        rescue = "none")
    expect_identical(result$data, x)
    expect_true(all(is.na(result$cutoffs)))
    expect_identical(result$groups$A$absent, rownames(x))
    expect_identical(result$profiles$A$metadata$profile_type, "absent_only")
    expect_identical(nrow(result$profiles$A$grid), 0L)
    expect_identical(result$cutoff_diagnostics$A$source, "not_needed")
    expect_silent(render_policy_profile(result, "A"))
    expect_output(print(result), "2 absent")
})

test_that("globally missing data have an empty but inspectable no-rescue result", {
    x <- matrix(NA_real_, 2, 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2")))
    result <- classify_missingness(x, c("A", "B"), rescue = "none")
    expect_identical(dim(result$data), c(0L, 2L))
    expect_identical(result$feature_status$drop_reason, rep("all_missing", 2))
    expect_identical(names(result$groups), c("A", "B"))
    expect_identical(names(result$profiles), c("A", "B"))
    expect_identical(nrow(result$classifications), 0L)
    expect_identical(result$profiles$A$metadata$profile_type, "empty")
    expect_identical(summary(result)$features[["dropped"]], 2L)
    expect_silent(render_policy_profile(result, "A"))
})

test_that("legacy rescue exposes original counts and synthetic completeness", {
    x <- rbind(on_off = c(NA, 12), anchor = c(8, 14))
    colnames(x) <- c("a", "b")
    result <- classify_missingness(x, c("A", "B"))
    block <- classification_row(result, "on_off", "A")
    expect_identical(block$state, "complete")
    expect_identical(block$observed_count, 1L)
    expect_identical(block$original_observed_count, 0L)
    expect_identical(block$original_missing_count, 1L)
    expect_true(is.na(block$original_mean_intensity))
    expect_identical(summary(result)$seeded_complete, 1L)
    expect_output(print(result), "synthetic completeness")
    expect_identical(result$policies$rescue, "minimum")
    plain <- classify_missingness(x, c("A", "B"), rescue = "none")
    expect_identical(classification_row(plain, "on_off", "A")$state, "absent")
})

test_that("no-rescue mode preserves representation, permutation, and random state", {
    fixture <- observation_policy_fixture()
    with_preserved_random_state({
        set.seed(312)
        before <- .Random.seed
        result <- classify_missingness(fixture$x, fixture$group,
            cutoffs = c(A = 12), rescue = "none")
        expect_identical(.Random.seed, before)
        reversed <- classify_missingness(fixture$x[5:1, 6:1],
            fixture$group[6:1], cutoffs = c(A = 12), rescue = "none")
        expect_identical(reversed$data, result$data[4:1, 6:1])
        expect_identical(reversed$cutoffs, result$cutoffs)
        se <- rich_summarized_experiment(fixture$x, fixture$group)
        fit_se <- classify_missingness(se, group_col = "condition",
            assay = "intensity", cutoffs = c(A = 12), rescue = "none")
        expect_identical(SummarizedExperiment::assay(fit_se$data, "intensity"),
            result$data)
        expect_identical(fit_se$classifications, result$classifications)
        expect_identical(.Random.seed, before)
    })
})

test_that("observation and ambiguity policies use explicit names", {
    fixture <- observation_policy_fixture()
    for (value in list(NA_character_, "", "unknown", c("none", "minimum"))) {
        expect_error(classify_missingness(fixture$x, fixture$group,
            rescue = value), "rescue")
    }
    for (value in list(NA_character_, "", "unknown", c("warn", "error"))) {
        expect_error(classify_missingness(fixture$x, fixture$group,
            ambiguity = value), "ambiguity")
    }
})

test_that("absence does not rescue unsupported whole-feature combinations", {
    x <- rbind(
        absent_mnar = c(NA, NA, NA, 8, NA, NA),
        absent_insufficient = c(NA, NA, NA, 15, NA, NA),
        retained = c(NA, NA, NA, 15, 16, 17),
        anchor = c(8, 9, 10, 15, 16, 17)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))
    result <- classify_missingness(x, rep(c("A", "B"), each = 3),
        cutoffs = c(B = 12), rescue = "none")
    expect_identical(rownames(result$data), c("retained", "anchor"))
    expect_identical(result$feature_status$drop_reason[1:2],
        c("MNAR_or_absent_all_conditions", "insufficient:B"))
    expect_identical(result$groups$A$absent, "retained")
    expect_false(any(c("absent_mnar", "absent_insufficient") %in%
        unlist(result$groups, use.names = FALSE)))
    expect_true("MNAR_or_absent_all_conditions" %in% summary(result)$drops$reason)
    expect_output(print(summary(result)), "MNAR or absent in every condition")
})
