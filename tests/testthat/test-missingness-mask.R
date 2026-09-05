test_that("masks partition remaining missing cells without selecting observations", {
    fit <- classify_fixture()
    before <- serialize(fit, NULL)
    mnar <- matrix(FALSE, nrow(fit$data), ncol(fit$data),
        dimnames = dimnames(fit$data))
    mnar[c("on_off", "sparse_mnar"), 1:4] <-
        is.na(fit$data[c("on_off", "sparse_mnar"), 1:4])
    mar <- mnar & FALSE
    mar["mar_both", c("s3", "s6")] <- TRUE
    expect_identical(missingness_mask(fit, "MNAR"), mnar)
    expect_identical(missingness_mask(fit, "MAR"), mar)
    expect_identical(missingness_mask(fit, "absent"), mnar & FALSE)
    expect_identical(missingness_mask(fit), is.na(fit$data))
    expect_identical(mnar | mar, is.na(fit$data))
    expect_false(any(mnar & mar))
    expect_false(any(missingness_mask(fit)[!is.na(fit$data)]))
    expect_false(anyNA(missingness_mask(fit)))
    expect_identical(missingness_mask(fit, character()), mnar & FALSE)
    expect_identical(missingness_mask(fit, c("MAR", "MAR")), mar)
    expect_identical(serialize(fit, NULL), before)
})

test_that("no-rescue masks keep absence separate from MAR and MNAR", {
    fixture <- classification_fixture()
    fit <- classify_missingness(fixture$x, fixture$group,
        cutoffs = fixture$cutoffs, rescue = "none")
    absent <- matrix(FALSE, nrow(fit$data), ncol(fit$data),
        dimnames = dimnames(fit$data))
    absent["on_off", 1:4] <- TRUE
    expect_identical(missingness_mask(fit, "absent"), absent)
    mechanisms <- missingness_mask(fit, c("MAR", "MNAR"))
    expect_false(any(mechanisms & absent))
    expect_identical(mechanisms | absent, is.na(fit$data))
    expect_identical(rownames(absent), rownames(fit$data))
})

test_that("masks follow reordered inputs and classification tables", {
    fixture <- classification_fixture()
    original <- classify_fixture()
    rows <- rev(seq_len(nrow(fixture$x)))
    columns <- c(8L, 1L, 6L, 3L, 4L, 5L, 2L, 7L)
    fit <- classify_missingness(fixture$x[rows, columns],
        fixture$group[columns], cutoffs = fixture$cutoffs)
    fit$classifications <- fit$classifications[
        rev(seq_len(nrow(fit$classifications))), ]
    for (state in c("MNAR", "MAR", "absent")) {
        expect_identical(missingness_mask(fit, state),
            missingness_mask(original, state)[
                rownames(fit$data), colnames(fit$data), drop = FALSE])
    }
})

test_that("experiment masks use the recorded assay, including after serialization", {
    fixture <- classification_fixture()
    auxiliary <- fixture$x
    auxiliary[] <- 100
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(auxiliary = auxiliary, intensity = fixture$x),
        colData = data.frame(condition = fixture$group,
            row.names = colnames(fixture$x)),
        metadata = list(source = "mask fixture")
    )
    chosen <- "intensity"
    fit <- classify_missingness(se, group_col = "condition", assay = chosen,
        cutoffs = fixture$cutoffs)
    before <- serialize(fit, NULL)
    expect_identical(fit$assay_index, 2L)
    expect_identical(missingness_mask(fit), missingness_mask(classify_fixture()))
    expect_identical(missingness_mask(unserialize(before)),
        missingness_mask(fit))
    expect_identical(serialize(fit, NULL), before)
    fit$assay_index <- NULL
    expect_error(missingness_mask(fit), "recorded assay")
    fit$data <- fit$data[, , drop = FALSE]
    SummarizedExperiment::assays(fit$data) <-
        SummarizedExperiment::assays(fit$data)[2L]
    expect_identical(missingness_mask(fit), missingness_mask(classify_fixture()))
})

test_that("empty outputs and one-sample rescue retain matrix dimensions", {
    x <- rbind(absent = rep(NA_real_, 3), also_absent = rep(NA_real_, 3))
    colnames(x) <- c("a", "b", "c")
    fit <- classify_missingness(x, c("A", "B", "B"), rescue = "none")
    expect_identical(missingness_mask(fit), is.na(fit$data))
    expect_identical(dim(missingness_mask(fit)), c(0L, 3L))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(intensity = x),
        colData = data.frame(condition = c("A", "B", "B"),
            row.names = colnames(x))
    )
    empty_se <- classify_missingness(se, group_col = "condition",
        rescue = "none")
    expect_identical(missingness_mask(empty_se), missingness_mask(fit))
    x <- rbind(on_off = c(NA, 14, 15), anchor = c(8, 8, 8))
    colnames(x) <- c("a", "b", "c")
    for (rescue in c("minimum", "none")) {
        fit <- classify_missingness(x, c("A", "B", "B"), rescue = rescue)
        expect_identical(missingness_mask(fit), is.na(fit$data))
        expect_identical(unname(missingness_mask(fit)["on_off", "a"]),
            rescue == "none")
    }
})

test_that("invalid selections and inconsistent result mappings are rejected", {
    fit <- classify_fixture()
    expect_error(missingness_mask(fit$data), "imputefinder_result")
    for (states in list(NULL, NA_character_, "M", "complete", "insufficient",
        1, TRUE, factor("MNAR"), matrix("MNAR"))) {
        expect_error(missingness_mask(fit, states), "`states`")
    }
    invalid <- fit
    invalid$groups_by_sample <- rev(invalid$groups_by_sample)
    expect_error(missingness_mask(invalid), "output sample order")
    invalid <- fit
    invalid$classifications <- invalid$classifications[-1L, ]
    expect_error(missingness_mask(invalid), "one classification")
    invalid <- fit
    invalid$classifications <- rbind(invalid$classifications,
        invalid$classifications[1L, ])
    expect_error(missingness_mask(invalid), "one classification")
    invalid <- fit
    invalid$classifications$state[1L] <- NA_character_
    expect_error(missingness_mask(invalid), "classification")
    invalid <- fit
    invalid$classifications$feature <- NULL
    expect_error(missingness_mask(invalid), "classifications")
    invalid <- fit
    invalid$classifications$state[1L] <- "unknown"
    expect_error(missingness_mask(invalid), "state vocabulary")
    fit$assay_index <- NULL
    expect_identical(missingness_mask(fit), is.na(fit$data))
})
