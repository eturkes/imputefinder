normative_fixture <- function() {
    x <- rbind(
        on_off = c(NA, NA, NA, NA, 15, 16, 15, 17),
        mar_both = c(14, 15, NA, 16, 14, NA, 15, 16),
        sparse_mar = c(20, NA, NA, NA, 20, 21, 22, 23),
        sparse_mnar = c(8, NA, NA, NA, 14, 15, 16, 17),
        all_mnar = c(8, NA, NA, NA, 9, NA, NA, NA),
        globally_absent = rep(NA, 8),
        complete_low = c(8, 8, 8, 8, 9, 9, 9, 9)
    )
    colnames(x) <- paste0("s", seq_len(ncol(x)))

    list(
        x = x,
        group = rep(c("A", "B"), each = 4),
        cutoffs = c(A = 12, B = 12)
    )
}

classify_normative_fixture <- function() {
    fixture <- normative_fixture()
    classify_missingness(
        x = fixture$x,
        group = fixture$group,
        cutoffs = fixture$cutoffs
    )
}

classification_row <- function(result, feature, condition) {
    result$classifications[
        result$classifications$feature == feature &
            result$classifications$condition == condition,
        ,
        drop = FALSE
    ]
}

matrix_input_fixture <- function() {
    x <- matrix(
        c(0, -2, NA, 3, 4, 5),
        nrow = 2,
        dimnames = list(c("protein_b", "protein_a"), c("s2", "s1", "s3"))
    )

    list(x = x, group = c("B", "A", "B"))
}

prepare_matrix_input <- function(x, group, group_col = NULL, assay = NULL) {
    imputefinder:::.prepare_matrix_input(
        x = x,
        group = group,
        group_col = group_col,
        assay = assay
    )
}

summarized_experiment_input_fixture <- function(multiple_assays = FALSE) {
    fixture <- matrix_input_fixture()
    assays <- list(intensity = fixture$x)
    if (multiple_assays) {
        assays$auxiliary <- fixture$x + 100
    }

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        rowData = data.frame(
            annotation = c("second", "first"),
            row.names = rownames(fixture$x)
        ),
        colData = data.frame(
            batch = c(2L, 1L, 2L),
            condition = factor(fixture$group),
            row.names = colnames(fixture$x)
        ),
        metadata = list(source = "adapter fixture")
    )

    list(se = se, x = fixture$x, group = fixture$group)
}

prepare_input <- function(x, group = NULL, group_col = NULL, assay = NULL) {
    imputefinder:::.prepare_input(
        x = x,
        group = group,
        group_col = group_col,
        assay = assay
    )
}

restore_output_data <- function(data, prepared, original) {
    imputefinder:::.restore_output_data(
        data = data,
        prepared = prepared,
        original = original
    )
}

automatic_cutoff_statistics <- function(
    n = 800L,
    right = 12,
    left = right - 0.6,
    height = 0.72,
    mar_rate = 0.05,
    center = 13,
    scale = 2.2,
    score_offset = 0L,
    missing = NULL,
    condition = "A"
) {
    index <- seq_len(n)
    mean_intensity <- stats::qnorm(
        (index - 0.5) / n,
        mean = center,
        sd = scale
    )
    if (is.null(missing)) {
        score <- (((index + score_offset) * 104729) %% 1000003) / 1000003
        ramp <- pmin(pmax((right - mean_intensity) / (right - left), 0), 1)
        missing_probability <- 1 - (1 - mar_rate) * (1 - height * ramp)
        missing <- score < missing_probability
    }

    data.frame(
        feature = sprintf("feature_%05d", index),
        condition = rep(condition, n),
        sample_count = rep(8L, n),
        observed_count = ifelse(missing, 7L, 8L),
        missing_count = ifelse(missing, 1L, 0L),
        missing_fraction = ifelse(missing, 1 / 8, 0),
        mean_intensity = mean_intensity,
        seeded = rep(FALSE, n),
        stringsAsFactors = FALSE
    )
}

automatic_cutoff_profile <- function(...) {
    statistics <- automatic_cutoff_statistics(...)
    imputefinder:::.condition_missingness_profile(statistics)
}

automatic_cutoff_matrix_fixture <- function() {
    a <- automatic_cutoff_statistics(
        n = 1200L,
        right = 11,
        left = 10.2,
        height = 0.68,
        center = 12.5,
        scale = 2.4
    )
    b <- automatic_cutoff_statistics(
        n = 1200L,
        right = 14,
        left = 13.2,
        height = 0.68,
        center = 12.5,
        scale = 2.4,
        score_offset = 37L,
        condition = "B"
    )
    condition_block <- function(statistics) {
        block <- matrix(
            rep(statistics$mean_intensity, 4L),
            nrow = nrow(statistics)
        )
        missing <- statistics$missing_count > 0L
        block[cbind(which(missing), rep(4L, sum(missing)))] <- NA_real_
        block
    }
    x <- cbind(condition_block(a), condition_block(b))
    rownames(x) <- a$feature
    colnames(x) <- paste0("sample_", seq_len(ncol(x)))

    list(
        x = x,
        group = rep(c("A", "B"), each = 4L),
        expected = c(A = 10.698082335694377, B = 13.815875319932115)
    )
}

rich_summarized_experiment <- function(x, group) {
    reference <- matrix(
        seq_len(length(x)),
        nrow = nrow(x),
        dimnames = dimnames(x)
    )
    feature_index <- match(rownames(x), sort(rownames(x), method = "radix"))
    sample_index <- match(colnames(x), sort(colnames(x), method = "radix"))

    SummarizedExperiment::SummarizedExperiment(
        assays = list(reference = reference, intensity = x),
        rowData = data.frame(
            feature_index = feature_index,
            annotation = paste0("annotation_", feature_index),
            row.names = rownames(x)
        ),
        colData = data.frame(
            sample_index = sample_index,
            condition = factor(group, levels = rev(sort(unique(group)))),
            row.names = colnames(x)
        ),
        metadata = list(
            source = "M6b representation fixture",
            nested = list(preserved = TRUE)
        )
    )
}

post_seed_statistics <- function(x, group, seed = 1L) {
    prepared <- prepare_matrix_input(x, group)
    rescued <- imputefinder:::.seed_missing_conditions(
        prepared$data,
        prepared$groups_by_sample,
        seed
    )

    list(
        rescued = rescued,
        statistics = imputefinder:::.feature_condition_statistics(
            rescued$data,
            prepared$groups_by_sample,
            rescued$seed_log
        )
    )
}

missingness_profiles <- function(statistics) {
    imputefinder:::.build_missingness_profiles(statistics)
}
