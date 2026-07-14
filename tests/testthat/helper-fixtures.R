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
