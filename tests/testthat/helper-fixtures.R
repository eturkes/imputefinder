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
