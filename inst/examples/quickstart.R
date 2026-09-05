# Run from a writable folder. Re-running replaces these example output files.
# All inputs are synthetic. Manual cutoffs explain the rules, not a fitted model.
if (!requireNamespace("imputefinder", quietly = TRUE)) {
    stop("Install imputefinder before running this example.", call. = FALSE)
}

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
condition <- rep(c("A", "B"), each = 4)
fit <- imputefinder::classify_missingness(
    x, group = condition, cutoffs = c(A = 12, B = 12),
    rescue = "minimum", seed = 1L
)

# Include globally absent features, which do not enter classification.
states <- expand.grid(feature = rownames(x), column = c("A", "B"),
    stringsAsFactors = FALSE)
key <- paste(states$feature, states$column, sep = ":")
classified_key <- paste(fit$classifications$feature,
    fit$classifications$condition, sep = ":")
matched <- match(key, classified_key)
states$state <- fit$classifications$state[matched]
states$state[is.na(states$state)] <- "Not classified"
states$label <- states$state
rescued <- !is.na(matched) & fit$classifications$seeded[matched]
states$label[which(rescued)] <- paste0(states$label[which(rescued)], "\n(rescued)")
decisions <- data.frame(
    feature = fit$feature_status$feature,
    column = "Decision",
    state = ifelse(fit$feature_status$retained, "Keep", "Drop"),
    stringsAsFactors = FALSE
)
decisions$label <- decisions$state
state_table <- rbind(states, decisions)
state_table$feature <- factor(state_table$feature, levels = rev(rownames(x)))
state_table$column <- factor(state_table$column, levels = c("A", "B", "Decision"))
state_plot <- ggplot2::ggplot(state_table,
    ggplot2::aes(x = column, y = feature, fill = state)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.8) +
    ggplot2::scale_fill_manual(values = c(
        complete = "#D7EAE7", MAR = "#C9E1F0", MNAR = "#F7D8B0",
        insufficient = "#EBC9C9", "Not classified" = "#EEEEEE",
        Keep = "#D7EAE7", Drop = "#EBC9C9")) +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
        title = "Separate condition states, one retention decision",
        subtitle = "Synthetic log2 intensities; four samples per condition; cutoffs A = B = 12",
        x = NULL, y = NULL,
        caption = paste(
            "Rescue inserts one artificial value (8) for on_off in condition A.",
            "MAR and MNAR are working labels, not established causes. No imputation is performed.",
            sep = "\n")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        legend.position = "none",
        axis.text.x = ggplot2::element_text(face = "bold"),
        plot.title.position = "plot",
        plot.caption = ggplot2::element_text(hjust = 0)
    )

output_dir <- "imputefinder-results"
if (!dir.exists(output_dir) && !dir.create(output_dir)) {
    stop("Cannot create imputefinder-results in the working folder.", call. = FALSE)
}
utils::write.csv(fit$classifications,
    file.path(output_dir, "classifications.csv"), row.names = FALSE)
utils::write.csv(fit$feature_status,
    file.path(output_dir, "feature-status.csv"), row.names = FALSE)
utils::write.csv(data.frame(feature = rownames(fit$data), fit$data,
    check.names = FALSE), file.path(output_dir, "retained-intensities.csv"),
    row.names = FALSE)
utils::write.csv(fit$seed_log, file.path(output_dir, "seed-log.csv"),
    row.names = FALSE)
saveRDS(fit, file.path(output_dir, "result.rds"))
ggplot2::ggsave(file.path(output_dir, "classifications.png"), state_plot,
    width = 9, height = 5.8, units = "in", dpi = 180, bg = "white")
print(fit)
print(fit$feature_status, row.names = FALSE)
message("Example files: ", normalizePath(output_dir))
if (interactive()) {
    print(state_plot)
}
