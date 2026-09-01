#!/usr/bin/env Rscript

# M10 release-candidate performance gate. Run after installing the current
# source tree, from the repository root:
# Rscript --vanilla dev/performance-validation.R

.performance_specification <- list(
    feature_count = 10000L,
    condition_count = 5L,
    samples_per_condition = 10L,
    timing_repetitions = 3L,
    maximum_median_seconds = 10,
    maximum_total_allocation_ratio = 100,
    maximum_largest_allocation_ratio = 2,
    maximum_result_size_ratio = 12,
    maximum_peak_rss_mib = 768,
    maximum_automatic_cutoff_error = 1
)

.performance_condition_block <- function(
    means,
    samples,
    right_boundary,
    score_offset
) {
    index <- seq_along(means)
    score <- (((index + score_offset) * 104729) %% 1000003) / 1000003
    ramp <- pmin(
        pmax((right_boundary - means) / 0.8, 0),
        1
    )
    missing_probability <- 1 - (1 - 0.05) * (1 - 0.68 * ramp)
    missing <- score < missing_probability
    block <- matrix(rep(means, samples), nrow = length(means))
    missing_sample <- 1L + ((index[missing] + score_offset) %% samples)
    block[cbind(index[missing], missing_sample)] <- NA_real_
    block
}

.performance_fixture <- function(specification) {
    n <- specification$feature_count
    samples <- specification$samples_per_condition
    index <- seq_len(n)
    features <- sprintf("protein_%05d", index)
    means <- stats::qnorm(
        (index - 0.5) / n,
        mean = 12.5,
        sd = 2.4
    )
    conditions <- LETTERS[seq_len(specification$condition_count)]
    boundaries <- stats::setNames(
        seq(11, by = 0.5, length.out = length(conditions)),
        conditions
    )
    blocks <- lapply(
        seq_along(conditions),
        function(condition_index) {
            .performance_condition_block(
                means,
                samples,
                boundaries[[condition_index]],
                score_offset = 37L * (condition_index - 1L)
            )
        }
    )
    x <- do.call(cbind, blocks)
    rownames(x) <- features
    colnames(x) <- unlist(
        lapply(
            conditions,
            function(condition) {
                sprintf("%s_%02d", condition, seq_len(samples))
            }
        ),
        use.names = FALSE
    )
    groups <- stats::setNames(
        rep(conditions, each = samples),
        colnames(x)
    )

    off_count <- 40L
    for (condition_index in seq_along(conditions)) {
        rows <- seq.int(
            from = 1L + (condition_index - 1L) * off_count,
            length.out = off_count
        )
        columns <- groups == conditions[[condition_index]]
        x[rows, columns] <- NA_real_
    }

    list(x = x, groups = groups, cutoffs = boundaries)
}

.profiled_call <- function(code) {
    path <- tempfile("imputefinder-rprofmem-", fileext = ".log")
    on.exit(
        {
            try(utils::Rprofmem(NULL), silent = TRUE)
            unlink(path)
        },
        add = TRUE
    )
    utils::Rprofmem(path, threshold = 0)
    timing <- system.time(value <- force(code))
    utils::Rprofmem(NULL)
    lines <- readLines(path, warn = FALSE)
    allocation <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
    allocation <- allocation[is.finite(allocation)]

    list(
        value = value,
        elapsed_seconds = unname(timing[["elapsed"]]),
        allocation_count = length(allocation),
        allocation_bytes = sum(allocation),
        largest_allocation_bytes = max(allocation, 0)
    )
}

.retained_changes_are_logged <- function(input, result) {
    output <- result$data
    original <- input[rownames(output), colnames(output), drop = FALSE]
    changed <- (is.na(original) != is.na(output)) |
        (!is.na(original) & !is.na(output) & original != output)
    expected <- matrix(
        FALSE,
        nrow(output),
        ncol(output),
        dimnames = dimnames(output)
    )
    log <- result$seed_log[
        result$seed_log$feature %in% rownames(output),
        ,
        drop = FALSE
    ]
    indices <- cbind(
        match(log$feature, rownames(output)),
        match(log$sample, colnames(output))
    )
    expected[indices] <- TRUE

    identical(changed, expected) &&
        identical(output[!is.na(original)], original[!is.na(original)]) &&
        all(is.na(log$old_value)) &&
        identical(unname(output[indices]), unname(log$inserted_value))
}

.timed_classification <- function(fixture, cutoffs, seed) {
    timing <- system.time({
        result <- classify_missingness(
            fixture$x,
            fixture$groups,
            cutoffs = cutoffs,
            seed = seed
        )
    })
    stopifnot(.retained_changes_are_logged(fixture$x, result))
    unname(timing[["elapsed"]])
}

.benchmark_path <- function(name, fixture, cutoffs, specification) {
    timings <- vapply(
        seq_len(specification$timing_repetitions),
        function(repetition) {
            gc()
            .timed_classification(
                fixture,
                cutoffs,
                seed = 100L + repetition
            )
        },
        numeric(1L)
    )
    gc()
    profiled <- .profiled_call(classify_missingness(
        fixture$x,
        fixture$groups,
        cutoffs = cutoffs,
        seed = 101L
    ))
    result <- profiled$value
    stopifnot(.retained_changes_are_logged(fixture$x, result))
    input_bytes <- as.numeric(object.size(fixture$x))

    data.frame(
        path = name,
        median_seconds = stats::median(timings),
        minimum_seconds = min(timings),
        maximum_seconds = max(timings),
        profiled_seconds = profiled$elapsed_seconds,
        allocation_count = profiled$allocation_count,
        allocation_mib = profiled$allocation_bytes / 1024^2,
        allocation_ratio = profiled$allocation_bytes / input_bytes,
        largest_allocation_mib = profiled$largest_allocation_bytes / 1024^2,
        largest_allocation_ratio =
            profiled$largest_allocation_bytes / input_bytes,
        result_mib = as.numeric(object.size(result)) / 1024^2,
        result_size_ratio = as.numeric(object.size(result)) / input_bytes,
        retained_features = nrow(result$data),
        seed_insertions = nrow(result$seed_log),
        maximum_cutoff_error = if (identical(name, "automatic")) {
            max(abs(result$cutoffs - fixture$cutoffs))
        } else {
            0
        },
        stringsAsFactors = FALSE
    )
}

.peak_rss_mib <- function() {
    status <- "/proc/self/status"
    if (!file.exists(status)) {
        return(NA_real_)
    }
    line <- grep("^VmHWM:", readLines(status, warn = FALSE), value = TRUE)
    if (length(line) != 1L) {
        return(NA_real_)
    }
    kib <- suppressWarnings(as.numeric(gsub("[^0-9]", "", line)))
    kib / 1024
}

.assess_performance <- function(results, peak_rss, specification) {
    checks <- c(
        median_seconds = all(
            results$median_seconds <= specification$maximum_median_seconds
        ),
        total_allocation = all(
            results$allocation_ratio <=
                specification$maximum_total_allocation_ratio
        ),
        largest_allocation = all(
            results$largest_allocation_ratio <=
                specification$maximum_largest_allocation_ratio
        ),
        result_size = all(
            results$result_size_ratio <=
                specification$maximum_result_size_ratio
        ),
        automatic_cutoff = all(
            results$maximum_cutoff_error <=
                specification$maximum_automatic_cutoff_error
        ),
        peak_rss = is.na(peak_rss) ||
            peak_rss <= specification$maximum_peak_rss_mib
    )
    if (!all(checks)) {
        stop(
            "Performance gate failures: ",
            paste(names(checks)[!checks], collapse = ", "),
            call. = FALSE
        )
    }
    checks
}

if (!requireNamespace("imputefinder", quietly = TRUE)) {
    stop("Install the current imputefinder source tree first.", call. = FALSE)
}
library(imputefinder)

specification <- .performance_specification
fixture <- .performance_fixture(specification)
original <- fixture$x
manual <- .benchmark_path("manual", fixture, fixture$cutoffs, specification)
automatic <- .benchmark_path("automatic", fixture, NULL, specification)
results <- rbind(manual, automatic)
peak_rss <- .peak_rss_mib()
checks <- .assess_performance(results, peak_rss, specification)
stopifnot(identical(fixture$x, original))

print(results, row.names = FALSE, digits = 6)
cat(sprintf("peak_rss_mib: %.3f\n", peak_rss))
cat("gates: ", paste(names(checks), checks, sep = "=", collapse = "; "), "\n", sep = "")
