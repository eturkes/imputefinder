#!/usr/bin/env Rscript

# M11f stable/sidecar differential + overhead gate. Run after installing the
# current source tree, from the repository root:
# Rscript --vanilla dev/sidecar-differential-validation.R

.sidecar_specification <- list(
    timing_repetitions = 3L,
    maximum_median_seconds = 30,
    maximum_elapsed_overhead_ratio = 8,
    maximum_total_allocation_ratio = 250,
    maximum_allocation_overhead_ratio = 4,
    maximum_largest_allocation_ratio = 2,
    maximum_result_size_ratio = 8,
    maximum_sidecar_to_classic_size_ratio = 1.5
)

.sidecar_timed_call <- function(fixture, design, base_fit, seed) {
    gc()
    timing <- system.time({
        analysis <- analyze_missingness(
            fixture$x,
            design,
            base_fit = base_fit,
            seed = seed
        )
    })
    stopifnot(identical(analysis$classic, base_fit))
    unname(timing[["elapsed"]])
}

.sidecar_comparator_gate <- function(reference) {
    compare <- getFromNamespace(".compare_classic_fits", "imputefinder")
    tolerance <- getFromNamespace(
        ".classic_compatibility_tolerance",
        "imputefinder"
    )()[["absolute"]]
    exact <- reference
    exact$classifications$state[[1L]] <- "MAR"
    canonical <- reference
    canonical$profiles[[1L]]$metadata$warnings <- "fixture warning"
    within <- reference
    within$profiles[[1L]]$grid$intensity[[1L]] <-
        within$profiles[[1L]]$grid$intensity[[1L]] + tolerance / 4
    beyond <- reference
    beyond$profiles[[1L]]$grid$intensity[[1L]] <-
        beyond$profiles[[1L]]$grid$intensity[[1L]] + tolerance * 100
    reports <- lapply(
        list(exact, canonical, within, beyond),
        compare,
        reference = reference
    )

    c(
        exact = identical(reports[[1L]]$exact, "classifications"),
        canonical = identical(reports[[2L]]$canonical, "profiles"),
        tolerance_within = isTRUE(reports[[3L]]$compatible),
        tolerance_beyond = identical(reports[[4L]]$tolerance, "profiles")
    )
}

.assess_sidecar <- function(
    result,
    comparator,
    original,
    specification
) {
    checks <- c(
        exact_classic = isTRUE(result$exact_classic),
        exact_input = identical(original, result$input),
        exact_mask_size = isTRUE(result$exact_mask_size),
        comparator,
        median_seconds =
            result$median_seconds <= specification$maximum_median_seconds,
        elapsed_overhead =
            result$elapsed_overhead_ratio <=
                specification$maximum_elapsed_overhead_ratio,
        total_allocation =
            result$allocation_ratio <=
                specification$maximum_total_allocation_ratio,
        allocation_overhead =
            result$allocation_overhead_ratio <=
                specification$maximum_allocation_overhead_ratio,
        largest_allocation =
            result$largest_allocation_ratio <=
                specification$maximum_largest_allocation_ratio,
        result_size =
            result$result_size_ratio <=
                specification$maximum_result_size_ratio,
        sidecar_size =
            result$sidecar_to_classic_size_ratio <=
                specification$maximum_sidecar_to_classic_size_ratio
    )
    if (!all(checks)) {
        stop(
            "Sidecar differential gate failures: ",
            paste(names(checks)[!checks], collapse = ", "),
            call. = FALSE
        )
    }
    checks
}

source("dev/performance-validation.R", local = TRUE)

sidecar_specification <- .sidecar_specification
sidecar_original <- fixture$x
sidecar_design <- missingness_design(
    data.frame(
        condition = unname(fixture$groups),
        row.names = colnames(fixture$x)
    ),
    condition = "condition"
)
sidecar_fit <- classify_missingness(
    fixture$x,
    fixture$groups,
    cutoffs = fixture$cutoffs,
    seed = 101L
)
sidecar_timings <- vapply(
    seq_len(sidecar_specification$timing_repetitions),
    function(repetition) {
        .sidecar_timed_call(fixture, sidecar_design, sidecar_fit, 101L)
    },
    numeric(1L)
)
gc()
sidecar_profiled <- .profiled_call(analyze_missingness(
    fixture$x,
    sidecar_design,
    base_fit = sidecar_fit,
    seed = 101L
))
sidecar <- sidecar_profiled$value
input_bytes <- as.numeric(object.size(fixture$x))
classic_bytes <- as.numeric(object.size(sidecar_fit))
stable_allocation_bytes <- manual$allocation_mib[[1L]] * 1024^2
sidecar_result <- list(
    input = fixture$x,
    exact_classic = identical(sidecar$classic, sidecar_fit),
    exact_mask_size = length(sidecar$input$original_mask$bytes) ==
        ceiling(length(fixture$x) / 8),
    median_seconds = stats::median(sidecar_timings),
    minimum_seconds = min(sidecar_timings),
    maximum_seconds = max(sidecar_timings),
    elapsed_overhead_ratio = stats::median(sidecar_timings) /
        manual$median_seconds[[1L]],
    profiled_seconds = sidecar_profiled$elapsed_seconds,
    allocation_mib = sidecar_profiled$allocation_bytes / 1024^2,
    allocation_ratio = sidecar_profiled$allocation_bytes / input_bytes,
    allocation_overhead_ratio = sidecar_profiled$allocation_bytes /
        stable_allocation_bytes,
    largest_allocation_mib =
        sidecar_profiled$largest_allocation_bytes / 1024^2,
    largest_allocation_ratio =
        sidecar_profiled$largest_allocation_bytes / input_bytes,
    result_mib = as.numeric(object.size(sidecar)) / 1024^2,
    result_size_ratio = as.numeric(object.size(sidecar)) / input_bytes,
    sidecar_to_classic_size_ratio =
        as.numeric(object.size(sidecar)) / classic_bytes
)
comparator_checks <- .sidecar_comparator_gate(sidecar_fit)
sidecar_checks <- .assess_sidecar(
    sidecar_result,
    comparator_checks,
    sidecar_original,
    sidecar_specification
)

print(as.data.frame(sidecar_result[!names(sidecar_result) %in% "input"]))
cat(
    "gates: ",
    paste(names(sidecar_checks), sidecar_checks, sep = "=", collapse = "; "),
    "\n",
    sep = ""
)
