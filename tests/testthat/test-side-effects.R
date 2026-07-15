.capture_process_state <- function() {
    list(
        working_directory = getwd(),
        options = options(),
        rng_kind = RNGkind(),
        random_seed = get(
            ".Random.seed",
            envir = globalenv(),
            inherits = FALSE
        ),
        current_device = grDevices::dev.cur(),
        device_list = grDevices::dev.list(),
        graphics_parameters = graphics::par(no.readonly = TRUE)
    )
}

.expect_only_logged_output_changes <- function(input, result) {
    output <- result$data
    original <- input[
        rownames(output),
        colnames(output),
        drop = FALSE
    ]
    observed <- !is.na(original)
    expect_identical(output[observed], original[observed])

    changed <- (is.na(original) != is.na(output)) |
        (!is.na(original) & !is.na(output) & original != output)
    expected <- matrix(
        FALSE,
        nrow = nrow(output),
        ncol = ncol(output),
        dimnames = dimnames(output)
    )
    retained_log <- result$seed_log[
        result$seed_log$feature %in% rownames(output),
        ,
        drop = FALSE
    ]
    indices <- cbind(
        match(retained_log$feature, rownames(output)),
        match(retained_log$sample, colnames(output))
    )
    expected[indices] <- TRUE

    expect_identical(changed, expected)
    expect_true(all(is.na(retained_log$old_value)))
    expect_identical(
        unname(output[indices]),
        unname(retained_log$inserted_value)
    )
}

test_that("public classification leaves process-global state unchanged", {
    simulation <- scientific_routine_fixture("automatic_cliff")

    with_preserved_random_state({
        RNGkind("Knuth-TAOCP-2002", "Box-Muller", "Rejection")
        set.seed(90210L)
        grDevices::pdf(NULL)
        device <- grDevices::dev.cur()
        on.exit(grDevices::dev.off(device), add = TRUE)
        graphics::par(mar = c(3.1, 3.2, 3.3, 3.4))
        before <- .capture_process_state()

        manual <- classify_missingness(
            simulation$data,
            simulation$groups,
            cutoffs = c(A = 12, B = 12),
            seed = 17L
        )
        automatic <- classify_missingness(
            simulation$data,
            simulation$groups,
            seed = 17L
        )
        stored_plot <- plot_missingness(automatic, "A")
        after <- .capture_process_state()

        expect_identical(after, before)
        expect_s3_class(manual, "imputefinder_result")
        expect_s3_class(automatic, "imputefinder_result")
        expect_s3_class(stored_plot, "ggplot")
    })
})

test_that("final matrix changes are exactly the retained logged seeds", {
    simulation <- scientific_routine_fixture("automatic_cliff")
    original <- simulation$data
    manual <- classify_missingness(
        simulation$data,
        simulation$groups,
        cutoffs = c(A = 12, B = 12),
        seed = 29L
    )
    automatic <- classify_missingness(
        simulation$data,
        simulation$groups,
        seed = 29L
    )

    expect_identical(simulation$data, original)
    .expect_only_logged_output_changes(original, manual)
    .expect_only_logged_output_changes(original, automatic)
})
