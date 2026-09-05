.run_quickstart_example <- function() {
    script <- system.file("examples", "quickstart.R", package = "imputefinder",
        mustWork = TRUE)
    work <- tempfile("imputefinder-example-")
    dir.create(work)
    previous <- setwd(work)
    on.exit({
        setwd(previous)
        unlink(work, recursive = TRUE)
    }, add = TRUE)
    example <- new.env(parent = globalenv())
    capture.output(suppressMessages(sys.source(script, envir = example)))
    files <- list.files(example$output_dir, full.names = TRUE)
    checksums <- tools::md5sum(files)
    capture.output(suppressMessages(sys.source(script, envir = example)))
    expect_identical(tools::md5sum(files), checksums)
    expect_identical(list.files(work), "imputefinder-results")
    tables <- lapply(files[grepl("\\.csv$", files)], utils::read.csv,
        check.names = FALSE)
    names(tables) <- basename(files[grepl("\\.csv$", files)])
    list(
        example = example, files = basename(files), tables = tables,
        saved = readRDS(file.path(example$output_dir, "result.rds")),
        png = readBin(file.path(example$output_dir, "classifications.png"),
            what = "raw", n = 24L)
    )
}

test_that("the installed quickstart exports classifications and rescue provenance", {
    run <- .run_quickstart_example()
    fit <- run$example$fit
    expect_setequal(run$files, c("classifications.csv", "feature-status.csv",
        "retained-intensities.csv", "seed-log.csv", "result.rds",
        "classifications.png"))
    expect_identical(rownames(fit$data),
        c("on_off", "mar_both", "sparse_mnar", "complete_low"))
    expect_identical(fit$feature_status$retained,
        c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
    expect_identical(fit$feature_status$drop_reason,
        c(NA_character_, NA_character_, "insufficient:A", NA_character_,
          "MNAR_all_conditions", "all_missing", NA_character_))
    expect_equal(run$tables[["feature-status.csv"]], fit$feature_status)
    expect_equal(run$tables[["classifications.csv"]], fit$classifications)
    seed_log <- run$tables[["seed-log.csv"]]
    expect_true(all(is.na(seed_log$old_value)))
    # CSV cannot preserve a numeric type for an entirely missing column.
    seed_log$old_value <- as.numeric(seed_log$old_value)
    expect_equal(seed_log, fit$seed_log)
    expect_identical(run$tables[["retained-intensities.csv"]]$feature,
        rownames(fit$data))
    expect_identical(names(run$tables[["retained-intensities.csv"]])[-1L],
        colnames(fit$data))
    expect_equal(unname(as.matrix(run$tables[["retained-intensities.csv"]][-1L])),
        unname(fit$data))
    expect_identical(run$saved, fit)
    expect_equal(fit$cutoffs, c(A = 12, B = 12))
    expect_equal(nrow(fit$seed_log), 1L)
    expect_identical(fit$seed_log$feature, "on_off")
    expect_equal(fit$seed_log$inserted_value, 8)
    block <- subset(fit$classifications, feature == "on_off" & condition == "A")
    expect_identical(block$state, "MNAR")
    expect_equal(block$original_observed_count, 0L)
    expect_equal(block$observed_count, 1L)
    original <- run$example$x[rownames(fit$data), , drop = FALSE]
    expect_equal(sum(is.na(original) & !is.na(fit$data)), 1L)
    expect_identical(fit$data[!is.na(original)], original[!is.na(original)])
    expect_equal(nrow(run$example$state_plot$data), 21L)
    expect_equal(sum(run$example$state_table$label == "MNAR\n(rescued)"), 1L)
    expect_equal(sum(run$example$state_table$state == "Not classified"), 2L)
    expect_identical(run$png[1:8], as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)))
    expect_equal(sum(as.integer(run$png[17:20]) * 256^(3:0)), 1620)
    expect_equal(sum(as.integer(run$png[21:24]) * 256^(3:0)), 1044)
})
