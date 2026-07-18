analysis_schema_fixture <- function(classic = NULL) {
    fixture <- normative_fixture()
    design <- missingness_design(
        data.frame(
            condition = fixture$group,
            row.names = colnames(fixture$x)
        ),
        condition = "condition"
    )
    input <- imputefinder:::.new_sidecar_input(fixture$x)
    if (is.null(classic)) {
        classic <- classify_normative_fixture()
    }

    imputefinder:::.new_imputefinder_analysis(
        classic = classic,
        design = design,
        input = input,
        call = quote(analyze_missingness(x, design))
    )
}

test_that("matrix fingerprints use a frozen canonical SHA-256 protocol", {
    x <- matrix(
        c(1, NA, -0, -2),
        nrow = 2L,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    original <- x

    fingerprint <- imputefinder:::.matrix_fingerprint(x)

    expect_identical(
        fingerprint,
        list(
            algorithm = "sha256",
            canonicalization = "matrix_core_be_v1",
            value = paste0(
                "a09c1f5501caee9b7d9a6dcce7322e56",
                "c2028a19c103c1f3bba7764d6338a24a"
            )
        )
    )
    expect_match(fingerprint$value, "^[0-9a-f]{64}$")
    expect_identical(imputefinder:::.matrix_fingerprint(x), fingerprint)
    expect_identical(
        imputefinder:::.matrix_fingerprint(
            unserialize(serialize(x, NULL, version = 3L))
        ),
        fingerprint
    )
    expect_identical(x, original)

    wrong_mask <- imputefinder:::.pack_original_mask(x)
    wrong_mask$bytes[[1L]] <- xor(
        wrong_mask$bytes[[1L]],
        as.raw(1L)
    )
    expect_error(
        imputefinder:::.matrix_fingerprint(x, wrong_mask),
        class = "imputefinder_fingerprint_schema_error"
    )
})

test_that("fingerprints cover ordered names, storage, values, and mask", {
    x <- matrix(
        c(1, NA, 3, 4),
        nrow = 2L,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    reference <- imputefinder:::.matrix_fingerprint(x)$value
    candidates <- list(
        value = { y <- x; y[[1L]] <- 2; y },
        mask = { y <- x; y[[2L]] <- 2; y },
        feature_name = { y <- x; rownames(y)[[1L]] <- "other"; y },
        sample_name = { y <- x; colnames(y)[[1L]] <- "other"; y },
        feature_order = x[2:1, , drop = FALSE],
        sample_order = x[, 2:1, drop = FALSE],
        storage = matrix(
            as.integer(x),
            nrow = 2L,
            dimnames = dimnames(x)
        )
    )

    candidate_hashes <- vapply(
        candidates,
        function(candidate) {
            imputefinder:::.matrix_fingerprint(candidate)$value
        },
        character(1L)
    )
    expect_false(any(candidate_hashes == reference))
    expect_length(unique(candidate_hashes), length(candidate_hashes))
})

test_that("text identifiers canonicalize to UTF-8 and bytes stay exact", {
    utf8_name <- "caf\u00e9"
    latin1_name <- iconv(utf8_name, from = "UTF-8", to = "latin1")
    Encoding(latin1_name) <- "latin1"
    bytes_name <- latin1_name
    Encoding(bytes_name) <- "bytes"
    x <- matrix(
        1:2,
        nrow = 1L,
        dimnames = list("feature", c("s1", "s2"))
    )
    utf8 <- latin1 <- bytes <- x
    colnames(utf8)[[1L]] <- utf8_name
    colnames(latin1)[[1L]] <- latin1_name
    colnames(bytes)[[1L]] <- bytes_name

    expect_identical(
        imputefinder:::.matrix_fingerprint(utf8),
        imputefinder:::.matrix_fingerprint(latin1)
    )
    expect_false(identical(
        imputefinder:::.matrix_fingerprint(utf8),
        imputefinder:::.matrix_fingerprint(bytes)
    ))
})

test_that("original masks are compact exact column-major bitsets", {
    x <- matrix(
        c(NA, 2, 3, 4, NA, 6, 7, 8, NA),
        nrow = 3L,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:3))
    )
    original <- x

    packed <- imputefinder:::.pack_original_mask(x)
    unpacked <- imputefinder:::.unpack_original_mask(
        packed,
        dimensions = dim(x),
        feature_names = rownames(x),
        sample_names = colnames(x)
    )

    expect_identical(
        names(packed),
        c("encoding", "cells", "bytes")
    )
    expect_identical(packed$encoding, "base_packbits_lsb0_v1")
    expect_identical(packed$cells, 9)
    expect_type(packed$bytes, "raw")
    expect_length(packed$bytes, 2L)
    expect_identical(unpacked, is.na(x))
    expect_identical(x, original)

    bits <- as.logical(rawToBits(packed$bytes))
    expect_false(any(bits[10:16]))
})

test_that("mask validation rejects malformed lengths and padding", {
    x <- matrix(
        c(NA, 2, 3, 4, 5, 6, 7, 8, 9),
        nrow = 3L,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:3))
    )
    packed <- imputefinder:::.pack_original_mask(x)

    extra_byte <- packed
    extra_byte$bytes <- c(extra_byte$bytes, as.raw(0L))
    expect_error(
        imputefinder:::.unpack_original_mask(extra_byte, dim(x)),
        class = "imputefinder_mask_schema_error"
    )

    nonzero_padding <- packed
    nonzero_padding$bytes[[2L]] <- as.raw(128L)
    expect_error(
        imputefinder:::.unpack_original_mask(nonzero_padding, dim(x)),
        class = "imputefinder_mask_schema_error"
    )

    wrong_encoding <- packed
    wrong_encoding$encoding <- "base_packbits_msb0_v2"
    error <- tryCatch(
        imputefinder:::.unpack_original_mask(wrong_encoding, dim(x)),
        imputefinder_mask_lifecycle_error = identity
    )
    expect_s3_class(error, "imputefinder_mask_lifecycle_error")
    expect_identical(error$expected_encoding, "base_packbits_lsb0_v1")
    expect_identical(error$actual_encoding, "base_packbits_msb0_v2")
})

test_that("sidecar input records identity without retaining numeric values", {
    fixture <- normative_fixture()
    acquisition <- stats::setNames(
        rep("DIA", ncol(fixture$x)),
        colnames(fixture$x)
    )

    input <- imputefinder:::.new_sidecar_input(
        fixture$x,
        representation = "matrix",
        acquisition = acquisition
    )

    expect_identical(
        names(input),
        c(
            "dimensions",
            "feature_names",
            "sample_names",
            "storage",
            "representation",
            "assay",
            "fingerprint",
            "original_mask",
            "scale",
            "acquisition"
        )
    )
    expect_identical(
        input$dimensions,
        c(features = 7L, samples = 8L)
    )
    expect_identical(input$feature_names, rownames(fixture$x))
    expect_identical(input$sample_names, colnames(fixture$x))
    expect_identical(input$storage, "double")
    expect_identical(input$representation, "matrix")
    expect_identical(input$assay, NA_character_)
    expect_identical(input$scale, "unnormalised_log2")
    expect_identical(input$acquisition, acquisition)
    expect_false(any(vapply(input, is.matrix, logical(1L))))
    expect_lt(length(input$original_mask$bytes), length(fixture$x))

    se_input <- imputefinder:::.new_sidecar_input(
        fixture$x,
        representation = "SummarizedExperiment",
        assay = "intensity"
    )
    expect_identical(se_input$representation, "SummarizedExperiment")
    expect_identical(se_input$assay, "intensity")
    expect_invisible(imputefinder:::.validate_sidecar_input(se_input))

    noncanonical <- input
    noncanonical$acquisition <- factor(noncanonical$acquisition)
    names(noncanonical$acquisition) <- input$sample_names
    expect_error(
        imputefinder:::.validate_sidecar_input(noncanonical),
        class = "imputefinder_input_schema_error"
    )
})

test_that("sidecar input verification reports exact mismatch families", {
    fixture <- normative_fixture()
    input <- imputefinder:::.new_sidecar_input(fixture$x)

    expect_invisible(
        imputefinder:::.verify_sidecar_input(input, fixture$x)
    )

    changed_value <- fixture$x
    changed_value[[2L]] <- 10
    value_error <- tryCatch(
        imputefinder:::.verify_sidecar_input(input, changed_value),
        imputefinder_input_mismatch_error = identity
    )
    expect_s3_class(value_error, "imputefinder_input_mismatch_error")
    expect_identical(value_error$differences, "fingerprint")
    expect_identical(
        value_error$expected_fingerprint,
        input$fingerprint$value
    )
    expect_match(value_error$actual_fingerprint, "^[0-9a-f]{64}$")

    changed_mask <- fixture$x
    changed_mask[[2L]] <- NA_real_
    mask_error <- tryCatch(
        imputefinder:::.verify_sidecar_input(input, changed_mask),
        imputefinder_input_mismatch_error = identity
    )
    expect_identical(
        mask_error$differences,
        c("original_mask", "fingerprint")
    )

    reordered <- fixture$x[nrow(fixture$x):1, , drop = FALSE]
    order_error <- tryCatch(
        imputefinder:::.verify_sidecar_input(input, reordered),
        imputefinder_input_mismatch_error = identity
    )
    expect_identical(
        order_error$differences,
        c("feature_names", "original_mask", "fingerprint")
    )

    wrong_storage <- input
    wrong_storage$storage <- "integer"
    storage_error <- tryCatch(
        imputefinder:::.verify_sidecar_input(wrong_storage, fixture$x),
        imputefinder_input_mismatch_error = identity
    )
    expect_identical(storage_error$differences, "storage")
})

test_that("fingerprint protocol changes require an explicit schema path", {
    fixture <- normative_fixture()
    input <- imputefinder:::.new_sidecar_input(fixture$x)
    input$fingerprint$canonicalization <- "matrix_core_be_v2"

    error <- tryCatch(
        imputefinder:::.verify_sidecar_input(input, fixture$x),
        imputefinder_fingerprint_lifecycle_error = identity
    )

    expect_s3_class(error, "imputefinder_fingerprint_lifecycle_error")
    expect_identical(error$expected_algorithm, "sha256")
    expect_identical(error$actual_algorithm, "sha256")
    expect_identical(
        error$expected_canonicalization,
        "matrix_core_be_v1"
    )
    expect_identical(
        error$actual_canonicalization,
        "matrix_core_be_v2"
    )
})

test_that("the minimum analysis schema round-trips exactly", {
    analysis <- analysis_schema_fixture()
    round_trip <- unserialize(serialize(analysis, NULL, version = 3L))

    expect_s3_class(analysis, "imputefinder_analysis")
    expect_identical(
        names(analysis),
        c(
            "classic",
            "spec",
            "design",
            "input",
            "sentinel",
            "stability",
            "provenance"
        )
    )
    expect_s3_class(analysis$classic, "imputefinder_result")
    expect_identical(
        names(analysis$design),
        c("declared", "estimability", "unavailable_roles")
    )
    expect_s3_class(analysis$design$declared, "missingness_design")
    expect_null(analysis$design$estimability)
    expect_identical(
        analysis$design$unavailable_roles,
        c("nuisance", "block", "acquisition")
    )
    expect_null(analysis$sentinel)
    expect_null(analysis$stability)
    expect_identical(
        names(analysis$provenance),
        c(
            "call",
            "seeds",
            "hashes",
            "warnings",
            "failures",
            "assumptions",
            "training_scope"
        )
    )
    expect_identical(round_trip, analysis)
    expect_invisible(
        imputefinder:::.validate_imputefinder_analysis(round_trip)
    )
})

test_that("analysis lifecycle rejects silent upgrades", {
    analysis <- analysis_schema_fixture()

    wrong_schema <- analysis
    wrong_schema$spec$schema <- "imputefinder_analysis_v2"
    error <- tryCatch(
        imputefinder:::.validate_imputefinder_analysis(wrong_schema),
        imputefinder_analysis_lifecycle_error = identity
    )
    expect_s3_class(error, "imputefinder_analysis_lifecycle_error")
    expect_identical(error$expected_schema, "imputefinder_analysis_v1")
    expect_identical(error$actual_schema, "imputefinder_analysis_v2")

    wrong_lifecycle <- analysis
    wrong_lifecycle$spec$lifecycle <- "stable"
    expect_error(
        imputefinder:::.validate_imputefinder_analysis(wrong_lifecycle),
        class = "imputefinder_analysis_lifecycle_error"
    )

    malformed <- analysis
    malformed$input$original_mask$bytes <- raw()
    expect_error(
        imputefinder:::.validate_imputefinder_analysis(malformed),
        class = "imputefinder_mask_schema_error"
    )

    malformed_provenance <- analysis
    malformed_provenance$provenance$assumptions$acquisition <- "declared"
    expect_error(
        imputefinder:::.validate_imputefinder_analysis(
            malformed_provenance
        ),
        class = "imputefinder_analysis_schema_error"
    )
})

test_that("classic failures are a portable alternative to full results", {
    failure <- imputefinder:::.new_classic_failure(
        simpleError("automatic cutoff evidence is ambiguous"),
        stage = "cutoff",
        call = quote(classify_missingness(x, group))
    )
    analysis <- analysis_schema_fixture(classic = failure)
    round_trip <- unserialize(serialize(analysis, NULL, version = 3L))

    expect_s3_class(failure, "imputefinder_classic_failure")
    expect_identical(
        names(failure),
        c("stage", "class", "message", "call")
    )
    expect_identical(failure$stage, "cutoff")
    expect_true(all(c("simpleError", "error", "condition") %in% failure$class))
    expect_identical(
        failure$message,
        "automatic cutoff evidence is ambiguous"
    )
    expect_identical(round_trip, analysis)
    expect_invisible(
        imputefinder:::.validate_imputefinder_analysis(round_trip)
    )

    partial_result <- classify_normative_fixture()
    partial_result$profiles <- NULL
    expect_error(
        analysis_schema_fixture(classic = partial_result),
        class = "imputefinder_analysis_schema_error"
    )
})

test_that("analysis-schema work leaves rules_v1 byte-exact", {
    fixture <- normative_fixture()
    before <- classify_normative_fixture()
    before_bytes <- serialize(before, NULL, version = 3L)
    classifier_formals <- formals(classify_missingness)

    analysis_schema_fixture()
    imputefinder:::.matrix_fingerprint(fixture$x)
    imputefinder:::.pack_original_mask(fixture$x)
    after <- classify_normative_fixture()

    expect_identical(formals(classify_missingness), classifier_formals)
    expect_identical(serialize(after, NULL, version = 3L), before_bytes)
    expect_identical(names(attributes(after)), c("names", "class"))
})
