.INPUT_FINGERPRINT_ALGORITHM <- "sha256"
.INPUT_FINGERPRINT_CANONICALIZATION <- "matrix_core_be_v1"
.ORIGINAL_MASK_ENCODING <- "base_packbits_lsb0_v1"
.SIDECAR_SCALE_DECLARATION <- "unnormalised_log2"

.SIDECAR_INPUT_FIELDS <- c(
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

.canonical_sidecar_names <- function(x) {
    output <- unname(x)
    text <- Encoding(output) != "bytes"
    output[text] <- enc2utf8(output[text])
    if (anyNA(output)) {
        .abort_sidecar(
            "Matrix identifiers must convert losslessly to UTF-8 or bytes.",
            "imputefinder_fingerprint_schema_error",
            field = "names"
        )
    }
    output
}

.fingerprint_be_i32 <- function(x) {
    writeBin(as.integer(x), raw(), size = 4L, endian = "big")
}

.encode_fingerprint_name <- function(value) {
    if (identical(Encoding(value), "bytes")) {
        tag <- as.raw(2L)
    } else {
        value <- enc2utf8(value)
        if (is.na(value)) {
            .abort_sidecar(
                "A matrix identifier could not be encoded for hashing.",
                "imputefinder_fingerprint_schema_error",
                field = "names"
            )
        }
        tag <- as.raw(1L)
    }

    bytes <- charToRaw(value)
    if (length(bytes) > .Machine$integer.max) {
        .abort_sidecar(
            "A matrix identifier exceeds the fingerprint protocol limit.",
            "imputefinder_fingerprint_schema_error",
            field = "names"
        )
    }
    c(tag, .fingerprint_be_i32(length(bytes)), bytes)
}

.encode_fingerprint_names <- function(values) {
    unlist(
        lapply(unname(values), .encode_fingerprint_name),
        use.names = FALSE
    )
}

.pack_original_mask <- function(data) {
    .validate_matrix_data(data)
    mask <- as.vector(is.na(data))
    padding <- as.integer((-length(mask)) %% 8)

    list(
        encoding = .ORIGINAL_MASK_ENCODING,
        cells = as.double(length(mask)),
        bytes = packBits(
            c(mask, rep(FALSE, padding)),
            type = "raw"
        )
    )
}

.validate_mask_dimensions <- function(dimensions) {
    valid <- is.numeric(dimensions) &&
        length(dimensions) == 2L &&
        !anyNA(dimensions) &&
        all(is.finite(dimensions)) &&
        all(dimensions >= 1) &&
        all(dimensions == floor(dimensions)) &&
        all(dimensions <= .Machine$integer.max)
    if (!valid) {
        .abort_sidecar(
            "Original-mask dimensions must be two positive integers.",
            "imputefinder_mask_schema_error",
            field = "dimensions"
        )
    }

    as.integer(unname(dimensions))
}

.mask_encoding_value <- function(mask) {
    if (is.list(mask) && "encoding" %in% names(mask)) {
        .sidecar_scalar_character(mask$encoding)
    } else {
        NA_character_
    }
}

.validate_original_mask_shape <- function(mask, dimensions) {
    valid <- is.list(mask) &&
        identical(names(mask), c("encoding", "cells", "bytes")) &&
        is.numeric(mask$cells) &&
        length(mask$cells) == 1L &&
        is.finite(mask$cells) &&
        is.raw(mask$bytes)
    cell_count <- prod(as.double(dimensions))
    expected_bytes <- ceiling(cell_count / 8)
    valid <- valid &&
        identical(as.double(mask$cells), cell_count) &&
        identical(as.double(length(mask$bytes)), expected_bytes)
    if (!valid) {
        .abort_sidecar(
            "Stored original-mask bytes do not match their dimensions.",
            "imputefinder_mask_schema_error",
            field = "original_mask"
        )
    }

    cell_count
}

.validate_original_mask <- function(mask, dimensions) {
    dimensions <- .validate_mask_dimensions(dimensions)
    actual_encoding <- .mask_encoding_value(mask)
    if (!identical(actual_encoding, .ORIGINAL_MASK_ENCODING)) {
        .abort_sidecar(
            paste0(
                "Unsupported original-mask encoding; reconstruct the ",
                "analysis with the current `analyze_missingness()`."
            ),
            "imputefinder_mask_lifecycle_error",
            expected_encoding = .ORIGINAL_MASK_ENCODING,
            actual_encoding = actual_encoding
        )
    }

    cell_count <- .validate_original_mask_shape(mask, dimensions)
    bits <- as.logical(rawToBits(mask$bytes))
    if (length(bits) > cell_count &&
        any(bits[seq.int(cell_count + 1, length(bits))])) {
        .abort_sidecar(
            "Stored original-mask padding bits must be zero.",
            "imputefinder_mask_schema_error",
            field = "padding"
        )
    }

    bits[seq_len(cell_count)]
}

.unpack_original_mask <- function(
    mask,
    dimensions,
    feature_names = NULL,
    sample_names = NULL
) {
    dimensions <- .validate_mask_dimensions(dimensions)
    bits <- .validate_original_mask(mask, dimensions)
    if (!is.null(feature_names)) {
        .validate_axis_names(
            feature_names,
            dimensions[[1L]],
            "feature",
            "Stored input"
        )
    }
    if (!is.null(sample_names)) {
        .validate_axis_names(
            sample_names,
            dimensions[[2L]],
            "sample",
            "Stored input"
        )
    }

    matrix(
        bits,
        nrow = dimensions[[1L]],
        ncol = dimensions[[2L]],
        dimnames = list(feature_names, sample_names)
    )
}

.canonical_matrix_bytes <- function(data, original_mask) {
    .validate_matrix_data(data)
    stored_missing <- .validate_original_mask(original_mask, dim(data))
    if (!identical(stored_missing, as.vector(is.na(data)))) {
        .abort_sidecar(
            "Original-mask bytes do not describe the matrix missing cells.",
            "imputefinder_fingerprint_schema_error",
            field = "original_mask"
        )
    }
    values <- as.vector(data)
    values[is.na(values)] <- if (is.integer(values)) 0L else 0
    storage_tag <- if (is.integer(values)) as.raw(1L) else as.raw(2L)
    value_size <- if (is.integer(values)) 4L else 8L

    c(
        charToRaw("imputefinder:matrix-core:be:v1"),
        as.raw(0L),
        storage_tag,
        .fingerprint_be_i32(dim(data)),
        .encode_fingerprint_names(rownames(data)),
        .encode_fingerprint_names(colnames(data)),
        original_mask$bytes,
        writeBin(
            values,
            raw(),
            size = value_size,
            endian = "big"
        )
    )
}

.matrix_fingerprint <- function(data, original_mask = NULL) {
    .validate_matrix_data(data)
    if (is.null(original_mask)) {
        original_mask <- .pack_original_mask(data)
    }
    bytes <- .canonical_matrix_bytes(data, original_mask)

    list(
        algorithm = .INPUT_FINGERPRINT_ALGORITHM,
        canonicalization = .INPUT_FINGERPRINT_CANONICALIZATION,
        value = unname(tools::sha256sum(bytes = bytes))
    )
}

.fingerprint_protocol_values <- function(fingerprint) {
    if (!is.list(fingerprint)) {
        return(c(algorithm = NA_character_, canonicalization = NA_character_))
    }
    c(
        algorithm = .sidecar_scalar_character(fingerprint$algorithm),
        canonicalization = .sidecar_scalar_character(
            fingerprint$canonicalization
        )
    )
}

.validate_matrix_fingerprint <- function(fingerprint) {
    actual <- .fingerprint_protocol_values(fingerprint)
    expected <- c(
        algorithm = .INPUT_FINGERPRINT_ALGORITHM,
        canonicalization = .INPUT_FINGERPRINT_CANONICALIZATION
    )
    if (!identical(actual, expected)) {
        .abort_sidecar(
            paste0(
                "Unsupported fingerprint protocol; reconstruct the ",
                "analysis with the current `analyze_missingness()`."
            ),
            "imputefinder_fingerprint_lifecycle_error",
            expected_algorithm = unname(expected[["algorithm"]]),
            actual_algorithm = unname(actual[["algorithm"]]),
            expected_canonicalization = unname(
                expected[["canonicalization"]]
            ),
            actual_canonicalization = unname(
                actual[["canonicalization"]]
            )
        )
    }

    valid <- identical(
        names(fingerprint),
        c("algorithm", "canonicalization", "value")
    ) &&
        is.character(fingerprint$value) &&
        length(fingerprint$value) == 1L &&
        !is.na(fingerprint$value) &&
        grepl("^[0-9a-f]{64}$", fingerprint$value)
    if (!valid) {
        .abort_sidecar(
            "Stored fingerprint metadata is malformed.",
            "imputefinder_fingerprint_schema_error",
            field = "fingerprint"
        )
    }

    invisible(fingerprint)
}

.normalise_sidecar_acquisition <- function(acquisition, sample_names) {
    if (is.null(acquisition)) {
        return(NULL)
    }
    valid <- is.atomic(acquisition) &&
        is.null(dim(acquisition)) &&
        length(acquisition) == length(sample_names) &&
        !anyNA(acquisition) &&
        identical(names(acquisition), sample_names)
    values <- if (valid) as.character(acquisition) else character()
    valid <- valid && !anyNA(values) && all(nzchar(values))
    if (!valid) {
        .abort_sidecar(
            paste0(
                "Acquisition declarations must be complete and named in ",
                "input sample order."
            ),
            "imputefinder_input_schema_error",
            field = "acquisition"
        )
    }

    stats::setNames(values, sample_names)
}

.validate_sidecar_representation <- function(representation, assay) {
    valid_representation <- is.character(representation) &&
        length(representation) == 1L &&
        !is.na(representation) &&
        representation %in% c("matrix", "SummarizedExperiment")
    valid_assay <- is.character(assay) &&
        length(assay) == 1L &&
        (is.na(assay) || nzchar(assay))
    matrix_assay <- identical(representation, "matrix") && is.na(assay)
    se_assay <- identical(representation, "SummarizedExperiment") &&
        valid_assay
    if (!valid_representation || !valid_assay ||
        (!matrix_assay && !se_assay)) {
        .abort_sidecar(
            "Stored representation and assay provenance is malformed.",
            "imputefinder_input_schema_error",
            field = "representation"
        )
    }

    invisible(representation)
}

.new_sidecar_input <- function(
    data,
    representation = "matrix",
    assay = NA_character_,
    acquisition = NULL
) {
    .validate_matrix_data(data)
    .validate_sidecar_representation(representation, assay)
    feature_names <- .canonical_sidecar_names(rownames(data))
    sample_names <- .canonical_sidecar_names(colnames(data))
    acquisition <- .normalise_sidecar_acquisition(
        acquisition,
        sample_names
    )
    original_mask <- .pack_original_mask(data)

    list(
        dimensions = c(
            features = as.integer(nrow(data)),
            samples = as.integer(ncol(data))
        ),
        feature_names = feature_names,
        sample_names = sample_names,
        storage = typeof(data),
        representation = representation,
        assay = assay,
        fingerprint = .matrix_fingerprint(data, original_mask),
        original_mask = original_mask,
        scale = .SIDECAR_SCALE_DECLARATION,
        acquisition = acquisition
    )
}

.validate_sidecar_input_shape <- function(input) {
    valid <- is.list(input) &&
        identical(names(input), .SIDECAR_INPUT_FIELDS) &&
        identical(
            names(input$dimensions),
            c("features", "samples")
        ) &&
        is.character(input$feature_names) &&
        is.character(input$sample_names) &&
        is.character(input$storage) &&
        length(input$storage) == 1L &&
        input$storage %in% c("integer", "double") &&
        identical(input$scale, .SIDECAR_SCALE_DECLARATION)
    if (!valid) {
        .abort_sidecar(
            "Stored input does not satisfy the sidecar input schema.",
            "imputefinder_input_schema_error",
            field = "input"
        )
    }

    invisible(input)
}

.validate_sidecar_input <- function(input) {
    .validate_sidecar_input_shape(input)
    dimensions <- .validate_mask_dimensions(input$dimensions)
    canonical_dimensions <- c(
        features = dimensions[[1L]],
        samples = dimensions[[2L]]
    )
    canonical_features <- .canonical_sidecar_names(input$feature_names)
    canonical_samples <- .canonical_sidecar_names(input$sample_names)
    canonical_acquisition <- .normalise_sidecar_acquisition(
        input$acquisition,
        input$sample_names
    )
    canonical <- identical(input$dimensions, canonical_dimensions) &&
        identical(input$feature_names, canonical_features) &&
        identical(input$sample_names, canonical_samples) &&
        identical(input$acquisition, canonical_acquisition)
    if (!canonical) {
        .abort_sidecar(
            "Stored input identifiers or declarations are not canonical.",
            "imputefinder_input_schema_error",
            field = "input"
        )
    }
    .validate_axis_names(
        input$feature_names,
        dimensions[[1L]],
        "feature",
        "Stored input"
    )
    .validate_axis_names(
        input$sample_names,
        dimensions[[2L]],
        "sample",
        "Stored input"
    )
    .validate_sidecar_representation(input$representation, input$assay)
    .validate_matrix_fingerprint(input$fingerprint)
    .validate_original_mask(input$original_mask, dimensions)

    invisible(input)
}

.input_mismatch_differences <- function(input, data, fingerprint, mask) {
    dimensions <- c(
        features = as.integer(nrow(data)),
        samples = as.integer(ncol(data))
    )
    differences <- character()
    if (!identical(input$dimensions, dimensions)) {
        differences <- c(differences, "dimensions")
    }
    if (!identical(
        input$feature_names,
        .canonical_sidecar_names(rownames(data))
    )) {
        differences <- c(differences, "feature_names")
    }
    if (!identical(
        input$sample_names,
        .canonical_sidecar_names(colnames(data))
    )) {
        differences <- c(differences, "sample_names")
    }
    if (!identical(input$original_mask, mask)) {
        differences <- c(differences, "original_mask")
    }
    if (!identical(input$storage, typeof(data))) {
        differences <- c(differences, "storage")
    }
    if (!identical(input$fingerprint$value, fingerprint$value)) {
        differences <- c(differences, "fingerprint")
    }

    differences
}

.verify_sidecar_input <- function(input, data) {
    .validate_sidecar_input(input)
    .validate_matrix_data(data)
    mask <- .pack_original_mask(data)
    fingerprint <- .matrix_fingerprint(data, mask)
    differences <- .input_mismatch_differences(
        input,
        data,
        fingerprint,
        mask
    )
    if (length(differences) > 0L) {
        .abort_sidecar(
            paste0(
                "Supplied input does not match the analysis fingerprint: ",
                paste(differences, collapse = ", "),
                "."
            ),
            "imputefinder_input_mismatch_error",
            differences = differences,
            expected_fingerprint = input$fingerprint$value,
            actual_fingerprint = fingerprint$value
        )
    }

    invisible(data)
}
