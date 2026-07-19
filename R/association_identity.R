.ASSOCIATION_PROTOCOL_ID <- "m13_a_association_protocol_v3"
.ASSOCIATION_PROTOCOL_HASH <- "254e9cb03d3f981bbbcf621b8e13985c8abb43916ab9b79901b32bac13b55c59"
.ASSOCIATION_CONTRACT_HASH <- "04e824321ac68c73633277064397f1a3fa96ba73ff342fbfe575df7fc6b1f4a6"
.ASSOCIATION_INPUT_HASH <- "association_mask_design_be_v1"
.ASSOCIATION_HASH_ROLES <- c(
    "condition", "nuisance", "block", "acquisition"
)

.abort_association <- function(
    message,
    subclass = "imputefinder_association_error",
    ...
) {
    .abort_sidecar(message, subclass, ...)
}

.association_be_i32 <- function(x) {
    valid <- is.numeric(x) && !anyNA(x) && all(x >= 0) &&
        all(x <= .Machine$integer.max) && all(x == floor(x))
    if (!valid) {
        .abort_association(
            "Association identity lengths must fit signed 32-bit integers.",
            "imputefinder_association_identity_error"
        )
    }
    writeBin(as.integer(x), raw(), size = 4L, endian = "big")
}

.association_canonical_text <- function(x) {
    output <- unname(as.character(x))
    text <- Encoding(output) != "bytes"
    output[text] <- enc2utf8(output[text])
    if (anyNA(output)) {
        .abort_association(
            "Association identity text must encode losslessly.",
            "imputefinder_association_identity_error"
        )
    }
    output
}

.association_encode_text <- function(value) {
    value <- .association_canonical_text(value)
    if (length(value) != 1L) {
        .abort_association(
            "Association identity text records must be scalar.",
            "imputefinder_association_identity_error"
        )
    }
    tag <- if (identical(Encoding(value), "bytes")) {
        as.raw(2L)
    } else {
        as.raw(1L)
    }
    bytes <- charToRaw(value)
    c(tag, .association_be_i32(length(bytes)), bytes)
}

.association_encode_text_vector <- function(values) {
    values <- .association_canonical_text(values)
    c(
        .association_be_i32(length(values)),
        unlist(lapply(values, .association_encode_text), use.names = FALSE)
    )
}

.association_raw_hex <- function(bytes) {
    paste(sprintf("%02x", as.integer(bytes)), collapse = "")
}

.association_hash_roles <- function(roles, sample_data) {
    valid <- is.list(roles) &&
        setequal(names(roles), .ASSOCIATION_HASH_ROLES)
    if (!valid) {
        .abort_association(
            "Association identity roles are malformed.",
            "imputefinder_association_identity_error"
        )
    }
    roles <- roles[.ASSOCIATION_HASH_ROLES]
    roles <- lapply(seq_along(roles), function(index) {
        value <- .association_canonical_text(roles[[index]])
        if (identical(.ASSOCIATION_HASH_ROLES[[index]], "nuisance")) {
            value <- sort(value, method = "radix")
        }
        value
    })
    names(roles) <- .ASSOCIATION_HASH_ROLES
    cardinality <- lengths(roles)
    declared <- unname(unlist(roles, use.names = FALSE))
    valid <- cardinality[[1L]] == 1L &&
        all(cardinality[c(3L, 4L)] <= 1L) &&
        !anyDuplicated(declared) &&
        all(declared %in% names(sample_data))
    if (!valid) {
        .abort_association(
            "Association identity roles have invalid cardinality or columns.",
            "imputefinder_association_identity_error"
        )
    }
    roles
}

.association_hash_interactions <- function(interactions, declared) {
    if (is.null(interactions)) {
        interactions <- list()
    }
    if (!is.list(interactions)) {
        .abort_association(
            "Association identity interactions must be a list.",
            "imputefinder_association_identity_error"
        )
    }
    output <- lapply(interactions, function(components) {
        components <- sort(
            .association_canonical_text(components),
            method = "radix"
        )
        valid <- length(components) >= 2L &&
            !anyDuplicated(components) && all(components %in% declared)
        if (!valid) {
            .abort_association(
                "Association identity interaction components are malformed.",
                "imputefinder_association_identity_error"
            )
        }
        components
    })
    if (length(output)) {
        encoded <- lapply(output, .association_encode_text_vector)
        keys <- vapply(encoded, .association_raw_hex, character(1L))
        if (anyDuplicated(keys)) {
            .abort_association(
                "Association identity interactions must be unique.",
                "imputefinder_association_identity_error"
            )
        }
        output <- output[order(keys, method = "radix")]
    }
    output
}

.association_encode_roles <- function(roles) {
    c(
        .association_be_i32(length(roles)),
        unlist(lapply(names(roles), function(role) {
            c(
                .association_encode_text(role),
                .association_encode_text_vector(roles[[role]])
            )
        }), use.names = FALSE)
    )
}

.association_encode_interactions <- function(interactions) {
    c(
        .association_be_i32(length(interactions)),
        unlist(
            lapply(interactions, .association_encode_text_vector),
            use.names = FALSE
        )
    )
}

.association_encode_design_column <- function(column, values) {
    numeric <- !is.factor(values) &&
        (is.integer(values) || is.double(values))
    if (numeric) {
        values <- as.double(values)
        values[values == 0] <- 0
        if (anyNA(values) || any(!is.finite(values))) {
            .abort_association(
                "Numeric association identity values must be finite.",
                "imputefinder_association_identity_error"
            )
        }
        payload <- c(
            as.raw(1L),
            .association_be_i32(length(values)),
            writeBin(values, raw(), size = 8L, endian = "big")
        )
    } else {
        values <- .association_canonical_text(values)
        if (!length(values) || any(!nzchar(values))) {
            .abort_association(
                "Categorical association identity values must be nonempty.",
                "imputefinder_association_identity_error"
            )
        }
        payload <- c(as.raw(2L), .association_encode_text_vector(values))
    }
    c(.association_encode_text(column), payload)
}

.association_mask_design_bytes <- function(
    missing_mask,
    sample_data,
    roles,
    interactions = list()
) {
    valid_mask <- is.matrix(missing_mask) && is.logical(missing_mask) &&
        !anyNA(missing_mask) && nrow(missing_mask) > 0L &&
        ncol(missing_mask) > 0L && !is.null(rownames(missing_mask)) &&
        !is.null(colnames(missing_mask)) &&
        !anyNA(rownames(missing_mask)) && !anyNA(colnames(missing_mask)) &&
        all(nzchar(rownames(missing_mask))) &&
        all(nzchar(colnames(missing_mask))) &&
        !anyDuplicated(rownames(missing_mask)) &&
        !anyDuplicated(colnames(missing_mask))
    valid_data <- is.data.frame(sample_data) &&
        nrow(sample_data) == ncol(missing_mask) &&
        !is.null(rownames(sample_data)) && !anyNA(rownames(sample_data)) &&
        !anyDuplicated(rownames(sample_data)) &&
        setequal(rownames(sample_data), colnames(missing_mask))
    if (!valid_mask || !valid_data) {
        .abort_association(
            "Association mask/design identity input is malformed.",
            "imputefinder_association_identity_error"
        )
    }

    feature_names <- .association_canonical_text(rownames(missing_mask))
    sample_names <- .association_canonical_text(colnames(missing_mask))
    data_names <- .association_canonical_text(rownames(sample_data))
    column_names <- .association_canonical_text(names(sample_data))
    valid_names <- !anyDuplicated(feature_names) &&
        !anyDuplicated(sample_names) && !anyDuplicated(data_names) &&
        !anyDuplicated(column_names) && setequal(sample_names, data_names)
    if (!valid_names) {
        .abort_association(
            "Canonical association identifiers must remain unique and aligned.",
            "imputefinder_association_identity_error"
        )
    }
    rownames(sample_data) <- data_names
    names(sample_data) <- column_names

    feature_order <- order(feature_names, method = "radix")
    sample_order <- order(sample_names, method = "radix")
    feature_names <- feature_names[feature_order]
    sample_names <- sample_names[sample_order]
    missing_mask <- missing_mask[
        feature_order,
        sample_order,
        drop = FALSE
    ]
    sample_data <- sample_data[
        match(sample_names, data_names),
        ,
        drop = FALSE
    ]
    rownames(sample_data) <- sample_names

    roles <- .association_hash_roles(roles, sample_data)
    declared <- unname(unlist(roles, use.names = FALSE))
    interactions <- .association_hash_interactions(interactions, declared)
    mask_bits <- as.vector(t(missing_mask))
    padding <- as.integer((-length(mask_bits)) %% 8L)
    mask_bytes <- packBits(
        c(mask_bits, rep(FALSE, padding)),
        type = "raw"
    )
    design_bytes <- unlist(lapply(declared, function(column) {
        .association_encode_design_column(column, sample_data[[column]])
    }), use.names = FALSE)

    c(
        charToRaw("imputefinder:association-mask-design:be:v1"),
        as.raw(0L),
        as.raw(1L), .association_be_i32(dim(missing_mask)),
        as.raw(2L), .association_encode_text_vector(feature_names),
        as.raw(3L), .association_encode_text_vector(sample_names),
        as.raw(4L), .association_be_i32(length(mask_bits)),
        .association_be_i32(length(mask_bytes)), mask_bytes,
        as.raw(5L), .association_encode_roles(roles),
        as.raw(6L), .association_encode_interactions(interactions),
        as.raw(7L), .association_be_i32(length(declared)), design_bytes
    )
}

.association_mask_design_sha256 <- function(
    missing_mask,
    sample_data,
    roles,
    interactions = list()
) {
    unname(tools::sha256sum(bytes = .association_mask_design_bytes(
        missing_mask,
        sample_data,
        roles,
        interactions
    )))
}

.association_input_sha256 <- function(data, design) {
    .validate_matrix_data(data)
    design <- .align_missingness_design(design, colnames(data))
    .association_mask_design_sha256(
        is.na(data),
        design$sample_data,
        design$roles,
        design$interactions
    )
}

.association_stratum_id <- function(acquisition = NULL) {
    if (is.null(acquisition)) {
        return("all")
    }
    valid <- is.character(acquisition) && length(acquisition) == 1L &&
        !is.na(acquisition) && nzchar(acquisition)
    if (!valid) {
        .abort_association(
            "Acquisition stratum labels must be nonempty scalars.",
            "imputefinder_association_identity_error"
        )
    }
    bytes <- c(
        charToRaw("imputefinder:association-stratum:v1"),
        as.raw(0L),
        .association_encode_text(acquisition)
    )
    paste0("s_", unname(tools::sha256sum(bytes = bytes)))
}

.association_hypothesis_id <- function(
    stratum,
    coefficient,
    label,
    term_id
) {
    values <- c(stratum, coefficient, label, term_id)
    valid <- is.character(values) && length(values) == 4L &&
        !anyNA(values) && all(nzchar(values))
    if (!valid) {
        .abort_association(
            "Association hypothesis identity fields must be nonempty scalars.",
            "imputefinder_association_identity_error"
        )
    }
    bytes <- c(
        charToRaw("imputefinder:association-hypothesis:v1"),
        as.raw(0L),
        .association_encode_text_vector(values)
    )
    paste0("a_", unname(tools::sha256sum(bytes = bytes)))
}
