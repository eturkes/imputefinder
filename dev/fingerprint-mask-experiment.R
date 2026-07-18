.experiment_be_i32 <- function(x) {
    writeBin(as.integer(x), raw(), size = 4L, endian = "big")
}

.experiment_encode_name <- function(value) {
    if (identical(Encoding(value), "bytes")) {
        tag <- as.raw(2L)
    } else {
        value <- enc2utf8(value)
        stopifnot(length(value) == 1L, !is.na(value))
        tag <- as.raw(1L)
    }
    bytes <- charToRaw(value)
    c(tag, .experiment_be_i32(length(bytes)), bytes)
}

.experiment_encode_names <- function(values) {
    unlist(
        lapply(unname(values), .experiment_encode_name),
        use.names = FALSE
    )
}

.experiment_pack_mask <- function(x) {
    mask <- as.vector(is.na(x))
    padding <- as.integer((-length(mask)) %% 8)
    packBits(c(mask, rep(FALSE, padding)), type = "raw")
}

.experiment_canonical_bytes <- function(x) {
    mask <- .experiment_pack_mask(x)
    values <- as.vector(x)
    values[is.na(values)] <- if (is.integer(values)) 0L else 0
    storage_tag <- if (is.integer(values)) as.raw(1L) else as.raw(2L)
    value_size <- if (is.integer(values)) 4L else 8L

    c(
        charToRaw("imputefinder:matrix-core:be:v1"),
        as.raw(0L),
        storage_tag,
        .experiment_be_i32(dim(x)),
        .experiment_encode_names(rownames(x)),
        .experiment_encode_names(colnames(x)),
        mask,
        writeBin(
            values,
            raw(),
            size = value_size,
            endian = "big"
        )
    )
}

.experiment_hex <- function(x) {
    paste(format(x), collapse = "")
}

.experiment_hashers <- function(bytes) {
    hashers <- list(
        tools = function() tools::sha256sum(bytes = bytes)
    )
    if (requireNamespace("digest", quietly = TRUE)) {
        hashers$digest <- function() {
            digest::digest(bytes, algo = "sha256", serialize = FALSE)
        }
    }
    if (requireNamespace("openssl", quietly = TRUE)) {
        hashers$openssl <- function() {
            .experiment_hex(openssl::sha256(bytes))
        }
    }
    hashers
}

.experiment_median_ms <- function(fun, iterations = 20L, repeats = 7L) {
    elapsed <- replicate(
        repeats,
        system.time(for (index in seq_len(iterations)) fun())[["elapsed"]]
    )
    1000 * median(elapsed) / iterations
}

.experiment_matrix <- function() {
    set.seed(17L)
    x <- matrix(
        stats::rnorm(500000L),
        nrow = 10000L,
        dimnames = list(
            sprintf("protein_%05d", seq_len(10000L)),
            sprintf("sample_%02d", seq_len(50L))
        )
    )
    x[sample.int(length(x), 100000L)] <- NA_real_
    x
}

.experiment_hash_benchmark <- function(x) {
    bytes <- .experiment_canonical_bytes(x)
    hashers <- .experiment_hashers(bytes)
    hashes <- vapply(hashers, function(hasher) hasher(), character(1L))
    stopifnot(length(unique(hashes)) == 1L)

    data.frame(
        candidate = names(hashers),
        median_ms = vapply(
            hashers,
            function(hasher) .experiment_median_ms(hasher),
            numeric(1L)
        ),
        payload_bytes = length(bytes),
        hash = unname(hashes),
        row.names = NULL
    )
}

.experiment_mask_benchmark <- function(cell_count = 500000L) {
    score <- ((seq_len(cell_count) * 104729) %% 1000003) / 1000003
    rows <- lapply(
        c(0.01, 0.2, 0.5),
        function(rate) {
            mask <- score < rate
            padding <- as.integer((-length(mask)) %% 8)
            padded <- c(mask, rep(FALSE, padding))
            packed <- packBits(padded, type = "raw")
            unpack <- function() {
                as.logical(rawToBits(packed))[seq_along(mask)]
            }
            stopifnot(identical(mask, unpack()))

            representations <- list(
                logical = mask,
                raw = as.raw(mask),
                packed = packed,
                sparse_indices = which(mask)
            )
            sizes <- vapply(
                representations,
                function(value) {
                    length(serialize(value, NULL, version = 3L))
                },
                integer(1L)
            )
            data.frame(
                missing_rate = rate,
                missing_cells = sum(mask),
                logical_bytes = sizes[["logical"]],
                raw_bytes = sizes[["raw"]],
                packed_bytes = sizes[["packed"]],
                sparse_index_bytes = sizes[["sparse_indices"]],
                pack_median_ms = .experiment_median_ms(
                    function() packBits(padded, type = "raw"),
                    iterations = 50L
                ),
                unpack_median_ms = .experiment_median_ms(
                    unpack,
                    iterations = 50L
                )
            )
        }
    )
    do.call(rbind, rows)
}

x <- .experiment_matrix()
hash_benchmark <- .experiment_hash_benchmark(x)
mask_benchmark <- .experiment_mask_benchmark()

known <- matrix(
    c(1, NA, -0, -2),
    nrow = 2L,
    dimnames = list(c("f1", "f2"), c("s1", "s2"))
)
known_hash <- tools::sha256sum(
    bytes = .experiment_canonical_bytes(known)
)
stopifnot(
    identical(
        known_hash,
        "a09c1f5501caee9b7d9a6dcce7322e56c2028a19c103c1f3bba7764d6338a24a"
    )
)

utf8_name <- "caf\u00e9"
latin1_name <- iconv(utf8_name, from = "UTF-8", to = "latin1")
Encoding(latin1_name) <- "latin1"
utf8_matrix <- latin1_matrix <- known
rownames(utf8_matrix)[[1L]] <- utf8_name
rownames(latin1_matrix)[[1L]] <- latin1_name
stopifnot(identical(
    tools::sha256sum(bytes = .experiment_canonical_bytes(utf8_matrix)),
    tools::sha256sum(bytes = .experiment_canonical_bytes(latin1_matrix))
))

print(hash_benchmark, row.names = FALSE)
print(mask_benchmark, row.names = FALSE)
cat("known fixture SHA-256:", known_hash, "\n")
cat("R:", R.version.string, "\n")
