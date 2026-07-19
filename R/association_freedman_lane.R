.ASSOCIATION_FREEDMAN_LANE_CANDIDATE <- "a_fraction_freedman_lane"
.ASSOCIATION_EXACT_PERMUTATION_MAX <- 100000
.ASSOCIATION_MONTE_CARLO_DRAWS <- 9999L
.ASSOCIATION_MAP_RETRY_LIMIT <- 1000000L

.abort_association_permutation <- function(message) {
    .abort_association(
        message,
        "imputefinder_association_permutation_error"
    )
}

.association_permutation_map_sha256 <- function(map) {
    valid <- is.integer(map) && length(map) > 0L && !anyNA(map) &&
        identical(sort(map, method = "radix"), seq_along(map))
    if (!valid) {
        .abort_association_permutation(
            "Association permutation maps must be exact integer bijections."
        )
    }
    bytes <- c(
        charToRaw("imputefinder:association-permutation-map:v1"),
        as.raw(0L),
        .association_be_i32(length(map)),
        .association_be_i32(map)
    )
    unname(tools::sha256sum(bytes = bytes))
}

.association_seed_manifest <- function(input_sha256, instances, draw_ids) {
    valid <- .association_sha256(input_sha256) &&
        is.data.frame(instances) && identical(
            names(instances),
            c("candidate_id", "acquisition", "hypothesis_id")
        ) && nrow(instances) > 0L &&
        all(vapply(instances, is.character, logical(1L))) &&
        !anyNA(instances) &&
        all(instances$candidate_id == .ASSOCIATION_FREEDMAN_LANE_CANDIDATE) &&
        all(nzchar(instances$acquisition)) &&
        all(grepl("^a_[0-9a-f]{64}$", instances$hypothesis_id)) &&
        !anyDuplicated(instances) &&
        is.integer(draw_ids) && length(draw_ids) > 0L &&
        !anyNA(draw_ids) && all(draw_ids > 0L) &&
        !anyDuplicated(draw_ids)
    if (!valid) {
        .abort_association_permutation(
            "Association seed-manifest input is malformed."
        )
    }
    order_index <- do.call(
        order,
        c(unname(instances), list(method = "radix"))
    )
    instances <- instances[order_index, , drop = FALSE]
    row.names(instances) <- NULL
    grid <- instances[
        rep(seq_len(nrow(instances)), each = length(draw_ids)),
        ,
        drop = FALSE
    ]
    grid$draw_id <- rep(
        sort(draw_ids, method = "radix"),
        nrow(instances)
    )
    row.names(grid) <- NULL
    used <- new.env(hash = TRUE, parent = emptyenv())
    seed <- seed_nonce <- integer(nrow(grid))
    for (index in seq_len(nrow(grid))) {
        candidate_nonce <- 0L
        repeat {
            key <- c(
                .ASSOCIATION_PROTOCOL_ID,
                input_sha256,
                grid$candidate_id[[index]],
                grid$acquisition[[index]],
                grid$hypothesis_id[[index]],
                as.character(grid$draw_id[[index]]),
                as.character(candidate_nonce)
            )
            digest <- unname(tools::sha256sum(
                bytes = c(
                    charToRaw(
                        "imputefinder:association-permutation-seed:v2"
                    ),
                    as.raw(0L),
                    .association_encode_text_vector(key)
                )
            ))
            candidate_seed <- as.integer(
                strtoi(substr(digest, 1L, 7L), 16L) + 1L
            )
            seed_key <- as.character(candidate_seed)
            if (!exists(seed_key, envir = used, inherits = FALSE)) {
                assign(seed_key, TRUE, envir = used)
                break
            }
            candidate_nonce <- candidate_nonce + 1L
        }
        seed[[index]] <- candidate_seed
        seed_nonce[[index]] <- candidate_nonce
    }
    data.frame(
        protocol_id = .ASSOCIATION_PROTOCOL_ID,
        candidate_id = grid$candidate_id,
        acquisition = grid$acquisition,
        hypothesis_id = grid$hypothesis_id,
        draw_id = as.integer(grid$draw_id),
        seed = seed,
        seed_nonce = seed_nonce,
        stringsAsFactors = FALSE
    )
}

.association_canonical_permutation_groups <- function(groups, n) {
    valid <- is.list(groups) && length(groups) > 0L &&
        all(vapply(groups, is.integer, logical(1L))) &&
        all(lengths(groups) > 0L)
    if (!valid) {
        .abort_association_permutation(
            "Independent Monte Carlo groups are malformed."
        )
    }
    groups <- lapply(groups, sort, method = "radix")
    first <- vapply(groups, `[[`, integer(1L), 1L)
    groups <- groups[order(first, method = "radix")]
    if (!identical(
        sort(unlist(groups, use.names = FALSE), method = "radix"),
        seq_len(n)
    )) {
        .abort_association_permutation(
            "Independent Monte Carlo groups must partition rows."
        )
    }
    groups
}

.association_canonical_permutation_swaps <- function(swaps, n) {
    identity <- seq_len(n)
    valid <- is.list(swaps) && length(swaps) > 0L &&
        all(vapply(swaps, function(map) {
            is.integer(map) && length(map) == n && !anyNA(map) &&
                identical(sort(map, method = "radix"), identity) &&
                identical(map[map], identity) &&
                !identical(map, identity)
        }, logical(1L)))
    if (!valid) {
        .abort_association_permutation(
            "Blocked Monte Carlo swaps are malformed."
        )
    }
    moved <- lapply(swaps, function(map) which(map != identity))
    if (anyDuplicated(unlist(moved, use.names = FALSE))) {
        .abort_association_permutation(
            "Blocked Monte Carlo swaps must move disjoint rows."
        )
    }
    swaps
}

.association_restore_rng <- function(kind, had_seed, seed) {
    do.call(RNGkind, as.list(kind))
    if (had_seed) {
        assign(".Random.seed", seed, envir = .GlobalEnv)
    } else if (exists(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
    )) {
        rm(".Random.seed", envir = .GlobalEnv)
    }
}

.association_monte_carlo_maps <- function(
    seed_manifest,
    design,
    n,
    groups = NULL,
    swaps = NULL
) {
    seed_fields <- c(
        "protocol_id", "candidate_id", "acquisition", "hypothesis_id",
        "draw_id", "seed", "seed_nonce"
    )
    valid <- is.data.frame(seed_manifest) &&
        identical(names(seed_manifest), seed_fields) &&
        nrow(seed_manifest) > 0L &&
        all(vapply(seed_manifest[c(
            "protocol_id", "candidate_id", "acquisition", "hypothesis_id"
        )], is.character, logical(1L))) &&
        !anyNA(seed_manifest) &&
        all(seed_manifest$protocol_id == .ASSOCIATION_PROTOCOL_ID) &&
        is.integer(seed_manifest$draw_id) &&
        identical(
            seed_manifest$draw_id,
            sort(seed_manifest$draw_id, method = "radix")
        ) && all(seed_manifest$draw_id > 0L) &&
        !anyDuplicated(seed_manifest$draw_id) &&
        length(unique(seed_manifest$candidate_id)) == 1L &&
        identical(
            unique(seed_manifest$candidate_id),
            .ASSOCIATION_FREEDMAN_LANE_CANDIDATE
        ) && length(unique(seed_manifest$acquisition)) == 1L &&
        nzchar(unique(seed_manifest$acquisition)) &&
        length(unique(seed_manifest$hypothesis_id)) == 1L &&
        grepl(
            "^a_[0-9a-f]{64}$",
            unique(seed_manifest$hypothesis_id)
        ) &&
        is.integer(seed_manifest$seed) && all(seed_manifest$seed > 0L) &&
        !anyDuplicated(seed_manifest$seed) &&
        is.integer(seed_manifest$seed_nonce) &&
        all(seed_manifest$seed_nonce >= 0L) &&
        is.integer(n) && length(n) == 1L && !is.na(n) && n > 0L &&
        .association_scalar_character(design) &&
        design %in% c("independent", "blocked")
    if (!valid) {
        .abort_association_permutation(
            "Monte Carlo seed rows are malformed."
        )
    }
    if (identical(design, "independent")) {
        groups <- .association_canonical_permutation_groups(groups, n)
        if (sum(lgamma(lengths(groups) + 1)) <=
            log(.ASSOCIATION_EXACT_PERMUTATION_MAX)) {
            .abort_association_permutation(
                "Monte Carlo permutation requires more than 100000 maps."
            )
        }
    } else {
        swaps <- .association_canonical_permutation_swaps(swaps, n)
        if (length(swaps) * log(2) <=
            log(.ASSOCIATION_EXACT_PERMUTATION_MAX)) {
            .abort_association_permutation(
                "Monte Carlo permutation requires more than 100000 maps."
            )
        }
    }

    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) {
        get(".Random.seed", envir = .GlobalEnv)
    } else {
        NULL
    }
    on.exit(
        .association_restore_rng(old_kind, had_seed, old_seed),
        add = TRUE
    )
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")

    identity <- seq_len(n)
    maps <- vector("list", nrow(seed_manifest))
    map_retry <- integer(nrow(seed_manifest))
    map_sha256 <- character(nrow(seed_manifest))
    accepted <- new.env(hash = TRUE, parent = emptyenv())
    for (index in seq_len(nrow(seed_manifest))) {
        set.seed(seed_manifest$seed[[index]])
        retry <- 0L
        repeat {
            map <- identity
            if (identical(design, "independent")) {
                for (group in groups) {
                    map[group] <- group[sample.int(length(group))]
                }
            } else {
                selected <- sample.int(
                    2L,
                    length(swaps),
                    replace = TRUE
                ) - 1L
                for (swap_index in which(selected == 1L)) {
                    map <- swaps[[swap_index]][map]
                }
            }
            hash <- .association_permutation_map_sha256(map)
            duplicate <- exists(hash, envir = accepted, inherits = FALSE)
            if (duplicate && !identical(
                get(hash, envir = accepted),
                map
            )) {
                .abort_association_permutation(
                    "Association permutation-map SHA-256 collision."
                )
            }
            if (!identical(map, identity) && !duplicate) {
                assign(hash, map, envir = accepted)
                break
            }
            retry <- retry + 1L
            if (retry > .ASSOCIATION_MAP_RETRY_LIMIT) {
                .abort_association(
                    "Association Monte Carlo map retry limit exceeded.",
                    "imputefinder_association_permutation_retry_error"
                )
            }
        }
        maps[[index]] <- map
        map_retry[[index]] <- retry
        map_sha256[[index]] <- hash
    }
    names(maps) <- as.character(seed_manifest$draw_id)
    manifest <- cbind(
        seed_manifest,
        map_retry = map_retry,
        map_sha256 = map_sha256,
        stringsAsFactors = FALSE
    )
    list(manifest = manifest, maps = maps)
}

.association_permutation_tolerance <- function(z0) {
    if (!is.matrix(z0) || !is.numeric(z0) || any(!is.finite(z0))) {
        return(NA_real_)
    }
    maximum <- if (length(z0)) max(abs(z0)) else 0
    sqrt(.Machine$double.eps) * max(1, maximum)
}

.association_constrained_null <- function(z, contrast, response) {
    valid <- is.matrix(z) && is.numeric(z) && nrow(z) > 0L &&
        ncol(z) > 0L && all(is.finite(z)) &&
        is.double(contrast) && length(contrast) == ncol(z) &&
        all(is.finite(contrast)) &&
        is.double(response) && length(response) == nrow(z) &&
        all(is.finite(response))
    magnitude <- if (valid) sum(contrast^2) else NA_real_
    if (!valid || !is.finite(magnitude) || magnitude <= 0) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    rank <- ncol(z)
    projector <- diag(rank) - tcrossprod(contrast) / magnitude
    projector <- (projector + t(projector)) / 2
    dimnames(projector) <- list(colnames(z), colnames(z))
    basis <- tryCatch(
        .design_null_basis(projector, rank - 1L),
        error = function(error) NULL
    )
    if (is.null(basis)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    z0 <- z %*% basis
    if (ncol(z0)) {
        decomposition <- tryCatch(
            .design_svd(z0),
            error = function(error) NULL
        )
        if (is.null(decomposition) ||
            decomposition$rank != ncol(z0)) {
            return(.association_robust_failure(
                "association_numerical_failure"
            ))
        }
        left <- decomposition$decomposition$u[
            ,
            seq_len(decomposition$rank),
            drop = FALSE
        ]
        fitted <- as.vector(left %*% crossprod(left, response))
    } else {
        fitted <- rep(0, nrow(z))
    }
    residual <- response - fitted
    tolerance <- .association_permutation_tolerance(z0)
    if (any(!is.finite(c(fitted, residual, tolerance)))) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    list(
        ok = TRUE,
        basis = basis,
        z0 = z0,
        fitted = fitted,
        residual = residual,
        tolerance = tolerance
    )
}

.association_lexicographic_permutations <- function(values) {
    values <- sort(as.integer(values), method = "radix")
    output <- list(values)
    if (length(values) < 2L) {
        return(output)
    }
    repeat {
        rises <- which(values[-length(values)] < values[-1L])
        if (!length(rises)) {
            break
        }
        pivot <- max(rises)
        swap <- max(which(values > values[[pivot]]))
        held <- values[[pivot]]
        values[[pivot]] <- values[[swap]]
        values[[swap]] <- held
        suffix <- seq.int(pivot + 1L, length(values))
        values[suffix] <- rev(values[suffix])
        output[[length(output) + 1L]] <- values
    }
    output
}

.association_exact_group_maps <- function(groups, n) {
    identity <- seq_len(n)
    maps <- list(identity)
    for (group in groups) {
        permutations <- .association_lexicographic_permutations(group)
        expanded <- vector("list", length(maps) * length(permutations))
        cursor <- 0L
        for (map in maps) {
            for (permutation in permutations) {
                cursor <- cursor + 1L
                expanded[[cursor]] <- map
                expanded[[cursor]][group] <- permutation
            }
        }
        maps <- expanded
    }
    matrix_maps <- do.call(rbind, maps)
    order_arguments <- lapply(seq_len(n), function(column) {
        matrix_maps[, column]
    })
    ordered <- do.call(order, c(order_arguments, list(method = "radix")))
    maps[ordered]
}

.association_group_components <- function(adjacency) {
    n <- nrow(adjacency)
    remaining <- rep(TRUE, n)
    output <- list()
    while (any(remaining)) {
        start <- which(remaining)[[1L]]
        component <- start
        frontier <- start
        remaining[[start]] <- FALSE
        while (length(frontier)) {
            neighbours <- which(colSums(
                adjacency[frontier, , drop = FALSE]
            ) > 0L & remaining)
            if (!length(neighbours)) {
                break
            }
            component <- c(component, neighbours)
            remaining[neighbours] <- FALSE
            frontier <- neighbours
        }
        output[[length(output) + 1L]] <- sort(
            as.integer(component),
            method = "radix"
        )
    }
    output
}

.association_transformation_summary <- function(log_count, count) {
    if (!is.double(log_count) || length(log_count) != 1L ||
        !is.finite(log_count) || !is.double(count) || length(count) != 1L ||
        is.na(count) || count <= 0) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    allowable <- if (is.finite(count)) count else .Machine$double.xmax
    mode <- if (log_count < log(20)) {
        "low_resolution"
    } else if (log_count <= log(.ASSOCIATION_EXACT_PERMUTATION_MAX)) {
        "exact"
    } else {
        "monte_carlo"
    }
    list(
        ok = TRUE,
        allowable_transformations = as.double(allowable),
        mode = mode
    )
}

.association_integer_product_double <- function(factors) {
    valid <- is.integer(factors) && !anyNA(factors) && all(factors >= 1L)
    if (!valid) {
        return(NA_real_)
    }
    base <- 2^20
    digits <- 1
    for (factor in factors) {
        carry <- 0
        for (index in seq_along(digits)) {
            total <- digits[[index]] * factor + carry
            carry <- floor(total / base)
            digits[[index]] <- total - carry * base
        }
        while (carry > 0) {
            quotient <- floor(carry / base)
            digits[[length(digits) + 1L]] <- carry - quotient * base
            carry <- quotient
        }
        highest_bits <- floor(log2(digits[[length(digits)]])) + 1L
        if ((length(digits) - 1L) * 20L + highest_bits > 1024L) {
            return(.Machine$double.xmax)
        }
    }
    bits <- unlist(lapply(digits, function(digit) {
        as.integer(intToBits(as.integer(digit))[seq_len(20L)])
    }), use.names = FALSE)
    highest <- max(which(bits != 0L))
    bits <- bits[seq_len(highest)]
    bit_count <- length(bits)
    if (bit_count > 1024L) {
        return(.Machine$double.xmax)
    }
    shift <- max(0L, bit_count - 53L)
    quotient_bits <- bits[seq.int(shift + 1L, bit_count)]
    quotient <- 0
    for (bit in rev(quotient_bits)) {
        quotient <- quotient * 2 + bit
    }
    if (bit_count == 1024L &&
        quotient == 2^53 - 1 && shift > 0L &&
        any(bits[seq_len(shift)] != 0L)) {
        return(.Machine$double.xmax)
    }
    if (shift > 0L && bits[[shift]] == 1L) {
        above_half <- shift > 1L && any(
            bits[seq_len(shift - 1L)] != 0L
        )
        if (above_half || quotient %% 2 == 1) {
            quotient <- quotient + 1
            if (quotient == 2^53) {
                quotient <- 2^52
                shift <- shift + 1L
            }
        }
    }
    output <- quotient * 2^shift
    if (is.finite(output)) as.double(output) else .Machine$double.xmax
}

.association_factorial_product_double <- function(sizes) {
    if (!is.integer(sizes) || !length(sizes) || anyNA(sizes) ||
        any(sizes < 1L)) {
        return(NA_real_)
    }
    sizes <- sizes[sizes > 1L]
    if (!length(sizes)) {
        return(1.0)
    }
    if (any(sizes >= 171L) ||
        sum(lgamma(sizes + 1)) > log(.Machine$double.xmax) + 1e-8) {
        return(.Machine$double.xmax)
    }
    factors <- unlist(lapply(sizes, function(size) {
        if (size < 2L) integer() else seq.int(2L, size)
    }), use.names = FALSE)
    .association_integer_product_double(as.integer(factors))
}

.association_independent_permutations <- function(z0, enumerate = TRUE) {
    tolerance <- .association_permutation_tolerance(z0)
    if (!is.logical(enumerate) || length(enumerate) != 1L ||
        is.na(enumerate) || !is.finite(tolerance) || nrow(z0) < 1L) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    n <- nrow(z0)
    adjacency <- matrix(TRUE, nrow = n, ncol = n)
    if (ncol(z0)) {
        for (row in seq_len(n)) {
            adjacency[row, ] <- apply(
                abs(sweep(z0, 2L, z0[row, ], `-`)),
                1L,
                max
            ) <= tolerance
        }
    }
    groups <- .association_group_components(adjacency)
    transitive <- all(vapply(groups, function(group) {
        all(adjacency[group, group, drop = FALSE])
    }, logical(1L)))
    if (!transitive) {
        return(.association_robust_failure(
            "association_incompatible_permutation"
        ))
    }
    groups <- .association_canonical_permutation_groups(groups, n)
    summary <- .association_transformation_summary(
        sum(lgamma(lengths(groups) + 1)),
        .association_factorial_product_double(lengths(groups))
    )
    if (!summary$ok) {
        return(summary)
    }
    allowable <- summary$allowable_transformations
    mode <- summary$mode
    maps <- if (enumerate && mode %in% c("low_resolution", "exact")) {
        .association_exact_group_maps(groups, n)
    } else {
        NULL
    }
    list(
        ok = TRUE,
        design = "independent",
        groups = groups,
        swaps = NULL,
        allowable_transformations = allowable,
        mode = mode,
        maps = maps,
        tolerance = tolerance
    )
}

.association_row_signatures <- function(data, columns) {
    if (!length(columns)) {
        return(rep("00000000", nrow(data)))
    }
    vapply(seq_len(nrow(data)), function(index) {
        bytes <- unlist(lapply(columns, function(column) {
            .association_encode_design_column(column, data[[column]][index])
        }), use.names = FALSE)
        .association_raw_hex(bytes)
    }, character(1L))
}

.association_block_swap <- function(
    block_rows,
    condition,
    reference,
    target,
    signatures,
    n
) {
    reference_rows <- block_rows[condition[block_rows] == reference]
    target_rows <- block_rows[condition[block_rows] == target]
    if (!length(reference_rows) || !length(target_rows)) {
        return(NULL)
    }
    reference_groups <- split(reference_rows, signatures[reference_rows])
    target_groups <- split(target_rows, signatures[target_rows])
    reference_names <- sort(names(reference_groups), method = "radix")
    target_names <- sort(names(target_groups), method = "radix")
    if (!identical(reference_names, target_names) ||
        any(lengths(reference_groups[reference_names]) !=
            lengths(target_groups[reference_names]))) {
        return(NULL)
    }
    map <- seq_len(n)
    for (signature in reference_names) {
        low <- sort(reference_groups[[signature]], method = "radix")
        high <- sort(target_groups[[signature]], method = "radix")
        map[low] <- high
        map[high] <- low
    }
    as.integer(map)
}

.association_exact_block_maps <- function(swaps, n) {
    count <- length(swaps)
    total <- 2^count
    identity <- seq_len(n)
    lapply(0:(total - 1), function(mask) {
        map <- identity
        selected <- as.logical(intToBits(as.integer(mask))[seq_len(count)])
        for (index in which(selected)) {
            map <- swaps[[index]][map]
        }
        as.integer(map)
    })
}

.association_map_preserves_null <- function(map, z0, tolerance) {
    identical(sort(map, method = "radix"), seq_len(nrow(z0))) &&
        if (ncol(z0)) {
            max(abs(z0[map, , drop = FALSE] - z0)) <= tolerance
        } else {
            TRUE
        }
}

.association_blocked_permutations <- function(
    stratum,
    hypothesis,
    z0,
    enumerate = TRUE
) {
    if (!is.logical(enumerate) || length(enumerate) != 1L ||
        is.na(enumerate)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    design <- stratum$design
    condition_column <- design$roles$condition
    components <- hypothesis$components[[1L]]
    encodings <- hypothesis$component_encodings[[1L]]
    pivot <- sort(
        intersect(components, condition_column),
        method = "radix"
    )
    if (!length(pivot) ||
        !identical(unname(encodings[pivot[[1L]]]), "treatment")) {
        return(.association_robust_failure(
            "association_incompatible_permutation"
        ))
    }
    specifications <- .association_hypothesis_components(
        stratum$core,
        hypothesis
    )
    specification <- specifications[[pivot[[1L]]]]
    if (is.null(specification) || is.na(specification$reference) ||
        is.na(specification$target)) {
        return(.association_robust_failure(
            "association_incompatible_permutation"
        ))
    }
    data <- .design_canonical_data(design)
    x <- stratum$core$model$matrix
    if (!identical(rownames(data), rownames(x)) ||
        !is.matrix(z0) || nrow(z0) != nrow(x)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    condition <- .canonical_design_text(data[[pivot[[1L]]]])
    block <- .canonical_design_text(data[[design$roles$block]])
    signatures <- .association_row_signatures(
        data,
        design$roles$nuisance
    )
    blocks <- sort(unique(block), method = "radix")
    swaps <- lapply(blocks, function(value) {
        .association_block_swap(
            which(block == value),
            condition,
            specification$reference,
            specification$target,
            signatures,
            nrow(x)
        )
    })
    swaps <- swaps[!vapply(swaps, is.null, logical(1L))]
    if (!length(swaps)) {
        return(.association_robust_failure(
            "association_incompatible_permutation"
        ))
    }
    swaps <- .association_canonical_permutation_swaps(swaps, nrow(x))
    tolerance <- .association_permutation_tolerance(z0)
    invariant <- is.finite(tolerance) && all(vapply(swaps, function(map) {
        .association_map_preserves_null(map, z0, tolerance)
    }, logical(1L)))
    coefficient <- match(
        hypothesis$coefficient,
        colnames(stratum$core$model$matrix)
    )
    touched <- !is.na(coefficient) && any(vapply(swaps, function(map) {
        any(x[map, coefficient] != x[, coefficient])
    }, logical(1L)))
    if (!invariant || !touched) {
        return(.association_robust_failure(
            "association_incompatible_permutation"
        ))
    }
    summary <- .association_transformation_summary(
        length(swaps) * log(2),
        2^length(swaps)
    )
    if (!summary$ok) {
        return(summary)
    }
    allowable <- summary$allowable_transformations
    mode <- summary$mode
    maps <- if (enumerate && mode %in% c("low_resolution", "exact")) {
        .association_exact_block_maps(swaps, nrow(x))
    } else {
        NULL
    }
    list(
        ok = TRUE,
        design = "blocked",
        groups = NULL,
        swaps = swaps,
        allowable_transformations = allowable,
        mode = mode,
        maps = maps,
        tolerance = tolerance
    )
}

.association_freedman_lane_restriction <- function(
    fit,
    hypothesis,
    stratum
) {
    if (!is.list(fit) || !isTRUE(fit$ok)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    columns <- colnames(stratum$core$model$matrix)
    column <- match(hypothesis$coefficient, columns)
    if (is.na(column)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    axis <- numeric(length(columns))
    axis[[column]] <- 1
    contrast <- as.double(crossprod(fit$basis, axis))
    restriction <- .association_constrained_null(
        fit$z,
        contrast,
        as.double(stratum$response)
    )
    if (restriction$ok) {
        restriction$contrast <- contrast
    }
    restriction
}

.association_freedman_lane_records <- function(preparation) {
    hypotheses <- preparation$hypotheses
    records <- vector("list", nrow(hypotheses))
    names(records) <- hypotheses$hypothesis
    for (stratum in preparation$strata) {
        selected <- which(hypotheses$stratum == stratum$stratum)
        if (!length(selected)) {
            next
        }
        support <- preparation$support[selected, , drop = FALSE]
        residual_df <- as.integer(
            length(stratum$samples) - stratum$core$rank$rank
        )
        degenerate <- length(unique(unname(stratum$response))) == 1L
        fit <- NULL
        for (position in seq_along(selected)) {
            index <- selected[[position]]
            hypothesis <- hypotheses[index, , drop = FALSE]
            code <- if (!hypothesis$estimable[[1L]]) {
                "association_nonestimable"
            } else if (!support$eligible[[position]]) {
                support$code[[position]]
            } else if (residual_df < 3L) {
                "association_low_residual_df"
            } else if (degenerate) {
                "association_degenerate_response"
            } else {
                NULL
            }
            record <- list(
                stratum = stratum,
                hypothesis = hypothesis,
                code = code,
                fit = NULL,
                observed = NULL,
                restriction = NULL,
                permutations = NULL,
                maps = NULL
            )
            if (!is.null(code)) {
                records[[index]] <- record
                next
            }
            if (is.null(fit)) {
                fit <- .association_robust_algebra(stratum)
            }
            result <- if (fit$ok) {
                .association_robust_contrast(
                    fit,
                    hypothesis$coefficient,
                    colnames(stratum$core$model$matrix)
                )
            } else {
                fit
            }
            if (!result$ok) {
                record$code <- result$code
                records[[index]] <- record
                next
            }
            restriction <- .association_freedman_lane_restriction(
                fit,
                hypothesis,
                stratum
            )
            if (!restriction$ok) {
                record$code <- "association_numerical_failure"
                records[[index]] <- record
                next
            }
            permutations <- if (identical(
                support$design[[position]],
                "independent"
            )) {
                .association_independent_permutations(
                    restriction$z0,
                    enumerate = FALSE
                )
            } else {
                .association_blocked_permutations(
                    stratum,
                    hypothesis,
                    restriction$z0,
                    enumerate = FALSE
                )
            }
            if (!permutations$ok) {
                record$code <- permutations$code
            } else if (identical(
                permutations$mode,
                "low_resolution"
            )) {
                record$code <- "association_low_permutation_resolution"
            }
            record$fit <- fit
            record$observed <- result
            record$restriction <- restriction
            record$permutations <- permutations
            records[[index]] <- record
        }
    }
    records
}

.association_freedman_lane_instances <- function(records) {
    selected <- which(vapply(records, function(record) {
        is.null(record$code) &&
            identical(record$permutations$mode, "monte_carlo")
    }, logical(1L)))
    if (!length(selected)) {
        return(data.frame(
            candidate_id = character(),
            acquisition = character(),
            hypothesis_id = character(),
            stringsAsFactors = FALSE
        ))
    }
    data.frame(
        candidate_id = rep(
            .ASSOCIATION_FREEDMAN_LANE_CANDIDATE,
            length(selected)
        ),
        acquisition = vapply(records[selected], function(record) {
            record$stratum$acquisition
        }, character(1L)),
        hypothesis_id = vapply(records[selected], function(record) {
            record$hypothesis$hypothesis
        }, character(1L)),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

.association_freedman_lane_seeds <- function(records, input_sha256) {
    instances <- .association_freedman_lane_instances(records)
    if (!nrow(instances)) {
        return(.empty_association_seed_manifest()[
            .ASSOCIATION_SEED_MANIFEST_FIELDS[seq_len(7L)]
        ])
    }
    .association_seed_manifest(
        input_sha256,
        instances,
        seq_len(.ASSOCIATION_MONTE_CARLO_DRAWS)
    )
}

.association_freedman_lane_generate_maps <- function(record, seeds) {
    rows <- seeds$acquisition == record$stratum$acquisition &
        seeds$hypothesis_id == record$hypothesis$hypothesis
    generated <- tryCatch(
        .association_monte_carlo_maps(
            seeds[rows, , drop = FALSE],
            record$permutations$design,
            as.integer(length(record$stratum$response)),
            groups = record$permutations$groups,
            swaps = record$permutations$swaps
        ),
        imputefinder_association_permutation_retry_error = function(error) {
            NULL
        }
    )
    if (is.null(generated)) {
        return(list(
            ok = FALSE,
            code = "association_numerical_failure"
        ))
    }
    c(list(ok = TRUE), generated)
}

.association_bind_freedman_lane_seed_manifests <- function(rows) {
    rows <- rows[!vapply(rows, is.null, logical(1L))]
    if (!length(rows)) {
        return(.empty_association_seed_manifest())
    }
    output <- do.call(rbind, unname(rows))
    ordered <- do.call(order, c(
        unname(output[c(
            "candidate_id", "acquisition", "hypothesis_id", "draw_id"
        )]),
        list(method = "radix")
    ))
    output <- output[ordered, , drop = FALSE]
    row.names(output) <- NULL
    output[.ASSOCIATION_SEED_MANIFEST_FIELDS]
}

.association_freedman_lane_exact_maps <- function(record) {
    if (identical(record$permutations$design, "independent")) {
        .association_exact_group_maps(
            record$permutations$groups,
            length(record$stratum$response)
        )
    } else {
        .association_exact_block_maps(
            record$permutations$swaps,
            length(record$stratum$response)
        )
    }
}

.association_byte_identical_numeric <- function(left, right) {
    is.double(left) && is.double(right) &&
        identical(left, right, num.eq = FALSE)
}

.association_freedman_lane_evaluate <- function(record) {
    maps <- if (identical(record$permutations$mode, "exact")) {
        if (is.null(record$permutations$maps)) {
            .association_freedman_lane_exact_maps(record)
        } else {
            record$permutations$maps
        }
    } else {
        record$maps
    }
    expected <- if (identical(record$permutations$mode, "exact")) {
        as.integer(record$permutations$allowable_transformations)
    } else {
        .ASSOCIATION_MONTE_CARLO_DRAWS
    }
    if (!is.list(maps) || length(maps) != expected) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    observed_wald <- (record$observed$effect /
        record$observed$standard_error)^2
    if (!is.finite(observed_wald)) {
        return(.association_robust_failure("association_numerical_failure"))
    }
    exceedance <- 0L
    columns <- colnames(record$stratum$core$model$matrix)
    for (map in maps) {
        invariant <- .association_map_preserves_null(
            map,
            record$restriction$z0,
            record$restriction$tolerance
        )
        if (!invariant) {
            return(.association_robust_failure(
                "association_numerical_failure"
            ))
        }
        transformed_residual <- record$restriction$residual[map]
        transformed_response <- if (.association_byte_identical_numeric(
            transformed_residual,
            record$restriction$residual
        )) {
            as.double(record$stratum$response)
        } else {
            record$restriction$fitted + transformed_residual
        }
        fit <- .association_robust_refit(
            record$fit,
            as.double(transformed_response)
        )
        transformed <- if (fit$ok) {
            .association_robust_contrast(
                fit,
                record$hypothesis$coefficient,
                columns
            )
        } else {
            fit
        }
        if (!transformed$ok) {
            return(.association_robust_failure(
                "association_numerical_failure"
            ))
        }
        transformed_wald <- (transformed$effect /
            transformed$standard_error)^2
        if (!is.finite(transformed_wald)) {
            return(.association_robust_failure(
                "association_numerical_failure"
            ))
        }
        if (transformed_wald >= observed_wald) {
            exceedance <- exceedance + 1L
        }
    }
    exact <- identical(record$permutations$mode, "exact")
    raw_p <- if (exact) {
        exceedance / expected
    } else {
        (1 + exceedance) / 10000
    }
    list(
        ok = TRUE,
        result = record$observed,
        statistic = as.double(observed_wald),
        raw_p = as.double(raw_p),
        permutation_count = if (exact) expected else 10000L,
        allowable_transformations =
            record$permutations$allowable_transformations,
        evaluated_transformations = as.integer(expected),
        exceedance_count = exceedance,
        permutation_mode = if (exact) "exact" else "monte_carlo"
    )
}

.new_association_freedman_lane_outcome <- function(result, hypothesis) {
    observed <- result$result
    diagnostics <- observed$diagnostics
    diagnostics$allowable_transformations <-
        result$allowable_transformations
    diagnostics$evaluated_transformations <-
        result$evaluated_transformations
    diagnostics$exceedance_count <- result$exceedance_count
    diagnostics$permutation_mode <- result$permutation_mode
    structure(
        list(
            status = "available",
            quantity = hypothesis$hypothesis,
            candidate = .ASSOCIATION_FREEDMAN_LANE_CANDIDATE,
            stratum = hypothesis$stratum,
            hypothesis = hypothesis$hypothesis,
            effect = observed$effect,
            standard_error = observed$standard_error,
            conf_low = observed$conf_low,
            conf_high = observed$conf_high,
            statistic = result$statistic,
            reference_df = observed$reference_df,
            raw_p = result$raw_p,
            adjusted_p = NA_real_,
            flag = NA,
            log_odds = NA_real_,
            log_odds_standard_error = NA_real_,
            log_odds_conf_low = NA_real_,
            log_odds_conf_high = NA_real_,
            permutation_count = result$permutation_count,
            diagnostics = diagnostics
        ),
        class = "imputefinder_association"
    )
}

.association_freedman_lane_run_record <- function(record, seeds) {
    unavailable <- function(code) {
        list(
            outcome = .new_association_candidate_unavailable(
                record$hypothesis$hypothesis,
                code
            ),
            seed_manifest = NULL
        )
    }
    if (!is.null(record$code)) {
        return(unavailable(record$code))
    }
    generated <- NULL
    if (identical(record$permutations$mode, "monte_carlo")) {
        generated <- .association_freedman_lane_generate_maps(record, seeds)
        if (!generated$ok) {
            return(unavailable(generated$code))
        }
        record$maps <- generated$maps
    }
    result <- .association_freedman_lane_evaluate(record)
    if (!result$ok) {
        return(unavailable(result$code))
    }
    list(
        outcome = .new_association_freedman_lane_outcome(
            result,
            record$hypothesis
        ),
        seed_manifest = if (is.null(generated)) {
            NULL
        } else {
            generated$manifest
        }
    )
}

.run_association_freedman_lane <- function(preparation) {
    .validate_association_preparation(preparation)
    hypotheses <- preparation$hypotheses
    records <- .association_freedman_lane_records(preparation)
    seeds <- .association_freedman_lane_seeds(
        records,
        preparation$input_sha256
    )
    evaluated <- lapply(
        records,
        .association_freedman_lane_run_record,
        seeds = seeds
    )
    outcomes <- lapply(evaluated, `[[`, "outcome")
    seed_manifests <- lapply(evaluated, `[[`, "seed_manifest")
    names(outcomes) <- hypotheses$hypothesis
    available <- unname(vapply(outcomes, function(outcome) {
        identical(class(outcome), "imputefinder_association")
    }, logical(1L)))
    for (stratum in unique(hypotheses$stratum)) {
        selected <- which(hypotheses$stratum == stratum & available)
        if (!length(selected)) {
            next
        }
        raw <- vapply(outcomes[selected], `[[`, numeric(1L), "raw_p")
        adjusted <- stats::p.adjust(raw, method = "holm")
        for (position in seq_along(selected)) {
            index <- selected[[position]]
            outcomes[[index]]$adjusted_p <- unname(adjusted[[position]])
            outcomes[[index]]$flag <-
                outcomes[[index]]$adjusted_p <= 0.05
        }
    }
    hypotheses$family_member <- available
    multiplicity <- .association_candidate_multiplicity(
        preparation$response,
        hypotheses,
        available
    )
    warnings <- sort(unique(unname(vapply(
        outcomes[!available],
        `[[`,
        character(1L),
        "code"
    ))), method = "radix")
    diagnostics <- list(
        input_sha256 = preparation$input_sha256,
        strata = .association_candidate_stratum_diagnostics(
            preparation,
            hypotheses,
            available
        ),
        seed_manifest = .association_bind_freedman_lane_seed_manifests(
            seed_manifests
        ),
        warnings = warnings
    )
    artifact <- structure(
        list(
            schema = .ASSOCIATION_CANDIDATE_ARTIFACT_SCHEMA,
            protocol = preparation$protocol,
            candidate = .ASSOCIATION_FREEDMAN_LANE_CANDIDATE,
            input_sha256 = preparation$input_sha256,
            response = preparation$response,
            hypotheses = hypotheses,
            support = preparation$support,
            outcomes = outcomes,
            multiplicity = multiplicity,
            diagnostics = diagnostics
        ),
        class = "imputefinder_association_candidate_artifact"
    )
    .validate_association_candidate_artifact(artifact, preparation)
    artifact
}
