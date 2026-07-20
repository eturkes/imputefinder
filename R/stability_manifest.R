.STABILITY_MANIFEST_SCHEMA <- "robustness_perturbation_manifest_v1"
.STABILITY_PROTOCOL_ID <- "m12_b_perturbation_protocol_v2"
.STABILITY_PROTOCOL_HASH <- paste0(
    "20912c4dfadeebb0c1da7100cb54dfa27",
    "f202ea39c3e0d991fb0ee38bab383c2"
)
.STABILITY_SEED_STREAM <- "m12c_b_perturbation_seed_stream_v2"
.STABILITY_INPUT_HASH <- "stability_matrix_design_clusters_be_v1"
.STABILITY_RNG <- "Mersenne-Twister/Inversion/Rejection"
.STABILITY_BOOTSTRAP_DRAWS <- 999L
.STABILITY_PERTURBATION_IDS <- c(
    sampling_bootstrap = "b_sampling_unit_bootstrap",
    sampling_leave_one_out = "b_sampling_leave_one_unit_out",
    estimator_bootstrap = "b_estimator_feature_bootstrap",
    assumption = "b_assumption_policy_panel"
)
.STABILITY_PERTURBATION_FAMILIES <- c(
    "sampling", "sampling", "estimator", "assumption"
)
.STABILITY_REPLACEMENT <- c(
    "with_replacement", "without_replacement", "with_replacement", "none"
)
.STABILITY_CUTOFF_POLICIES <- c(
    "baseline_fixed_v1", "baseline_fixed_v1",
    "automatic_reestimated_v1", "named_policy_v1"
)
.STABILITY_WEIGHTING <- c(
    "integer_multiplicity_v1", "binary_omission_v1",
    "integer_multiplicity_v1", "unit_weight_v1"
)
.STABILITY_FAILURE_POLICIES <- c(
    "unavailable_unconditional_denominator_v1",
    "unavailable_per_unit_v1",
    "cutoff_failure_by_class_v1",
    "unavailable_per_policy_v1"
)
.STABILITY_POLICY_IDS <- c(
    "b_policy_v1_baseline",
    "b_policy_v1_reestimated",
    "b_policy_v1_cutoff_sweep",
    "b_policy_v1_median",
    "b_policy_v1_trimmed",
    "b_policy_v1_half_or_more"
)
.STABILITY_DRAW_FIELDS <- c(
    "perturbation_id", "family", "instance_id", "draw_id", "seed",
    "nonce", "source_kind", "source_set", "weights", "policy_id"
)

.abort_stability_manifest <- function(message, ...) {
    .abort_sidecar(
        message,
        "imputefinder_stability_manifest_error",
        ...
    )
}

.stability_be_i32 <- function(x) {
    valid <- is.numeric(x) && !anyNA(x) && all(is.finite(x)) &&
        all(x >= 0) && all(x <= .Machine$integer.max) &&
        all(x == floor(x))
    if (!valid) {
        .abort_stability_manifest(
            "Stability identity lengths must fit signed 32-bit integers.",
            field = "identity"
        )
    }
    writeBin(as.integer(x), raw(), size = 4L, endian = "big")
}

.stability_encode_text <- function(value) {
    value <- .canonical_sidecar_names(as.character(value))
    if (length(value) != 1L) {
        .abort_stability_manifest(
            "Stability identity text records must be scalar.",
            field = "identity"
        )
    }
    tag <- if (identical(Encoding(value), "bytes")) {
        as.raw(2L)
    } else {
        as.raw(1L)
    }
    bytes <- charToRaw(value)
    c(tag, .stability_be_i32(length(bytes)), bytes)
}

.stability_encode_text_vector <- function(values) {
    values <- .canonical_sidecar_names(as.character(values))
    c(
        .stability_be_i32(length(values)),
        unlist(lapply(values, .stability_encode_text), use.names = FALSE)
    )
}

.stability_encode_roles <- function(roles) {
    roles <- roles[.MISSINGNESS_DESIGN_ROLES]
    c(
        .stability_be_i32(length(roles)),
        unlist(lapply(names(roles), function(role) {
            c(
                .stability_encode_text(role),
                .stability_encode_text_vector(roles[[role]])
            )
        }), use.names = FALSE)
    )
}

.stability_encode_interactions <- function(interactions) {
    c(
        .stability_be_i32(length(interactions)),
        unlist(
            lapply(interactions, .stability_encode_text_vector),
            use.names = FALSE
        )
    )
}

.stability_encode_design_column <- function(column, values) {
    numeric <- !is.factor(values) &&
        (is.integer(values) || is.double(values))
    if (numeric) {
        values <- as.double(values)
        values[values == 0] <- 0
        payload <- c(
            as.raw(1L),
            .stability_be_i32(length(values)),
            writeBin(values, raw(), size = 8L, endian = "big")
        )
    } else {
        payload <- c(
            as.raw(2L),
            .stability_encode_text_vector(values)
        )
    }
    c(.stability_encode_text(column), payload)
}

.stability_design_bytes <- function(design) {
    .validate_missingness_design(design)
    data <- .design_canonical_data(design)
    roles <- design$roles[.MISSINGNESS_DESIGN_ROLES]
    declared <- unname(unlist(roles, use.names = FALSE))
    columns <- unlist(lapply(declared, function(column) {
        .stability_encode_design_column(column, data[[column]])
    }), use.names = FALSE)

    c(
        charToRaw("imputefinder:stability-design:be:v1"),
        as.raw(0L),
        as.raw(1L), .stability_encode_text_vector(rownames(data)),
        as.raw(2L), .stability_encode_roles(roles),
        as.raw(3L), .stability_encode_interactions(design$interactions),
        as.raw(4L), .stability_be_i32(length(declared)), columns
    )
}

.stability_protocol <- function() {
    c(
        id = .STABILITY_PROTOCOL_ID,
        hash = .STABILITY_PROTOCOL_HASH,
        seed_stream = .STABILITY_SEED_STREAM,
        input_hash = .STABILITY_INPUT_HASH,
        rng = .STABILITY_RNG
    )
}

.stability_safe_tokens <- function(x) {
    is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x)) &&
        !anyDuplicated(x) &&
        !any(grepl("|", x, fixed = TRUE)) &&
        !any(grepl("\r", x, fixed = TRUE)) &&
        !any(grepl("\n", x, fixed = TRUE))
}

.stability_source_ids <- function(x) {
    is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x)) &&
        !anyDuplicated(x) &&
        identical(x, sort(x, method = "radix"))
}

.validate_stability_seed_request <- function(
    input_sha256,
    perturbation_id,
    instance_ids,
    draw_ids
) {
    valid <- is.character(input_sha256) && length(input_sha256) == 1L &&
        !is.na(input_sha256) && grepl("^[0-9a-f]{64}$", input_sha256) &&
        is.character(perturbation_id) && length(perturbation_id) == 1L &&
        !is.na(perturbation_id) &&
        perturbation_id %in% .STABILITY_PERTURBATION_IDS[c(
            "sampling_bootstrap", "estimator_bootstrap"
        )] &&
        .stability_safe_tokens(instance_ids) &&
        is.integer(draw_ids) && length(draw_ids) > 0L &&
        !anyNA(draw_ids) && all(draw_ids > 0L) && !anyDuplicated(draw_ids)
    if (!valid) {
        .abort_stability_manifest(
            "Stability seed-manifest input is malformed.",
            field = "seeds"
        )
    }
    invisible(TRUE)
}

.stability_allocate_seed <- function(
    input_sha256,
    perturbation_id,
    instance_id,
    draw_id,
    used
) {
    nonce <- 0L
    repeat {
        key <- paste(
            .STABILITY_PROTOCOL_ID, input_sha256, perturbation_id,
            instance_id, draw_id, nonce,
            sep = "|"
        )
        digest <- unname(tools::sha256sum(
            bytes = charToRaw(enc2utf8(key))
        ))
        seed <- as.integer(strtoi(substr(digest, 1L, 7L), 16L) + 1L)
        seed_key <- as.character(seed)
        if (!exists(seed_key, envir = used, inherits = FALSE)) {
            assign(seed_key, TRUE, envir = used)
            return(c(seed = seed, nonce = nonce))
        }
        nonce <- nonce + 1L
    }
}

.stability_seed_manifest <- function(
    input_sha256,
    perturbation_id,
    instance_ids,
    draw_ids
) {
    .validate_stability_seed_request(
        input_sha256,
        perturbation_id,
        instance_ids,
        draw_ids
    )

    grid <- expand.grid(
        instance_id = sort(unname(instance_ids), method = "radix"),
        draw_id = sort(unname(draw_ids), method = "radix"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    used <- new.env(hash = TRUE, parent = emptyenv())
    allocated <- vapply(seq_len(nrow(grid)), function(index) {
        .stability_allocate_seed(
            input_sha256,
            perturbation_id,
            grid$instance_id[[index]],
            grid$draw_id[[index]],
            used
        )
    }, numeric(2L))

    data.frame(
        protocol_id = rep(.STABILITY_PROTOCOL_ID, nrow(grid)),
        perturbation_id = rep(perturbation_id, nrow(grid)),
        instance_id = grid$instance_id,
        draw_id = as.integer(grid$draw_id),
        seed = as.integer(allocated["seed", ]),
        nonce = as.integer(allocated["nonce", ]),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

.stability_prefixed_hash <- function(prefix, tag, values) {
    bytes <- c(
        charToRaw(tag),
        as.raw(0L),
        .stability_encode_text_vector(values)
    )
    paste0(prefix, unname(tools::sha256sum(bytes = bytes)))
}

.stability_cluster_labels <- function(feature_names, feature_clusters) {
    if (is.null(feature_clusters)) {
        return(sort(feature_names, method = "radix"))
    }
    valid <- is.atomic(feature_clusters) && is.null(dim(feature_clusters)) &&
        length(feature_clusters) == length(feature_names) &&
        !is.null(names(feature_clusters)) && !anyNA(names(feature_clusters)) &&
        all(nzchar(names(feature_clusters))) &&
        !anyDuplicated(names(feature_clusters))
    cluster_names <- if (valid) {
        .canonical_sidecar_names(names(feature_clusters))
    } else {
        character()
    }
    valid <- valid && !anyDuplicated(cluster_names) &&
        setequal(cluster_names, feature_names)
    labels <- if (valid) {
        .canonical_sidecar_names(as.character(feature_clusters))
    } else {
        character()
    }
    if (!valid || length(labels) != length(feature_names) ||
        anyNA(labels) || any(!nzchar(labels))) {
        .abort_stability_manifest(
            paste0(
                "`feature_clusters` must be a complete named atomic ",
                "vector aligned exactly to the input features."
            ),
            field = "feature_clusters"
        )
    }
    names(labels) <- cluster_names
    unname(labels[match(sort(feature_names, method = "radix"), cluster_names)])
}

.stability_cluster_ids <- function(labels) {
    unique_labels <- sort(unique(labels), method = "radix")
    stats::setNames(vapply(unique_labels, function(label) {
        .stability_prefixed_hash(
            "fc_",
            "imputefinder:stability-feature-cluster:v1",
            label
        )
    }, character(1L)), unique_labels)
}

.normalise_stability_feature_clusters <- function(data, feature_clusters) {
    feature_names <- .canonical_sidecar_names(rownames(data))
    ordered_features <- sort(feature_names, method = "radix")
    labels <- .stability_cluster_labels(feature_names, feature_clusters)
    cluster_ids <- .stability_cluster_ids(labels)

    data.frame(
        feature = ordered_features,
        cluster = unname(cluster_ids[labels]),
        cluster_label = labels,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

.stability_instance_id <- function(role, conditions, acquisitions) {
    values <- c(
        paste0("role=", role),
        paste0("condition=", conditions),
        paste0("acquisition=", acquisitions)
    )
    .stability_prefixed_hash(
        "si_",
        "imputefinder:stability-sampling-instance:v1",
        values
    )
}

.new_stability_sample_sources <- function(design) {
    core <- .new_design_estimability(design)
    source <- core$units$sample[c("sample", "condition", "unit")]
    acquisition_role <- design$roles$acquisition
    if (length(acquisition_role)) {
        canonical_data <- .design_canonical_data(design)
        source$acquisition <- .canonical_design_text(
            canonical_data[[acquisition_role]]
        )
    } else {
        source$acquisition <- rep(NA_character_, nrow(source))
    }
    source <- source[c("sample", "condition", "acquisition", "unit")]

    units <- sort(unique(source$unit), method = "radix")
    unit_instance <- stats::setNames(vapply(units, function(unit) {
        selected <- source$unit == unit
        conditions <- sort(unique(source$condition[selected]), method = "radix")
        acquisitions <- if (length(acquisition_role)) {
            sort(unique(source$acquisition[selected]), method = "radix")
        } else {
            character()
        }
        .stability_instance_id(
            core$units$resampling_role,
            conditions,
            acquisitions
        )
    }, character(1L)), units)
    source$instance <- unname(unit_instance[source$unit])
    row.names(source) <- NULL
    source
}

.stability_input_bytes <- function(data, design, feature_sources) {
    feature_names <- .canonical_sidecar_names(rownames(data))
    sample_names <- .canonical_sidecar_names(colnames(data))
    rownames(data) <- feature_names
    colnames(data) <- sample_names
    data <- data[
        order(feature_names, method = "radix"),
        order(sample_names, method = "radix"),
        drop = FALSE
    ]
    design <- .align_missingness_design(design, sample_names)
    matrix_bytes <- .canonical_matrix_bytes(data, .pack_original_mask(data))
    design_bytes <- .stability_design_bytes(design)
    cluster_bytes <- c(
        .stability_encode_text_vector(feature_sources$feature),
        .stability_encode_text_vector(feature_sources$cluster_label)
    )

    c(
        charToRaw("imputefinder:stability-input:be:v1"),
        as.raw(0L),
        as.raw(1L), .stability_be_i32(length(matrix_bytes)), matrix_bytes,
        as.raw(2L), .stability_be_i32(length(design_bytes)), design_bytes,
        as.raw(3L), .stability_be_i32(length(cluster_bytes)), cluster_bytes
    )
}

.stability_input_sha256 <- function(data, design, feature_sources) {
    unname(tools::sha256sum(bytes = .stability_input_bytes(
        data,
        design,
        feature_sources
    )))
}

.stability_draw_frame <- function(
    perturbation_id,
    family,
    instance_id,
    draw_id,
    seed,
    nonce,
    source_kind,
    source_set,
    weights,
    policy_id
) {
    count <- length(draw_id)
    output <- data.frame(
        perturbation_id = rep(perturbation_id, count),
        family = rep(family, count),
        instance_id = unname(instance_id),
        draw_id = as.integer(draw_id),
        seed = as.integer(seed),
        nonce = as.integer(nonce),
        source_kind = rep(source_kind, count),
        source_set = unname(source_set),
        policy_id = unname(policy_id),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    output$weights <- I(weights)
    output[.STABILITY_DRAW_FIELDS]
}

.stability_with_preserved_rng_kind <- function(code) {
    caller_kind <- RNGkind()
    on.exit(do.call(RNGkind, as.list(caller_kind)), add = TRUE)
    force(code)
}

.stability_seeded_weights <- function(seed, source_count) {
    withr::with_preserve_seed(.stability_with_preserved_rng_kind(
        withr::with_seed(
            seed,
            tabulate(
                sample.int(source_count, source_count, replace = TRUE),
                nbins = source_count
            ),
            .rng_kind = "Mersenne-Twister",
            .rng_normal_kind = "Inversion",
            .rng_sample_kind = "Rejection"
        )
    ))
}

.stability_bootstrap_draws <- function(
    seeds,
    source_map,
    family,
    source_kind
) {
    valid <- is.list(source_map) && .stability_safe_tokens(names(source_map)) &&
        setequal(unique(seeds$instance_id), names(source_map)) &&
        all(vapply(source_map, .stability_source_ids, logical(1L)))
    if (!valid) {
        .abort_stability_manifest(
            "Stability bootstrap sources are malformed.",
            field = "sources"
        )
    }

    weights <- vector("list", nrow(seeds))
    for (index in seq_len(nrow(seeds))) {
        source <- source_map[[seeds$instance_id[[index]]]]
        weights[[index]] <- .stability_seeded_weights(
            seeds$seed[[index]],
            length(source)
        )
    }

    .stability_draw_frame(
        perturbation_id = unique(seeds$perturbation_id),
        family = family,
        instance_id = seeds$instance_id,
        draw_id = seeds$draw_id,
        seed = seeds$seed,
        nonce = seeds$nonce,
        source_kind = source_kind,
        source_set = seeds$instance_id,
        weights = weights,
        policy_id = rep(NA_character_, nrow(seeds))
    )
}

.stability_sampling_source_map <- function(sample_sources) {
    instances <- sort(unique(sample_sources$instance), method = "radix")
    output <- lapply(instances, function(instance) {
        sort(
            unique(sample_sources$unit[sample_sources$instance == instance]),
            method = "radix"
        )
    })
    names(output) <- instances
    output
}

.stability_leave_one_out_draws <- function(sample_sources) {
    units <- unique(sample_sources[c("unit", "instance")])
    units <- units[do.call(order, c(
        unname(units[c("instance", "unit")]),
        list(method = "radix")
    )), , drop = FALSE]
    row.names(units) <- NULL
    all_units <- sort(units$unit, method = "radix")
    weights <- lapply(units$unit, function(omitted) {
        as.integer(all_units != omitted)
    })

    .stability_draw_frame(
        perturbation_id = .STABILITY_PERTURBATION_IDS[[
            "sampling_leave_one_out"
        ]],
        family = "sampling",
        instance_id = units$instance,
        draw_id = seq_len(nrow(units)),
        seed = rep(NA_integer_, nrow(units)),
        nonce = rep(NA_integer_, nrow(units)),
        source_kind = "sampling_unit",
        source_set = rep("all_sampling_units", nrow(units)),
        weights = weights,
        policy_id = rep(NA_character_, nrow(units))
    )
}

.stability_estimator_instance <- function(feature_sources) {
    clusters <- sort(unique(feature_sources$cluster), method = "radix")
    .stability_prefixed_hash(
        "ei_",
        "imputefinder:stability-estimator-instance:v1",
        clusters
    )
}

.stability_assumption_draws <- function(feature_sources) {
    features <- feature_sources$feature
    count <- length(.STABILITY_POLICY_IDS)
    .stability_draw_frame(
        perturbation_id = .STABILITY_PERTURBATION_IDS[["assumption"]],
        family = "assumption",
        instance_id = rep("policy_all", count),
        draw_id = seq_len(count),
        seed = rep(NA_integer_, count),
        nonce = rep(NA_integer_, count),
        source_kind = "feature_condition",
        source_set = rep("all_feature_conditions", count),
        weights = rep(list(rep(1L, length(features))), count),
        policy_id = .STABILITY_POLICY_IDS
    )
}

.new_stability_perturbations <- function(sample_sources, resampling_role) {
    unit_count <- length(unique(sample_sources$unit))
    instance_count <- length(unique(sample_sources$instance))
    if (!is.character(resampling_role) || length(resampling_role) != 1L ||
        !resampling_role %in% c("sample", "block")) {
        .abort_stability_manifest(
            "Stability resampling role is malformed.",
            field = "resampling_unit"
        )
    }
    data.frame(
        perturbation_id = unname(.STABILITY_PERTURBATION_IDS),
        family = .STABILITY_PERTURBATION_FAMILIES,
        draw_count = c(
            .STABILITY_BOOTSTRAP_DRAWS,
            0L,
            .STABILITY_BOOTSTRAP_DRAWS,
            length(.STABILITY_POLICY_IDS)
        ),
        realized_draw_count = c(
            .STABILITY_BOOTSTRAP_DRAWS,
            unit_count,
            .STABILITY_BOOTSTRAP_DRAWS,
            length(.STABILITY_POLICY_IDS)
        ),
        instance_count = c(instance_count, instance_count, 1L, 1L),
        replacement = .STABILITY_REPLACEMENT,
        resampling_unit = c(
            resampling_role, resampling_role,
            "feature_cluster", "feature_condition"
        ),
        cutoff_policy = .STABILITY_CUTOFF_POLICIES,
        weighting = .STABILITY_WEIGHTING,
        failure_policy = .STABILITY_FAILURE_POLICIES,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

.new_stability_policies <- function() {
    data.frame(
        policy_id = .STABILITY_POLICY_IDS,
        cutoff = c(
            "fixed_baseline",
            "reestimated_automatic",
            "type8_q025_q25_q50_q75_q975_successful_estimator",
            "fixed_baseline",
            "fixed_baseline",
            "fixed_baseline"
        ),
        summary = c(
            "arithmetic_observed_mean",
            "arithmetic_observed_mean",
            "arithmetic_observed_mean",
            "observed_median",
            "observed_trimmed_mean_0.20",
            "arithmetic_observed_mean"
        ),
        majority = c(
            rep("strictly_more_than_half", 5L),
            "at_least_half"
        ),
        policy_version = rep("rules_v1_sensitivity_v1", 6L),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

.validate_new_stability_manifest <- function(manifest) {
    valid <- is.list(manifest) &&
        identical(
            names(manifest),
            c(
                "schema", "protocol", "input_sha256", "sources",
                "perturbations", "policies", "seeds", "draws"
            )
        ) &&
        identical(manifest$schema, .STABILITY_MANIFEST_SCHEMA) &&
        identical(manifest$protocol, .stability_protocol()) &&
        is.character(manifest$input_sha256) &&
        length(manifest$input_sha256) == 1L &&
        grepl("^[0-9a-f]{64}$", manifest$input_sha256) &&
        identical(names(manifest$sources), c("sample", "feature")) &&
        is.data.frame(manifest$sources$sample) &&
        is.data.frame(manifest$sources$feature) &&
        is.data.frame(manifest$perturbations) &&
        is.data.frame(manifest$policies) &&
        is.data.frame(manifest$seeds) &&
        is.data.frame(manifest$draws) &&
        identical(names(manifest$draws), .STABILITY_DRAW_FIELDS)
    if (!valid) {
        .abort_stability_manifest(
            "Constructed stability perturbation manifest is malformed.",
            field = "manifest"
        )
    }
    invisible(manifest)
}

.new_stability_stream_plan <- function(
    input_sha256,
    sample_sources,
    feature_sources
) {
    sampling_sources <- .stability_sampling_source_map(sample_sources)
    sampling_seeds <- .stability_seed_manifest(
        input_sha256,
        .STABILITY_PERTURBATION_IDS[["sampling_bootstrap"]],
        names(sampling_sources),
        seq_len(.STABILITY_BOOTSTRAP_DRAWS)
    )
    estimator_instance <- .stability_estimator_instance(feature_sources)
    estimator_sources <- list(sort(
        unique(feature_sources$cluster),
        method = "radix"
    ))
    names(estimator_sources) <- estimator_instance
    estimator_seeds <- .stability_seed_manifest(
        input_sha256,
        .STABILITY_PERTURBATION_IDS[["estimator_bootstrap"]],
        estimator_instance,
        seq_len(.STABILITY_BOOTSTRAP_DRAWS)
    )
    seeds <- rbind(sampling_seeds, estimator_seeds)
    row.names(seeds) <- NULL
    list(
        sampling_sources = sampling_sources,
        sampling_seeds = sampling_seeds,
        estimator_sources = estimator_sources,
        estimator_seeds = estimator_seeds,
        seeds = seeds
    )
}

.new_stability_draws <- function(plan, sample_sources, feature_sources) {
    draws <- rbind(
        .stability_bootstrap_draws(
            plan$sampling_seeds,
            plan$sampling_sources,
            family = "sampling",
            source_kind = "sampling_unit"
        ),
        .stability_leave_one_out_draws(sample_sources),
        .stability_bootstrap_draws(
            plan$estimator_seeds,
            plan$estimator_sources,
            family = "estimator",
            source_kind = "feature_cluster"
        ),
        .stability_assumption_draws(feature_sources)
    )
    row.names(draws) <- NULL
    draws
}

.new_stability_perturbation_manifest <- function(
    data,
    design,
    feature_clusters = NULL
) {
    .validate_matrix_data(data)
    design <- .align_missingness_design(design, colnames(data))
    feature_sources <- .normalise_stability_feature_clusters(
        data,
        feature_clusters
    )
    sample_sources <- .new_stability_sample_sources(design)
    input_sha256 <- .stability_input_sha256(
        data,
        design,
        feature_sources
    )
    plan <- .new_stability_stream_plan(
        input_sha256,
        sample_sources,
        feature_sources
    )

    manifest <- structure(
        list(
            schema = .STABILITY_MANIFEST_SCHEMA,
            protocol = .stability_protocol(),
            input_sha256 = input_sha256,
            sources = list(
                sample = sample_sources,
                feature = feature_sources
            ),
            perturbations = .new_stability_perturbations(
                sample_sources,
                if (length(design$roles$block)) "block" else "sample"
            ),
            policies = .new_stability_policies(),
            seeds = plan$seeds,
            draws = .new_stability_draws(
                plan,
                sample_sources,
                feature_sources
            )
        ),
        class = "imputefinder_stability_manifest"
    )
    .validate_new_stability_manifest(manifest)
    manifest
}

.validate_stability_perturbation_manifest <- function(
    manifest,
    data,
    design,
    feature_clusters = NULL
) {
    expected <- .new_stability_perturbation_manifest(
        data,
        design,
        feature_clusters
    )
    if (!identical(manifest, expected)) {
        .abort_stability_manifest(
            paste0(
                "Stored stability perturbation manifest is malformed or ",
                "does not match the supplied input, design, and grouping."
            ),
            field = "manifest"
        )
    }
    invisible(manifest)
}
