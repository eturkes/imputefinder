stability_manifest_fixture <- function() {
    metadata <- data.frame(
        condition = c("B", "A", "B", "A", "B", "A"),
        row.names = paste0("s", c(6L, 1L, 5L, 2L, 4L, 3L)),
        stringsAsFactors = FALSE
    )
    design <- missingness_design(metadata, condition = "condition")
    x <- matrix(
        c(
            10, 11, NA, 13, 14, 15,
            20, NA, 22, 23, 24, 25,
            30, 31, 32, NA, 34, 35,
            40, 41, 42, 43, NA, 45
        ),
        nrow = 4L,
        byrow = TRUE,
        dimnames = list(
            c("f4", "f1", "f3", "f2"),
            rownames(metadata)
        )
    )
    feature_clusters <- c(
        f1 = "module_1",
        f2 = "module_1",
        f3 = "module_2",
        f4 = "module_3"
    )

    list(
        x = x,
        metadata = metadata,
        design = design,
        feature_clusters = feature_clusters
    )
}

stability_block_fixture <- function() {
    metadata <- expand.grid(
        technical = c("t1", "t2"),
        condition = c("A", "B"),
        subject = paste0("p", seq_len(3L)),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    metadata$run <- seq_len(nrow(metadata))
    rownames(metadata) <- sprintf("sample_%02d", seq_len(nrow(metadata)))
    design <- missingness_design(
        metadata,
        condition = "condition",
        nuisance = c("technical", "run"),
        block = "subject"
    )
    x <- matrix(
        seq_len(3L * nrow(metadata)),
        nrow = 3L,
        dimnames = list(paste0("feature_", 1:3), rownames(metadata))
    )

    list(x = x, metadata = metadata, design = design)
}

test_that("stability seed manifests port the frozen collision-safe stream", {
    seeds <- imputefinder:::.stability_seed_manifest(
        paste(rep("0", 64L), collapse = ""),
        "b_sampling_unit_bootstrap",
        c("condition_b", "condition_a"),
        as.integer(c(3, 1, 2))
    )

    expect_identical(
        names(seeds),
        c(
            "protocol_id", "perturbation_id", "instance_id", "draw_id",
            "seed", "nonce"
        )
    )
    expect_identical(
        seeds$instance_id,
        rep(c("condition_a", "condition_b"), 3L)
    )
    expect_identical(seeds$draw_id, rep(1:3, each = 2L))
    expect_identical(
        seeds$seed,
        c(
            209599270L, 21962589L, 111781779L,
            59443635L, 71248539L, 244028140L
        )
    )
    expect_identical(seeds$nonce, rep(0L, 6L))
    expect_identical(anyDuplicated(seeds$seed), 0L)

    expect_error(
        imputefinder:::.stability_seed_manifest(
            paste(rep("0", 64L), collapse = ""),
            "b_sampling_unit_bootstrap",
            "unsafe|instance",
            1L
        ),
        class = "imputefinder_stability_manifest_error"
    )
})

test_that("independent manifests retain weighted named perturbations", {
    fixture <- stability_manifest_fixture()
    manifest <- imputefinder:::.new_stability_perturbation_manifest(
        fixture$x,
        fixture$design,
        fixture$feature_clusters
    )

    expect_s3_class(manifest, "imputefinder_stability_manifest")
    expect_identical(
        names(manifest),
        c(
            "schema", "protocol", "input_sha256", "sources",
            "perturbations", "policies", "seeds", "draws"
        )
    )
    expect_identical(manifest$schema, "robustness_perturbation_manifest_v1")
    expect_identical(
        manifest$protocol,
        c(
            id = "m12_b_perturbation_protocol_v2",
            hash = paste0(
                "20912c4dfadeebb0c1da7100cb54dfa27",
                "f202ea39c3e0d991fb0ee38bab383c2"
            ),
            seed_stream = "m12c_b_perturbation_seed_stream_v2",
            input_hash = "stability_matrix_design_clusters_be_v1",
            rng = "Mersenne-Twister/Inversion/Rejection"
        )
    )
    expect_match(manifest$input_sha256, "^[0-9a-f]{64}$")
    expect_identical(
        manifest$input_sha256,
        "fb03be68befb24470d22bea68b1020ee1b297425d7e2d5f57ffb7685efffba0d"
    )
    expect_identical(
        manifest$seeds$seed[1:4],
        c(196830613L, 167208946L, 146704639L, 88452635L)
    )
    expect_identical(names(manifest$sources), c("sample", "feature"))
    expect_identical(
        names(manifest$sources$sample),
        c("sample", "condition", "acquisition", "unit", "instance")
    )
    expect_identical(manifest$sources$sample$sample, paste0("s", 1:6))
    expect_identical(
        names(manifest$sources$feature),
        c("feature", "cluster", "cluster_label")
    )
    expect_identical(manifest$sources$feature$feature, paste0("f", 1:4))
    expect_identical(
        manifest$sources$feature$cluster[
            manifest$sources$feature$feature %in% c("f1", "f2")
        ],
        rep(manifest$sources$feature$cluster[[1L]], 2L)
    )

    expect_identical(
        manifest$perturbations$perturbation_id,
        c(
            "b_sampling_unit_bootstrap",
            "b_sampling_leave_one_unit_out",
            "b_estimator_feature_bootstrap",
            "b_assumption_policy_panel"
        )
    )
    expect_identical(
        manifest$perturbations$draw_count,
        c(999L, 0L, 999L, 6L)
    )
    expect_identical(
        manifest$perturbations$realized_draw_count,
        c(999L, 6L, 999L, 6L)
    )
    expect_identical(
        manifest$perturbations$replacement,
        c("with_replacement", "without_replacement", "with_replacement", "none")
    )
    expect_identical(
        manifest$perturbations$resampling_unit,
        c("sample", "sample", "feature_cluster", "feature_condition")
    )
    expect_identical(
        manifest$perturbations$cutoff_policy,
        c(
            "baseline_fixed_v1", "baseline_fixed_v1",
            "automatic_reestimated_v1", "named_policy_v1"
        )
    )
    expect_identical(
        manifest$perturbations$weighting,
        c(
            "integer_multiplicity_v1", "binary_omission_v1",
            "integer_multiplicity_v1", "unit_weight_v1"
        )
    )
    expect_identical(
        manifest$perturbations$failure_policy,
        c(
            "unavailable_unconditional_denominator_v1",
            "unavailable_per_unit_v1",
            "cutoff_failure_by_class_v1",
            "unavailable_per_policy_v1"
        )
    )
    expect_identical(
        names(manifest$policies),
        c("policy_id", "cutoff", "summary", "majority", "policy_version")
    )
    expect_identical(
        manifest$policies$policy_id,
        imputefinder:::.STABILITY_POLICY_IDS
    )
    expect_identical(
        manifest$policies$majority,
        c(rep("strictly_more_than_half", 5L), "at_least_half")
    )

    sample_source <- manifest$sources$sample
    expect_length(unique(sample_source$instance), 2L)
    expect_true(all(vapply(
        split(sample_source$condition, sample_source$instance),
        function(value) length(unique(value)) == 1L,
        logical(1L)
    )))

    sampling <- manifest$draws[
        manifest$draws$perturbation_id == "b_sampling_unit_bootstrap",
        ,
        drop = FALSE
    ]
    expect_identical(nrow(sampling), 2L * 999L)
    expect_true(all(vapply(seq_len(nrow(sampling)), function(index) {
        length(sampling$weights[[index]]) == 3L &&
            sum(sampling$weights[[index]]) == 3L
    }, logical(1L))))
    expect_identical(sampling$source_set, sampling$instance_id)
    expect_identical(
        unclass(sampling$weights[1:4]),
        list(c(1L, 0L, 2L), c(1L, 0L, 2L), c(2L, 0L, 1L), 2:0)
    )

    leave_one_out <- manifest$draws[
        manifest$draws$perturbation_id ==
            "b_sampling_leave_one_unit_out",
        ,
        drop = FALSE
    ]
    expect_identical(nrow(leave_one_out), 6L)
    expect_true(all(vapply(leave_one_out$weights, function(weight) {
        identical(sort(weight), c(0L, rep(1L, 5L)))
    }, logical(1L))))

    estimator <- manifest$draws[
        manifest$draws$perturbation_id ==
            "b_estimator_feature_bootstrap",
        ,
        drop = FALSE
    ]
    expect_identical(nrow(estimator), 999L)
    expect_true(all(vapply(estimator$weights, sum, integer(1L)) == 3L))

    assumption <- manifest$draws[
        manifest$draws$perturbation_id == "b_assumption_policy_panel",
        ,
        drop = FALSE
    ]
    expect_identical(
        assumption$policy_id,
        c(
            "b_policy_v1_baseline", "b_policy_v1_reestimated",
            "b_policy_v1_cutoff_sweep", "b_policy_v1_median",
            "b_policy_v1_trimmed", "b_policy_v1_half_or_more"
        )
    )
    expect_true(all(vapply(assumption$weights, function(weight) {
        all(weight == 1L)
    }, logical(1L))))
    expect_invisible(imputefinder:::.validate_stability_perturbation_manifest(
        manifest,
        fixture$x,
        fixture$design,
        fixture$feature_clusters
    ))
    expect_identical(
        unserialize(serialize(manifest, NULL, version = 3L)),
        manifest
    )
})

test_that("blocked manifests keep subjects and siblings inseparable", {
    fixture <- stability_block_fixture()
    manifest <- imputefinder:::.new_stability_perturbation_manifest(
        fixture$x,
        fixture$design
    )
    samples <- manifest$sources$sample

    expect_identical(unique(samples$unit), paste0("p", 1:3))
    expect_length(unique(samples$instance), 1L)
    expect_true(all(vapply(split(samples$sample, samples$unit), length,
                             integer(1L)) == 4L))
    expect_true(all(vapply(split(samples$condition, samples$unit), function(x) {
        identical(sort(unique(x)), c("A", "B"))
    }, logical(1L))))

    sampling <- manifest$draws[
        manifest$draws$perturbation_id == "b_sampling_unit_bootstrap",
        ,
        drop = FALSE
    ]
    expect_identical(nrow(sampling), 999L)
    expect_identical(sampling$source_set, sampling$instance_id)
    expect_true(all(vapply(sampling$weights, length, integer(1L)) == 3L))
    expect_true(all(vapply(sampling$weights, sum, integer(1L)) == 3L))

    leave_one_out <- manifest$draws[
        manifest$draws$perturbation_id ==
            "b_sampling_leave_one_unit_out",
        ,
        drop = FALSE
    ]
    expect_identical(nrow(leave_one_out), 3L)
    expect_true(all(vapply(leave_one_out$weights, function(weight) {
        sum(weight == 0L) == 1L && sum(weight) == 2L
    }, logical(1L))))

    sample_order <- rev(colnames(fixture$x))
    reencoded <- fixture$metadata[sample_order, , drop = FALSE]
    reencoded$condition <- factor(
        reencoded$condition,
        levels = rev(unique(reencoded$condition))
    )
    reencoded$technical <- factor(
        reencoded$technical,
        levels = rev(unique(reencoded$technical))
    )
    reencoded$subject <- factor(
        reencoded$subject,
        levels = rev(unique(reencoded$subject))
    )
    reencoded$run <- as.double(reencoded$run)
    recoded_design <- missingness_design(
        reencoded,
        condition = "condition",
        nuisance = c("run", "technical"),
        block = "subject"
    )
    recoded <- imputefinder:::.new_stability_perturbation_manifest(
        fixture$x[rev(rownames(fixture$x)), sample_order, drop = FALSE],
        recoded_design
    )
    expect_identical(recoded, manifest)
})

test_that("source labels remain unrestricted while stream tokens stay safe", {
    sample_names <- c("sample|one", "sample\ntwo")
    feature_names <- c("feature|one", "feature\ntwo")
    metadata <- data.frame(
        condition = c("A|left", "B\nright"),
        row.names = sample_names,
        stringsAsFactors = FALSE
    )
    design <- missingness_design(metadata, condition = "condition")
    x <- matrix(
        seq_len(4L),
        nrow = 2L,
        dimnames = list(feature_names, sample_names)
    )
    clusters <- stats::setNames(c("module|x", "module\ny"), feature_names)

    manifest <- imputefinder:::.new_stability_perturbation_manifest(
        x,
        design,
        clusters
    )

    expect_identical(
        sort(manifest$sources$sample$sample, method = "radix"),
        sort(sample_names, method = "radix")
    )
    expect_identical(
        sort(manifest$sources$feature$feature, method = "radix"),
        sort(feature_names, method = "radix")
    )
    expect_false(any(grepl("|", manifest$sources$sample$instance,
                           fixed = TRUE)))
    expect_false(any(grepl("\n", manifest$sources$sample$instance,
                           fixed = TRUE)))
})

test_that("manifest identity and draws are named-order invariant", {
    fixture <- stability_manifest_fixture()
    baseline <- imputefinder:::.new_stability_perturbation_manifest(
        fixture$x,
        fixture$design,
        fixture$feature_clusters
    )
    feature_order <- rev(rownames(fixture$x))
    sample_order <- rev(colnames(fixture$x))
    permuted_x <- fixture$x[feature_order, sample_order, drop = FALSE]
    permuted_metadata <- fixture$metadata[sample_order, , drop = FALSE]
    permuted_design <- missingness_design(
        permuted_metadata,
        condition = "condition"
    )
    permuted_clusters <- factor(
        fixture$feature_clusters[feature_order],
        levels = rev(unique(fixture$feature_clusters))
    )
    names(permuted_clusters) <- feature_order
    permuted <- imputefinder:::.new_stability_perturbation_manifest(
        permuted_x,
        permuted_design,
        permuted_clusters
    )

    expect_identical(permuted, baseline)

    changed_x <- fixture$x
    changed_x[!is.na(changed_x)][[1L]] <-
        changed_x[!is.na(changed_x)][[1L]] + 0.25
    changed_input <- imputefinder:::.new_stability_perturbation_manifest(
        changed_x,
        fixture$design,
        fixture$feature_clusters
    )
    expect_false(identical(changed_input$input_sha256, baseline$input_sha256))
    expect_false(identical(changed_input$seeds$seed, baseline$seeds$seed))

    changed_clusters <- fixture$feature_clusters
    changed_clusters[["f3"]] <- "module_1"
    regrouped <- imputefinder:::.new_stability_perturbation_manifest(
        fixture$x,
        fixture$design,
        changed_clusters
    )
    expect_false(identical(regrouped$input_sha256, baseline$input_sha256))
    expect_false(identical(regrouped$sources$feature, baseline$sources$feature))

    changed_metadata <- fixture$metadata
    changed_metadata["s6", "condition"] <- "A"
    changed_design <- missingness_design(
        changed_metadata,
        condition = "condition"
    )
    redesigned <- imputefinder:::.new_stability_perturbation_manifest(
        fixture$x,
        changed_design,
        fixture$feature_clusters
    )
    expect_false(identical(redesigned$input_sha256, baseline$input_sha256))
    expect_false(identical(redesigned$sources$sample, baseline$sources$sample))
})

test_that("manifest construction is isolated from fit-only and caller state", {
    fixture <- stability_manifest_fixture()
    original_x <- fixture$x
    original_design <- fixture$design
    original_clusters <- fixture$feature_clusters

    with_preserved_random_state({
        RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
        set.seed(48129L)
        caller_kind <- RNGkind()
        caller_seed <- .Random.seed
        manifest <- imputefinder:::.new_stability_perturbation_manifest(
            fixture$x,
            fixture$design,
            fixture$feature_clusters
        )
        expect_s3_class(manifest, "imputefinder_stability_manifest")
        expect_identical(RNGkind(), caller_kind)
        expect_identical(.Random.seed, caller_seed)
    })
    with_preserved_random_state({
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
        }
        manifest <- imputefinder:::.new_stability_perturbation_manifest(
            fixture$x,
            fixture$design,
            fixture$feature_clusters
        )
        expect_s3_class(manifest, "imputefinder_stability_manifest")
        expect_false(exists(
            ".Random.seed",
            envir = globalenv(),
            inherits = FALSE
        ))
    })

    expect_identical(fixture$x, original_x)
    expect_identical(fixture$design, original_design)
    expect_identical(fixture$feature_clusters, original_clusters)

    fit <- classify_missingness(
        fixture$x,
        fixture$metadata$condition,
        cutoffs = c(A = 12, B = 12)
    )
    unavailable <- imputefinder:::.fit_only_artifact(
        fit,
        "stability_perturbation_manifest"
    )
    expect_s3_class(unavailable, "imputefinder_unavailable")
    expect_identical(unavailable$code, "original_input_required")
    expect_identical(unavailable$requires, c("x", "missingness_design"))

    expect_error(
        imputefinder:::.new_stability_perturbation_manifest(
            fixture$x,
            fixture$design,
            fixture$feature_clusters[-1L]
        ),
        class = "imputefinder_stability_manifest_error"
    )

    manifest <- imputefinder:::.new_stability_perturbation_manifest(
        fixture$x,
        fixture$design,
        fixture$feature_clusters
    )
    tampered <- manifest
    tampered$draws$weights[[1L]][[1L]] <-
        tampered$draws$weights[[1L]][[1L]] + 1L
    expect_error(
        imputefinder:::.validate_stability_perturbation_manifest(
            tampered,
            fixture$x,
            fixture$design,
            fixture$feature_clusters
        ),
        class = "imputefinder_stability_manifest_error"
    )
})
