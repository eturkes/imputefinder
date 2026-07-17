test_that("stable scientific identifiers are explicit and versioned", {
    identifiers <- imputefinder:::.stable_scientific_identifiers()

    expect_identical(
        identifiers,
        c(
            method = "rules_v1",
            profile = "count_weighted_density_v1",
            cutoff = "manual_or_derivative_boundary_v1",
            rescue = "condition_minimum_single_cell_v1",
            policy = "four_state_reconciliation_v1"
        )
    )
    expect_true(all(grepl("_v[1-9][0-9]*$", identifiers)))
})

test_that("the minimum sidecar spec is deterministic and self-identifying", {
    spec <- imputefinder:::.new_sidecar_spec()

    expect_identical(
        names(spec),
        c("schema", "lifecycle", "modules", "scientific", "software")
    )
    expect_identical(spec$schema, "imputefinder_analysis_v1")
    expect_identical(spec$lifecycle, "experimental")
    expect_identical(
        spec$modules,
        c(
            sentinel = "design_confounding_sentinel_v1",
            stability = "robustness_certificate_v1"
        )
    )
    expect_identical(
        spec$scientific,
        imputefinder:::.stable_scientific_identifiers()
    )
    expect_identical(
        spec$software,
        c(
            package = "imputefinder",
            version = unname(as.character(
                getNamespaceVersion("imputefinder")
            ))
        )
    )
    expect_identical(imputefinder:::.new_sidecar_spec(), spec)
})

test_that("sidecar metadata leaves the stable result exactly unchanged", {
    before <- classify_normative_fixture()
    before_bytes <- serialize(before, NULL, version = 3L)

    imputefinder:::.new_sidecar_spec()
    after <- classify_normative_fixture()

    expect_identical(serialize(after, NULL, version = 3L), before_bytes)
    expect_identical(class(after), "imputefinder_result")
    expect_identical(
        names(after),
        c(
            "data",
            "classifications",
            "groups",
            "feature_status",
            "cutoffs",
            "cutoff_diagnostics",
            "profiles",
            "seed_log",
            "groups_by_sample",
            "call"
        )
    )
    expect_identical(names(attributes(after)), c("names", "class"))
})
