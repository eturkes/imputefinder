sidecar_differential_fixture <- function() {
    fixture <- normative_fixture()
    fit <- classify_missingness(
        fixture$x,
        fixture$group,
        cutoffs = fixture$cutoffs,
        seed = 19L
    )

    list(fixture = fixture, fit = fit)
}

test_that("the compact differential fixture freezes all comparator categories", {
    candidate <- sidecar_differential_fixture()
    reference <- candidate$fit
    tolerance <- unname(
        imputefinder:::.classic_compatibility_tolerance()[["absolute"]]
    )

    exact <- reference
    exact$classifications$state[[1L]] <- "MAR"
    canonical <- reference
    canonical$profiles$A$metadata$warnings <- "fixture warning"
    within <- reference
    within$profiles$A$grid$intensity[[1L]] <-
        within$profiles$A$grid$intensity[[1L]] + tolerance / 4
    beyond <- reference
    beyond$profiles$A$grid$intensity[[1L]] <-
        beyond$profiles$A$grid$intensity[[1L]] + tolerance * 100

    reports <- lapply(
        list(exact, canonical, within, beyond),
        imputefinder:::.compare_classic_fits,
        reference = reference
    )
    observed <- lapply(
        reports,
        function(report) {
            report[c("compatible", "exact", "canonical", "tolerance")]
        }
    )

    expect_identical(
        observed,
        list(
            list(
                compatible = FALSE,
                exact = "classifications",
                canonical = character(),
                tolerance = character()
            ),
            list(
                compatible = FALSE,
                exact = character(),
                canonical = "profiles",
                tolerance = character()
            ),
            list(
                compatible = TRUE,
                exact = character(),
                canonical = character(),
                tolerance = character()
            ),
            list(
                compatible = FALSE,
                exact = character(),
                canonical = character(),
                tolerance = "profiles"
            )
        )
    )
})

test_that("the exact fixture isolates sidecar work from rules_v1", {
    candidate <- sidecar_differential_fixture()
    design <- missingness_design(
        data.frame(
            condition = candidate$fixture$group,
            row.names = colnames(candidate$fixture$x)
        ),
        condition = "condition"
    )
    before <- serialize(candidate$fit, NULL, version = 3L)

    analysis <- analyze_missingness(
        candidate$fixture$x,
        design,
        base_fit = candidate$fit,
        seed = 19L,
        modules = c("sentinel", "stability")
    )

    expect_identical(serialize(analysis$classic, NULL, version = 3L), before)
    expect_identical(
        serialize(candidate$fit, NULL, version = 3L),
        before
    )
    expect_identical(names(attributes(candidate$fit)), c("names", "class"))
})
