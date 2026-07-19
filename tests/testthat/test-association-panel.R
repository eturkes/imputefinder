association_panel_detection_matrix <- function(counts, total) {
    samples <- sprintf("sample_%02d", seq_along(counts))
    output <- matrix(
        NA_real_,
        nrow = total,
        ncol = length(counts),
        dimnames = list(sprintf("feature_%03d", seq_len(total)), samples)
    )
    for (index in seq_along(counts)) {
        positions <- (seq_len(counts[[index]]) + index - 2L) %% total + 1L
        output[positions, index] <- positions
    }
    output
}

association_panel_fixture <- function(candidate = "quasibinomial") {
    data <- association_panel_detection_matrix(
        c(2L, 4L, 5L, 7L, 8L, 10L, 6L, 9L, 11L, 12L, 14L, 16L),
        20L
    )
    metadata <- data.frame(
        condition = rep(c("A", "B"), each = 6L),
        row.names = colnames(data),
        stringsAsFactors = FALSE
    )
    design <- missingness_design(metadata, condition = "condition")
    preparation <- imputefinder:::.new_association_preparation(data, design)
    artifact <- switch(
        candidate,
        robust = imputefinder:::.run_association_ols_hc3_cr2(preparation),
        quasibinomial = imputefinder:::.run_association_quasibinomial(
            preparation
        )
    )
    panel <- structure(
        list(
            schema = "sentinel_association_v1",
            protocol = list(
                id = "m13_a_association_protocol_v3",
                hash = "5f484f614d72d9560d96e2b8f75b71cc1027950958ea1459be93eece595069f5",
                contract_hash = "85d121ebd5ddd589be0f389553a3b7ffa6076a12a8fcbcb476a8969f7d987ae4",
                gate_registry_hash = "f2fd8d0171003e0ae9b31757acf2d43acf7477965d517f1757f4d895c6abe731",
                implementation_hash = strrep("1", 64L),
                candidate_evidence_hash = strrep("2", 64L),
                winner_state = "winner_locked",
                m12_contract_hash = "45d1cda936b21a42d6735d900917734398eecb3ff1c136cab486cf90ee2b21e5",
                m12_gate_registry_hash = "70d101b890de687be8973c9ff263caa454b61505f99470e5ec6061cb344db92f"
            ),
            candidate = artifact$candidate,
            methods = c(
                response = "a_sample_detection_fraction_v3",
                encoding = "canonical_treatment_contrasts_v1",
                rank = "svd_relative_v1",
                candidate = artifact$candidate,
                multiplicity = "Holm"
            ),
            response = artifact$response,
            hypotheses = artifact$hypotheses,
            support = artifact$support,
            outcomes = artifact$outcomes,
            multiplicity = artifact$multiplicity,
            diagnostics = artifact$diagnostics
        ),
        class = "imputefinder_association_panel"
    )
    list(
        data = data,
        design = design,
        preparation = preparation,
        artifact = artifact,
        panel = panel
    )
}

test_that("frozen public-panel shape validates without a materialization seam", {
    fixture <- association_panel_fixture()
    panel <- fixture$panel

    expect_identical(class(panel), "imputefinder_association_panel")
    expect_identical(
        names(panel),
        c(
            "schema", "protocol", "candidate", "methods", "response",
            "hypotheses", "support", "outcomes", "multiplicity",
            "diagnostics"
        )
    )
    expect_identical(panel$schema, "sentinel_association_v1")
    expect_false("input_sha256" %in% names(panel))
    expect_invisible(imputefinder:::.validate_association_panel_shape(
        panel,
        fixture$preparation
    ))

    malformed <- list(
        { x <- panel; x$protocol$candidate_evidence_hash <- "detached"; x },
        { x <- panel; x$methods[["candidate"]] <- "a_fraction_freedman_lane"; x },
        { x <- panel; x$outcomes[[1L]]$quantity <- "wrong"; x },
        fixture$artifact
    )
    for (candidate in malformed) {
        expect_error(
            imputefinder:::.validate_association_panel_shape(
                candidate,
                fixture$preparation
            ),
            class = "imputefinder_association_panel_error"
        )
    }
})

test_that("panel validation replays robust numerical provenance", {
    fixture <- association_panel_fixture("robust")
    altered <- fixture$panel
    outcome <- altered$outcomes[[1L]]
    outcome$effect <- outcome$effect + 0.01
    outcome$statistic <- outcome$effect / outcome$standard_error
    outcome$raw_p <- 2 * stats::pt(-abs(outcome$statistic), outcome$reference_df)
    outcome$adjusted_p <- outcome$raw_p
    outcome$flag <- outcome$adjusted_p <= 0.05
    critical <- stats::qt(0.975, outcome$reference_df)
    outcome$conf_low <- outcome$effect - critical * outcome$standard_error
    outcome$conf_high <- outcome$effect + critical * outcome$standard_error
    altered$outcomes[[1L]] <- outcome

    expect_error(
        imputefinder:::.validate_association_panel_shape(
            altered,
            fixture$preparation
        ),
        class = "imputefinder_association_panel_error"
    )
})

test_that("pre-winner sentinel lifecycle rejects an association slot", {
    fixture <- association_panel_fixture()
    input <- imputefinder:::.new_sidecar_input(fixture$data)
    record <- imputefinder:::.new_sidecar_design_record(fixture$design, input)
    condition <- fixture$design$sample_data[[fixture$design$roles$condition]]
    groups <- stats::setNames(
        condition,
        rownames(fixture$design$sample_data)
    )
    sentinel <- imputefinder:::.new_sidecar_sentinel(
        fixture$data,
        groups,
        fixture$design
    )
    sentinel$association <- fixture$panel

    expect_error(
        imputefinder:::.validate_sidecar_sentinel(sentinel, input, record),
        class = "imputefinder_analysis_schema_error"
    )
})
