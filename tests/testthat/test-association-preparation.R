association_factorial_fixture <- function() {
    metadata <- expand.grid(
        replicate = seq_len(4L),
        time = c(0, 1),
        run = c(0L, 1L),
        condition = c("A", "B"),
        acquisition = c("DDA", "DIA"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    rownames(metadata) <- sprintf("sample_%03d", seq_len(nrow(metadata)))
    design <- missingness_design(
        metadata,
        condition = "condition",
        nuisance = c("time", "run"),
        acquisition = "acquisition",
        interactions = list(c("time", "condition", "run"))
    )
    x <- matrix(
        1,
        nrow = 4L,
        ncol = nrow(metadata),
        dimnames = list(
            c("always", "run_one", "dda_only", "absent"),
            rownames(metadata)
        )
    )
    x["run_one", metadata$run == 0L] <- NA_real_
    x["dda_only", metadata$acquisition == "DIA"] <- NA_real_
    x["absent", ] <- NA_real_
    sample_order <- c(seq(64L, 2L, by = -2L), seq(63L, 1L, by = -2L))

    list(
        x = x[c("dda_only", "absent", "always", "run_one"), sample_order,
            drop = FALSE
        ],
        metadata = metadata,
        design = design
    )
}

association_support_fixture <- function(unit_count, blocked = FALSE) {
    if (blocked) {
        metadata <- expand.grid(
            technical = c("t1", "t2"),
            condition = c("A", "B"),
            subject = sprintf("p%02d", seq_len(unit_count)),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
    } else {
        metadata <- data.frame(
            condition = rep(c("A", "B"), each = unit_count),
            stringsAsFactors = FALSE
        )
    }
    rownames(metadata) <- sprintf("sample_%03d", seq_len(nrow(metadata)))
    design <- missingness_design(
        metadata,
        condition = "condition",
        block = if (blocked) "subject" else NULL
    )
    x <- matrix(
        1,
        nrow = 1L,
        ncol = nrow(metadata),
        dimnames = list("feature", rownames(metadata))
    )
    list(x = x, design = design)
}

association_hash_fixture <- function() {
    mask <- matrix(
        c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
        nrow = 2L,
        byrow = TRUE,
        dimnames = list(
            c("feature_b", "feature_a"),
            c("s4", "s2", "s1", "s3")
        )
    )
    sample_data <- data.frame(
        condition = factor(c("B", "A", "A", "B"), levels = c("B", "A")),
        batch = c("two", "one", "two", "one"),
        run = c(4L, 2L, 1L, 3L),
        acquisition = rep("DDA", 4L),
        row.names = colnames(mask),
        stringsAsFactors = FALSE
    )
    roles <- list(
        condition = "condition",
        nuisance = c("run", "batch"),
        block = character(),
        acquisition = "acquisition"
    )
    interactions <- list(
        c("condition", "run"),
        c("batch", "condition")
    )
    list(
        mask = mask,
        sample_data = sample_data,
        roles = roles,
        interactions = interactions
    )
}

test_that("association identities implement the frozen byte protocol", {
    fixture <- association_hash_fixture()
    expected <- "a8f005269eecbf507ffb5c0b4a1385f6ee60446e2c646201edbfaa6428386c0e"
    hash <- imputefinder:::.association_mask_design_sha256(
        fixture$mask,
        fixture$sample_data,
        fixture$roles,
        fixture$interactions
    )

    expect_identical(hash, expected)
    expect_identical(
        imputefinder:::.association_stratum_id("DDA"),
        paste0(
            "s_",
            "23419cf41eaaf424a742472631d1666f3abbb750c230854f095cc546565920d5"
        )
    )
    expect_identical(
        imputefinder:::.association_hypothesis_id(
            "all", "coef_0002", "condition[B]", "term_0002"
        ),
        paste0(
            "a_",
            "9bab8d1d33b4f63ad0c4d232a4373fe9e6c972f0510883601f876b6748ec3ef5"
        )
    )

    permuted <- fixture$mask[c(2L, 1L), c(3L, 1L, 4L, 2L), drop = FALSE]
    recoded <- fixture$sample_data[colnames(permuted), , drop = FALSE]
    recoded$condition <- as.character(recoded$condition)
    recoded$run <- as.double(recoded$run)
    roles <- fixture$roles[c("acquisition", "block", "nuisance", "condition")]
    roles$nuisance <- rev(roles$nuisance)
    expect_identical(
        imputefinder:::.association_mask_design_sha256(
            permuted,
            recoded,
            roles,
            lapply(rev(fixture$interactions), rev)
        ),
        expected
    )

    changed <- fixture$mask
    changed[[1L]] <- !changed[[1L]]
    expect_false(identical(
        imputefinder:::.association_mask_design_sha256(
            changed,
            fixture$sample_data,
            fixture$roles,
            fixture$interactions
        ),
        expected
    ))

    changed_design <- fixture$sample_data
    changed_design$run[[1L]] <- changed_design$run[[1L]] + 1L
    changed_roles <- fixture$roles
    changed_roles$nuisance <- c("run", "acquisition")
    changed_roles$acquisition <- "batch"
    expect_false(identical(
        imputefinder:::.association_mask_design_sha256(
            fixture$mask,
            changed_design,
            fixture$roles,
            fixture$interactions
        ),
        expected
    ))
    expect_false(identical(
        imputefinder:::.association_mask_design_sha256(
            fixture$mask,
            fixture$sample_data,
            changed_roles,
            fixture$interactions
        ),
        expected
    ))
    expect_false(identical(
        imputefinder:::.association_mask_design_sha256(
            fixture$mask,
            fixture$sample_data,
            fixture$roles,
            fixture$interactions[-1L]
        ),
        expected
    ))

    x <- matrix(
        1,
        nrow = nrow(fixture$mask),
        ncol = ncol(fixture$mask),
        dimnames = dimnames(fixture$mask)
    )
    x[fixture$mask] <- NA_real_
    design <- missingness_design(
        fixture$sample_data,
        condition = "condition",
        nuisance = c("run", "batch"),
        acquisition = "acquisition",
        interactions = fixture$interactions
    )
    expect_identical(
        imputefinder:::.association_input_sha256(x, design),
        expected
    )
    intensity_changed <- x
    intensity_changed[!is.na(intensity_changed)] <- seq_len(sum(!is.na(x)))
    expect_identical(
        imputefinder:::.association_input_sha256(intensity_changed, design),
        expected
    )
})

test_that("preparation rebuilds acquisition strata from the global response", {
    fixture <- association_factorial_fixture()
    preparation <- imputefinder:::.new_association_preparation(
        fixture$x,
        fixture$design
    )

    expect_s3_class(preparation, "imputefinder_association_preparation")
    expect_identical(
        names(preparation),
        c(
            "schema", "protocol", "input_sha256", "identity", "response",
            "hypotheses", "support", "strata"
        )
    )
    expect_identical(preparation$schema, "association_preparation_v1")
    expect_identical(
        preparation$protocol,
        c(
            id = "m13_a_association_protocol_v3",
            hash = "5f484f614d72d9560d96e2b8f75b71cc1027950958ea1459be93eece595069f5",
            contract_hash = "85d121ebd5ddd589be0f389553a3b7ffa6076a12a8fcbcb476a8969f7d987ae4",
            response = "a_sample_detection_fraction_v3",
            encoding = "canonical_treatment_contrasts_v1",
            rank = "svd_relative_v1",
            input_hash = "association_mask_design_be_v1"
        )
    )
    expect_identical(
        names(preparation$identity),
        c(
            "dimensions", "feature_names", "sample_names", "original_mask",
            "design"
        )
    )
    expect_identical(
        names(preparation$response),
        c(
            "stratum", "acquisition", "sample",
            "globally_observable_count", "detected_count",
            "detection_fraction"
        )
    )
    expect_identical(unique(preparation$response$acquisition), c("DDA", "DIA"))
    expect_identical(
        preparation$response$sample,
        unlist(lapply(c("DDA", "DIA"), function(acquisition) {
            sort(
                rownames(fixture$metadata)[
                    fixture$metadata$acquisition == acquisition
                ],
                method = "radix"
            )
        }), use.names = FALSE)
    )
    expect_identical(
        preparation$response$globally_observable_count,
        rep(3L, nrow(fixture$metadata))
    )
    response_rows <- match(
        preparation$response$sample,
        rownames(fixture$metadata)
    )
    expected_detected <- 1L +
        as.integer(fixture$metadata$run[response_rows] == 1L) +
        as.integer(
            fixture$metadata$acquisition[response_rows] == "DDA"
        )
    expect_identical(preparation$response$detected_count, expected_detected)
    expect_identical(
        preparation$response$detection_fraction,
        expected_detected / 3
    )

    expect_identical(
        names(preparation$hypotheses),
        c(
            "hypothesis", "stratum", "coefficient", "label", "term_id",
            "term", "kind", "components", "component_encodings", "role",
            "eligible", "estimable"
        )
    )
    expect_identical(
        preparation$hypotheses$label,
        rep(c("condition[B]", "run", "time", "condition[B]:run:time"), 2L)
    )
    expect_true(all(preparation$hypotheses$eligible))
    expect_true(all(preparation$hypotheses$estimable))
    interaction <- which(preparation$hypotheses$kind == "interaction")[[1L]]
    expect_identical(
        preparation$hypotheses$components[[interaction]],
        c("condition", "run", "time")
    )
    expect_identical(
        preparation$hypotheses$component_encodings[[interaction]],
        c(condition = "treatment", run = "numeric", time = "numeric")
    )

    numeric_support <- preparation$support$numeric_components[[interaction]]
    expect_identical(
        numeric_support,
        data.frame(
            component = c("run", "time"),
            median = c(0.5, 0.5),
            minimum = c(0, 0),
            maximum = c(1, 1),
            extrapolates_0_1 = c(FALSE, FALSE),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(preparation$support$side_count_low[interaction], 16L)
    expect_identical(preparation$support$side_count_high[interaction], 16L)
    expect_identical(preparation$support$cell_count_min[interaction], 4L)
    expect_true(all(preparation$support$eligible))
    expect_identical(names(preparation$strata), unique(preparation$response$stratum))
    expect_true(all(vapply(preparation$strata, function(stratum) {
        identical(rownames(stratum$core$model$matrix), stratum$samples) &&
            identical(names(stratum$response), stratum$samples) &&
            length(stratum$design$roles$acquisition) == 0L
    }, logical(1L))))
    expect_invisible(
        imputefinder:::.validate_association_preparation(preparation)
    )
})

test_that("preparation validation rejects detached identity and joins", {
    fixture <- association_factorial_fixture()
    preparation <- imputefinder:::.new_association_preparation(
        fixture$x,
        fixture$design
    )

    detached <- preparation
    first <- substr(detached$input_sha256, 1L, 1L)
    detached$input_sha256 <- paste0(
        if (identical(first, "0")) "1" else "0",
        substr(detached$input_sha256, 2L, 64L)
    )
    changed_mask <- preparation
    byte <- as.integer(changed_mask$identity$original_mask$bytes[[1L]])
    changed_mask$identity$original_mask$bytes[[1L]] <- as.raw(
        bitwXor(byte, 1L)
    )
    changed_design <- preparation
    changed_design$identity$design$sample_data$run[[1L]] <-
        changed_design$identity$design$sample_data$run[[1L]] + 1L

    omitted_sample <- preparation
    stratum <- omitted_sample$strata[[1L]]
    keep <- stratum$samples[-1L]
    design <- stratum$design
    reduced_design <- missingness_design(
        design$sample_data[keep, , drop = FALSE],
        condition = design$roles$condition,
        nuisance = design$roles$nuisance,
        block = design$roles$block,
        interactions = design$interactions
    )
    stratum$design <- reduced_design
    stratum$core <- imputefinder:::.new_design_estimability(reduced_design)
    stratum$samples <- rownames(stratum$core$model$matrix)
    stratum$response <- stratum$response[stratum$samples]
    omitted_sample$strata[[1L]] <- stratum
    selected <- omitted_sample$hypotheses$stratum == stratum$stratum
    hypotheses <- imputefinder:::.association_hypotheses(
        stratum$core,
        stratum$stratum
    )
    omitted_sample$hypotheses[selected, ] <- hypotheses
    omitted_sample$support[selected, ] <- imputefinder:::.association_support(
        stratum$design,
        stratum$core,
        hypotheses
    )

    relabeled <- preparation
    old <- names(relabeled$strata)[[1L]]
    replacement <- paste0("s_", paste(rep("0", 64L), collapse = ""))
    response_rows <- relabeled$response$stratum == old
    hypothesis_rows <- relabeled$hypotheses$stratum == old
    relabeled$response$stratum[response_rows] <- replacement
    relabeled$strata[[1L]]$stratum <- replacement
    names(relabeled$strata)[[1L]] <- replacement
    relabeled$hypotheses$stratum[hypothesis_rows] <- replacement
    relabeled$hypotheses$hypothesis[hypothesis_rows] <- vapply(
        which(hypothesis_rows),
        function(index) {
            imputefinder:::.association_hypothesis_id(
                replacement,
                relabeled$hypotheses$coefficient[[index]],
                relabeled$hypotheses$label[[index]],
                relabeled$hypotheses$term_id[[index]]
            )
        },
        character(1L)
    )
    relabeled$support$hypothesis[hypothesis_rows] <-
        relabeled$hypotheses$hypothesis[hypothesis_rows]

    missing_hypothesis <- preparation
    missing_hypothesis$hypotheses <-
        missing_hypothesis$hypotheses[-1L, , drop = FALSE]
    missing_hypothesis$support <-
        missing_hypothesis$support[-1L, , drop = FALSE]
    corrupted_support <- preparation
    corrupted_support$support$side_count_low[[1L]] <-
        corrupted_support$support$side_count_low[[1L]] + 1L
    empty_available <- preparation
    empty_available$hypotheses <-
        empty_available$hypotheses[FALSE, , drop = FALSE]
    empty_available$support <- empty_available$support[FALSE, , drop = FALSE]
    swapped_strata <- preparation
    swapped_strata$strata <- swapped_strata$strata[c(2L, 1L)]
    names(swapped_strata$strata) <- names(preparation$strata)

    for (corrupted in list(
        detached,
        changed_mask,
        changed_design,
        omitted_sample,
        relabeled,
        missing_hypothesis,
        corrupted_support,
        empty_available,
        swapped_strata
    )) {
        expect_error(
            imputefinder:::.validate_association_preparation(corrupted),
            class = "imputefinder_association_preparation_error"
        )
    }
})

test_that("support floors count biological units rather than siblings", {
    independent_pass <- association_support_fixture(4L)
    independent_fail <- association_support_fixture(3L)
    blocked_pass <- association_support_fixture(6L, blocked = TRUE)
    blocked_fail <- association_support_fixture(5L, blocked = TRUE)

    support <- lapply(
        list(independent_pass, independent_fail, blocked_pass, blocked_fail),
        function(fixture) {
            imputefinder:::.new_association_preparation(
                fixture$x,
                fixture$design
            )$support[1L, , drop = FALSE]
        }
    )
    expect_identical(support[[1L]]$side_count_low, 4L)
    expect_identical(support[[1L]]$side_count_high, 4L)
    expect_true(support[[1L]]$eligible)
    expect_identical(support[[2L]]$code, "association_low_independent_support")

    expect_identical(support[[3L]]$side_count_low, 6L)
    expect_identical(support[[3L]]$side_count_high, 6L)
    expect_identical(support[[3L]]$complete_block_count, 6L)
    expect_true(support[[3L]]$eligible)
    expect_identical(support[[4L]]$complete_block_count, 5L)
    expect_identical(support[[4L]]$code, "association_low_block_support")
})

test_that("categorical main support follows each encoded target level", {
    metadata <- data.frame(
        condition = c(rep("A", 4L), rep("B", 4L), "C"),
        row.names = sprintf("sample_%02d", seq_len(9L)),
        stringsAsFactors = FALSE
    )
    design <- missingness_design(metadata, condition = "condition")
    x <- matrix(
        1,
        nrow = 1L,
        ncol = nrow(metadata),
        dimnames = list("feature", rownames(metadata))
    )
    preparation <- imputefinder:::.new_association_preparation(x, design)

    expect_identical(
        preparation$hypotheses$label,
        c("condition[B]", "condition[C]")
    )
    expect_identical(
        preparation$support$side_definition,
        c("condition[A]/condition[B]", "condition[A]/condition[C]")
    )
    expect_identical(preparation$support$side_count_low, c(4L, 4L))
    expect_identical(preparation$support$side_count_high, c(4L, 1L))
    expect_identical(preparation$support$eligible, c(TRUE, FALSE))
    expect_identical(
        preparation$support$code,
        c(NA_character_, "association_low_independent_support")
    )
})

test_that("blocked numeric support weights unique positions and excludes ties", {
    positions <- expand.grid(
        dose = c(0.2, 0.5, 0.8),
        condition = c("A", "B"),
        subject = sprintf("p%02d", seq_len(6L)),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    siblings <- ifelse(positions$dose == 0.2, 5L, 1L)
    metadata <- positions[rep(seq_len(nrow(positions)), siblings), , drop = FALSE]
    rownames(metadata) <- sprintf("sample_%03d", seq_len(nrow(metadata)))
    design <- missingness_design(
        metadata,
        condition = "condition",
        nuisance = "dose",
        block = "subject"
    )
    x <- matrix(
        1,
        nrow = 1L,
        ncol = nrow(metadata),
        dimnames = list("feature", rownames(metadata))
    )
    preparation <- imputefinder:::.new_association_preparation(x, design)
    dose <- which(preparation$hypotheses$label == "dose")

    expect_length(dose, 1L)
    expect_identical(
        preparation$support$numeric_components[[dose]],
        data.frame(
            component = "dose",
            median = 0.5,
            minimum = 0.2,
            maximum = 0.8,
            extrapolates_0_1 = TRUE,
            stringsAsFactors = FALSE
        )
    )
    expect_identical(preparation$support$side_count_low[dose], 6L)
    expect_identical(preparation$support$side_count_high[dose], 6L)
    expect_identical(preparation$support$complete_block_count[dose], 6L)
    expect_true(preparation$support$eligible[dose])
})

test_that("blocked interaction support requires six complete blocks", {
    complete <- expand.grid(
        technical = c("t1", "t2"),
        batch = c("x", "y"),
        condition = c("A", "B"),
        subject = sprintf("p%02d", seq_len(6L)),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    incomplete <- expand.grid(
        technical = c("t1", "t2"),
        replicate = seq_len(6L),
        cell = seq_len(4L),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    incomplete$batch <- c("x", "y", "x", "y")[incomplete$cell]
    incomplete$condition <- c("A", "A", "B", "B")[incomplete$cell]
    incomplete$subject <- sprintf(
        "cell_%d_p%02d",
        incomplete$cell,
        incomplete$replicate
    )

    prepare <- function(metadata) {
        rownames(metadata) <- sprintf("sample_%03d", seq_len(nrow(metadata)))
        design <- missingness_design(
            metadata,
            condition = "condition",
            nuisance = "batch",
            block = "subject",
            interactions = list(c("condition", "batch"))
        )
        x <- matrix(
            1,
            nrow = 1L,
            ncol = nrow(metadata),
            dimnames = list("feature", rownames(metadata))
        )
        preparation <- imputefinder:::.new_association_preparation(x, design)
        interaction <- which(preparation$hypotheses$kind == "interaction")
        preparation$support[interaction, , drop = FALSE]
    }
    pass <- prepare(complete)
    fail <- prepare(incomplete)

    expect_identical(pass$cell_count_min, 6L)
    expect_identical(pass$complete_block_count, 6L)
    expect_true(pass$eligible)
    expect_identical(fail$cell_count_min, 6L)
    expect_identical(fail$complete_block_count, 0L)
    expect_false(fail$eligible)
    expect_identical(fail$code, "association_low_interaction_support")
})

test_that("eligible aliased axes remain while block and acquisition terms do not", {
    metadata <- expand.grid(
        condition = c("A", "B"),
        subject = sprintf("p%02d", seq_len(4L)),
        acquisition = c("DDA", "DIA"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    metadata$batch <- ifelse(metadata$condition == "A", "x", "y")
    rownames(metadata) <- sprintf("sample_%03d", seq_len(nrow(metadata)))
    design <- missingness_design(
        metadata,
        condition = "condition",
        nuisance = "batch",
        block = "subject",
        acquisition = "acquisition",
        interactions = list(
            c("condition", "batch"),
            c("condition", "subject"),
            c("condition", "acquisition")
        )
    )
    x <- matrix(
        1,
        nrow = 1L,
        ncol = nrow(metadata),
        dimnames = list("feature", rownames(metadata))
    )
    preparation <- imputefinder:::.new_association_preparation(x, design)

    expect_identical(
        preparation$hypotheses$term,
        rep(c("condition", "batch", "batch:condition"), 2L)
    )
    expect_true(all(!preparation$hypotheses$estimable))
    expect_false(any(vapply(
        preparation$hypotheses$components,
        function(components) any(components %in% c("subject", "acquisition")),
        logical(1L)
    )))
    expect_true(all(vapply(
        preparation$strata,
        function(stratum) {
            "condition:subject" %in% stratum$core$model$terms$term &&
                !"acquisition" %in% unlist(
                    stratum$core$model$terms$components,
                    use.names = FALSE
                )
        },
        logical(1L)
    )))
})

test_that("preparation abstains before fitting on empty endpoints", {
    metadata <- data.frame(
        condition = c("A", "B"),
        row.names = c("s1", "s2"),
        stringsAsFactors = FALSE
    )
    design <- missingness_design(metadata, condition = "condition")
    absent <- matrix(
        NA_real_,
        nrow = 2L,
        ncol = 2L,
        dimnames = list(c("f1", "f2"), rownames(metadata))
    )
    unavailable <- imputefinder:::.new_association_preparation(absent, design)
    expect_identical(class(unavailable), "imputefinder_unavailable")
    expect_identical(names(unavailable), imputefinder:::.UNAVAILABLE_FIELDS)
    expect_identical(
        unavailable,
        structure(
            list(
                status = "unavailable",
                quantity = "association",
                code = "association_no_observable_features",
                message = "No feature is observed in the supplied input.",
                requires = "one globally observable feature"
            ),
            class = "imputefinder_unavailable"
        )
    )
    expect_invisible(imputefinder:::.validate_unavailable(unavailable))

    one_level <- metadata
    one_level$condition <- "A"
    design <- missingness_design(one_level, condition = "condition")
    observed <- absent
    observed[] <- 1
    unavailable <- imputefinder:::.new_association_preparation(observed, design)
    expect_identical(class(unavailable), "imputefinder_unavailable")
    expect_identical(names(unavailable), imputefinder:::.UNAVAILABLE_FIELDS)
    expect_identical(
        unavailable,
        structure(
            list(
                status = "unavailable",
                quantity = "association",
                code = "association_no_testable_hypotheses",
                message = paste0(
                    "The rebuilt design has no eligible association ",
                    "coefficient."
                ),
                requires = paste0(
                    "one eligible condition, nuisance, or interaction ",
                    "coefficient"
                )
            ),
            class = "imputefinder_unavailable"
        )
    )
    expect_invisible(imputefinder:::.validate_unavailable(unavailable))
})
