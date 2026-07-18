design_metadata_fixture <- function() {
    data.frame(
        condition = factor(
            c("treated", "control", "treated", "control"),
            levels = c("treated", "control")
        ),
        batch = factor(c("b2", "b1", "b1", "b2")),
        run_order = c(4, 1, 3, 2),
        subject = c(2L, 1L, 2L, 1L),
        acquisition_mode = c("DIA", "DIA", "DIA", "DIA"),
        ignored = letters[1:4],
        row.names = c("s4", "s1", "s3", "s2")
    )
}

test_that("missingness_design freezes the minimum public typed design", {
    sample_data <- design_metadata_fixture()
    original <- sample_data

    design <- missingness_design(
        sample_data,
        condition = "condition",
        nuisance = c("run_order", "batch"),
        block = "subject",
        acquisition = "acquisition_mode",
        interactions = list(
            c("run_order", "condition"),
            c("batch", "condition")
        )
    )

    expect_identical(
        names(formals(missingness_design)),
        c(
            "sample_data",
            "condition",
            "nuisance",
            "block",
            "acquisition",
            "interactions"
        )
    )
    expect_true("missingness_design" %in% getNamespaceExports("imputefinder"))
    expect_s3_class(design, "missingness_design")
    expect_identical(
        names(design),
        c("spec", "sample_data", "roles", "interactions")
    )
    expect_identical(
        design$spec,
        list(
            schema = "missingness_design_v1",
            lifecycle = "experimental"
        )
    )
    expect_identical(
        design$roles,
        list(
            condition = "condition",
            nuisance = c("batch", "run_order"),
            block = "subject",
            acquisition = "acquisition_mode"
        )
    )
    expect_identical(
        names(design$sample_data),
        c(
            "condition",
            "batch",
            "run_order",
            "subject",
            "acquisition_mode"
        )
    )
    expect_identical(rownames(design$sample_data), rownames(sample_data))
    expect_identical(
        design$sample_data$condition,
        as.character(sample_data$condition)
    )
    expect_identical(
        design$sample_data$subject,
        as.character(sample_data$subject)
    )
    expect_identical(
        design$sample_data$acquisition_mode,
        sample_data$acquisition_mode
    )
    expect_identical(design$sample_data$batch, sample_data$batch)
    expect_identical(
        design$interactions,
        list(
            c("batch", "condition"),
            c("condition", "run_order")
        )
    )
    expect_identical(sample_data, original)
})

test_that("optional roles are explicit and structural support is deferred", {
    sample_data <- data.frame(
        condition = "only_condition",
        row.names = "sample_1"
    )

    design <- missingness_design(sample_data, condition = "condition")

    expect_identical(
        design$roles,
        list(
            condition = "condition",
            nuisance = character(),
            block = character(),
            acquisition = character()
        )
    )
    expect_identical(design$interactions, list())
    expect_identical(design$sample_data, sample_data)
    expect_invisible(imputefinder:::.validate_missingness_design(design))
})

test_that("typed designs align exactly without mutating their source", {
    design <- missingness_design(
        design_metadata_fixture(),
        condition = "condition",
        nuisance = "batch",
        block = "subject"
    )
    original <- design

    aligned <- imputefinder:::.align_missingness_design(
        design,
        c("s1", "s2", "s3", "s4")
    )

    expect_s3_class(aligned, "missingness_design")
    expect_identical(
        rownames(aligned$sample_data),
        c("s1", "s2", "s3", "s4")
    )
    expect_identical(
        aligned$sample_data$condition,
        c("control", "control", "treated", "treated")
    )
    expect_identical(aligned$spec, design$spec)
    expect_identical(aligned$roles, design$roles)
    expect_identical(design, original)
})

test_that("alignment failures carry portable typed context", {
    design <- missingness_design(
        design_metadata_fixture(),
        condition = "condition"
    )

    error <- tryCatch(
        imputefinder:::.align_missingness_design(
            design,
            c("s1", "s2", "s3", "other")
        ),
        imputefinder_design_alignment_error = identity
    )

    expect_s3_class(error, "imputefinder_design_alignment_error")
    expect_s3_class(error, "imputefinder_design_error")
    expect_identical(error$missing_samples, "s4")
    expect_identical(error$unexpected_samples, "other")
    expect_match(
        conditionMessage(error),
        "must match the analysis sample names exactly",
        fixed = TRUE
    )

    for (sample_names in list(
        c("s1", "s1", "s3", "s4"),
        c("s1", "", "s3", "s4"),
        c("s1", NA_character_, "s3", "s4")
    )) {
        expect_error(
            imputefinder:::.align_missingness_design(design, sample_names),
            class = "imputefinder_design_alignment_error"
        )
    }
})

test_that("sample metadata schema failures are typed and strict", {
    sample_data <- design_metadata_fixture()

    invalid_inputs <- list(
        as.matrix(sample_data),
        unname(sample_data),
        sample_data[FALSE, , drop = FALSE]
    )
    for (invalid in invalid_inputs) {
        expect_error(
            missingness_design(invalid, condition = "condition"),
            class = "imputefinder_design_schema_error"
        )
    }

    default_names <- data.frame(condition = c("A", "B"))
    expect_error(
        missingness_design(default_names, condition = "condition"),
        "explicit, non-empty, unique sample row names",
        fixed = TRUE,
        class = "imputefinder_design_schema_error"
    )

    duplicated_rows <- sample_data
    attr(duplicated_rows, "row.names") <- c("s4", "s1", "s1", "s2")
    expect_error(
        missingness_design(duplicated_rows, condition = "condition"),
        class = "imputefinder_design_schema_error"
    )

    duplicated_columns <- cbind(sample_data["condition"], sample_data["batch"])
    names(duplicated_columns) <- c("condition", "condition")
    expect_error(
        missingness_design(duplicated_columns, condition = "condition"),
        class = "imputefinder_design_schema_error"
    )
})

test_that("declared roles reject ambiguity and incomplete values", {
    sample_data <- design_metadata_fixture()

    for (condition in list(NULL, character(), c("condition", "batch"), NA, "")) {
        expect_error(
            missingness_design(sample_data, condition = condition),
            class = "imputefinder_design_role_error"
        )
    }
    expect_error(
        missingness_design(sample_data, condition = "absent"),
        "must name sample metadata columns exactly",
        fixed = TRUE,
        class = "imputefinder_design_role_error"
    )
    expect_error(
        missingness_design(
            sample_data,
            condition = "condition",
            nuisance = c("batch", "batch")
        ),
        class = "imputefinder_design_role_error"
    )
    expect_error(
        missingness_design(
            sample_data,
            condition = "condition",
            nuisance = "condition"
        ),
        "Each metadata column must have exactly one design role",
        fixed = TRUE,
        class = "imputefinder_design_role_error"
    )

    missing_value <- sample_data
    missing_value$batch[1] <- NA
    expect_error(
        missingness_design(
            missing_value,
            condition = "condition",
            nuisance = "batch"
        ),
        "Declared role columns must not contain missing or empty values",
        fixed = TRUE,
        class = "imputefinder_design_role_error"
    )

    non_finite <- sample_data
    non_finite$run_order[1] <- Inf
    expect_error(
        missingness_design(
            non_finite,
            condition = "condition",
            nuisance = "run_order"
        ),
        class = "imputefinder_design_role_error"
    )

    list_role <- sample_data
    list_role$batch <- I(as.list(list_role$batch))
    expect_error(
        missingness_design(
            list_role,
            condition = "condition",
            nuisance = "batch"
        ),
        class = "imputefinder_design_role_error"
    )
})

test_that("interactions are explicit canonical role sets", {
    sample_data <- design_metadata_fixture()

    design <- missingness_design(
        sample_data,
        condition = "condition",
        nuisance = c("run_order", "batch"),
        block = "subject",
        interactions = list(
            c("run_order", "condition"),
            c("condition", "batch", "run_order")
        )
    )
    reordered <- missingness_design(
        sample_data[, rev(names(sample_data)), drop = FALSE],
        condition = "condition",
        nuisance = c("batch", "run_order"),
        block = "subject",
        interactions = rev(list(
            c("run_order", "condition"),
            c("run_order", "batch", "condition")
        ))
    )

    expect_identical(design, reordered)
    expect_identical(
        design$interactions,
        list(
            c("batch", "condition", "run_order"),
            c("condition", "run_order")
        )
    )

    invalid_interactions <- list(
        "condition:batch",
        list("condition"),
        list(c("condition", "absent")),
        list(c("condition", "condition")),
        list(c("condition", "batch"), c("batch", "condition"))
    )
    for (interactions in invalid_interactions) {
        expect_error(
            missingness_design(
                sample_data,
                condition = "condition",
                nuisance = c("batch", "run_order"),
                interactions = interactions
            ),
            class = "imputefinder_design_role_error"
        )
    }
})

test_that("design lifecycle rejects silent schema upgrades", {
    design <- missingness_design(
        design_metadata_fixture(),
        condition = "condition"
    )
    round_trip <- unserialize(serialize(design, NULL, version = 3L))

    expect_identical(round_trip, design)
    expect_invisible(
        imputefinder:::.validate_missingness_design(round_trip)
    )

    wrong_schema <- design
    wrong_schema$spec$schema <- "missingness_design_v2"
    error <- tryCatch(
        imputefinder:::.validate_missingness_design(wrong_schema),
        imputefinder_design_lifecycle_error = identity
    )
    expect_s3_class(error, "imputefinder_design_lifecycle_error")
    expect_identical(error$expected_schema, "missingness_design_v1")
    expect_identical(error$actual_schema, "missingness_design_v2")

    wrong_lifecycle <- design
    wrong_lifecycle$spec$lifecycle <- "stable"
    expect_error(
        imputefinder:::.validate_missingness_design(wrong_lifecycle),
        class = "imputefinder_design_lifecycle_error"
    )

    malformed <- design
    malformed$sample_data$condition[1] <- NA_character_
    expect_error(
        imputefinder:::.validate_missingness_design(malformed),
        class = "imputefinder_design_schema_error"
    )
})

test_that("typed design construction leaves rules_v1 exactly isolated", {
    fixture <- normative_fixture()
    before <- classify_normative_fixture()
    before_bytes <- serialize(before, NULL, version = 3L)
    classifier_formals <- formals(classify_missingness)

    sample_data <- data.frame(
        condition = fixture$group,
        row.names = colnames(fixture$x)
    )
    design <- missingness_design(sample_data, condition = "condition")
    capture.output(print(design))
    after <- classify_normative_fixture()

    expect_identical(formals(classify_missingness), classifier_formals)
    expect_identical(serialize(after, NULL, version = 3L), before_bytes)
    expect_identical(names(attributes(after)), c("names", "class"))
})
