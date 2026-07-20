.UNAVAILABLE_FIELDS <- c(
    "status",
    "quantity",
    "code",
    "message",
    "requires"
)
.FIT_ONLY_ARTIFACTS <- c(
    "dropped_row_cells",
    "original_mask",
    "pre_rescue_evidence",
    "stability_perturbation_manifest"
)

.new_unavailable <- function(quantity, code, message, requires) {
    unavailable <- structure(
        list(
            status = "unavailable",
            quantity = quantity,
            code = code,
            message = message,
            requires = requires
        ),
        class = "imputefinder_unavailable"
    )
    .validate_unavailable(unavailable)
    unavailable
}

.validate_unavailable <- function(unavailable) {
    valid <- is.list(unavailable) &&
        identical(class(unavailable), "imputefinder_unavailable") &&
        identical(names(unavailable), .UNAVAILABLE_FIELDS) &&
        identical(unavailable$status, "unavailable") &&
        !is.na(.sidecar_scalar_character(unavailable$quantity)) &&
        !is.na(.sidecar_scalar_character(unavailable$code)) &&
        !is.na(.sidecar_scalar_character(unavailable$message)) &&
        is.character(unavailable$requires) &&
        length(unavailable$requires) > 0L &&
        !anyNA(unavailable$requires) &&
        all(nzchar(unavailable$requires)) &&
        !anyDuplicated(unavailable$requires)
    if (!valid) {
        .abort_sidecar(
            "Stored unavailable result is malformed.",
            "imputefinder_unavailable_schema_error",
            field = "unavailable"
        )
    }

    invisible(unavailable)
}

.normalise_fit_only_artifact <- function(artifact) {
    valid <- is.character(artifact) &&
        length(artifact) == 1L &&
        !is.na(artifact) &&
        artifact %in% .FIT_ONLY_ARTIFACTS
    if (!valid) {
        .abort_sidecar(
            "Unknown fit-only artifact request.",
            "imputefinder_unavailable_request_error",
            artifact = artifact,
            supported = .FIT_ONLY_ARTIFACTS
        )
    }

    unname(artifact)
}

.fit_only_artifact <- function(base_fit, artifact) {
    .validate_base_fit_result(base_fit)
    artifact <- .normalise_fit_only_artifact(artifact)
    .new_unavailable(
        quantity = artifact,
        code = "original_input_required",
        message = paste0(
            "A classic fit cannot reconstruct this artifact; supply the ",
            "original `x` and typed design to `analyze_missingness()`."
        ),
        requires = c("x", "missingness_design")
    )
}
