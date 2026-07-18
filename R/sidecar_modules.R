.sidecar_module_problem <- function(modules) {
    if (!is.character(modules)) {
        return("type")
    }
    if (anyNA(modules) || any(!nzchar(modules))) {
        return("missing_or_empty")
    }
    if (anyDuplicated(modules)) {
        return("duplicate")
    }
    if (!all(modules %in% .SIDECAR_MODULE_NAMES)) {
        return("unknown")
    }
    NULL
}

.normalise_sidecar_modules <- function(modules) {
    problem <- .sidecar_module_problem(modules)
    if (!is.null(problem)) {
        .abort_sidecar(
            paste0(
                "`modules` must be a unique character subset of: ",
                paste(.SIDECAR_MODULE_NAMES, collapse = ", "),
                "."
            ),
            "imputefinder_analysis_module_error",
            problem = problem,
            supported = .SIDECAR_MODULE_NAMES
        )
    }

    .SIDECAR_MODULE_NAMES[.SIDECAR_MODULE_NAMES %in% unname(modules)]
}

.sidecar_module_selected <- function(modules, module) {
    module %in% modules
}

.run_pre_rescue_sentinel <- function(modules, data, groups_by_sample) {
    if (!.sidecar_module_selected(modules, "sentinel")) {
        return(NULL)
    }
    .new_sidecar_sentinel(data, groups_by_sample)
}

.pending_sidecar_module <- function(module) {
    identifier <- unname(.SIDECAR_MODULE_IDENTIFIERS[[module]])
    .new_unavailable(
        quantity = module,
        code = "module_pending_validation",
        message = paste0(
            "The requested module remains unavailable until its frozen ",
            "validation and implementation milestones pass."
        ),
        requires = identifier
    )
}

.run_pending_stability <- function(modules) {
    if (!.sidecar_module_selected(modules, "stability")) {
        return(NULL)
    }
    .pending_sidecar_module("stability")
}

.validate_sidecar_stability <- function(stability) {
    if (is.null(stability)) {
        return(invisible(stability))
    }
    if (inherits(stability, "imputefinder_unavailable")) {
        .validate_unavailable(stability)
        valid <- identical(stability$quantity, "stability") &&
            identical(stability$code, "module_pending_validation") &&
            identical(
                stability$requires,
                unname(.SIDECAR_MODULE_IDENTIFIERS[["stability"]])
            )
        if (valid) {
            return(invisible(stability))
        }
    }
    .abort_sidecar(
        "Stored stability output does not satisfy its module schema.",
        "imputefinder_analysis_schema_error",
        field = "stability"
    )
}
