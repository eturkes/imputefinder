.abort_sidecar <- function(message, subclass, ...) {
    error <- structure(
        c(
            list(message = message, call = NULL),
            list(...)
        ),
        class = c(
            subclass,
            "imputefinder_analysis_error",
            "error",
            "condition"
        )
    )
    stop(error)
}

.sidecar_scalar_character <- function(x, default = NA_character_) {
    if (is.character(x) &&
        length(x) == 1L &&
        !is.na(x) &&
        nzchar(x)) {
        unname(x)
    } else {
        default
    }
}
