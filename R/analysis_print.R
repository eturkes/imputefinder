.analysis_module_status <- function(value) {
    if (is.null(value)) {
        return("not requested")
    }
    if (inherits(value, "imputefinder_unavailable")) {
        return(paste0("unavailable (", value$code, ")"))
    }
    "available"
}

#' @export
print.imputefinder_analysis <- function(x, ...) {
    .validate_imputefinder_analysis(x)
    classic <- if (inherits(x$classic, "imputefinder_classic_failure")) {
        paste0("failure at ", x$classic$stage)
    } else {
        "available"
    }
    dimensions <- x$input$dimensions

    cat("<imputefinder_analysis>\n")
    cat(sprintf(
        "Schema: %s (%s)\n",
        x$spec$schema,
        x$spec$lifecycle
    ))
    cat(sprintf(
        "Input: %d features x %d samples; classic: %s\n",
        dimensions[["features"]],
        dimensions[["samples"]],
        classic
    ))
    cat(sprintf(
        "Modules: sentinel: %s; stability: %s\n",
        .analysis_module_status(x$sentinel),
        .analysis_module_status(x$stability)
    ))

    invisible(x)
}
