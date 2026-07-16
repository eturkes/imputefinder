#!/usr/bin/env Rscript

# Fatal Bioconductor submission gates. The source tree must be a clean Git
# clone; the tarball must come from a vignette-bearing R CMD build of it.
#
# Rscript --vanilla dev/bioc-gates.R <git-clone> <source-tarball>

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
    stop(
        "Usage: bioc-gates.R <git-clone> <source-tarball>",
        call. = FALSE
    )
}

git_clone <- normalizePath(arguments[[1L]], mustWork = TRUE)
source_tarball <- normalizePath(arguments[[2L]], mustWork = TRUE)

.require_clean_gate <- function(label, result) {
    counts <- vapply(
        c("error", "warning", "note"),
        result$getNum,
        integer(1L)
    )
    message(
        sprintf(
            "%s: %d errors, %d warnings, %d notes",
            label,
            counts[["error"]],
            counts[["warning"]],
            counts[["note"]]
        )
    )
    if (counts[["error"]] > 0L || counts[["warning"]] > 0L) {
        stop(label, " failed", call. = FALSE)
    }
    invisible(result)
}

.require_clean_gate(
    "BiocCheckGitClone",
    BiocCheck::BiocCheckGitClone(git_clone)
)
.require_clean_gate(
    "BiocCheck(new-package = TRUE)",
    BiocCheck::BiocCheck(source_tarball, `new-package` = TRUE)
)
