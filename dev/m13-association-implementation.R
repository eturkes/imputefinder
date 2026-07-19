.m13a_implementation_fail <- function(message) {
    stop(message, call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (!identical(args, "--verify")) {
    .m13a_implementation_fail(
        "usage: Rscript --vanilla dev/m13-association-implementation.R --verify"
    )
}
if (!file.exists("DESCRIPTION") ||
    !identical(
        read.dcf("DESCRIPTION", fields = "Package")[[1L]],
        "imputefinder"
    )) {
    .m13a_implementation_fail(
        "M13 association implementation verification requires the repository root."
    )
}

pkgload::load_all(".", quiet = TRUE)
source_before <- imputefinder:::.new_association_source_manifest(".")
source_before_sha256 <-
    imputefinder:::.association_source_manifest_sha256(source_before)

results <- testthat::test_local(
    ".",
    filter = "^association-",
    reporter = "silent",
    load_package = "source"
)
summary <- as.data.frame(results)
required_summary <- c(
    "nb", "failed", "skipped", "error", "warning", "passed"
)
summary_valid <- is.data.frame(summary) && nrow(summary) > 0L &&
    all(required_summary %in% names(summary)) &&
    all(vapply(
        summary[c("nb", "failed", "warning", "passed")],
        is.numeric,
        logical(1L)
    )) && is.logical(summary$skipped) && is.logical(summary$error) &&
    !anyNA(summary[required_summary])
valid <- summary_valid &&
    all(summary$failed == 0L) &&
    all(!summary$skipped) &&
    all(!summary$error) &&
    all(summary$warning == 0L) &&
    all(summary$nb == summary$passed)
if (!valid) {
    .m13a_implementation_fail(
        "M13 association focused implementation verification failed."
    )
}
source_after <- imputefinder:::.new_association_source_manifest(".")
source_after_sha256 <-
    imputefinder:::.association_source_manifest_sha256(source_after)
if (!identical(source_before_sha256, source_after_sha256)) {
    .m13a_implementation_fail(
        "M13 association sources changed during focused verification."
    )
}
expectation_count <- sum(summary$passed)
if (!is.numeric(expectation_count) || length(expectation_count) != 1L ||
    is.na(expectation_count) || expectation_count < 1L ||
    expectation_count > .Machine$integer.max ||
    expectation_count != floor(expectation_count)) {
    .m13a_implementation_fail(
        "M13 association expectation count is malformed."
    )
}
namespace_sha256 <- imputefinder:::.association_namespace_sha256()
if (!is.character(namespace_sha256) || length(namespace_sha256) != 1L ||
    is.na(namespace_sha256) ||
    !grepl("^[0-9a-f]{64}$", namespace_sha256)) {
    .m13a_implementation_fail(
        "M13 association namespace fingerprint is malformed."
    )
}
cat(
    "m13a_test_evidence_v1\tpassed\t",
    as.integer(expectation_count),
    "\t",
    namespace_sha256,
    "\t",
    source_after_sha256,
    "\n",
    sep = ""
)
