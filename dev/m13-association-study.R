#!/usr/bin/env Rscript

# Sealed M13 Track A candidate-study harness. Run from repository root only.
# It allocates synthetic replicates 1-32 after the implementation manifest
# locks. Development replicates, external artifacts, and public materialization
# have no command path here.

.m13a_study_fail <- function(message) {
    stop(message, call. = FALSE)
}

.m13a_study_root <- ".agent/m13a-candidate-study"
.m13a_command_lock <- ".agent/m13a-candidate-study.command-lock"
.m13a_manifest_path <- file.path(
    .m13a_study_root,
    "implementation-manifest.rds"
)
.m13a_evidence_path <- file.path(
    .m13a_study_root,
    "candidate-evidence.rds"
)

.m13a_require_root <- function() {
    valid <- file.exists("DESCRIPTION") && identical(
        read.dcf("DESCRIPTION", fields = "Package")[[1L]],
        "imputefinder"
    ) && dir.exists("R") && dir.exists("dev")
    if (!valid) {
        .m13a_study_fail(
            "M13 association study commands require the repository root."
        )
    }
    invisible(TRUE)
}

.m13a_atomic_save_rds <- function(value, path) {
    if (file.exists(path) || dir.exists(path)) {
        .m13a_study_fail(paste0(
            "M13 association artifact already exists: ",
            path,
            "."
        ))
    }
    directory <- dirname(path)
    if (!dir.exists(directory) && !dir.create(
        directory,
        recursive = TRUE,
        showWarnings = FALSE
    )) {
        .m13a_study_fail(paste0(
            "Cannot create M13 association artifact directory: ",
            directory,
            "."
        ))
    }
    normalized_root <- normalizePath(
        .m13a_study_root,
        winslash = "/",
        mustWork = TRUE
    )
    normalized_directory <- normalizePath(
        directory,
        winslash = "/",
        mustWork = TRUE
    )
    if (!identical(normalized_directory, normalized_root) && !startsWith(
        normalized_directory,
        paste0(normalized_root, "/")
    )) {
        .m13a_study_fail(
            "M13 association artifact directory resolves outside its study root."
        )
    }
    temporary <- tempfile(
        ".m13a-partial-",
        tmpdir = normalized_directory
    )
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    saveRDS(value, temporary, version = 3L, compress = "xz")
    if (!file.link(temporary, path)) {
        .m13a_study_fail(paste0(
            "Cannot atomically lock M13 association artifact: ",
            path,
            "."
        ))
    }
    unlink(temporary, force = TRUE)
    invisible(path)
}

.m13a_with_command_lock <- function(command, action) {
    if (identical(command, "--status")) {
        return(action())
    }
    parent <- dirname(.m13a_command_lock)
    if (!dir.exists(parent) && !dir.create(
        parent,
        recursive = TRUE,
        showWarnings = FALSE
    )) {
        .m13a_study_fail("Cannot create the project-local agent directory.")
    }
    repository_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
    normalized_parent <- normalizePath(
        parent,
        winslash = "/",
        mustWork = TRUE
    )
    if (!startsWith(normalized_parent, paste0(repository_root, "/"))) {
        .m13a_study_fail(
            "The M13 association command-lock parent is outside the repository."
        )
    }
    if (!dir.create(.m13a_command_lock, showWarnings = FALSE)) {
        .m13a_study_fail(
            "Another M13 association study command holds the exclusive lock."
        )
    }
    on.exit(unlink(.m13a_command_lock, recursive = TRUE, force = TRUE))
    action()
}

.m13a_read_manifest <- function() {
    if (!file.exists(.m13a_manifest_path)) {
        .m13a_study_fail(
            "Lock the M13 association implementation before allocation."
        )
    }
    tryCatch(
        readRDS(.m13a_manifest_path),
        error = function(error) {
            .m13a_study_fail(
                "The M13 association implementation manifest is unreadable."
            )
        }
    )
}

.m13a_lock_implementation <- function() {
    if (file.exists(.m13a_study_root) || dir.exists(.m13a_study_root)) {
        .m13a_study_fail(
            "The M13 association study root already exists; lock is immutable."
        )
    }
    manifest <- imputefinder:::.new_association_implementation_manifest(".")
    parent <- dirname(.m13a_study_root)
    if (!dir.exists(parent) && !dir.create(
        parent,
        recursive = TRUE,
        showWarnings = FALSE
    )) {
        .m13a_study_fail("Cannot create the project-local agent directory.")
    }
    repository_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
    parent <- normalizePath(parent, winslash = "/", mustWork = TRUE)
    if (!startsWith(parent, paste0(repository_root, "/"))) {
        .m13a_study_fail(
            "The project-local agent directory resolves outside the repository."
        )
    }
    temporary <- tempfile("m13a-study-lock-", tmpdir = parent)
    if (!dir.create(temporary, showWarnings = FALSE)) {
        .m13a_study_fail("Cannot create the temporary M13 study lock.")
    }
    complete <- FALSE
    on.exit({
        if (!complete && dir.exists(temporary)) {
            unlink(temporary, recursive = TRUE, force = TRUE)
        }
    }, add = TRUE)
    saveRDS(
        manifest,
        file.path(temporary, basename(.m13a_manifest_path)),
        version = 3L,
        compress = "xz"
    )
    if (!file.rename(temporary, .m13a_study_root)) {
        .m13a_study_fail(
            "Cannot atomically install the M13 association study lock."
        )
    }
    complete <- TRUE
    cat("implementation_manifest: ", manifest$manifest_hash, "\n", sep = "")
    cat("state: locked; candidate allocation remains unopened\n")
    invisible(manifest)
}

.m13a_protocol_context <- function() {
    manifest <- .m13a_read_manifest()
    resolver <- imputefinder:::.new_association_protocol_study_resolver(
        manifest,
        "."
    )
    list(manifest = manifest, resolver = resolver)
}

.m13a_validate_or_write_input <- function(
    resolver,
    runtime,
    scenario_id,
    replicate
) {
    input <- imputefinder:::.association_protocol_input(
        runtime,
        scenario_id,
        replicate
    )
    path <- imputefinder:::.association_protocol_artifact_path(
        resolver,
        "input",
        scenario_id,
        replicate,
        must_exist = FALSE
    )
    if (file.exists(path)) {
        stored <- imputefinder:::.association_read_protocol_input(
            resolver,
            scenario_id,
            replicate
        )
        imputefinder:::.validate_association_study_input(
            stored,
            scenario_id,
            replicate
        )
        if (!identical(stored, input)) {
            .m13a_study_fail(
                "Stored M13 association input differs from generator replay."
            )
        }
        return(stored)
    }
    .m13a_atomic_save_rds(input, path)
    input
}

.m13a_validate_or_write_result <- function(
    resolver,
    candidate_id,
    scenario_id,
    replicate,
    input_sha256,
    preparation
) {
    path <- imputefinder:::.association_protocol_artifact_path(
        resolver,
        "result",
        scenario_id,
        replicate,
        candidate_id,
        must_exist = FALSE
    )
    if (file.exists(path)) {
        record <- imputefinder:::.association_read_protocol_result(
            resolver,
            candidate_id,
            scenario_id,
            replicate
        )
        imputefinder:::.validate_association_study_result_record(
            record,
            candidate_id,
            scenario_id,
            replicate,
            input_sha256
        )
        replay <- imputefinder:::.association_execute_protocol_candidate(
            candidate_id,
            preparation
        )
        if (!identical(record$result, replay)) {
            .m13a_study_fail(
                "Stored M13 association result differs from engine replay."
            )
        }
        return(record)
    }
    started <- unname(proc.time()[["elapsed"]])
    result <- imputefinder:::.association_execute_protocol_candidate(
        candidate_id,
        preparation
    )
    elapsed <- as.double(max(
        0,
        unname(proc.time()[["elapsed"]]) - started
    ))
    record <- imputefinder:::.new_association_study_result(
        candidate_id,
        scenario_id,
        replicate,
        input_sha256,
        result,
        elapsed
    )
    .m13a_atomic_save_rds(record, path)
    record
}

.m13a_run_candidates <- function() {
    if (file.exists(.m13a_evidence_path)) {
        .m13a_study_fail(
            "Candidate evidence is already locked; allocation cannot resume."
        )
    }
    context <- .m13a_protocol_context()
    resolver <- context$resolver
    runtime <- imputefinder:::.association_protocol_runtime(resolver)
    completed <- 0L
    total <- length(imputefinder:::.ASSOCIATION_STUDY_SCENARIOS) *
        length(imputefinder:::.ASSOCIATION_STUDY_REPLICATES) *
        length(imputefinder:::.ASSOCIATION_CANDIDATES)
    for (scenario_id in imputefinder:::.ASSOCIATION_STUDY_SCENARIOS) {
        for (replicate in imputefinder:::.ASSOCIATION_STUDY_REPLICATES) {
            input <- .m13a_validate_or_write_input(
                resolver,
                runtime,
                scenario_id,
                replicate
            )
            input_sha256 <- imputefinder:::.association_study_input_sha256(
                input
            )
            preparation <- imputefinder:::.new_association_preparation(
                input$data,
                input$design
            )
            for (candidate_id in imputefinder:::.ASSOCIATION_CANDIDATES) {
                .m13a_validate_or_write_result(
                    resolver,
                    candidate_id,
                    scenario_id,
                    replicate,
                    input_sha256,
                    preparation
                )
                completed <- completed + 1L
                cat(
                    "candidate_run: ", completed, "/", total, " ",
                    candidate_id, " ", scenario_id, " ", replicate, "\n",
                    sep = ""
                )
            }
        }
    }
    imputefinder:::.validate_association_implementation_manifest(
        context$manifest,
        "."
    )
    cat("state: candidate artifacts complete; evidence unresolved\n")
    invisible(TRUE)
}

.m13a_resolve_evidence <- function() {
    if (file.exists(.m13a_evidence_path)) {
        .m13a_study_fail("Candidate evidence is already locked.")
    }
    context <- .m13a_protocol_context()
    evidence <- imputefinder:::.authorize_association_candidate_evidence(
        context$resolver,
        context$manifest,
        "."
    )
    .m13a_atomic_save_rds(evidence, .m13a_evidence_path)
    cat("candidate_evidence: ", evidence$evidence_hash, "\n", sep = "")
    cat("state: ", evidence$state, "\n", sep = "")
    cat("selected_candidate: ", evidence$selected_candidate, "\n", sep = "")
    invisible(evidence)
}

.m13a_verify_evidence <- function() {
    if (!file.exists(.m13a_evidence_path)) {
        .m13a_study_fail("Candidate evidence has not been resolved.")
    }
    context <- .m13a_protocol_context()
    evidence <- tryCatch(
        readRDS(.m13a_evidence_path),
        error = function(error) {
            .m13a_study_fail("Candidate evidence is unreadable.")
        }
    )
    imputefinder:::.validate_association_candidate_evidence(
        evidence,
        context$resolver,
        context$manifest,
        "."
    )
    cat("candidate_evidence: ", evidence$evidence_hash, "\n", sep = "")
    cat("state: verified_", evidence$state, "\n", sep = "")
    invisible(evidence)
}

.m13a_status <- function() {
    manifest <- file.exists(.m13a_manifest_path)
    evidence <- file.exists(.m13a_evidence_path)
    input_count <- if (dir.exists(file.path(.m13a_study_root, "inputs"))) {
        length(list.files(
            file.path(.m13a_study_root, "inputs"),
            pattern = "^replicate-[0-9]{3}\\.rds$",
            recursive = TRUE
        ))
    } else {
        0L
    }
    result_count <- if (dir.exists(file.path(.m13a_study_root, "results"))) {
        length(list.files(
            file.path(.m13a_study_root, "results"),
            pattern = "^replicate-[0-9]{3}\\.rds$",
            recursive = TRUE
        ))
    } else {
        0L
    }
    cat("implementation_manifest_present: ", manifest, "\n", sep = "")
    cat("input_artifacts: ", input_count, "/416\n", sep = "")
    cat("result_artifacts: ", result_count, "/1248\n", sep = "")
    cat("candidate_evidence_present: ", evidence, "\n", sep = "")
    invisible(TRUE)
}

.m13a_require_root()
args <- commandArgs(trailingOnly = TRUE)
allowed <- c(
    "--lock-implementation", "--run-candidates", "--resolve-evidence",
    "--verify-evidence", "--status"
)
if (length(args) != 1L || !args %in% allowed) {
    .m13a_study_fail(paste0(
        "usage: Rscript --vanilla dev/m13-association-study.R ",
        "--lock-implementation|--run-candidates|--resolve-evidence|",
        "--verify-evidence|--status"
    ))
}

pkgload::load_all(".", quiet = TRUE)
action <- switch(
    args,
    "--lock-implementation" = .m13a_lock_implementation,
    "--run-candidates" = .m13a_run_candidates,
    "--resolve-evidence" = .m13a_resolve_evidence,
    "--verify-evidence" = .m13a_verify_evidence,
    "--status" = .m13a_status
)
.m13a_with_command_lock(args, action)
