.ASSOCIATION_IMPLEMENTATION_MANIFEST_SCHEMA <-
    "m13a_implementation_manifest_v1"
.ASSOCIATION_IMPLEMENTATION_MANIFEST_FIELDS <- c(
    "schema", "contract_hash", "protocol_hash", "gate_registry_hash",
    "effective_manifest_hash", "source_files", "environment", "tests",
    "state", "manifest_hash"
)
.ASSOCIATION_IMPLEMENTATION_SOURCE_FIELDS <- c("path", "role", "sha256")
.ASSOCIATION_IMPLEMENTATION_SOURCE_ROLES <- c(
    "contract", "engine", "evaluator", "schema", "test", "harness"
)
.ASSOCIATION_IMPLEMENTATION_ENVIRONMENT_FIELDS <- c(
    "name", "value", "required"
)
.ASSOCIATION_IMPLEMENTATION_TEST_FIELDS <- c(
    "command", "expectation_count", "status", "evidence_sha256"
)
.ASSOCIATION_IMPLEMENTATION_TEST_COMMANDS <- c(
    "Rscript --vanilla dev/m13-association-contract.R --verify",
    "Rscript --vanilla dev/m13-association-implementation.R --verify"
)
.ASSOCIATION_IMPLEMENTATION_TEST_SCRIPTS <- c(
    "dev/m13-association-contract.R",
    "dev/m13-association-implementation.R"
)

# The semantic resolver, immutable artifact inventory, and study harness are
# complete. Public construction still runs both registered verification rails.
.ASSOCIATION_IMPLEMENTATION_SEAL_READY <- TRUE

.abort_association_manifest <- function(message) {
    .abort_association(
        message,
        "imputefinder_association_manifest_error"
    )
}

.association_hash_escape <- function(x) {
    x <- enc2utf8(x)
    x <- gsub("%", "%25", x, fixed = TRUE)
    x <- gsub("\t", "%09", x, fixed = TRUE)
    x <- gsub("\r", "%0D", x, fixed = TRUE)
    gsub("\n", "%0A", x, fixed = TRUE)
}

.association_hash_canonical_value <- function(x) {
    if (is.integer(x)) {
        output <- as.character(x)
    } else if (is.numeric(x)) {
        output <- vapply(x, function(value) {
            if (is.na(value)) {
                return("<NA>")
            }
            format(
                value,
                digits = 17L,
                scientific = FALSE,
                trim = TRUE
            )
        }, character(1L))
    } else if (is.logical(x)) {
        output <- ifelse(
            is.na(x),
            "<NA>",
            ifelse(x, "TRUE", "FALSE")
        )
    } else {
        output <- x
    }
    output[is.na(output)] <- "<NA>"
    .association_hash_escape(output)
}

.association_hash_frame <- function(frame) {
    if (!is.data.frame(frame)) {
        .abort_association_manifest(
            "Association evidence components must be data frames."
        )
    }
    if (nrow(frame)) {
        frame <- frame[order(frame[[1L]], method = "radix"), , drop = FALSE]
    }
    header <- paste(
        .association_hash_escape(names(frame)),
        collapse = "\t"
    )
    rows <- if (nrow(frame)) {
        vapply(seq_len(nrow(frame)), function(index) {
            paste(vapply(frame, function(column) {
                .association_hash_canonical_value(column[[index]])
            }, character(1L)), collapse = "\t")
        }, character(1L))
    } else {
        character()
    }
    canonical <- paste(c(header, rows), collapse = "\n")
    unname(tools::sha256sum(bytes = charToRaw(enc2utf8(canonical))))
}

.association_evidence_bundle_hash <- function(components) {
    valid <- is.list(components) && length(components) > 0L &&
        !is.null(names(components)) && !anyNA(names(components)) &&
        all(nzchar(names(components))) && !anyDuplicated(names(components)) &&
        all(vapply(components, is.data.frame, logical(1L)))
    if (!valid) {
        .abort_association_manifest(
            "Association evidence-hash components are malformed."
        )
    }
    hashes <- vapply(
        components,
        .association_hash_frame,
        character(1L)
    )
    aggregate <- paste(names(hashes), hashes, sep = "\t", collapse = "\n")
    unname(tools::sha256sum(bytes = charToRaw(enc2utf8(aggregate))))
}

.association_output_bytes <- function(output) {
    if (is.raw(output)) {
        return(output)
    }
    if (is.character(output) && length(output) == 1L && !is.na(output)) {
        return(charToRaw(enc2utf8(output)))
    }
    .abort_association_manifest(
        "Association test output must be one text value or exact raw bytes."
    )
}

.association_output_sha256 <- function(output) {
    unname(tools::sha256sum(bytes = .association_output_bytes(output)))
}

.association_namespace_binding_hash <- function(value) {
    if (is.function(value)) {
        value <- utils::removeSource(value)
        attributes <- attributes(value)
        if (length(attributes)) {
            attributes[c("srcref", "srcfile", "wholeSrcref")] <- NULL
        }
        payload <- list(
            type = "closure",
            formals = formals(value),
            body = body(value),
            attributes = attributes
        )
    } else {
        type <- typeof(value)
        if (type %in% c(
            "environment", "externalptr", "weakref", "promise"
        )) {
            .abort_association_manifest(
                "The imputefinder namespace contains an unhashable binding."
            )
        }
        payload <- list(type = type, value = value)
    }
    unname(tools::sha256sum(bytes = serialize(
        payload,
        connection = NULL,
        version = 3L
    )))
}

.association_namespace_sha256 <- function() {
    namespace <- asNamespace("imputefinder")
    bindings <- sort(
        ls(namespace, all.names = TRUE),
        method = "radix"
    )
    bindings <- bindings[
        !startsWith(bindings, ".__") & bindings != ".packageName"
    ]
    if (!length(bindings) || anyDuplicated(bindings)) {
        .abort_association_manifest(
            "The live imputefinder namespace cannot be inventoried exactly."
        )
    }
    values <- lapply(
        bindings,
        get,
        envir = namespace,
        inherits = FALSE
    )
    frame <- data.frame(
        name = bindings,
        type = unname(vapply(values, typeof, character(1L))),
        sha256 = unname(vapply(
            values,
            .association_namespace_binding_hash,
            character(1L)
        )),
        stringsAsFactors = FALSE
    )
    .association_hash_frame(frame)
}

.association_namespace_version <- function(package) {
    version <- tryCatch(
        as.character(utils::packageVersion(package)),
        error = function(error) NA_character_
    )
    if (!.association_scalar_character(version)) {
        .abort_association_manifest(
            paste0("Cannot resolve loaded namespace version: ", package, ".")
        )
    }
    unname(version)
}

.association_namespace_is_base <- function(package) {
    priority <- tryCatch(
        utils::packageDescription(package, fields = "Priority"),
        error = function(error) NA_character_
    )
    identical(unname(priority), "base")
}

.association_implementation_environment <- function(
    namespace_sha256 = .association_namespace_sha256()
) {
    loaded <- loadedNamespaces()
    dependencies <- sort(
        setdiff(
            loaded[!vapply(
                loaded,
                .association_namespace_is_base,
                logical(1L)
            )],
            "imputefinder"
        ),
        method = "radix"
    )
    if (!"imputefinder" %in% loaded) {
        .abort_association_manifest(
            "The implementation seal requires the loaded imputefinder namespace."
        )
    }
    live_namespace_sha256 <- .association_namespace_sha256()
    if (!.association_sha256(namespace_sha256) ||
        !identical(namespace_sha256, live_namespace_sha256)) {
        .abort_association_manifest(
            paste0(
                "The live imputefinder namespace differs from the ",
                "source-loaded verification process."
            )
        )
    }
    names <- c(
        "R.platform",
        "R.version.string",
        "imputefinder.version",
        "imputefinder.namespace_sha256",
        paste0("dependency:", dependencies)
    )
    values <- c(
        R.version$platform,
        R.version.string,
        .association_namespace_version("imputefinder"),
        live_namespace_sha256,
        vapply(
            dependencies,
            .association_namespace_version,
            character(1L)
        )
    )
    environment <- data.frame(
        name = names,
        value = unname(values),
        required = rep(TRUE, length(names)),
        stringsAsFactors = FALSE
    )
    environment <- environment[
        order(environment$name, method = "radix"),
        ,
        drop = FALSE
    ]
    row.names(environment) <- NULL
    environment[.ASSOCIATION_IMPLEMENTATION_ENVIRONMENT_FIELDS]
}

.association_manifest_root <- function(root) {
    valid <- is.character(root) && length(root) == 1L && !is.na(root) &&
        nzchar(root) && dir.exists(root)
    if (!valid) {
        .abort_association_manifest(
            "Association manifest root must be one existing directory."
        )
    }
    normalizePath(root, winslash = "/", mustWork = TRUE)
}

.association_manifest_source_target <- function(path, root) {
    safe <- .association_scalar_character(path) &&
        !grepl("^/", path) && !grepl("\\\\", path) &&
        !grepl("(^|/)\\.{1,2}(/|$)", path) &&
        !grepl("//", path, fixed = TRUE) &&
        !grepl("[\t\r\n]", path)
    if (!safe) {
        .abort_association_manifest(
            "Association source paths must be safe canonical relative paths."
        )
    }
    target <- file.path(root, path)
    information <- file.info(target)
    if (nrow(information) != 1L || is.na(information$isdir) ||
        information$isdir) {
        .abort_association_manifest(
            paste0("Association source file is absent: ", path, ".")
        )
    }
    target <- normalizePath(target, winslash = "/", mustWork = TRUE)
    if (!startsWith(target, paste0(root, "/"))) {
        .abort_association_manifest(
            "Association source path resolves outside the manifest root."
        )
    }
    target
}

.association_relative_r_files <- function(root, directory) {
    target <- file.path(root, directory)
    if (!dir.exists(target)) {
        return(character())
    }
    files <- list.files(
        target,
        pattern = "\\.R$",
        recursive = TRUE,
        full.names = FALSE
    )
    paste0(directory, "/", files)
}

.association_source_role <- function(path) {
    if (path %in% c(
        "dev/m13-association-contract.R",
        "dev/m13-association-contract.md"
    )) {
        return("contract")
    }
    if (path %in% c(
        "dev/m13-association-implementation.R",
        "dev/m13-association-study.R"
    )) {
        return("harness")
    }
    if (grepl("^tests/", path)) {
        return("test")
    }
    if (grepl(
        "^R/association_(candidate|freedman_lane|preparation|quasibinomial|robust)\\.R$",
        path
    )) {
        return("engine")
    }
    if (grepl("^R/association_(evidence|manifest)\\.R$", path)) {
        return("evaluator")
    }
    "schema"
}

.association_required_source_spec <- function(root = ".") {
    root <- .association_manifest_root(root)
    paths <- c(
        "DESCRIPTION",
        "NAMESPACE",
        .association_relative_r_files(root, "R"),
        .association_relative_r_files(root, "dev"),
        "dev/m13-association-contract.md",
        "tests/testthat.R",
        .association_relative_r_files(root, "tests/testthat")
    )
    paths <- sort(unique(paths), method = "radix")
    paths <- paths[
        !startsWith(paths, "dev/") |
            grepl("^dev/m1[23].*\\.R$", paths) |
            paths == "dev/m13-association-contract.md"
    ]
    targets <- vapply(
        paths,
        .association_manifest_source_target,
        character(1L),
        root = root
    )
    if (length(targets) != length(paths)) {
        .abort_association_manifest(
            "Association source inventory could not be resolved exactly."
        )
    }
    output <- data.frame(
        path = paths,
        role = unname(vapply(paths, .association_source_role, character(1L))),
        stringsAsFactors = FALSE
    )
    if (!setequal(
        unique(output$role),
        .ASSOCIATION_IMPLEMENTATION_SOURCE_ROLES
    )) {
        .abort_association_manifest(
            "Association source inventory lacks a required implementation role."
        )
    }
    output
}

.new_association_source_manifest <- function(root = ".") {
    root <- .association_manifest_root(root)
    source_files <- .association_required_source_spec(root)
    targets <- vapply(
        source_files$path,
        .association_manifest_source_target,
        character(1L),
        root = root
    )
    output <- data.frame(
        path = source_files$path,
        role = source_files$role,
        sha256 = unname(tools::sha256sum(targets)),
        stringsAsFactors = FALSE
    )
    output[.ASSOCIATION_IMPLEMENTATION_SOURCE_FIELDS]
}

.association_source_manifest_sha256 <- function(source_files) {
    valid <- is.data.frame(source_files) &&
        identical(
            names(source_files),
            .ASSOCIATION_IMPLEMENTATION_SOURCE_FIELDS
        ) && is.character(source_files$path) &&
        is.character(source_files$role) && is.character(source_files$sha256) &&
        !anyNA(source_files) && !anyDuplicated(source_files$path) &&
        identical(
            source_files$path,
            sort(source_files$path, method = "radix")
        ) && all(vapply(
            source_files$sha256,
            .association_sha256,
            logical(1L)
        ))
    if (!valid) {
        .abort_association_manifest(
            "Association source inventory cannot be hashed exactly."
        )
    }
    .association_hash_frame(source_files)
}

.association_read_raw_file <- function(path) {
    size <- file.info(path)$size
    if (is.na(size) || size < 0 || size > .Machine$integer.max) {
        .abort_association_manifest(
            "Association harness output cannot be read exactly."
        )
    }
    readBin(path, what = "raw", n = as.integer(size))
}

.association_run_registered_test <- function(index, root) {
    if (!.association_integer_scalar(index, 1L) ||
        index > length(.ASSOCIATION_IMPLEMENTATION_TEST_COMMANDS)) {
        .abort_association_manifest(
            "Association test-command index is outside the frozen registry."
        )
    }
    stdout <- tempfile("m13a-stdout-")
    stderr <- tempfile("m13a-stderr-")
    previous <- getwd()
    on.exit({
        setwd(previous)
        unlink(c(stdout, stderr), force = TRUE)
    }, add = TRUE)
    setwd(root)
    status <- suppressWarnings(system2(
        file.path(R.home("bin"), "Rscript"),
        c(
            "--vanilla",
            .ASSOCIATION_IMPLEMENTATION_TEST_SCRIPTS[[index]],
            "--verify"
        ),
        stdout = stdout,
        stderr = stderr
    ))
    output <- .association_read_raw_file(stdout)
    errors <- .association_read_raw_file(stderr)
    list(
        command = .ASSOCIATION_IMPLEMENTATION_TEST_COMMANDS[[index]],
        status = as.integer(status),
        output = output,
        errors = errors
    )
}

.association_contract_expectation_count <- function(output) {
    text <- rawToChar(.association_output_bytes(output))
    lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
    tests <- lines[startsWith(lines, "self_tests: ")]
    valid <- length(tests) == 1L &&
        "contract_version: m13_association_contract_v4" %in% lines &&
        paste0(
            "state: frozen_unrun; candidate=1-32; ",
            "development=sealed; confirmation=sealed"
        ) %in% lines
    entries <- if (valid) {
        strsplit(sub("^self_tests: ", "", tests), "; ", fixed = FALSE)[[1L]]
    } else {
        character()
    }
    valid <- valid && length(entries) > 0L && all(grepl(
        "^[A-Za-z0-9_.]+=(TRUE)$",
        entries
    )) && !anyDuplicated(sub("=TRUE$", "", entries))
    if (!valid) {
        .abort_association_manifest(
            "Frozen contract verification output is malformed or failed."
        )
    }
    as.integer(length(entries))
}

.association_implementation_test_receipt <- function(output) {
    text <- rawToChar(.association_output_bytes(output))
    valid <- grepl(
        paste0(
            "^m13a_test_evidence_v1\\tpassed\\t[1-9][0-9]*\\t",
            "[0-9a-f]{64}\\t[0-9a-f]{64}\\n$"
        ),
        text
    )
    if (!valid) {
        .abort_association_manifest(
            "Focused implementation-test output is malformed or failed."
        )
    }
    fields <- strsplit(
        sub("\\n$", "", text),
        "\t",
        fixed = TRUE
    )[[1L]]
    receipt <- list(
        expectation_count = as.integer(fields[[3L]]),
        namespace_sha256 = fields[[4L]],
        source_sha256 = fields[[5L]]
    )
    if (!.association_integer_scalar(receipt$expectation_count, 1L) ||
        !.association_sha256(receipt$namespace_sha256) ||
        !.association_sha256(receipt$source_sha256)) {
        .abort_association_manifest(
            "Focused implementation-test receipt values are malformed."
        )
    }
    receipt
}

.association_test_receipt_from_runs <- function(runs) {
    valid <- is.list(runs) &&
        length(runs) == length(.ASSOCIATION_IMPLEMENTATION_TEST_COMMANDS) &&
        identical(
            vapply(runs, `[[`, character(1L), "command"),
            .ASSOCIATION_IMPLEMENTATION_TEST_COMMANDS
        ) && all(vapply(runs, function(run) {
            is.list(run) && identical(
                names(run),
                c("command", "status", "output", "errors")
            ) && .association_integer_scalar(run$status, 0L) &&
                identical(run$status, 0L) && is.raw(run$output) &&
                is.raw(run$errors) && !length(run$errors)
        }, logical(1L)))
    if (!valid) {
        .abort_association_manifest(
            "Registered association test execution failed or was malformed."
        )
    }
    implementation <- .association_implementation_test_receipt(
        runs[[2L]]$output
    )
    counts <- c(
        .association_contract_expectation_count(runs[[1L]]$output),
        implementation$expectation_count
    )
    evidence <- data.frame(
        command = .ASSOCIATION_IMPLEMENTATION_TEST_COMMANDS,
        expectation_count = counts,
        status = rep("passed", length(counts)),
        evidence_sha256 = unname(vapply(
            runs,
            function(run) .association_output_sha256(run$output),
            character(1L)
        )),
        stringsAsFactors = FALSE
    )
    list(
        tests = evidence[.ASSOCIATION_IMPLEMENTATION_TEST_FIELDS],
        namespace_sha256 = implementation$namespace_sha256,
        source_sha256 = implementation$source_sha256
    )
}

.association_test_evidence_from_runs <- function(runs) {
    .association_test_receipt_from_runs(runs)$tests
}

.run_association_test_receipt <- function(root = ".") {
    root <- .association_manifest_root(root)
    runs <- lapply(
        seq_along(.ASSOCIATION_IMPLEMENTATION_TEST_COMMANDS),
        .association_run_registered_test,
        root = root
    )
    .association_test_receipt_from_runs(runs)
}

.run_association_test_evidence <- function(root = ".") {
    .run_association_test_receipt(root)$tests
}

.association_implementation_manifest_components <- function(manifest) {
    scalar_fields <- c(
        "schema", "contract_hash", "protocol_hash", "gate_registry_hash",
        "effective_manifest_hash", "state"
    )
    data <- data.frame(
        position = match(
            scalar_fields,
            .ASSOCIATION_IMPLEMENTATION_MANIFEST_FIELDS
        ),
        field = scalar_fields,
        value = unname(vapply(
            scalar_fields,
            function(field) manifest[[field]],
            character(1L)
        )),
        type = rep("character", length(scalar_fields)),
        stringsAsFactors = FALSE
    )
    list(
        scalar_fields = data,
        source_files = manifest$source_files,
        environment = manifest$environment,
        tests = manifest$tests
    )
}

.association_implementation_manifest_hash <- function(manifest) {
    valid <- is.list(manifest) && identical(
        names(manifest),
        .ASSOCIATION_IMPLEMENTATION_MANIFEST_FIELDS
    ) && all(vapply(manifest[c(
        "schema", "contract_hash", "protocol_hash", "gate_registry_hash",
        "effective_manifest_hash", "state", "manifest_hash"
    )], is.character, logical(1L))) &&
        all(vapply(manifest[c(
            "source_files", "environment", "tests"
        )], is.data.frame, logical(1L)))
    if (!valid) {
        .abort_association_manifest(
            "Association implementation manifest cannot be canonicalized."
        )
    }
    .association_evidence_bundle_hash(
        .association_implementation_manifest_components(manifest)
    )
}

.validate_association_source_manifest <- function(source_files, root) {
    expected <- .association_required_source_spec(root)
    valid <- is.data.frame(source_files) &&
        identical(names(source_files), .ASSOCIATION_IMPLEMENTATION_SOURCE_FIELDS) &&
        identical(source_files[c("path", "role")], expected) &&
        is.character(source_files$sha256) && !anyNA(source_files$sha256) &&
        all(vapply(source_files$sha256, .association_sha256, logical(1L)))
    if (!valid) {
        .abort_association_manifest(
            "Stored association source inventory is malformed or incomplete."
        )
    }
    targets <- vapply(
        source_files$path,
        .association_manifest_source_target,
        character(1L),
        root = root
    )
    if (!identical(
        source_files$sha256,
        unname(tools::sha256sum(targets))
    )) {
        .abort_association_manifest(
            "Stored association source hashes differ from live file bytes."
        )
    }
    invisible(source_files)
}

.validate_association_implementation_environment <- function(
    environment,
    namespace_sha256
) {
    valid <- is.data.frame(environment) && nrow(environment) >= 3L &&
        identical(
            names(environment),
            .ASSOCIATION_IMPLEMENTATION_ENVIRONMENT_FIELDS
        ) && is.character(environment$name) &&
        is.character(environment$value) && is.logical(environment$required) &&
        !anyNA(environment) && all(nzchar(environment$name)) &&
        all(nzchar(environment$value)) && !anyDuplicated(environment$name) &&
        identical(
            environment$name,
            sort(environment$name, method = "radix")
        ) && all(environment$required)
    if (!valid || !identical(
        environment,
        .association_implementation_environment(namespace_sha256)
    )) {
        .abort_association_manifest(
            "Stored association environment differs from all live namespaces."
        )
    }
    invisible(environment)
}

.validate_association_test_evidence <- function(tests, observed_tests) {
    valid <- is.data.frame(tests) &&
        identical(names(tests), .ASSOCIATION_IMPLEMENTATION_TEST_FIELDS) &&
        identical(tests$command, .ASSOCIATION_IMPLEMENTATION_TEST_COMMANDS) &&
        is.integer(tests$expectation_count) &&
        all(tests$expectation_count > 0L) &&
        identical(tests$status, rep("passed", nrow(tests))) &&
        is.character(tests$evidence_sha256) &&
        all(vapply(
            tests$evidence_sha256,
            .association_sha256,
            logical(1L)
        )) && identical(tests, observed_tests)
    if (!valid) {
        .abort_association_manifest(
            "Stored association tests differ from registered passing executions."
        )
    }
    invisible(tests)
}

.validate_association_implementation_manifest_core <- function(
    manifest,
    root,
    observed_receipt
) {
    root <- .association_manifest_root(root)
    valid <- is.list(manifest) && identical(
        names(manifest),
        .ASSOCIATION_IMPLEMENTATION_MANIFEST_FIELDS
    ) && identical(
        manifest$schema,
        .ASSOCIATION_IMPLEMENTATION_MANIFEST_SCHEMA
    ) && identical(manifest$contract_hash, .ASSOCIATION_CONTRACT_HASH) &&
        identical(manifest$protocol_hash, .ASSOCIATION_PROTOCOL_HASH) &&
        identical(
            manifest$gate_registry_hash,
            .ASSOCIATION_GATE_REGISTRY_HASH
        ) && identical(
            manifest$effective_manifest_hash,
            .ASSOCIATION_EFFECTIVE_MANIFEST_HASH
        ) && identical(manifest$state, "locked") &&
        .association_sha256(manifest$manifest_hash)
    if (!valid) {
        .abort_association_manifest(
            "Stored association implementation-manifest header is malformed."
        )
    }
    .validate_association_source_manifest(manifest$source_files, root)
    receipt_valid <- is.list(observed_receipt) && identical(
        names(observed_receipt),
        c("tests", "namespace_sha256", "source_sha256")
    ) && is.data.frame(observed_receipt$tests) &&
        .association_sha256(observed_receipt$namespace_sha256) &&
        .association_sha256(observed_receipt$source_sha256)
    if (!receipt_valid) {
        .abort_association_manifest(
            "Observed association implementation receipt is malformed."
        )
    }
    .validate_association_implementation_environment(
        manifest$environment,
        observed_receipt$namespace_sha256
    )
    .validate_association_test_evidence(
        manifest$tests,
        observed_receipt$tests
    )
    if (!identical(
        observed_receipt$source_sha256,
        .association_source_manifest_sha256(manifest$source_files)
    )) {
        .abort_association_manifest(
            "Tested association source snapshot differs from the manifest."
        )
    }
    expected_hash <- .association_implementation_manifest_hash(manifest)
    if (!identical(manifest$manifest_hash, expected_hash)) {
        .abort_association_manifest(
            "Stored association implementation self-hash is detached."
        )
    }
    invisible(manifest)
}

.validate_association_implementation_manifest <- function(
    manifest,
    root = "."
) {
    if (!isTRUE(.ASSOCIATION_IMPLEMENTATION_SEAL_READY)) {
        .abort_association_manifest(
            paste0(
                "Association implementation seal remains frozen until the ",
                "semantic candidate-evidence resolver is complete."
            )
        )
    }
    observed_receipt <- .run_association_test_receipt(root)
    .validate_association_implementation_manifest_core(
        manifest,
        root,
        observed_receipt
    )
}

.new_association_implementation_manifest <- function(root = ".") {
    if (!isTRUE(.ASSOCIATION_IMPLEMENTATION_SEAL_READY)) {
        .abort_association_manifest(
            paste0(
                "Association implementation seal remains frozen until the ",
                "semantic candidate-evidence resolver is complete."
            )
        )
    }
    root <- .association_manifest_root(root)
    observed_receipt <- .run_association_test_receipt(root)
    manifest <- list(
        schema = .ASSOCIATION_IMPLEMENTATION_MANIFEST_SCHEMA,
        contract_hash = .ASSOCIATION_CONTRACT_HASH,
        protocol_hash = .ASSOCIATION_PROTOCOL_HASH,
        gate_registry_hash = .ASSOCIATION_GATE_REGISTRY_HASH,
        effective_manifest_hash = .ASSOCIATION_EFFECTIVE_MANIFEST_HASH,
        source_files = .new_association_source_manifest(root),
        environment = .association_implementation_environment(
            observed_receipt$namespace_sha256
        ),
        tests = observed_receipt$tests,
        state = "locked",
        manifest_hash = NA_character_
    )
    manifest$manifest_hash <-
        .association_implementation_manifest_hash(manifest)
    .validate_association_implementation_manifest_core(
        manifest,
        root,
        observed_receipt
    )
    manifest
}
