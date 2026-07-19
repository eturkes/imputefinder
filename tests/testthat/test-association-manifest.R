association_manifest_root <- function() {
    root <- testthat::test_path("..", "..")
    if (!file.exists(file.path(
        root,
        "dev",
        "m13-association-contract.R"
    ))) {
        testthat::skip("repository-only implementation-seal rail")
    }
    normalizePath(root, mustWork = TRUE)
}

association_manifest_runs <- function(source_sha256 = strrep("a", 64L)) {
    commands <- c(
        "Rscript --vanilla dev/m13-association-contract.R --verify",
        "Rscript --vanilla dev/m13-association-implementation.R --verify"
    )
    namespace_sha256 <- imputefinder:::.association_namespace_sha256()
    list(
        list(
            command = commands[[1L]],
            status = 0L,
            output = charToRaw(paste0(
                "contract_version: m13_association_contract_v3\n",
                "state: frozen_unrun; candidate=1-32; ",
                "development=sealed; confirmation=sealed\n",
                "self_tests: alpha=TRUE; beta=TRUE\n"
            )),
            errors = raw()
        ),
        list(
            command = commands[[2L]],
            status = 0L,
            output = charToRaw(paste0(
                "m13a_test_evidence_v1\tpassed\t7\t",
                namespace_sha256,
                "\t",
                source_sha256,
                "\n"
            )),
            errors = raw()
        )
    )
}

association_manifest_fixture <- function() {
    root <- association_manifest_root()
    source_files <- imputefinder:::.new_association_source_manifest(root)
    source_sha256 <-
        imputefinder:::.association_source_manifest_sha256(source_files)
    receipt <- imputefinder:::.association_test_receipt_from_runs(
        association_manifest_runs(source_sha256)
    )
    manifest <- list(
        schema = "m13a_implementation_manifest_v1",
        contract_hash = "85d121ebd5ddd589be0f389553a3b7ffa6076a12a8fcbcb476a8969f7d987ae4",
        protocol_hash = "5f484f614d72d9560d96e2b8f75b71cc1027950958ea1459be93eece595069f5",
        gate_registry_hash = "f2fd8d0171003e0ae9b31757acf2d43acf7477965d517f1757f4d895c6abe731",
        effective_manifest_hash = "5dcf4f6408aed4ac7f9ccf07f4502805d34d2ea513ac68ac1c8e1259f23372f6",
        source_files = source_files,
        environment = imputefinder:::.association_implementation_environment(
            receipt$namespace_sha256
        ),
        tests = receipt$tests,
        state = "locked",
        manifest_hash = NA_character_
    )
    manifest$manifest_hash <-
        imputefinder:::.association_implementation_manifest_hash(manifest)
    list(manifest = manifest, receipt = receipt)
}

test_that("evidence hashing exactly ports the frozen escaped-frame protocol", {
    frame <- data.frame(
        key = c("b", "a%\t\r\n"),
        integer = c(2L, NA_integer_),
        double = c(-0, 1 / 3),
        logical = c(TRUE, NA),
        text = c(NA_character_, "é"),
        stringsAsFactors = FALSE
    )
    empty <- data.frame(
        key = character(),
        value = double(),
        stringsAsFactors = FALSE
    )

    expect_identical(
        imputefinder:::.association_hash_frame(frame),
        "373bfe1f37e4d7dc447bfac7fcc58889022454e509c133dcf77925353283e292"
    )
    expect_identical(
        imputefinder:::.association_evidence_bundle_hash(list(
            frame = frame,
            empty = empty
        )),
        "7d423b99dda20674d0d23e65c9005b409570fdb6d0560c9962e3e51a96d37e49"
    )
    expect_identical(
        imputefinder:::.association_hash_frame(frame[2:1, , drop = FALSE]),
        imputefinder:::.association_hash_frame(frame)
    )
})

test_that("test evidence derives status and counts from registered executions", {
    runs <- association_manifest_runs()
    evidence <- imputefinder:::.association_test_evidence_from_runs(runs)
    receipt <- imputefinder:::.association_test_receipt_from_runs(runs)

    expect_identical(evidence$command, vapply(
        runs,
        `[[`,
        character(1L),
        "command"
    ))
    expect_identical(evidence$expectation_count, c(2L, 7L))
    expect_identical(evidence$status, c("passed", "passed"))
    expect_true(all(grepl("^[0-9a-f]{64}$", evidence$evidence_sha256)))
    expect_identical(
        receipt$namespace_sha256,
        imputefinder:::.association_namespace_sha256()
    )
    expect_identical(receipt$source_sha256, strrep("a", 64L))

    failures <- list(
        { x <- runs; x[[1L]]$status <- 1L; x },
        { x <- runs; x[[1L]]$errors <- charToRaw("warning\n"); x },
        {
            x <- runs
            x[[1L]]$output <- charToRaw(paste0(
                "contract_version: m13_association_contract_v3\n",
                "state: frozen_unrun; candidate=1-32; ",
                "development=sealed; confirmation=sealed\n",
                "self_tests: alpha=FALSE\n"
            ))
            x
        },
        {
            x <- runs
            x[[2L]]$output <- charToRaw(
                paste0(
                    "m13a_test_evidence_v1\tfailed\t7\t",
                    imputefinder:::.association_namespace_sha256(),
                    "\t",
                    strrep("a", 64L),
                    "\n"
                )
            )
            x
        },
        { x <- runs; x[[1L]]$command <- "false"; x }
    )
    for (candidate in failures) {
        expect_error(
            imputefinder:::.association_test_evidence_from_runs(candidate),
            class = "imputefinder_association_manifest_error"
        )
    }
})

test_that("source inventory is canonical and complete across implementation roles", {
    root <- association_manifest_root()
    specification <- imputefinder:::.association_required_source_spec(root)
    manifest <- imputefinder:::.new_association_source_manifest(root)

    expect_identical(
        specification$path,
        sort(specification$path, method = "radix")
    )
    expect_setequal(
        specification$role,
        c("contract", "engine", "evaluator", "schema", "test", "harness")
    )
    expect_true(all(c(
        "R/association_candidate.R",
        "R/association_freedman_lane.R",
        "R/association_quasibinomial.R",
        "R/association_robust.R",
        "R/association_manifest.R",
        "dev/m13-association-contract.R",
        "dev/m13-association-implementation.R",
        "tests/testthat/test-association-manifest.R"
    ) %in% specification$path))
    expect_identical(manifest[c("path", "role")], specification)
    expect_true(all(grepl("^[0-9a-f]{64}$", manifest$sha256)))
})

test_that("manifest core rejects every detached live seal component", {
    fixture <- association_manifest_fixture()
    manifest <- fixture$manifest
    reseal <- function(candidate) {
        candidate$manifest_hash <-
            imputefinder:::.association_implementation_manifest_hash(candidate)
        candidate
    }

    expect_invisible(
        imputefinder:::.validate_association_implementation_manifest_core(
            manifest,
            association_manifest_root(),
            fixture$receipt
        )
    )
    corruptions <- list(
        { x <- manifest; x$manifest_hash <- strrep("0", 64L); x },
        {
            x <- manifest
            x$source_files$sha256[[1L]] <- strrep("0", 64L)
            reseal(x)
        },
        {
            x <- manifest
            x$source_files <- x$source_files[-1L, , drop = FALSE]
            reseal(x)
        },
        {
            x <- manifest
            x$source_files$role[[1L]] <- "engine"
            reseal(x)
        },
        {
            x <- manifest
            x$environment$value[[1L]] <- "drift"
            reseal(x)
        },
        {
            x <- manifest
            x$tests$expectation_count[[1L]] <- 999L
            reseal(x)
        }
    )
    for (candidate in corruptions) {
        expect_error(
            imputefinder:::.validate_association_implementation_manifest_core(
                candidate,
                association_manifest_root(),
                fixture$receipt
            ),
            class = "imputefinder_association_manifest_error"
        )
    }

    detached_receipt <- fixture$receipt
    detached_receipt$tests$evidence_sha256[[1L]] <- strrep("0", 64L)
    expect_error(
        imputefinder:::.validate_association_implementation_manifest_core(
            manifest,
            association_manifest_root(),
            detached_receipt
        ),
        class = "imputefinder_association_manifest_error"
    )

    detached_namespace <- fixture$receipt
    detached_namespace$namespace_sha256 <- strrep("0", 64L)
    expect_error(
        imputefinder:::.validate_association_implementation_manifest_core(
            manifest,
            association_manifest_root(),
            detached_namespace
        ),
        class = "imputefinder_association_manifest_error"
    )

    detached_source <- fixture$receipt
    detached_source$source_sha256 <- strrep("0", 64L)
    expect_error(
        imputefinder:::.validate_association_implementation_manifest_core(
            manifest,
            association_manifest_root(),
            detached_source
        ),
        class = "imputefinder_association_manifest_error"
    )
})

test_that("authorization remains frozen until semantic evidence is complete", {
    fixture <- association_manifest_fixture()

    expect_error(
        imputefinder:::.new_association_implementation_manifest(
            association_manifest_root()
        ),
        class = "imputefinder_association_manifest_error"
    )
    expect_error(
        imputefinder:::.validate_association_implementation_manifest(
            fixture$manifest,
            association_manifest_root()
        ),
        class = "imputefinder_association_manifest_error"
    )
})
