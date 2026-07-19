.ASSOCIATION_PANEL_SCHEMA <- "sentinel_association_v1"
.ASSOCIATION_PANEL_FIELDS <- c(
    "schema", "protocol", "candidate", "methods", "response",
    "hypotheses", "support", "outcomes", "multiplicity", "diagnostics"
)
.ASSOCIATION_PANEL_PROTOCOL_FIELDS <- c(
    "id", "hash", "contract_hash", "gate_registry_hash",
    "implementation_hash", "candidate_evidence_hash", "winner_state",
    "m12_contract_hash", "m12_gate_registry_hash"
)
.ASSOCIATION_PANEL_WINNER_STATES <- c(
    "winner_locked", "development_pass", "confirmation_pass"
)
.abort_association_panel <- function(message) {
    .abort_association(
        message,
        "imputefinder_association_panel_error"
    )
}

.association_panel_methods <- function(candidate) {
    c(
        response = .ASSOCIATION_RESPONSE_METHOD,
        encoding = unname(.DESIGN_ESTIMABILITY_METHODS[["encoding"]]),
        rank = unname(.DESIGN_ESTIMABILITY_METHODS[["rank"]]),
        candidate = candidate,
        multiplicity = "Holm"
    )
}

.association_panel_candidate_artifact <- function(panel) {
    structure(
        list(
            schema = .ASSOCIATION_CANDIDATE_ARTIFACT_SCHEMA,
            protocol = .association_protocol(),
            candidate = panel$candidate,
            input_sha256 = panel$diagnostics$input_sha256,
            response = panel$response,
            hypotheses = panel$hypotheses,
            support = panel$support,
            outcomes = panel$outcomes,
            multiplicity = panel$multiplicity,
            diagnostics = panel$diagnostics
        ),
        class = "imputefinder_association_candidate_artifact"
    )
}

.validate_association_panel_header <- function(panel) {
    expected_methods <- if (is.list(panel) &&
        .association_scalar_character(panel$candidate)) {
        .association_panel_methods(panel$candidate)
    } else {
        NULL
    }
    protocol <- if (is.list(panel) && is.list(panel$protocol)) {
        panel$protocol
    } else {
        NULL
    }
    valid <- is.list(panel) &&
        identical(class(panel), "imputefinder_association_panel") &&
        identical(names(panel), .ASSOCIATION_PANEL_FIELDS) &&
        identical(panel$schema, .ASSOCIATION_PANEL_SCHEMA) &&
        .association_scalar_character(panel$candidate) &&
        panel$candidate %in% .ASSOCIATION_CANDIDATES &&
        is.list(protocol) &&
        identical(names(protocol), .ASSOCIATION_PANEL_PROTOCOL_FIELDS) &&
        identical(protocol$id, .ASSOCIATION_PROTOCOL_ID) &&
        identical(protocol$hash, .ASSOCIATION_PROTOCOL_HASH) &&
        identical(protocol$contract_hash, .ASSOCIATION_CONTRACT_HASH) &&
        identical(
            protocol$gate_registry_hash,
            .ASSOCIATION_GATE_REGISTRY_HASH
        ) &&
        .association_sha256(protocol$implementation_hash) &&
        .association_sha256(protocol$candidate_evidence_hash) &&
        .association_scalar_character(protocol$winner_state) &&
        protocol$winner_state %in% .ASSOCIATION_PANEL_WINNER_STATES &&
        identical(
            protocol$m12_contract_hash,
            .ASSOCIATION_M12_CONTRACT_HASH
        ) &&
        identical(
            protocol$m12_gate_registry_hash,
            .ASSOCIATION_M12_GATE_REGISTRY_HASH
        ) &&
        identical(panel$methods, expected_methods)
    if (!valid) {
        .abort_association_panel(
            "Stored association panel header or winner binding is malformed."
        )
    }
    invisible(panel)
}

.validate_association_panel_shape <- function(panel, preparation) {
    if (missing(preparation) ||
        !inherits(
            preparation,
            "imputefinder_association_preparation"
        )) {
        .abort_association_panel(
            "Association panel validation requires its available preparation."
        )
    }
    .validate_association_panel_header(panel)
    artifact <- .association_panel_candidate_artifact(panel)
    error <- tryCatch(
        {
            .validate_association_candidate_artifact(
                artifact,
                preparation
            )
            NULL
        },
        imputefinder_association_artifact_error = identity
    )
    if (!is.null(error)) {
        .abort_association_panel(
            paste0(
                "Stored association panel payload is malformed: ",
                conditionMessage(error)
            )
        )
    }
    invisible(panel)
}
