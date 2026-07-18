.STABLE_SCIENTIFIC_IDENTIFIERS <- c(
    method = "rules_v1",
    profile = "count_weighted_density_v1",
    cutoff = "manual_or_derivative_boundary_v1",
    rescue = "condition_minimum_single_cell_v1",
    policy = "four_state_reconciliation_v1"
)

.SIDECAR_SCHEMA_IDENTIFIER <- "imputefinder_analysis_v1"
.SIDECAR_LIFECYCLE <- "experimental"
.SIDECAR_MODULE_IDENTIFIERS <- c(
    sentinel = "design_confounding_sentinel_v1",
    stability = "robustness_certificate_v1"
)
.SIDECAR_MODULE_NAMES <- names(.SIDECAR_MODULE_IDENTIFIERS)

.stable_scientific_identifiers <- function() {
    .STABLE_SCIENTIFIC_IDENTIFIERS
}

.sidecar_software_identity <- function() {
    c(
        package = "imputefinder",
        version = unname(as.character(
            getNamespaceVersion("imputefinder")
        ))
    )
}

.new_sidecar_spec <- function() {
    list(
        schema = .SIDECAR_SCHEMA_IDENTIFIER,
        lifecycle = .SIDECAR_LIFECYCLE,
        modules = .SIDECAR_MODULE_IDENTIFIERS,
        scientific = .stable_scientific_identifiers(),
        software = .sidecar_software_identity()
    )
}
