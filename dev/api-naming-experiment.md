# M11 API naming experiment

Date: 2026-07-18
Decision scope: experimental sidecar golden path + public/internal boundary

## Protocol

Criteria, ordered: call-site meaning → function/class distinction → consistency
with stable `classify_missingness()`/`plot_missingness()` → collision risk →
brevity. This is an engineering rubric, not user-research evidence.

Current primary guidance:

- [Bioconductor R code - naming](https://contributions.bioconductor.org/r-code.html#naming-of-packages-functions-and-classes): constructors are ordinarily plain functions; avoid confusing existing names.
- [Bioconductor namespaces](https://contributions.bioconductor.org/namespace.html#function-names): underscore/camel-case exports are valid; dot-prefixed functions communicate internal scope.
- [`stats::model.matrix()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/model.matrix.html): “design matrix” already means expanded model terms in core R. The metadata-role constructor must not imply that it creates this matrix.

Collision probes:

1. Exact searches restricted to current official CRAN/Bioconductor package
   documentation found none of the candidates below.
2. A dependency-complete local namespace scan found zero exact exports across
   145 installed packages:

   ```r
   candidates <- c(
       "missingness_design", "design_missingness", "imputefinder_design",
       "define_missingness_design", "analyze_missingness",
       "analyse_missingness", "missingness_analysis", "test_detectability",
       "detectability_contrast"
   )
   packages <- installed.packages()[, "Package"]
   lapply(candidates, function(candidate) {
       packages[vapply(packages, function(package) {
           candidate %in% tryCatch(
               getNamespaceExports(package),
               error = function(error) character()
           )
       }, logical(1))]
   })
   ```

Collision absence is a tie-breaker, not a uniqueness guarantee.

## Constructor result

| Candidate | Call-site reading | Decision |
|---|---|---|
| `missingness_design()` | declare typed metadata roles for missingness analysis | select |
| `design_missingness()` | sounds like an action that generates a missingness process | reject |
| `imputefinder_design()` | package prefix adds little at a namespaced call and weakens conceptual search | reject |
| `define_missingness_design()` | accurate but longer without distinguishing a second constructor action | reject |

Function and S3 class intentionally share `missingness_design`: the function is
a constructor, not a statistical operation. It stores metadata roles and
declared interactions; model-matrix construction, rank, aliasing, estimability,
and resampling-unit derivation remain later internal operations.

## Golden-path result

- `missingness_design()` - exported constructor; selected now.
- `analyze_missingness()` - selected future verb for running the sidecar.
  `missingness_analysis()` reads like another constructor/class;
  `analyse_missingness()` breaks existing US-English verb consistency.
- `test_detectability()` - selected future verb for one immutable contrast
  result. `detectability_contrast()` reads like a contrast constructor rather
  than the detection + observed-abundance analysis.

Only implemented functions are exported. Dot-prefixed schema constructors,
validators, alignment, fingerprints, compatibility comparators, and module
runners stay internal. Typed error classes are catchable API; their constructors
stay internal. M11b freezes design error families:
`imputefinder_design_schema_error`, `imputefinder_design_role_error`,
`imputefinder_design_alignment_error`, and
`imputefinder_design_lifecycle_error`, all inheriting from
`imputefinder_design_error`.
