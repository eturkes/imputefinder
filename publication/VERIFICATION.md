# Verification

| Gate | Observed result |
|---|---|
| `Rscript tools/check-package.R --bioc` | R CMD check = `Status: OK`; tests=729 passes,0 failures/warnings/skips; BiocCheck1.49.30=0 errors,0 warnings,4 advisory notes |
| `Rscript publication/check-results.R` | Pass: dimensions, unique units, diagnostic/actionable denominators, failures, uncertainty, solver validity, spike truth |
| `python3 publication/test-reproduction.py` |4/4 pass: valid inventory, actual orphan-file rejection, missing-file rejection, Windows-style manifest-ID normalization |
| `python3 publication/reproduce.py --verify` |43/43 declared numerical/image artifacts byte-identical after full regeneration from empty result/figure directories |
| Repeated `Rscript publication/render.R` | HTML byte-identical; Word ZIP CRC + document XML valid; three embedded HTML figures; final references present |
| Visual QA | Three main scientific figures; rendered article opening + controlled-mixture table legible |
| `python3 publication/tidy-generated.py` twice | Second pass changes0/3 generated diagnostic files; upstream metadata untouched |

HTML SHA256=`ddb3607be2ddbd74a1fa9db4f7c1fd2f178aeb730c589ad571f0ec4fd81af22d`.
HTML equality applies to the final source and evaluated rendering environment; Word/PDF container metadata is outside the byte-equality claim.

## Evidence inventory

- Simulation=720 datasets ×10 routes →7,200 outcomes; uncertainty units=datasets.
- Downstream=90 datasets ×8 pipelines →720 outcomes; proDA cap100 convergence=90/90; default20 stopping=89/90.
- Observational source=PXD020717;2,645 filtered groups,17 donors;100 feature bootstraps per rescue policy;20 rescue seeds.
- Measured-cell masking=932 complete groups ×2 masks ×30 seeds ×5 pipelines →300 outcomes; complete-case selection explicit.
- Controlled source=iPRG2015 reprocessed in PXD028038;3,262 groups,12 technical runs,6 known spikes;84 pipeline/contrast/normalization outcomes +504 spike estimates; proDA convergence=12/12.
- Versions=R4.6.1; imputefinder0.99.8; detector3; limma3.99.0; proDA1.27.0; imputeLCMD2.1. Complete loaded versions=`results/*session-info.txt`.

Reproduction invocation uses caller-supplied absolute `R_LIBS_USER` paths and `PUBLICATION_WORKERS=4`; packages are not installed or mutated by the gate.
All source-table snapshots and extraction manifests remain checksum-pinned.

BiocCheck notes concern code/document formatting and admin-only subscription verification. These are advisory; they are not suppressed.
No hosted CI, journal acceptance, archival release DOI, or completed submission is claimed.
