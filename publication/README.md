# Article and evidence

Read [the rendered article](manuscript.html) or edit [its source](manuscript.Rmd).
An [editable Word article](manuscript.docx) supports journal preparation.
The [protocol](PROTOCOL.md) defines the estimands and review checks.
The [audit](AUDIT.md) records their adjudication.

## Reproduce

Use R 4.6 or later, Python 3, and Pandoc.
Install the package dependencies from `DESCRIPTION` in a compatible Bioconductor library.
Also install `limma`, `proDA`, `imputeLCMD`, `rmarkdown`, `knitr`, and `scales`.
The session records list the evaluated versions.
Set `R_LIBS_USER` to your library paths before running:
On Unix, set `PUBLICATION_WORKERS=4` to parallelise independent downstream datasets.
The default uses one worker; Windows also uses one worker.

```sh
python3 publication/reproduce.py
python3 publication/reproduce.py --verify
Rscript tools/check-package.R
Rscript tools/check-package.R --bioc
```

The first command clears the script-owned `results/` and `figures/` directories, then regenerates their contents and the article.
Keep handwritten files outside those directories.
The second command reruns every analysis and requires byte-identical numerical and image evidence.
Both commands require an exact output inventory; stale or missing artifacts cannot satisfy the gate.
It uses one numerical thread to limit machine-dependent variation.
HTML, Word, and PDF presentation files are regenerated but excluded from byte-equality claims.
The scripts source the checked-out package code, not an independently installed package copy.
Both public protein tables have checksum-pinned compressed snapshots for offline regeneration.

## Evidence map

| Files | Contents |
|---|---|
| `data/*manifest.json` | Source identity, extraction rules, and hashes |
| `results/simulation*.csv` | 720 datasets × ten operational routes; diagnostic and actionable denominators |
| `results/downstream*.csv` | 90 datasets × eight pipelines; availability, retention, reconstruction, effects, FDP, power, paired contrasts |
| `results/real-*.csv` | Full public-table decisions, bootstraps, run aggregation, and rescue sensitivity |
| `results/holdout*.csv` | Complete-feature measured-value masking; selection limitations apply |
| `results/spikein*.csv` | Controlled iPRG mixtures; all six spikes remain in every contrast denominator |
| `results/*session-info.txt` | Evaluated R and package versions |
| `figures/` | Vector SVG/PDF and PNG research figures |
| `artifact-manifest.json` | Deterministic evidence hashes |

MAR/MNAR labels are operational routes.
The article makes no general mechanism-identification, superiority, or FDR-control claim.
The iPRG contrasts share four mixtures and are not independent biological experiments.
