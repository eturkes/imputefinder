# ImputeFinder

- Scope = condition-aware operational missingness routing; MAR/MNAR labels establish neither mechanism nor biological absence.
- API = `R/`; contract tests = `tests/testthat/`; package gate = `Rscript tools/check-package.R [--bioc]`.
- Publication = `publication/PROTOCOL.md` fixes estimands + ten audit rows; `publication/reproduce.py` regenerates pinned public data, simulations, downstream/holdout evidence, figures + article. `--verify` reruns then compares deterministic evidence hashes.
- Dependencies = R/Bioconductor from `DESCRIPTION`; publication additionally requires `limma`, `imputeLCMD`, `proDA`, `rmarkdown`, `knitr`, `scales`, Pandoc, Python standard library. Resolve `R_LIBS_USER` to absolute paths before package checks.
- Article claims bind to generated tables; cache validity = content hash, not file existence. Report automatic failures, ambiguity, retention, conditional denominators + single-imputation inference limits.
- No superiority/novelty claim without direct evidence against applicable prior art: DEP, DAPAR/ProStaR, msImpute, intensity/missingness sector routing. Prior hybrid strategy = Fowler et al.; current package automatic diagnostics require separate evaluation.
