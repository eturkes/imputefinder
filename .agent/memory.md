# Project memory

- Local R dependencies → `.agent/R-library` (gitignored). Prefix commands with `R_LIBS_USER="$PWD/.agent/R-library"`.
- Runtime target → ordinary numeric matrices + `SummarizedExperiment`; one pure matrix core. `PLAN.md` owns scientific semantics.
