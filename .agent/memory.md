# Project memory

- Local R dependencies → `.agent/R-library` (gitignored). Prefix commands with `R_LIBS_USER="$PWD/.agent/R-library"`.
- Full local checks → prepend `.agent/TinyTeX/bin/x86_64-linux` + `.agent/html-tidy/root/usr/bin` to `PATH`; prepend `.agent/html-tidy/root/usr/lib/x86_64-linux-gnu` to `LD_LIBRARY_PATH`.
- Documentation → roxygen2 8 records `Config/roxygen2/version`; current BiocCheck still detects generated Rd via `RoxygenNote`. Retain both accurate fields after regeneration until BiocCheck catches up.
- Runtime target → ordinary numeric matrices + `SummarizedExperiment`; one pure matrix core. `PLAN.md` owns scientific semantics.
