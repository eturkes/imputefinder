script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
setwd(dirname(dirname(normalizePath(script))))
read <- function(name) read.csv(file.path("publication/results",name))
s <- read("simulation.csv")
stopifnot(nrow(s)==7200L,!anyDuplicated(s[c("scenario","nrep","seed","method")]),
    all(table(s$method)==720L),all(is.finite(s$available)),all(s$available %in% 0:1),
    all(s$actionable_coverage<=s$routing_coverage+1e-12),
    all(s$actionable_correct<=s$correct_all_missing+1e-12),
    all(s$actionable_correct<=s$actionable_coverage+1e-12),
    all(s$correct_all_missing<=s$routing_coverage+1e-12),
    all(s$assigned_cells<=s$missing_cells),
    all(s$available[s$reason!=""]==0),
    all(s$actionable_coverage[s$available==0]==0))
d <- read("downstream.csv")
stopifnot(nrow(d)==720L,!anyDuplicated(d[c("scenario","seed","method")]),
    all(!is.na(d$proda_converged[d$method=="proDA" & d$available==1])),
    all(table(d$method)==90L),all(d$power[d$available==0]==0),
    all(is.na(d$FDP[d$available==0])),all(is.na(d$missing_RMSE[d$method %in% c("observed_only","proDA")])) )
h <- read("holdout.csv")
stopifnot(nrow(h)==300L,!anyDuplicated(h[c("mask","seed","method")]),
    all(h$reconstructed_fraction<=1),all(h$reconstructed_fraction[h$available==0]==0))
b <- read("real-bootstrap.csv")
stopifnot(nrow(b)==400L,!anyDuplicated(b[c("rescue","seed","condition")]))
r <- read("real-rescue-seeds.csv")
stopifnot(nrow(r)==340L,all(is.finite(r$median_shift_vs_original)))
i <- read("spikein.csv");ie <- read("spikein-effects.csv")
probe <- read("proda-solver-probe.csv")
stopifnot(nrow(probe)==6L,all(probe$iterations<=probe$cap),all(is.finite(probe$error)))
stopifnot(nrow(i)==84L,nrow(ie)==504L,all(i$spikes_total==6L),
    all(!is.na(i$proda_converged[i$method=="proDA" & i$available==1])),
    all(i$spikes_detected<=6L),all(abs(i$spike_power-i$spikes_detected/6)<1e-12),
    all(ie$expected_log2FC!=0),!anyDuplicated(ie[c("contrast","normalization","method","feature")]),
    all(i$spike_power[i$available==0]==0))
stopifnot(all(d$available[d$method=="proDA" & d$proda_converged==0]==0),
    all(i$available[i$method=="proDA" & i$proda_converged==0]==0),
    all(d$proda_iterations[d$method=="proDA"]<=100),
    all(i$proda_iterations[i$method=="proDA"]<=100),
    all(d$proda_default_converged[d$method=="proDA"]==as.numeric(d$proda_converged[d$method=="proDA"]==1 & d$proda_iterations[d$method=="proDA"]<=20)))
for(file in c("simulation-summary.csv","downstream-summary.csv","holdout-summary.csv")) {
    summary <- read(file)
    stopifnot(all(summary$n<=summary$total),all(summary$n>=0),
        all(is.na(summary$estimate[summary$n==0])))
    availability <- subset(summary,metric=="available")
    stopifnot(all(is.finite(availability$lower)),all(is.finite(availability$upper)),
        all(availability$lower>=-1e-12),all(availability$upper<=1+1e-12))
}
cat("PASS: dimensions, unique independent units, causal/actionable denominators, failures, uncertainty, spike truth.\n")
