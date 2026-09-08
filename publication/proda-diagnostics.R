script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
setwd(dirname(dirname(normalizePath(script))))
Sys.setenv(IMPUTEFINDER_HELPERS_ONLY="1")
source("publication/benchmark.R")
Sys.unsetenv("IMPUTEFINDER_HELPERS_ONLY")
probe <- function(scenario) {
    sim <- simulate_data(scenario,6L,1L)
    x <- sim$x[rowSums(!is.na(sim$x))>0,,drop=FALSE]
    design <- model.matrix(~factor(sim$groups,levels=c("A","B")))
    colnames(design) <- c("Intercept","BvsA")
    records <- list()
    for(cap in c(20L,100L)) {
        set.seed(3401L)
        elapsed <- system.time(fit <- proDA::proDA(x,design=design,max_iter=cap,verbose=FALSE))[["elapsed"]]
        convergence <- proDA::convergence(fit)
        records[[as.character(cap)]] <- data.frame(scenario=scenario,seed=1,cap=cap,
            converged=convergence$successful,iterations=convergence$iterations,error=convergence$error)
        message(scenario," cap=",cap," converged=",convergence$successful,
            " iterations=",convergence$iterations," elapsed_seconds=",elapsed)
    }
    do.call(rbind,records)
}
records <- parallel::mclapply(c("MCAR","mixed","heteroscedastic_batch"),probe,
    mc.cores=if(.Platform$OS.type=="windows")1L else 3L,mc.set.seed=FALSE)
stopifnot(!any(vapply(records,inherits,logical(1),"try-error")))
write_csv(do.call(rbind,records),"proda-solver-probe.csv")
