script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
setwd(dirname(dirname(normalizePath(script))))
Sys.setenv(IMPUTEFINDER_HELPERS_ONLY="1")
source("publication/benchmark.R")
Sys.unsetenv("IMPUTEFINDER_HELPERS_ONLY")
x <- as.matrix(read.csv("publication/data/pxd020717-donors.csv",row.names=1,check.names=FALSE))
groups <- ifelse(startsWith(colnames(x),"AD"),"AD","control")
records <- list(); diagnostic <- list(); k <- 0L
for(rescue in c("none","minimum")) {
    auto <- run_fit(x,groups,rescue)
    diagnostic[[rescue]] <- if(inherits(auto,"error")) list(reason=conditionMessage(auto),diagnostic=auto$diagnostic) else auto$cutoff_diagnostics
    # Fixed empirical alternative when no unique detector output exists; not validated thresholds.
    cut <- if(!inherits(auto,"error")) auto$cutoffs else sapply(sort(unique(groups)),function(g) {
        m <- rowMeans(x[,groups==g,drop=FALSE],na.rm=TRUE); unname(quantile(m,.25,na.rm=TRUE))
    })
    for(delta in c(NA_real_,-1,0,1)) {
        fit <- if(is.na(delta)) auto else run_fit(x,groups,rescue,cut+delta)
        k <- k+1L
        available <- !inherits(fit,"error")
        base <- data.frame(rescue=rescue,cutoff_delta=delta,mode=if(is.na(delta)) "automatic" else "sensitivity",
            anchor=if(inherits(auto,"error")) "empirical_q25" else "automatic",available=available,
            retained=if(available) nrow(fit$data) else NA_integer_,
            rescued=if(available) nrow(fit$seed_log) else NA_integer_,
            reason=if(available) "" else conditionMessage(fit))
        if(available) {
            write_csv(fit$classifications,paste0("real-states-",rescue,"-",if(is.na(delta)) "auto" else delta,".csv"))
            count <- table(factor(fit$classifications$state,levels=c("complete","MAR","MNAR","absent","insufficient")))
            base <- cbind(base,as.data.frame(as.list(setNames(as.integer(count),names(count)))))
            base$cutoff_AD <- fit$cutoffs[["AD"]];base$cutoff_control <- fit$cutoffs[["control"]]
            if(is.na(delta)) {
                for(g in sort(unique(groups))) {
                    p <- plot_missingness(fit,g)
                    ggplot2::ggsave(paste0("publication/figures/profile-",rescue,"-",g,".svg"),p,width=6,height=4,device=grDevices::svg)
                }
            }
        } else {
            for(s in c("complete","MAR","MNAR","absent","insufficient","cutoff_AD","cutoff_control")) base[[s]] <- NA_real_
        }
        records[[k]] <- base
    }
}
write_csv(do.call(rbind,records),"real-summary.csv")
writeLines(capture.output(dput(diagnostic)),"publication/results/real-diagnostics.R")
bootstrap <- list(); k <- 0L
for(rescue in c("none","minimum")) for(seed in 1:100) {
    set.seed(seed+9000L)
    sample_x <- x[sample(nrow(x),replace=TRUE),,drop=FALSE]
    rownames(sample_x) <- paste0("bootstrap",seq_len(nrow(x)))
    fit <- run_fit(sample_x,groups,rescue)
    good <- !inherits(fit,"error")
    for(g in sort(unique(groups))) {
        k <- k+1L
        bootstrap[[k]] <- data.frame(rescue=rescue,seed=seed,condition=g,available=good,
            cutoff=if(good) fit$cutoffs[[g]] else NA_real_,
            ambiguous=if(good) isTRUE(fit$cutoff_diagnostics[[g]]$quality$derivative$credible_lobe_count>1) else NA,
            reason=if(good) "" else conditionMessage(fit))
    }
}
write_csv(do.call(rbind,bootstrap),"real-bootstrap.csv")
quality <- data.frame(features=nrow(x),donors=ncol(x),AD=sum(groups=="AD"),control=sum(groups=="control"),
    complete_features=sum(rowSums(is.na(x))==0),missing_fraction=mean(is.na(x)),
    globally_absent=sum(rowSums(!is.na(x))==0))
write_csv(quality,"real-quality.csv")
run_x <- as.matrix(read.csv("publication/data/pxd020717-runs.csv",row.names=1,check.names=FALSE))
run_groups <- ifelse(startsWith(colnames(run_x),"AD"),"AD","control")
run_records <- lapply(c("none","minimum"),function(rescue) {
    fit <- run_fit(run_x,run_groups,rescue)
    good <- !inherits(fit,"error")
    data.frame(rescue=rescue,available=good,retained=if(good)nrow(fit$data) else NA_integer_,
        cutoff_AD=if(good)fit$cutoffs[["AD"]] else NA_real_,
        cutoff_control=if(good)fit$cutoffs[["control"]] else NA_real_,
        reason=if(good)"" else conditionMessage(fit))
})
write_csv(do.call(rbind,run_records),"real-run-sensitivity.csv")
seed_records <- list();reference <- NULL
for(seed in 1:20) {
    fit <- run_fit(x,groups,"minimum",seed=seed)
    stopifnot(!inherits(fit,"error"))
    if(is.null(reference)) reference <- fit
    stopifnot(identical(fit$classifications,reference$classifications),
        identical(fit$feature_status,reference$feature_status),identical(fit$cutoffs,reference$cutoffs))
    original <- x[rownames(fit$data),,drop=FALSE]
    stopifnot(identical(fit$data[!is.na(original)],original[!is.na(original)]))
    median_rescue <- apply(fit$data,2,median,na.rm=TRUE)
    median_without <- apply(original,2,median,na.rm=TRUE)
    median_reference <- apply(reference$data,2,median,na.rm=TRUE)
    seed_records[[seed]] <- data.frame(seed=seed,sample=colnames(x),
        median_shift_vs_original=median_rescue-median_without,
        median_shift_vs_seed1=median_rescue-median_reference,
        inserted_in_retained=colSums(is.na(original)&!is.na(fit$data)))
}
write_csv(do.call(rbind,seed_records),"real-rescue-seeds.csv")
message("Real data complete.")
