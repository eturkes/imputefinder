script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
setwd(dirname(dirname(normalizePath(script))))
Sys.setenv(IMPUTEFINDER_DOWNSTREAM_HELPERS_ONLY="1")
source("publication/downstream.R")
Sys.unsetenv("IMPUTEFINDER_DOWNSTREAM_HELPERS_ONLY")
full <- as.matrix(read.csv("publication/data/pxd020717-donors.csv",row.names=1,check.names=FALSE))
complete <- full[rowSums(is.na(full))==0,,drop=FALSE]
groups <- ifelse(startsWith(colnames(complete),"AD"),"AD","control")
records <- list();k <- 0L
for(mask in c("MCAR","low_tail")) for(seed in 1:30) {
    set.seed(seed+70000L)
    chance <- if(mask=="MCAR") .2 else .6*plogis((quantile(complete,.25)-complete)/.5)
    missing <- matrix(runif(length(complete)),nrow(complete))<chance
    x <- complete;x[missing] <- NA_real_
    for(method in setdiff(methods,c("manual_12","observed_only","proDA"))) {
        warnings <- character()
        ans <- tryCatch(withCallingHandlers(imputed_pipeline(x,groups,method,seed+800L),
            warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
        good <- !inherits(ans,"error")
        k <- k+1L
        metrics <- c(available=as.numeric(good),retained_fraction=0,reconstructed_fraction=0,missing_RMSE=NA_real_,mean_bias=NA_real_)
        if(good) {
            idx <- match(rownames(ans),rownames(x)); held <- missing[idx,,drop=FALSE]
            error <- ans-complete[idx,,drop=FALSE]
            metrics <- c(available=1,retained_fraction=nrow(ans)/nrow(x),
                reconstructed_fraction=sum(held)/sum(missing),missing_RMSE=sqrt(mean(error[held]^2)),mean_bias=mean(error[held]))
        }
        records[[k]] <- cbind(data.frame(mask=mask,seed=seed,method=method,heldout_cells=sum(missing),
            reason=if(good) "" else conditionMessage(ans),warnings=paste(unique(warnings),collapse=" | ")),
            as.data.frame(as.list(metrics)))
    }
}
d <- do.call(rbind,records)
write_csv(d,"holdout.csv")
write_csv(summarise_metrics(d,c("mask","method"),c("available","retained_fraction","reconstructed_fraction","missing_RMSE","mean_bias")),"holdout-summary.csv")
message("Held-out complete subset: ",nrow(complete)," features; ",nrow(d)," pipelines.")
