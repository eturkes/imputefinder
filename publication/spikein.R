script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
setwd(dirname(dirname(normalizePath(script))))
Sys.setenv(IMPUTEFINDER_DOWNSTREAM_HELPERS_ONLY="1")
source("publication/downstream.R")
Sys.unsetenv("IMPUTEFINDER_DOWNSTREAM_HELPERS_ONLY")
x <- as.matrix(read.csv("publication/data/iprg2015-intensities.csv",row.names=1,check.names=FALSE))
samples <- read.csv("publication/data/iprg2015-samples.csv")
truth <- read.csv("publication/data/iprg2015-truth.csv")
spikes <- unique(truth$feature)
stopifnot(length(spikes)==6L,all(spikes %in% rownames(x)),identical(colnames(x),samples$sample))
contrasts <- combn(1:4,2,simplify=FALSE)
records <- list();details <- list();offset_records <- list();k <- 0L;j <- 0L
for(pair in contrasts) {
    select <- samples$condition %in% pair
    input <- x[,select,drop=FALSE]
    groups <- ifelse(samples$condition[select]==pair[1],"A","B")
    amount_A <- truth[truth$condition==pair[1],];amount_B <- truth[truth$condition==pair[2],]
    expected <- setNames(log2(amount_B$amount_fmol[match(amount_A$feature,amount_B$feature)]/amount_A$amount_fmol),amount_A$feature)
    stopifnot(all(expected!=0))
    medians <- apply(input,2,median,na.rm=TRUE)
    offsets <- medians-median(medians)
    name <- paste0(pair[2],"_vs_",pair[1])
    offset_records[[name]] <- data.frame(contrast=name,sample=colnames(input),offset=offsets)
    for(normalization in c("observed_median","none")) for(method in setdiff(methods,"manual_12")) {
        warnings <- character()
        convergence <- list(successful=NA,iterations=NA_integer_,error=NA_real_)
        ans <- tryCatch(withCallingHandlers({
            working <- if(method=="proDA") input[rowSums(!is.na(input))>0,,drop=FALSE] else imputed_pipeline(input,groups,method,900L)
            if(normalization=="observed_median") working <- sweep(working,2,offsets,"-")
            design <- model.matrix(~factor(groups,levels=c("A","B")))
            colnames(design) <- c("Intercept","BvsA")
            if(method=="proDA") {
                set.seed(3000L+pair[1]*10L+pair[2])
                fit <- proDA::proDA(working,design=design,max_iter=100L,verbose=FALSE)
                convergence <- proDA::convergence(fit)
                if(!isTRUE(convergence$successful)) stop("proDA failed numerical convergence within 100 iterations.")
                result <- proDA::test_diff(fit,contrast="BvsA",sort_by=NULL)
                result <- result[match(rownames(working),result$name),]
                effect <- result$diff;pvalue <- result$pval
            } else {
                fit <- limma::eBayes(limma::lmFit(working,design),trend=FALSE)
                effect <- fit$coefficients[,2];pvalue <- fit$p.value[,2]
            }
            allp <- setNames(rep(1,nrow(input)),rownames(input))
            allp[rownames(working)] <- ifelse(is.finite(pvalue),pvalue,1)
            significant <- names(allp)[p.adjust(allp,"BH")<.05]
            spike_estimate <- setNames(effect,rownames(working))[spikes]
            spike_error <- spike_estimate-expected[spikes]
            n_estimated <- sum(is.finite(spike_error))
            fp <- sum(!significant %in% spikes);tp <- sum(spikes %in% significant)
            list(summary=data.frame(available=1,retained=nrow(working),tested=sum(is.finite(pvalue)),
                spikes_detected=tp,spikes_total=6,spike_power=tp/6,background_discoveries=fp,
                realized_FDP=if(tp+fp)fp/(tp+fp) else 0,spike_effect_RMSE=sqrt(safe_mean(spike_error^2)),
                spikes_estimated=n_estimated,proda_converged=as.numeric(convergence$successful),
                proda_iterations=convergence$iterations,proda_error=convergence$error,
                proda_default_converged=if(is.na(convergence$successful))NA_real_ else as.numeric(convergence$successful && convergence$iterations<=20)),
                detail=data.frame(feature=spikes,expected_log2FC=expected[spikes],estimated_log2FC=spike_estimate,
                    discovered=spikes %in% significant))
        },warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
        good <- !inherits(ans,"error")
        if(!good) ans <- list(summary=data.frame(available=0,retained=0,tested=0,spikes_detected=0,spikes_total=6,
            spike_power=0,background_discoveries=NA_integer_,realized_FDP=NA_real_,spike_effect_RMSE=NA_real_,spikes_estimated=0,
            proda_converged=as.numeric(convergence$successful),proda_iterations=convergence$iterations,proda_error=convergence$error,
            proda_default_converged=if(is.na(convergence$successful))NA_real_ else as.numeric(convergence$successful && convergence$iterations<=20)),
            detail=data.frame(feature=spikes,expected_log2FC=expected[spikes],estimated_log2FC=NA_real_,discovered=FALSE),
            reason=conditionMessage(ans))
        k <- k+1L;j <- j+1L
        records[[k]] <- cbind(data.frame(contrast=name,normalization=normalization,method=method,
            reason=if(good)"" else ans$reason,warnings=paste(unique(warnings),collapse=" | ")),ans$summary)
        details[[j]] <- cbind(data.frame(contrast=name,normalization=normalization,method=method),ans$detail)
    }
}
write_csv(do.call(rbind,records),"spikein.csv")
write_csv(do.call(rbind,details),"spikein-effects.csv")
write_csv(do.call(rbind,offset_records),"spikein-normalization.csv")
message("iPRG2015: ",k," method-contrast-normalization outcomes; no independent-contrast confidence intervals.")
