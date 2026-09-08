script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
setwd(dirname(dirname(normalizePath(script))))
Sys.setenv(IMPUTEFINDER_HELPERS_ONLY="1")
source("publication/benchmark.R")
Sys.unsetenv("IMPUTEFINDER_HELPERS_ONLY")
stopifnot(requireNamespace("limma"), requireNamespace("imputeLCMD"),requireNamespace("proDA"))
methods <- c("observed_only","proDA","uniform_KNN","uniform_MinProb","DEP_absence",
    "automatic_minimum","automatic_none","manual_12")

imputed_pipeline <- function(x,groups,method,seed,cutoffs=NULL) {
    keep <- rowSums(!is.na(x))>0
    working <- x[keep,,drop=FALSE]
    state <- baseline(working,groups,"DEP_absence")$state
    if(method %in% c("automatic_minimum","automatic_none","manual_12")) {
        if(method=="manual_12" && is.null(cutoffs)) cutoffs <- setNames(c(12,12),sort(unique(groups)))
        fit <- run_fit(x,groups,if(method=="automatic_minimum") "minimum" else "none",cutoffs)
        if(inherits(fit,"error")) stop(fit)
        working <- fit$data
        state <- state_matrix(fit,working,groups)
    }
    if(nrow(working)<11L) stop("Fewer than 11 retained features for kNN comparison.")
    if(method=="observed_only") return(working)
    for(g in sort(unique(groups))) {
        columns <- groups==g; block <- working[,columns,drop=FALSE]
        condition_seed <- seed + 10000L * match(g,sort(unique(groups)))
        set.seed(condition_seed+1000L)
        if(method!="uniform_MinProb") {
            invisible(capture.output(knn <- imputeLCMD::impute.wrapper.KNN(block,K=10)))
        }
        set.seed(condition_seed+2000L)
        if(method!="uniform_KNN") {
            invisible(capture.output(low <- imputeLCMD::impute.MinProb(block,q=.01,tune.sigma=1)))
        }
        candidate <- if(method=="uniform_MinProb") low else knn
        if(!method %in% c("uniform_MinProb","uniform_KNN")) {
            low_rows <- state[,g] %in% c("MNAR","absent")
            candidate[low_rows,] <- low[low_rows,,drop=FALSE]
        }
        stopifnot(identical(candidate[!is.na(block)],block[!is.na(block)]))
        working[,columns] <- candidate
    }
    if(anyNA(working) || any(!is.finite(working))) stop("Imputation left non-finite values.")
    working
}

evaluate_downstream <- function(sim, method, seed) {
    warnings <- character()
    convergence <- list(successful=NA,iterations=NA_integer_,error=NA_real_)
    ans <- tryCatch(withCallingHandlers({
        working <- if(method=="proDA") sim$x[rowSums(!is.na(sim$x))>0,,drop=FALSE] else imputed_pipeline(sim$x,sim$groups,method,seed)
        indices <- match(rownames(working),rownames(sim$x))
        design <- model.matrix(~factor(sim$groups,levels=c("A","B")))
        colnames(design) <- c("Intercept","BvsA")
        if(method=="proDA") {
            set.seed(seed+3000L)
            fit <- proDA::proDA(working,design=design,max_iter=100L,verbose=FALSE)
            convergence <- proDA::convergence(fit)
            if(!isTRUE(convergence$successful)) stop("proDA failed numerical convergence within 100 iterations.")
            result <- proDA::test_diff(fit,contrast="BvsA",sort_by=NULL)
            result <- result[match(rownames(working),result$name),]
            estimate <- result$diff; pvalue <- result$pval
        } else {
            fit <- limma::eBayes(limma::lmFit(working,design),trend=FALSE)
            estimate <- fit$coefficients[,2];pvalue <- fit$p.value[,2]
        }
        p <- rep(1,nrow(sim$x)); p[indices] <- ifelse(is.finite(pvalue),pvalue,1)
        discoveries <- p.adjust(p,"BH")<.05
        changed <- sim$effect!=0
        targets <- is.na(sim$x[indices,,drop=FALSE])
        finite <- targets & is.finite(working)
        error <- working-sim$latent[indices,,drop=FALSE]
        testable <- is.finite(pvalue)
        c(available=1,retained_fraction=nrow(working)/nrow(sim$x),tested_fraction=sum(testable)/nrow(sim$x),
            reconstructed_missing_fraction=sum(finite)/sum(is.na(sim$x)),
            missing_RMSE=if(any(finite)) sqrt(mean(error[finite]^2)) else NA_real_,
            effect_RMSE=sqrt(safe_mean((estimate-sim$effect[indices])^2)),
            FDP=if(any(discoveries)) mean(!changed[discoveries]) else 0,
            power=sum(discoveries & changed)/sum(changed),discoveries=sum(discoveries),
            proda_converged=as.numeric(convergence$successful),proda_iterations=convergence$iterations,
            proda_error=convergence$error,proda_default_converged=if(is.na(convergence$successful))NA_real_ else as.numeric(convergence$successful && convergence$iterations<=20))
    },warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")}),error=function(e)e)
    if(inherits(ans,"error")) return(data.frame(available=0,retained_fraction=0,tested_fraction=0,
        reconstructed_missing_fraction=0,missing_RMSE=NA_real_,effect_RMSE=NA_real_,FDP=NA_real_,power=0,
        discoveries=0,proda_converged=as.numeric(convergence$successful),proda_iterations=convergence$iterations,proda_error=convergence$error,
        proda_default_converged=if(is.na(convergence$successful))NA_real_ else as.numeric(convergence$successful && convergence$iterations<=20),
        warnings=paste(unique(warnings),collapse=" | "),reason=conditionMessage(ans)))
    cbind(as.data.frame(as.list(ans)),warnings=paste(unique(warnings),collapse=" | "),reason="")
}

if(Sys.getenv("IMPUTEFINDER_DOWNSTREAM_HELPERS_ONLY")!="1") {
{
    output <- list(); k <- 0L
    workers <- as.integer(Sys.getenv("PUBLICATION_WORKERS","1"))
    stopifnot(length(workers)==1L,is.finite(workers),workers>=1L,workers<=8L)
    if(.Platform$OS.type=="windows") workers <- 1L
    for(scenario in c("MCAR","mixed","heteroscedastic_batch")) {
        message("Downstream: ",scenario)
        batches <- parallel::mclapply(1:30,function(seed) {
            sim <- simulate_data(scenario,6L,seed)
            do.call(rbind,lapply(methods,function(method)
                cbind(data.frame(scenario=scenario,seed=seed,method=method),evaluate_downstream(sim,method,seed+400L))))
        },mc.cores=workers,mc.set.seed=FALSE)
        if(any(vapply(batches,inherits,logical(1),"try-error"))) stop("A downstream worker failed.")
        k <- k+1L;output[[k]] <- do.call(rbind,batches)
    }
    write_csv(do.call(rbind,output),"downstream.csv")
}
d <- read.csv("publication/results/downstream.csv")
summary <- summarise_metrics(d,c("scenario","method"),c("available","retained_fraction","tested_fraction",
    "reconstructed_missing_fraction","missing_RMSE","effect_RMSE","FDP","power","discoveries"))
bounded <- !summary$metric %in% c("missing_RMSE","effect_RMSE","discoveries")
summary$lower[bounded] <- pmax(0,summary$lower[bounded]); summary$upper[bounded] <- pmin(1,summary$upper[bounded])
write_csv(summary,"downstream-summary.csv")
pairs <- list(); k <- 0L
reference <- d[d$method=="DEP_absence",]
for(method in setdiff(methods,"DEP_absence")) {
    candidate <- d[d$method==method,]
    candidate <- merge(candidate,reference,by=c("scenario","seed"),suffixes=c("_candidate","_reference"))
    delta <- candidate[c("scenario","seed")];delta$method <- method
    for(metric in c("available","retained_fraction","missing_RMSE","effect_RMSE","FDP","power")) {
        delta[[metric]] <- candidate[[paste0(metric,"_candidate")]]-candidate[[paste0(metric,"_reference")]]
    }
    k <- k+1L;pairs[[k]] <- delta
}
paired <- do.call(rbind,pairs)
write_csv(summarise_metrics(paired,c("scenario","method"),c("available","retained_fraction","missing_RMSE","effect_RMSE","FDP","power"),availability_binomial=FALSE),"downstream-paired-vs-DEP.csv")
writeLines(capture.output(sessionInfo()),"publication/results/downstream-session-info.txt")
}
