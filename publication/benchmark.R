# Reproduce from repository root: R_LIBS_USER=<library> Rscript publication/benchmark.R
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- dirname(dirname(normalizePath(script)))
setwd(root)
for (f in sort(list.files("R", "[.]R$", full.names = TRUE))) source(f)
dir.create("publication/results", showWarnings = FALSE)
options(digits = 16)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
write_csv <- function(x, name) write.csv(x, file.path("publication/results", name), row.names = FALSE, na = "NA")
safe_mean <- function(x) if (length(x) && any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
safe_div <- function(x, y) if (y > 0) x / y else NA_real_

simulate_data <- function(scenario, nrep, seed, n = 600L) {
    set.seed(100000L + seed + match(scenario, scenarios) * 1000L + nrep * 10000L)
    mu <- rnorm(n, 14, 2.5)
    effect <- rep(0, n)
    effect[sample(n, n / 5)] <- sample(c(-1.5, 1.5), n / 5, TRUE)
    groups <- rep(c("A", "B"), each = nrep)
    means <- cbind(A = mu, B = mu + effect)
    sd <- if (scenario == "heteroscedastic_batch") 0.2 + 0.12 * abs(mu - 16) else rep(0.4, n)
    latent <- means[, groups] + matrix(rnorm(n * 2L * nrep), n) * sd
    if (scenario == "heteroscedastic_batch") {
        shift <- rep(seq(-0.7, 0.7, length.out = nrep), 2L)
        latent <- sweep(latent, 2L, shift, "+")
    }
    covariate <- mu + rnorm(n)
    pm <- switch(scenario, MCAR = 0.25,
        covariate_MAR = matrix(rep(0.04 + 0.7 * plogis((12 - covariate)/0.8), 2*nrep), n),
        mixed = 0.12, heteroscedastic_batch = 0.1, 0.03)
    pn <- switch(scenario,
        smooth_MNAR = 0.85 * plogis((12-latent)/0.8),
        hard_censor = 0.95 * (latent < 12),
        dual_transition = 0.45*plogis((11-latent)/0.3) + 0.4*plogis((15-latent)/0.3),
        high_tail_MNAR = 0.8*plogis((latent-15)/0.5),
        heteroscedastic_batch = 0.85*plogis((12-latent)/0.8),
        mixed = 0.85*plogis((12-latent)/0.8), 0)
    cause <- matrix(0L, n, 2*nrep)
    cause[matrix(runif(length(latent)), n) < pm] <- 1L
    cause[cause == 0L & matrix(runif(length(latent)), n) < pn] <- 2L
    onoff <- sample(n, n / 20)
    latent[onoff, groups == "A"] <- 9 + matrix(rnorm(length(onoff)*nrep, sd=0.4), length(onoff))
    latent[onoff, groups == "B"] <- 16 + matrix(rnorm(length(onoff)*nrep, sd=0.4), length(onoff))
    cause[onoff, groups == "A"] <- 3L
    cause[onoff, groups == "B"] <- 0L
    effect[onoff] <- 7
    dimnames(latent) <- list(sprintf("p%04d", seq_len(n)), paste0(groups, rep(seq_len(nrep), 2)))
    dimnames(cause) <- dimnames(latent)
    x <- latent; x[cause > 0L] <- NA_real_
    list(x=x, latent=latent, cause=cause, groups=groups, effect=effect, onoff=onoff,
        observed_covariate=covariate)
}

run_fit <- function(x, groups, rescue="none", cutoffs=NULL, seed=1L) {
    tryCatch(classify_missingness(x, group=groups, rescue=rescue,
        cutoffs=cutoffs, ambiguity="lowest", seed=seed), error=function(e) e)
}

state_matrix <- function(fit, x, groups) {
    out <- matrix("unassigned", nrow(x), length(unique(groups)),
        dimnames=list(rownames(x), sort(unique(groups))))
    if (inherits(fit,"error")) return(out)
    cl <- fit$classifications[fit$classifications$feature %in% rownames(x),,drop=FALSE]
    out[cbind(match(cl$feature,rownames(out)), match(cl$condition,colnames(out)))] <- cl$state
    out
}

baseline <- function(x, groups, method) {
    cond <- sort(unique(groups))
    count <- sapply(cond, function(g) rowSums(!is.na(x[,groups==g,drop=FALSE])))
    means <- sapply(cond, function(g) rowMeans(x[,groups==g,drop=FALSE], na.rm=TRUE))
    threshold <- if (method=="pooled_q25") quantile(rowMeans(x,na.rm=TRUE),.25,na.rm=TRUE) else apply(means,2,quantile,.25,na.rm=TRUE)
    state <- matrix("MAR",nrow(x),length(cond),dimnames=list(rownames(x),cond))
    if (method=="pooled_q25") {
        pooled <- rowMeans(x,na.rm=TRUE)
        state[is.finite(pooled) & pooled<threshold,] <- "MNAR"
    }
    if (method=="condition_q25") state[sweep(means,2,threshold,"<") & is.finite(means)] <- "MNAR"
    if (method=="observed_majority") state[sweep(count,2,as.numeric(table(groups))/2,"<=")] <- "MNAR"
    if (method=="DEP_absence") state[apply(count==0,1,any),] <- "MNAR"
    if (method!="DEP_absence") state[count==0] <- "absent"
    list(state=state, retained=rowSums(count)>0, count=count)
}

score <- function(sim, state, retained, available=TRUE) {
    routed <- state[,sim$groups,drop=FALSE]
    eligible <- sim$cause %in% c(1L,2L)
    assigned <- routed %in% c("MAR","MNAR")
    correct <- (routed=="MAR" & sim$cause==1L) | (routed=="MNAR" & sim$cause==2L)
    actionable <- matrix(retained,nrow(sim$x),ncol(sim$x))
    c(available=as.numeric(available), routing_coverage=safe_div(sum(assigned & eligible),sum(eligible)),
        correct_all_missing=safe_div(sum(correct & eligible),sum(eligible)),
        actionable_coverage=safe_div(sum(assigned & eligible & actionable),sum(eligible)),
        actionable_correct=safe_div(sum(correct & eligible & actionable),sum(eligible)),
        selective_accuracy=safe_div(sum(correct & eligible),sum(assigned & eligible)),
        MAR_recall=safe_div(sum(routed=="MAR" & sim$cause==1L),sum(sim$cause==1L)),
        MNAR_recall=safe_div(sum(routed=="MNAR" & sim$cause==2L),sum(sim$cause==2L)),
        retained_fraction=mean(retained), onoff_retention=mean(retained[sim$onoff]),
        missing_cells=sum(eligible), assigned_cells=sum(assigned & eligible))
}

summarise_metrics <- function(data, keys, metrics, availability_binomial=TRUE) {
    sets <- split(seq_len(nrow(data)), interaction(data[keys], drop=TRUE, lex.order=TRUE))
    do.call(rbind,lapply(sets,function(idx) do.call(rbind,lapply(metrics,function(metric) {
        x <- data[[metric]][idx]; x <- x[is.finite(x)]; n <- length(x)
        estimate <- if(n) mean(x) else NA_real_
        margin <- if(n>1) qt(.975,n-1)*sd(x)/sqrt(n) else NA_real_
        lower <- estimate-margin; upper <- estimate+margin
        if(metric=="available" && n && availability_binomial) {
            z <- qnorm(.975); center <- (estimate+z^2/(2*n))/(1+z^2/n)
            width <- z*sqrt(estimate*(1-estimate)/n+z^2/(4*n^2))/(1+z^2/n)
            lower <- center-width; upper <- center+width
        }
        cbind(data[idx[1],keys,drop=FALSE],data.frame(metric=metric,n=n,total=length(idx),
            estimate=estimate,lower=lower,upper=upper))
    }))))
}

scenarios <- c("MCAR","covariate_MAR","smooth_MNAR","hard_censor",
    "dual_transition","high_tail_MNAR","heteroscedastic_batch","mixed")
if (Sys.getenv("IMPUTEFINDER_HELPERS_ONLY") != "1") {
{
    records <- list(); failures <- list(); k <- 0L
    for (scenario in scenarios) for (nrep in c(3L,6L,10L)) {
        message(scenario," n=",nrep)
        for (seed in 1:30) {
            sim <- simulate_data(scenario,nrep,seed)
            for (method in c("automatic_minimum","automatic_none","manual_11","manual_12","manual_13",
                "pooled_q25","condition_q25","observed_majority","all_MAR","DEP_absence")) {
                available <- TRUE; ambiguous <- FALSE; cut <- c(NA_real_,NA_real_); reason <- ""
                if (startsWith(method,"automatic") || startsWith(method,"manual")) {
                    cutoff <- if(startsWith(method,"manual")) setNames(rep(as.numeric(sub("manual_","",method)),2),c("A","B")) else NULL
                    fit <- run_fit(sim$x,sim$groups,if(method=="automatic_minimum") "minimum" else "none",cutoff)
                    available <- !inherits(fit,"error")
                    state <- state_matrix(fit,sim$x,sim$groups)
                    retained <- if(available) fit$feature_status$retained[match(rownames(sim$x),fit$feature_status$feature)] else rep(FALSE,nrow(sim$x))
                    if(available) {
                        cut <- unname(fit$cutoffs)
                        ambiguous <- any(vapply(fit$cutoff_diagnostics,function(d) isTRUE(d$quality$derivative$credible_lobe_count>1),logical(1)))
                    } else reason <- conditionMessage(fit)
                } else {
                    b <- baseline(sim$x,sim$groups,method); state <- b$state; retained <- b$retained
                }
                k <- k+1L
                records[[k]] <- cbind(data.frame(scenario=scenario,nrep=nrep,seed=seed,method=method,
                    ambiguous=ambiguous,cutoff_A=cut[1],cutoff_B=cut[2],reason=reason),as.data.frame(as.list(score(sim,state,retained,available))))
            }
        }
    }
    write_csv(do.call(rbind,records),"simulation.csv")
}
simulation <- read.csv("publication/results/simulation.csv")
metrics <- c("available","routing_coverage","correct_all_missing","actionable_coverage","actionable_correct","selective_accuracy","MAR_recall","MNAR_recall","retained_fraction","onoff_retention")
summary <- summarise_metrics(simulation,c("scenario","nrep","method"),metrics)
summary$lower <- pmax(0,summary$lower); summary$upper <- pmin(1,summary$upper)
write_csv(summary,"simulation-summary.csv")
writeLines(capture.output(sessionInfo()),"publication/results/session-info.txt")
message("Simulation complete: ",nrow(simulation)," method-datasets.")
}
