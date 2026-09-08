script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))
setwd(dirname(dirname(normalizePath(script))))
library(ggplot2)
theme_set(theme_minimal(base_size=11,base_family="sans") +
    theme(panel.grid.minor=element_blank(),legend.position="bottom",
        plot.title=element_text(face="bold"),strip.text=element_text(face="bold")))
methods <- c("automatic_minimum","automatic_none","manual_11","manual_12","manual_13",
    "pooled_q25","condition_q25","observed_majority","all_MAR","DEP_absence")
labels <- c("Automatic + rescue","Automatic, no rescue","Manual 11","Manual 12","Manual 13",
    "Pooled Q25","Condition Q25","Observed majority","All MAR","DEP absence rule")
label_method <- function(x) factor(x,levels=rev(methods),labels=rev(labels))
export <- function(p,name,width=10,height=7) {
    for(ext in c("svg","pdf","png")) {
        path <- paste0("publication/figures/",name,".",ext)
        if(ext=="svg") ggsave(path,p,width=width,height=height,device=grDevices::svg)
        else if(ext=="pdf") ggsave(path,p,width=width,height=height,device=grDevices::cairo_pdf)
        else ggsave(path,p,width=width,height=height,dpi=220)
    }
}
s <- read.csv("publication/results/simulation-summary.csv")
d <- subset(s,nrep==6 & metric %in% c("available","actionable_correct"))
d$scenario <- factor(d$scenario,levels=c("MCAR","covariate_MAR","smooth_MNAR","hard_censor",
    "dual_transition","high_tail_MNAR","heteroscedastic_batch","mixed"),
    labels=c("MCAR background","Covariate MAR","Smooth MNAR","Hard censoring",
        "Dual transition","High-tail MNAR","Variable noise + batch","Mixed"))
d$method <- label_method(d$method)
d$metric <- factor(d$metric,levels=c("available","actionable_correct"),
    labels=c("Pipeline availability","Correct actionable routes / original missing cells"))
d$label <- sprintf("%.0f",100*d$estimate)
p <- ggplot(d,aes(scenario,method,fill=estimate))+geom_tile(color="white")+
    geom_text(aes(label=label),size=2.7)+facet_wrap(~metric,ncol=1)+
    scale_fill_gradient(low="#f6eee7",high="#2e8e96",limits=c(0,1),labels=scales::label_percent())+
    labs(x=NULL,y=NULL,fill=NULL,title="Routing depends on the missingness scenario",
        subtitle="Six samples per condition; 30 independent datasets per cell. Labels are percentages.")+
    theme(axis.text.x=element_text(angle=25,hjust=1))
export(p,"simulation",11,9)

d <- read.csv("publication/results/downstream-summary.csv")
d <- subset(d,metric %in% c("FDP","power"))
d$scenario <- factor(d$scenario,levels=c("MCAR","mixed","heteroscedastic_batch"),
    labels=c("MCAR background","Mixed","Variable noise + batch"))
d$metric <- factor(d$metric,levels=c("FDP","power"),labels=c("FDP","Power"))
d$method <- factor(d$method,levels=rev(c("observed_only","proDA","uniform_KNN","uniform_MinProb","DEP_absence","manual_12","automatic_none","automatic_minimum")),
    labels=rev(c("Observed only","proDA","Uniform kNN","Uniform MinProb","DEP absence rule","Manual 12","Automatic, no rescue","Automatic + rescue")))
p <- ggplot(d,aes(estimate,method,color=method))+geom_pointrange(aes(xmin=lower,xmax=upper),orientation="y",na.rm=TRUE)+
    facet_grid(scenario~metric)+geom_vline(data=data.frame(metric="FDP"),aes(xintercept=.05),linetype=2,color="#a34831")+
    scale_x_continuous(labels=scales::label_percent(),limits=c(0,1))+
    labs(x="Mean across datasets (95% Monte Carlo interval)",y=NULL,title="Retention and imputation change discovery trade-offs",
        subtitle="FDP is conditional on successful pipelines; missing points indicate no completed analyses.")+
    theme(legend.position="none")
export(p,"downstream",11,9)

d <- read.csv("publication/results/real-bootstrap.csv")
primary <- subset(read.csv("publication/results/real-summary.csv"),mode=="automatic")
primary <- rbind(data.frame(rescue=primary$rescue,condition="AD",cutoff=primary$cutoff_AD),
    data.frame(rescue=primary$rescue,condition="control",cutoff=primary$cutoff_control))
p <- ggplot(d,aes(cutoff,fill=rescue))+geom_histogram(binwidth=.15,boundary=0,color="white")+
    geom_vline(data=primary,aes(xintercept=cutoff),color="#a34831",linewidth=.8)+
    facet_grid(condition~rescue)+scale_fill_manual(values=c(minimum="#2e8e96",none="#465d83"))+
    labs(x="Selected log2 cutoff",y="Feature-bootstrap draws",title="Real-data cutoffs are not uniquely determined",
        subtitle="100 feature resamples per rescue policy. All fits have multiple candidate transitions; red lines show full-table selections.")+
    theme(legend.position="none")
export(p,"real-cutoffs",10,7)
