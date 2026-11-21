source("scripts/common.R")
assert_private(c("miic","ppcor","Rcpp","scales","igraph","Biostrings"))
x<-readRDS("data/input_unpooled.rds");m<-read.csv("data/metadata_clustered.csv")
stopifnot(nrow(x)==6039,ncol(x)==26,!anyDuplicated(m$domain),!anyNA(x),
          identical(as.integer(x$selection_label),m$label))
for(cl in c("cluster80","cluster90","cluster95")) {
 stopifnot(all(vapply(split(nchar(m$cdr3),m[[cl]]),function(z)length(unique(z))==1,logical(1))))
}
pred<-read.csv("results/predictions.csv")
stopifnot(nrow(pred)==6039*9,all(is.finite(pred$probability)),all(pred$probability>=0 & pred$probability<=1))
for(md in c("heldout_cluster80","heldout_cluster90")) for(r in 1:3) {
 z<-pred[pred$mode==md & pred$repetition==r,];gr<-m[[sub("heldout_","",md)]][match(z$variant_id,m$variant_id)]
 stopifnot(all(vapply(split(z$fold,gr),function(v)length(unique(v))==1,logical(1))))
}
e<-readRDS("results/network_edges.rds")
for(z in e)stopifnot(!anyDuplicated(z$key),all(z$x %in% names(x)),all(z$y %in% names(x)),!any(z$x==z$y))
ns<-read.csv("results/network_summary.csv")
for(i in seq_len(nrow(ns))) {
 z<-e[[paste(ns$dataset[i],ns$method[i],sep="_")]]
 stopifnot(nrow(z)==ns$edges[i],sum(z$x=="selection_label"|z$y=="selection_label")==ns$label_degree[i])
}
stopifnot(nrow(read.csv("results/matched_controls.csv"))==90,
 nrow(read.csv("results/null_label_networks.csv"))==50,
 nrow(read.csv("results/bootstrap_run_sizes.csv"))==150,
 nrow(read.csv("results/subsampling_runs.csv"))==150,
 nrow(read.csv("results/encoding_sensitivity.csv"))==6)
# Every non-base namespace loaded by our analysis must resolve inside the private library.
suppressPackageStartupMessages(library(miic));suppressPackageStartupMessages(library(igraph));suppressPackageStartupMessages(library(Biostrings))
paths<-vapply(loadedNamespaces(),function(p) find.package(p),character(1))
basepkgs<-rownames(installed.packages(priority="base"))
stopifnot(all(startsWith(paths[!names(paths)%in%basepkgs],private)))
save_csv(data.frame(package=names(paths),path=paths),"environment/loaded_namespace_paths.csv")
cat("PASS: unique domains, input alignment, group separation in validation folds, graph tables, completed-run counts and package isolation.\n")
