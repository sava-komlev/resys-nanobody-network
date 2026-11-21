source("scripts/common.R");source("scripts/network_functions.R")
assert_private("miic")
x<-readRDS("data/input_unpooled.rds");m<-read.csv("data/metadata_clustered.csv")
groups<-split(seq_len(nrow(x)),m$cluster90)
B<-50L;allres<-list();fitsizes<-list()
for(mode in c("row_bootstrap","cluster_bootstrap","representative_bootstrap")) {
 for(b in seq_len(B)) {
  set.seed(30000+b+match(mode,c("row_bootstrap","cluster_bootstrap","representative_bootstrap"))*1000)
  if(mode=="row_bootstrap") ix<-sample.int(nrow(x),nrow(x),replace=TRUE)
  else {
    gg<-groups[sample.int(length(groups),length(groups),replace=TRUE)]
    if(mode=="cluster_bootstrap") ix<-unlist(gg,use.names=FALSE)
    else ix<-vapply(gg,function(v) v[sample.int(length(v),1)],integer(1))
  }
  dd<-pool_categories(x[ix,]);f<-fit_miic(dd,shuffles=10);ee<-miic_edges(f)
  allres[[paste(mode,b)]]<-data.frame(mode=mode,iteration=b,x=ee$x,y=ee$y,key=ee$key)
  fitsizes[[paste(mode,b)]]<-data.frame(mode=mode,iteration=b,n=nrow(dd),positives=sum(dd$selection_label=="1"),edges=nrow(ee))
  if(b%%10==0) {log_line(mode,b,"of",B);saveRDS(allres,"results/bootstrap_edges.rds")}
 }
}
raw<-do.call(rbind,allres);save_csv(raw,"results/bootstrap_edges.csv")
st<-aggregate(iteration~mode+key+x+y,raw,length);names(st)[names(st)=="iteration"]<-"retained";st$frequency<-st$retained/B
save_csv(st,"results/edge_stability.csv");save_csv(do.call(rbind,fitsizes),"results/bootstrap_run_sizes.csv")
# Label permutations on a representative-per-cluster dataset. Null preserves features.
reps<-readRDS("data/representatives.rds")$cluster90
dd<-pool_categories(x[reps,]);nulls<-list()
for(b in 1:50) {
 set.seed(50000+b);z<-dd;z$selection_label<-sample(z$selection_label)
 f<-fit_miic(z,shuffles=10);ee<-miic_edges(f)
 nulls[[b]]<-data.frame(iteration=b,label_degree=sum(ee$x=="selection_label"|ee$y=="selection_label"),edges=nrow(ee))
 if(b%%10==0)log_line("Label null",b,"of 50")
}
save_csv(do.call(rbind,nulls),"results/null_label_networks.csv")
# Sensitivity to rare-state pooling and length representation, same representative rows.
ss<-list();sse<-list()
for(threshold in c(.005,.01,.02)) for(length_mode in c("binned","continuous")) {
 dd<-x[reps,]
 if(length_mode=="continuous")dd$cdr3_length<-nchar(m$cdr3[reps])
 dd<-pool_categories(dd,threshold)
 if(length_mode=="continuous")dd$cdr3_length<-nchar(m$cdr3[reps])
 f<-fit_miic(dd,shuffles=10);ee<-miic_edges(f);nm<-paste(threshold,length_mode,sep="_")
 ss[[nm]]<-data.frame(threshold=threshold,length_mode=length_mode,edges=nrow(ee),
   label_degree=sum(ee$x=="selection_label"|ee$y=="selection_label"),
   label_neighbours=paste(sort(unique(c(ee$x[ee$y=="selection_label"],ee$y[ee$x=="selection_label"]))),collapse=";"))
 sse[[nm]]<-ee
}
save_csv(do.call(rbind,ss),"results/encoding_sensitivity.csv");saveRDS(sse,"results/encoding_sensitivity_edges.rds")
log_line("Robustness complete")
