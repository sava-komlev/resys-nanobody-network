# Stability without replacement, avoiding duplicated contingency-table records.
source("scripts/common.R");source("scripts/network_functions.R")
x<-readRDS("data/input_unpooled.rds");m<-read.csv("data/metadata_clustered.csv")
groups<-split(seq_len(nrow(x)),m$cluster90);B<-50L
edges<-list();runs<-list()
for(mode in c("row_subsample","cluster_subsample","representative_subsample")) for(b in 1:B) {
 set.seed(90000+1000*match(mode,c("row_subsample","cluster_subsample","representative_subsample"))+b)
 if(mode=="row_subsample")ix<-sample.int(nrow(x),floor(.8*nrow(x)))
 else {
   gg<-groups[sample.int(length(groups),floor(.8*length(groups)))]
   ix<-if(mode=="cluster_subsample")unlist(gg,use.names=FALSE)else vapply(gg,function(v) v[sample.int(length(v),1)],integer(1))
 }
 stopifnot(!anyDuplicated(ix),length(ix)>1000)
 f<-fit_miic(pool_categories(x[ix,]),shuffles=10);ee<-miic_edges(f)
 edges[[paste(mode,b)]]<-data.frame(mode=mode,iteration=b,x=ee$x,y=ee$y,key=ee$key)
 runs[[paste(mode,b)]]<-data.frame(mode=mode,iteration=b,n=length(ix),positives=sum(m$label[ix]),
   edges=nrow(ee),label_degree=sum(ee$x=="selection_label"|ee$y=="selection_label"))
 if(b%%10==0)log_line(mode,b,"of 50")
}
raw<-do.call(rbind,edges);save_csv(raw,"results/subsampling_edges.csv")
st<-aggregate(iteration~mode+key+x+y,raw,length);names(st)[names(st)=="iteration"]<-"retained";st$frequency<-st$retained/B
save_csv(st,"results/subsampling_stability.csv");save_csv(do.call(rbind,runs),"results/subsampling_runs.csv")
log_line("Subsampling complete")
