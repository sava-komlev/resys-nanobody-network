source("scripts/common.R");source("scripts/network_functions.R")
x<-readRDS("data/input_unpooled.rds");m<-read.csv("data/metadata_clustered.csv")
reps<-readRDS("data/representatives.rds")$cluster90
n<-length(reps);npos<-sum(m$label[reps]);groups<-split(seq_len(nrow(x)),m$cluster90)
res<-list();eres<-list()
for(mode in c("random_size_matched","random_size_and_class_matched","alternative_cluster_representatives")) for(b in 1:30) {
 set.seed(60000+1000*match(mode,c("random_size_matched","random_size_and_class_matched","alternative_cluster_representatives"))+b)
 if(mode=="random_size_matched")ix<-sample.int(nrow(x),n)
 if(mode=="random_size_and_class_matched")ix<-c(sample(which(m$label==1),npos),sample(which(m$label==0),n-npos))
 if(mode=="alternative_cluster_representatives")ix<-vapply(groups,function(v) v[sample.int(length(v),1)],integer(1))
 dd<-pool_categories(x[ix,]);ff<-fit_miic(dd,shuffles=10);ee<-miic_edges(ff)
 nn<-sort(unique(c(ee$x[ee$y=="selection_label"],ee$y[ee$x=="selection_label"])))
 res[[paste(mode,b)]]<-data.frame(mode=mode,iteration=b,n=length(ix),positives=sum(m$label[ix]),
   edges=nrow(ee),label_degree=length(nn),label_neighbours=paste(nn,collapse=";"))
 eres[[paste(mode,b)]]<-ee
 if(b%%10==0)log_line(mode,b,"of 30")
}
save_csv(do.call(rbind,res),"results/matched_controls.csv");saveRDS(eres,"results/matched_control_edges.rds")
# Descriptive cluster and CDR3-length composition, avoiding mutation-effect claims.
cs<-do.call(rbind,lapply(groups,function(ix) data.frame(cluster=m$cluster90[ix[1]],n=length(ix),
 positive=sum(m$label[ix]),fraction=mean(m$label[ix]),cdr3_length=nchar(m$cdr3[ix[1]]))))
cs<-cs[order(-cs$positive,-cs$n),];save_csv(cs,"results/cluster_composition.csv")
rr<-list()
for(nm in c("full","cluster90")) {
 ix<-if(nm=="full")seq_len(nrow(x))else reps
 d<-x[ix,];dd<-aggregate(as.integer(selection_label)~cdr3_length,d,mean)
 names(dd)[2]<-"positive_fraction";dd$n<-as.integer(table(factor(d$cdr3_length,levels=dd$cdr3_length)))
 dd$dataset<-nm;rr[[nm]]<-dd
}
save_csv(do.call(rbind,rr),"results/length_outcomes.csv")
log_line("Matched controls complete")
