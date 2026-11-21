# Secondary diagnostic, not the principal network analysis; fixed Laplace-smoothed
# categorical naive Bayes, no label-based feature selection or tuning.
source("scripts/common.R")
m<-read.csv("data/metadata_clustered.csv")
a<-read.csv("data/avida_imgt_H.csv",check.names=FALSE)
a<-a[match(m$variant_id,a$Id),]
positions<-grep("^[0-9]+[A-Z]*$",names(a),value=TRUE)
candidates<-positions[as.integer(sub("[A-Z]+$","",positions)) %in% c(27:38,56:65,105:117)]
len<-as.character(cut(nchar(m$cdr3),c(-Inf,12,15,18,21,Inf),labels=c("<=12","13-15","16-18","19-21",">=22")))
select_positions<-function(train) {
  fixed <- getOption("resys.fixed_positions")
  if (!is.null(fixed)) return(fixed)
  qc<-do.call(rbind,lapply(candidates,function(p) {
    v<-a[[p]][train];ff<-prop.table(table(v[v!="-"]));data.frame(p=p,coverage=mean(v!="-"),
     entropy=-sum(ff*log2(ff)),major=max(ff),states=sum(ff>=.01))
  }))
  qc<-qc[qc$coverage>=.95 & qc$major<=.95 & qc$states>=2,]
  qc<-qc[order(-qc$entropy,qc$p),];stopifnot(nrow(qc)>=24);qc$p[1:24]
}
predict_nb<-function(train,test,ys) {
 sel<-select_positions(train)
 lp<-matrix(rep(log(c(sum(ys[train]==0)+1,sum(ys[train]==1)+1)/(length(train)+2)),each=length(test)),ncol=2)
 for(p in c(sel,"cdr3_length")) {
   v<-if(p=="cdr3_length")len else a[[p]]
   ff<-prop.table(table(v[train]));keep<-names(ff)[ff>=.01 | names(ff)=="-"]
   if(p=="cdr3_length")keep<-c("<=12","13-15","16-18","19-21",">=22")
   vv<-ifelse(v %in% keep,v,"OTHER");lev<-unique(c(keep,"OTHER"))
   counts<-table(factor(vv[train],levels=lev),factor(ys[train],levels=0:1))
   probs<-sweep(counts+1,2,colSums(counts)+length(lev),"/")
   lp<-lp+log(probs[match(vv[test],rownames(probs)),,drop=FALSE])
 }
 list(p=plogis(lp[,2]-lp[,1]),selected=paste(sel,collapse=";"))
}
source("scripts/prediction_metrics.R")
predictions<-list();summaries<-list();selection<-list();size<-table(m$cluster90)
stopifnot(abs(metrics(c(0,1),c(.1,.9))$roc_auc-1)<1e-12,
          abs(metrics(c(0,1),c(.5,.5))$roc_auc-.5)<1e-12,
          abs(metrics(c(0,1),c(.5,.5))$average_precision-.5)<1e-12)
for(mode in c("random_rows","heldout_cluster90","heldout_cluster80")) for(repetition in 1:3) {
 set.seed(80000+repetition)
 if(mode=="random_rows")fold<-sample(rep(1:5,length.out=nrow(m)))
 else {
   gr<-m[[sub("heldout_","",mode)]];ug<-unique(gr)
   assignment<-sample(rep(1:5,length.out=length(ug)));fold<-assignment[match(gr,ug)]
   stopifnot(all(vapply(split(fold,gr),function(z) length(unique(z))==1,logical(1))))
 }
 pp<-numeric(nrow(m))
 for(k in 1:5) {
  train<-which(fold!=k);test<-which(fold==k)
  ans<-predict_nb(train,test,m$label);pp[test]<-ans$p
  selection[[paste(mode,repetition,k)]]<-data.frame(mode=mode,repetition=repetition,fold=k,selected=ans$selected)
 }
 for(weighting in c("per_sequence","equal_cluster90")) {
  w<-if(weighting=="per_sequence")rep(1,nrow(m))else 1/as.numeric(size[as.character(m$cluster90)])
  mm<-metrics(m$label,pp,w);mm$mode<-mode;mm$repetition<-repetition;mm$weighting<-weighting
  summaries[[paste(mode,repetition,weighting)]]<-mm
 }
 predictions[[paste(mode,repetition)]]<-data.frame(variant_id=m$variant_id,mode=mode,repetition=repetition,
   fold=fold,label=m$label,probability=pp)
 log_line(mode,"repeat",repetition,"AUC",round(metrics(m$label,pp)$roc_auc,3),"AP",round(metrics(m$label,pp)$average_precision,3))
}
save_csv(do.call(rbind,summaries),"results/prediction_metrics.csv")
save_csv(do.call(rbind,predictions),"results/predictions.csv")
save_csv(do.call(rbind,selection),"results/prediction_training_features.csv")
log_line("Prediction diagnostic complete")
