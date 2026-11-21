source("scripts/common.R")
assert_private("igraph")
ns<-read.csv("results/network_summary.csv");cl<-read.csv("results/cluster_summary.csv")
ctrl<-read.csv("results/matched_controls.csv");pr<-read.csv("results/prediction_metrics.csv")
blue<-"#246A91";orange<-"#D56A35";teal<-"#248A80";ink<-"#183340"
png("figures/network_summary.png",width=1800,height=1300,res=150)
par(mfrow=c(2,2),mar=c(5,5,4,1),col.axis=ink,col.lab=ink,col.main=ink)
dat<-rbind(full=c(6039,1311),representatives90=c(3785,216))
bp<-barplot(c(1311/6039,216/3785)*100,names.arg=c("All sequences","One per CDR3 group"),
 col=c(blue,teal),ylim=c(0,28),ylab="Positive labels (%)",main="A. Positive labels concentrate in sequence groups")
text(bp,c(1311/6039,216/3785)*100+1.5,c("21.7%\n1,311 / 6,039","5.7%\n216 / 3,785"),cex=.85)
mm<-matrix(c(17,0,1,1,2,4),nrow=2)
barplot(mm,beside=TRUE,col=c(blue,teal),names.arg=c("MIIC","ARACNE 0","ARACNE 0.1"),
 ylim=c(0,20),ylab="Features connected to selection label",main="B. Network conclusions change with sampling")
legend("topright",c("All sequences","One per 90% CDR3 group"),fill=c(blue,teal),bty="n",cex=.8)
vals<-list(ctrl$label_degree[ctrl$mode=="random_size_matched"],
 ctrl$label_degree[ctrl$mode=="random_size_and_class_matched"],
 ctrl$label_degree[ctrl$mode=="alternative_cluster_representatives"])
boxplot(vals,names=c("Size matched","Size + labels\nmatched","One per group"),col=c(blue,orange,teal),
 ylim=c(0,19),ylab="MIIC links to selection label",main="C. Counts alone do not explain the collapse")
mtext("30 runs per condition; 3,785 sequences each",side=3,line=.3,cex=.8)
pp<-subset(pr,weighting=="per_sequence")
av<-aggregate(cbind(roc_auc,average_precision)~mode,pp,mean)
av<-av[match(c("random_rows","heldout_cluster90","heldout_cluster80"),av$mode),]
matplot(1:3,av[c("roc_auc","average_precision")],type="b",pch=c(16,17),lty=1,col=c(blue,orange),
 xaxt="n",ylim=c(0,1),xlab="Train/test split",ylab="Prediction score",main="D. Holding out sequence groups reduces performance")
axis(1,at=1:3,labels=c("Random rows","90% groups","80% groups"))
abline(h=.5,col=blue,lty=3);abline(h=1311/6039,col=orange,lty=3)
legend("bottomleft",c("ROC AUC","Average precision"),col=c(blue,orange),pch=c(16,17),lty=1,bty="n",cex=.8)
mtext("Fixed naive Bayes; mean of 3 x 5-fold CV. Dotted: chance baselines.",side=3,line=.3,cex=.7)
dev.off()

all_edges<-readRDS("results/network_edges.rds");nodes<-names(readRDS("data/input_unpooled.rds"))
position<-suppressWarnings(as.integer(sub("IMGT_","",nodes)))
grp<-ifelse(is.na(position),ifelse(nodes=="selection_label","Outcome","Length"),
 ifelse(position<40,"CDR1",ifelse(position<70,"CDR2","CDR3")))
pal<-c(CDR1=blue,CDR2=teal,CDR3=orange,Length="#866B9C",Outcome=ink)
ordered<-c(unlist(lapply(c("CDR1","CDR2","CDR3"),function(z) { ix<-which(grp==z);ix[order(position[ix])] })),which(grp=="Length"))
xy<-matrix(0,length(nodes),2);angles<-seq(0,2*pi,length.out=length(ordered)+1)[seq_along(ordered)]+pi/2
xy[ordered,]<-cbind(cos(angles),sin(angles));labels<-sub("IMGT_","",nodes)
labels[nodes=="selection_label"]<-"Selection";labels[nodes=="cdr3_length"]<-"H3 len"
png("figures/miic_network_comparison.png",width=1900,height=1000,res=150)
par(mfrow=c(1,2),mar=c(2,1,4,1))
for(nm in c("full","cluster90")) {
 e<-all_edges[[paste0(nm,"_MIIC")]]
 g<-igraph::graph_from_data_frame(e[c("x","y")],directed=FALSE,vertices=nodes)
 edgecol<-ifelse(e$x=="selection_label"|e$y=="selection_label","#C45429AA","#97ADB766")
 plot(g,layout=xy,vertex.label=labels,vertex.color=pal[grp],vertex.frame.color=NA,
  vertex.label.color="white",vertex.label.cex=.7,vertex.size=ifelse(grp=="Outcome",22,ifelse(grp=="Length",17,12)),
  edge.color=edgecol,edge.width=ifelse(edgecol=="#C45429AA",1.6,.7),rescale=FALSE,
  xlim=c(-1.18,1.18),ylim=c(-1.18,1.18),asp=1,
  main=if(nm=="full")"All 6,039 sequences\n89 edges; 17 outcome links"else"One per 90% CDR3 group: 3,785 sequences\n22 edges; no outcome links")
 legend("bottomleft",legend=c("CDR1","CDR2","CDR3"),pt.bg=pal[c("CDR1","CDR2","CDR3")],pch=21,bty="n",cex=.8)
}
dev.off()

st<-read.csv("results/subsampling_stability.csv");features<-setdiff(nodes,"selection_label")
hm<-matrix(0,length(features),3,dimnames=list(features,c("row_subsample","cluster_subsample","representative_subsample")))
for(i in seq_len(nrow(st))) {
 r<-st[i,];if(r$x!="selection_label"&&r$y!="selection_label")next
 p<-if(r$x=="selection_label")r$y else r$x;hm[p,r$mode]<-r$frequency
}
png("figures/outcome_edge_stability.png",width=1100,height=1500,res=150)
par(mar=c(5,7,4,2))
image(1:3,1:nrow(hm),t(hm[nrow(hm):1,]),col=colorRampPalette(c("#F2F6F7",teal))(100),
 axes=FALSE,xlab="",ylab="",zlim=c(0,1),main="Selection-label edges under 80% subsampling\n50 runs per sampling scheme")
axis(1,1:3,c("Random\nrows","Whole CDR3\ngroups","One per\ngroup"),las=1,cex.axis=.8)
axis(2,1:nrow(hm),rev(gsub("cdr3_length","CDR3 length",sub("IMGT_","Position ",features))),las=1,cex.axis=.8)
for(i in 1:3)for(j in 1:nrow(hm))text(i,j,sprintf("%.0f%%",100*hm[nrow(hm)+1-j,i]),cex=.7,
 col=if(hm[nrow(hm)+1-j,i]>.55)"white"else ink)
dev.off()
save_csv(data.frame(feature=features,hm,check.names=FALSE),"results/outcome_subsampling_stability.csv")
log_line("Figures complete")
