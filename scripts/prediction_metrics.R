metrics<-function(y,p,w=rep(1,length(y))) {
 ord<-order(-p);y<-y[ord];p<-p[ord];w<-w[ord]
 ends<-which(!duplicated(p,fromLast=TRUE));tp<-cumsum(w*y)[ends];fp<-cumsum(w*(1-y))[ends]
 recall<-tp/sum(w*y);fpr<-fp/sum(w*(1-y));prec<-tp/(tp+fp)
 ap<-sum(diff(c(0,recall))*prec)
 auc<-sum(diff(c(0,fpr))*(c(0,recall[-length(recall)])+recall)/2)
 data.frame(roc_auc=auc,average_precision=ap,prevalence=sum(w*y)/sum(w),
   brier=sum(w*(p-y)^2)/sum(w))
}
