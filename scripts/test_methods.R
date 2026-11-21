source("scripts/common.R");source("scripts/network_functions.R")
# Statistical helper checks: independent exhaustive binary table and perfect identity.
a<-rep(c(1L,1L,2L,2L),100);b<-rep(c(1L,2L,1L,2L),100)
stopifnot(abs(mi_integer(a,b))<1e-12,abs(mi_integer(a,a)-1)<1e-12)
mm<-matrix(c(0,.8,.2,.8,0,.6,.2,.6,0),3,3,dimnames=list(letters[1:3],letters[1:3]))
ee<-data.frame(x=c("a","a","b"),y=c("b","c","c"))
ar<-aracne_edges(mm,ee,0)
stopifnot(setequal(ar$key,c("a--b","b--c")))
mm[,]<-1;diag(mm)<-0
stopifnot(nrow(aracne_edges(mm,ee,0))==3)
# MIIC should recover a strong noisy binary chain without the indirect endpoint edge.
set.seed(26);a<-sample(0:1,10000,TRUE);b<-bitwXor(a,rbinom(10000,1,.1));c<-bitwXor(b,rbinom(10000,1,.1))
dd<-data.frame(a=factor(a),b=factor(b),c=factor(c))
ff<-fit_miic(dd,shuffles=0)
stopifnot(setequal(miic_edges(ff)$key,c("a--b","b--c")))
cat("MI calculation, ARACNE pruning/ties and MIIC noisy-chain checks passed.\n")
