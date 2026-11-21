source("scripts/common.R");source("scripts/network_functions.R")
assert_private(c("miic","igraph"))
x<-readRDS("data/input_unpooled.rds"); reps<-readRDS("data/representatives.rds")
datasets<-c(list(full=seq_len(nrow(x))),reps)
allstats<-list();fits<-list();all_edges<-list()
for(nm in names(datasets)) {
  ix<-datasets[[nm]];dat<-pool_categories(x[ix,]);stopifnot(nrow(dat)>=1000,ncol(dat)>=20)
  log_line(nm,"MIIC starting",nrow(dat),"rows",sum(dat$selection_label=="1"),"positive")
  t0<-proc.time()[3];f<-fit_miic(dat);secs<-proc.time()[3]-t0
  e<-miic_edges(f);saveRDS(f,paste0("results/miic_",nm,".rds"));save_csv(f$summary,paste0("results/miic_summary_",nm,".csv"))
  fits[[nm]]<-f;all_edges[[paste(nm,"MIIC",sep="_")]]<-e
  allstats[[paste(nm,"MIIC")]]<-network_stats(e,names(dat),nm,"MIIC",secs)
  log_line(nm,"MIIC edges",nrow(e),"seconds",round(secs,2),"label neighbours",paste(allstats[[paste(nm,"MIIC")]]$label_neighbours))
  t0<-proc.time()[3];pw<-pairwise_network(dat,B=199);secs<-proc.time()[3]-t0
  save_csv(pw,paste0("results/pairwise_",nm,".csv"));sig<-pw[pw$q_BH<=.05,]
  all_edges[[paste(nm,"MI",sep="_")]]<-sig
  allstats[[paste(nm,"MI")]]<-network_stats(sig,names(dat),nm,"MI",secs)
  mat<-mi_matrix(dat)
  for(tau in c(0,.1)) {
    method<-paste0("ARACNE_",tau);e<-aracne_edges(mat,sig,tolerance=tau)
    all_edges[[paste(nm,method,sep="_")]]<-e
    allstats[[paste(nm,method)]]<-network_stats(e,names(dat),nm,method,secs)
  }
  save_csv(do.call(rbind,allstats),"results/network_summary.csv")
  saveRDS(all_edges,"results/network_edges.rds")
  log_line(nm,"complete")
}
saveRDS(fits,"results/miic_fits.rds")
ip<-installed.packages(lib.loc=private)
save_csv(ip[,c("Package","Version","Built"),drop=FALSE],"environment/package_versions.csv")
writeLines(capture.output(sessionInfo()),"environment/sessionInfo.txt")
jsonlite::write_json(list(R=as.character(getRversion()),platform=R.version$platform,
 packages=lapply(setNames(ip[,"Version"],ip[,"Package"]),function(v) list(Version=v))),
 "environment/versions.lock.json",pretty=TRUE,auto_unbox=TRUE)
print(do.call(rbind,allstats))
