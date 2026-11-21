source("scripts/common.R")
assert_private(c("Biostrings", "igraph"))
x <- read.csv("data/avida_network_input.csv", check.names = FALSE)
m <- read.csv("data/avida_variant_metadata.csv")
stopifnot(nrow(x) == getOption("resys.domains", 6039L), ncol(x) == 26, !anyDuplicated(m$domain), identical(x$selection_label,m$label))
unique_cdr3 <- unique(m$cdr3)
lens <- nchar(unique_cdr3)
thresholds <- c(.80, .90, .95)
edges <- setNames(lapply(thresholds, function(z) matrix(integer(),ncol=2)), as.character(thresholds))
log_line("Computing within-length CDR3 Hamming distances for",length(unique_cdr3),"unique strings")
for (L in sort(unique(lens))) {
  ids <- which(lens == L)
  if (length(ids) < 2) next
  dm <- as.matrix(Biostrings::stringDist(Biostrings::AAStringSet(unique_cdr3[ids]),method="hamming"))
  for (t in thresholds) {
    # Single-linkage clusters, same CDR3 length, at least t sequence identity.
    ij <- which(upper.tri(dm) & dm <= floor(L*(1-t)+1e-8), arr.ind=TRUE)
    if (nrow(ij)) edges[[as.character(t)]] <- rbind(edges[[as.character(t)]], cbind(ids[ij[,1]],ids[ij[,2]]))
  }
}
summaries <- list(); representatives <- list()
for (t in thresholds) {
  g <- igraph::make_empty_graph(length(unique_cdr3), directed=FALSE)
  ee <- edges[[as.character(t)]]
  if (nrow(ee)) g <- igraph::add_edges(g,as.vector(t(ee)))
  cluster <- igraph::components(g)$membership[match(m$cdr3,unique_cdr3)]
  cn <- paste0("cluster",round(t*100));m[[cn]] <- cluster
  sizes <- table(cluster)
  bycluster <- split(seq_len(nrow(m)),cluster)
  set.seed(20261030)
  # sample.int avoids R's sample(single_numeric_value) trap.
  reps <- vapply(bycluster,function(ix) ix[sample.int(length(ix),1)],integer(1))
  representatives[[cn]] <- reps
  summaries[[cn]] <- data.frame(identity=t, clusters=length(sizes), max_size=max(sizes),
     singleton_clusters=sum(sizes==1), median_size=median(sizes),
     mixed_label_clusters=sum(vapply(bycluster,function(ix) length(unique(m$label[ix]))>1,logical(1))),
     representatives=length(reps), positive_representatives=sum(m$label[reps]),
     representative_prevalence=mean(m$label[reps]))
}
save_csv(m,"data/metadata_clustered.csv")
save_csv(do.call(rbind,summaries),"results/cluster_summary.csv")
saveRDS(representatives,"data/representatives.rds")
print(do.call(rbind,summaries))
for (p in names(x)) {
 if(p=="cdr3_length") x[[p]] <- as.character(cut(x[[p]],breaks=c(-Inf,12,15,18,21,Inf),
        labels=c("<=12","13-15","16-18","19-21",">=22")))
 else x[[p]] <- as.character(x[[p]])
}
saveRDS(x,"data/input_unpooled.rds")
save_csv(x,"data/input_categorical.csv")
save_csv(data.frame(file=list.files("data",full.names=TRUE),md5=unname(tools::md5sum(list.files("data",full.names=TRUE)))),"results/input_checksums.csv")
log_line("Preparation complete")
