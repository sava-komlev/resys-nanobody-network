pool_categories <- function(dat, threshold=.01) {
 out <- dat
 for(p in names(out)) {
  v <- as.character(out[[p]])
  if(p != "selection_label" && p != "cdr3_length") {
    f <- prop.table(table(v)); rare <- names(f)[f < threshold & names(f) != "-"]
    v[v %in% rare] <- "OTHER"
  }
  out[[p]] <- factor(v)
 }
 out
}
mi_integer <- function(a,b,ka=max(a),kb=max(b)) {
 n <- length(a); counts <- matrix(tabulate(a+ka*(b-1L),nbins=ka*kb),ka,kb)
 expected <- outer(rowSums(counts),colSums(counts))/n
 nz <- counts>0
 sum(counts[nz]*log2(counts[nz]/expected[nz]))/n
}
mi_matrix <- function(dat) {
 v <- lapply(dat,function(z) as.integer(factor(z))); k <- vapply(v,max,integer(1))
 p <- length(v); out <- matrix(0,p,p,dimnames=list(names(dat),names(dat)))
 for(i in seq_len(p-1L)) for(j in (i+1L):p) out[i,j] <- out[j,i] <- mi_integer(v[[i]],v[[j]],k[i],k[j])
 out
}
pairwise_network <- function(dat,B=199) {
 v <- lapply(dat,function(z) as.integer(factor(z))); k <- vapply(v,max,integer(1))
 pairs <- combn(seq_along(v),2)
 out <- lapply(seq_len(ncol(pairs)),function(ii) {
  i <- pairs[1,ii]; j <- pairs[2,ii];obs <- mi_integer(v[[i]],v[[j]],k[i],k[j])
  null <- replicate(B,mi_integer(v[[i]],sample(v[[j]]),k[i],k[j]))
  data.frame(x=names(dat)[i],y=names(dat)[j],mi_bits=obs,
    mi_excess=obs-mean(null),p_perm=(1+sum(null>=obs-1e-12))/(B+1))
 })
 out <- do.call(rbind,out);out$q_BH <- p.adjust(out$p_perm,"BH");out$key <- edge_key(out$x,out$y)
 out
}
aracne_edges <- function(mi,significant,tolerance=0) {
 nodes <- rownames(mi); adjacency <- matrix(FALSE,nrow(mi),ncol(mi),dimnames=dimnames(mi))
 for(i in seq_len(nrow(significant))) adjacency[significant$x[i],significant$y[i]] <- adjacency[significant$y[i],significant$x[i]] <- TRUE
 remove <- matrix(FALSE,nrow(mi),ncol(mi)); triples <- combn(seq_along(nodes),3)
 for(i in seq_len(ncol(triples))) {
   a<-triples[1,i];b<-triples[2,i];c<-triples[3,i]
   if(!all(c(adjacency[a,b],adjacency[a,c],adjacency[b,c]))) next
   pairs<-rbind(c(a,b),c(a,c),c(b,c));vals<-c(mi[a,b],mi[a,c],mi[b,c])
   for(e in 1:3) if(vals[e] < (1-tolerance)*min(vals[-e])-1e-12)
     remove[pairs[e,1],pairs[e,2]] <- remove[pairs[e,2],pairs[e,1]] <- TRUE
 }
 ix <- which(upper.tri(adjacency)&adjacency&!remove,arr.ind=TRUE)
 out<-data.frame(x=nodes[ix[,1]],y=nodes[ix[,2]])
 out$key<-edge_key(out$x,out$y);out
}
fit_miic <- function(dat, shuffles=10) {
 miic::miic(dat,orientation=FALSE,latent="no",consistent="no",n_threads=1,
            n_shuffles=shuffles,conf_threshold=if(shuffles>0) .01 else 0,verbose=FALSE)
}
miic_edges <- function(fit) {
 out<-fit$summary[fit$summary$type=="P",c("x","y"),drop=FALSE]
 out$key<-edge_key(out$x,out$y);out
}
network_stats <- function(edges,nodes,name,method,seconds) {
 g<-igraph::graph_from_data_frame(edges[c("x","y")],directed=FALSE,vertices=nodes)
 neighbours <- if("selection_label" %in% igraph::V(g)$name) names(igraph::neighbors(g,"selection_label")) else character()
 data.frame(dataset=name,method=method,nodes=length(nodes),edges=nrow(edges),
   density=igraph::edge_density(g),components=igraph::components(g)$no,
   label_degree=length(neighbours),label_neighbours=paste(sort(neighbours),collapse=";"),seconds=seconds)
}
