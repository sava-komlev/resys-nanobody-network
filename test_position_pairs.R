# Run after both donors' prediction checks, from the project root.
source("scripts/common.R")
source("scripts/prediction_metrics.R")
# First two positions by IMGT number in each loop of the original 24-position panel.
pairs <- list(CDR1=c("28","29"), CDR2=c("56","57"), CDR3=c("105","106"))
models <- c("single_positions", names(pairs), "all_three_pairs")
# Compare one joint factor to two independent factors for the same letters.
# Add 50 pseudo-observations distributed according to the independent model.
prior_mass <- 50
encode <- function(v, train, length_feature) {
  freq <- prop.table(table(v[train]))
  keep <- names(freq)[freq >= .01 | names(freq) == "-"]
  if (length_feature) keep <- c("<=12","13-15","16-18","19-21",">=22")
  factor(ifelse(v %in% keep, v, "OTHER"), levels=unique(c(keep,"OTHER")))
}
records <- list()
for (donor in c("Puta","Christy")) {
  dir <- if (donor == "Puta") "." else "second_donor"
  read_table <- function(path) read.csv(file.path(dir,path),check.names=FALSE)
  meta <- read_table("data/metadata_clustered.csv")
  aln <- read_table("data/avida_imgt_H.csv")
  stopifnot(all(meta$variant_id %in% aln$Id))
  aln <- aln[match(meta$variant_id, aln$Id), ]
  aln$cdr3_length <- as.character(cut(nchar(meta$cdr3),c(-Inf,12,15,18,21,Inf),
      labels=c("<=12","13-15","16-18","19-21",">=22")))
  splits <- read_table("results/predictions.csv")
  panels <- read_table("results/prediction_training_features.csv")
  for (mode in c("heldout_cluster90","heldout_cluster80")) for (rep in 1:3) {
    old <- splits[splits$mode==mode & splits$repetition==rep, ]
    stopifnot(identical(meta$variant_id,old$variant_id), identical(meta$label,old$label))
    group <- meta[[sub("heldout_","",mode)]]
    scores <- matrix(NA_real_, nrow(meta), length(models),dimnames=list(NULL,models))
    for (fold in 1:5) {
      test <- which(old$fold==fold)
      train <- which(old$fold!=fold)
      stopifnot(!any(group[test] %in% group[train]))
      panel <- panels$selected[panels$mode==mode & panels$repetition==rep & panels$fold==fold]
      vars <- c(strsplit(panel,";",fixed=TRUE)[[1]],"cdr3_length")
      stopifnot(all(unlist(pairs) %in% vars))
      y <- factor(meta$label[train], levels=0:1)
      n <- as.numeric(table(y))
      base <- rep(log((n[2]+1)/(n[1]+1)), length(test))
      states <- probs <- list()
      for (v in vars) {
        a <- encode(aln[[v]], train, v=="cdr3_length")
        p <- sweep(table(a[train], y)+1, 2, n+nlevels(a), "/")
        base <- base + log(p[a[test],2]/p[a[test],1])
        states[[v]] <- a
        probs[[v]] <- p
      }
      scores[test,1] <- plogis(base)
      increments <- matrix(0,length(test),length(pairs),dimnames=list(NULL,names(pairs)))
      for (pair in names(pairs)) {
        vars2 <- pairs[[pair]]
        a <- states[[vars2[1]]]; b <- states[[vars2[2]]]
        for (label in 0:1) {
          k <- label+1
          ix <- train[meta$label[train]==label]
          independent <- outer(probs[[vars2[1]]][,k], probs[[vars2[2]]][,k])
          joint <- (as.matrix(table(a[ix],b[ix])) + prior_mass*independent)/(n[k]+prior_mass)
          ij <- cbind(as.integer(a[test]),as.integer(b[test]))
          increments[,pair] <- increments[,pair] + (2*label-1)*log(joint[ij]/independent[ij])
        }
        scores[test,pair] <- plogis(base+increments[,pair])
      }
      scores[test,"all_three_pairs"] <- plogis(base+rowSums(increments))
    }
    stopifnot(all(is.finite(scores)), max(abs(old$probability-scores[,1])) < 1e-10)
    for (model in models) for(weighting in c("sequence","group")) {
      w <- if(weighting=="group") 1/as.numeric(table(group)[as.character(group)]) else rep(1,length(group))
      records[[length(records)+1L]] <- data.frame(donor,mode,repetition=rep,model,weighting,
        metrics(meta$label,scores[,model],w))
    }
  }
}
save_csv(do.call(rbind,records),"results/position_pair_validation.csv")
