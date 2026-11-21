# Run in the project root after both within-donor analyses.
source("scripts/common.R")
panels <- lapply(c(".", "second_donor"), function(dir) {
  net <- read.csv(file.path(dir, "data/avida_network_input.csv"))
  meta <- read.csv(file.path(dir, "data/metadata_clustered.csv"))
  groups <- readRDS(file.path(dir, "data/representatives.rds"))
  stats <- read.csv(file.path(dir, "results/network_summary.csv"))
  stopifnot(!anyDuplicated(meta$domain), identical(net$selection_label, meta$label),
            nrow(net) == nrow(meta), nrow(stats) == 16)
  list(net = net, meta = meta, groups = groups, stats = stats)
})
names(panels) <- c("Puta", "Christy")
stopifnot(identical(names(panels$Puta$net), names(panels$Christy$net)))
controls <- read.csv("second_donor/results/matched_controls.csv")
stopifnot(nrow(controls) == 90L, all(table(controls$mode) == 30L))
matched <- subset(controls, mode == "random_size_and_class_matched")
ref <- panels$Christy$groups$cluster90
stopifnot(all(controls$n == length(ref)),
          all(matched$positives == sum(panels$Christy$meta$label[ref])))
prediction <- read.csv("second_donor/results/prediction_metrics.csv")
stopifnot(nrow(prediction) == 18L, all(is.finite(prediction$roc_auc)))
rows <- list()
for (donor in names(panels)) {
  p <- panels[[donor]]
  for (name in c("full", "cluster95", "cluster90", "cluster80")) {
    ix <- if (name == "full") seq_len(nrow(p$net)) else p$groups[[name]]
    if (name != "full") stopifnot(!anyDuplicated(p$meta[ix,name]),
        length(ix) == length(unique(p$meta[[name]])))
    z <- p$net[ix, ]
    pos <- z$selection_label == 1
    stat <- subset(p$stats, dataset == name & method == "MIIC")
    rows[[length(rows)+1L]] <- data.frame(donor, dataset=name, n=nrow(z), positive=sum(pos),
      links=stat$label_degree, features=stat$label_neighbours,
      r28_positive=mean(z$IMGT_28[pos]=="R"), r28_negative=mean(z$IMGT_28[!pos]=="R"),
      w109_positive=mean(z$IMGT_109[pos]=="W"), w109_negative=mean(z$IMGT_109[!pos]=="W"),
      length_positive=median(z$cdr3_length[pos]), length_negative=median(z$cdr3_length[!pos]))
  }
}
summary <- do.call(rbind, rows)
save_csv(summary, "results/donor_comparison.csv")
png("figures/donor_comparison.png", width=1200, height=720, res=140)
counts <- rbind(summary$links[summary$donor=="Puta"],summary$links[summary$donor=="Christy"])
bars <- barplot(counts, beside=TRUE, col=c("#245D79","#D87B3A"), border=NA,
  names.arg=c("All sequences","95% identity","90% identity","80% identity"),
  ylab="Features linked to selection by MIIC", ylim=c(0,max(counts)*1.22),
  main="Selection links across two alpacas", las=1)
text(bars, counts+0.5, counts)
legend("topright", names(panels), fill=c("#245D79","#D87B3A"), bty="n")
dev.off()
# Check whether exact sequences present in both animals drive the full-data result.
source("scripts/network_functions.R")
common_domains <- intersect(panels$Puta$meta$domain, panels$Christy$meta$domain)
input <- readRDS("second_donor/data/input_unpooled.rds")
keep <- !(panels$Christy$meta$domain %in% common_domains)
edges <- miic_edges(fit_miic(pool_categories(input[keep, ])))
sensitivity <- network_stats(edges, names(input), "Christy_without_shared_domains", "MIIC", NA_real_)
sensitivity$excluded <- sum(!keep)
sensitivity$n <- sum(keep)
save_csv(sensitivity, "results/donor_overlap_sensitivity.csv")
