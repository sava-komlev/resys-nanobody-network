# Run from the project root with sh run.sh test_domain_groups.R
source("scripts/common.R")
source("scripts/network_functions.R")
assert_private(c("Biostrings", "igraph", "miic"))

# Use the original equal-length Hamming / single-linkage definition, now on
# the extracted VHH domain. Identity thresholds describe different regions
# and therefore need not yield comparable numbers of groups.
domain_groups <- function(domains, cutoffs) {
  lengths <- nchar(domains)
  membership <- matrix(NA_integer_, length(domains), length(cutoffs))
  for (rows in split(seq_along(domains), lengths)) {
    if (length(rows) == 1L) {
      membership[rows, ] <- rows
      next
    }
    distances <- as.matrix(Biostrings::stringDist(
      Biostrings::AAStringSet(domains[rows]), method = "hamming"))
    for (j in seq_along(cutoffs)) {
      max_mismatches <- (lengths[rows[1]] * (100L - cutoffs[j])) %/% 100L
      graph <- igraph::graph_from_adjacency_matrix(distances <= max_mismatches,
                                                    mode = "undirected", diag = FALSE)
      component <- igraph::components(graph)$membership
      membership[rows, j] <- rows[match(component, component)]
    }
  }
  stopifnot(!anyNA(membership))
  membership
}

# Compare co-grouped pairs without enumerating the sequence pairs.
partition_overlap <- function(cdr, domain) {
  co_pairs <- function(z) sum(choose(as.numeric(table(z)), 2))
  n_cdr <- co_pairs(cdr)
  n_domain <- co_pairs(domain)
  shared <- co_pairs(paste(cdr, domain, sep = ":"))
  data.frame(cdr_groups = length(unique(cdr)), domain_groups = length(unique(domain)),
             cdr_pairs = n_cdr, domain_pairs = n_domain, shared_pairs = shared,
             pair_jaccard = shared / (n_cdr + n_domain - shared))
}

outcome_features <- function(edges) {
  sort(union(edges$x[edges$y == "selection_label"], edges$y[edges$x == "selection_label"]))
}
cutoffs <- c(90L, 95L, 98L)
runs <- list()
agreement <- list()
for (donor in c("Puta", "Christy")) {
  folder <- if (donor == "Puta") "." else "second_donor"
  metadata <- read.csv(file.path(folder, "data/metadata_clustered.csv"))
  x <- readRDS(file.path(folder, "data/input_unpooled.rds"))
  full <- outcome_features(readRDS(file.path(folder, "results/network_edges.rds"))$full_MIIC)
  stopifnot(!anyDuplicated(metadata$domain), !anyNA(metadata$domain),
            identical(as.character(metadata$label), x$selection_label))
  groups <- domain_groups(metadata$domain, cutoffs)
  for (identity in intersect(cutoffs, c(90L, 95L))) {
    comparison <- partition_overlap(metadata[[paste0("cluster", identity)]],
                                    groups[, match(identity, cutoffs)])
    agreement[[length(agreement) + 1L]] <- data.frame(donor, identity, comparison)
  }
  strata <- split(seq_len(nrow(x)), metadata$label)
  for (j in seq_along(cutoffs)) {
    families <- split(seq_len(nrow(x)), groups[, j])
    n <- length(families)
    if (j > 1L) stopifnot(all(vapply(families, function(ix)
      length(unique(groups[ix, j-1L])) == 1L, logical(1))))
    stopifnot(n >= 1000L)
    for (draw in seq_len(20L)) {
      set.seed(271000L + 10000L * match(donor, c("Puta", "Christy")) + 100L * j + draw)
      representatives <- vapply(families, function(ix) ix[sample.int(length(ix), 1L)], integer(1))
      positives <- sum(metadata$label[representatives])
      matched <- c(sample(strata[["1"]], positives), sample(strata[["0"]], n - positives))
      stopifnot(!anyDuplicated(groups[representatives, j]), !anyDuplicated(matched),
                sum(metadata$label[matched]) == positives)
      samples <- list(representatives = representatives, matched = matched)
      for (scheme in names(samples)) {
        edges <- miic_edges(fit_miic(pool_categories(x[samples[[scheme]], , drop = FALSE])))
        features <- outcome_features(edges)
        runs[[length(runs) + 1L]] <- data.frame(
          donor, identity = cutoffs[j], draw, scheme, n, positives,
          label_edges = length(features), shared_with_full = length(intersect(features, full)),
          features = paste(features, collapse = ";"))
      }
    }
    message(donor, " ", cutoffs[j], "%: ", n, " groups")
  }
}
runs <- do.call(rbind, runs)
stopifnot(nrow(runs) == 240L)
save_csv(runs, "results/domain_group_validation.csv")
save_csv(do.call(rbind, agreement), "results/domain_group_agreement.csv")

# Plot retained full-data neighbours. A new link at a different position does
# not count as a surviving original association.
runs <- read.csv("results/domain_group_validation.csv")
png("figures/domain_group_validation.png", width = 1750, height = 760, res = 150)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.9, 0.6), oma = c(3.4, 0.5, 1.7, 0.8))
colours <- c(representatives = "#146966", matched = "#A45331")
for (donor_name in c("Puta", "Christy")) {
  rows <- subset(runs, donor == donor_name)
  full_degree <- if (donor_name == "Puta") 17L else 3L
  plot(NA_real_, NA_real_, xlim = c(0.6, 3.4), ylim = c(-0.3, full_degree + 0.4),
       xlab = "Identity across the entire VHH sequence", ylab = "Original selection links retained",
       axes = FALSE, main = paste0(donor_name, " (", full_degree, " original links)"))
  ticks <- if (donor_name == "Puta") c(0, 4, 8, 12, 17) else 0:3
  abline(h = ticks, col = "#E2E8E7")
  axis(1, at = 1:3, labels = c("90%", "95%", "98%"), tick = FALSE)
  axis(2, at = ticks, las = 1)
  for (scheme_name in names(colours)) {
    sub <- rows[rows$scheme == scheme_name, ]
    values <- split(sub$shared_with_full, factor(sub$identity, levels = c(90, 95, 98)))
    xs <- 1:3 + if (scheme_name == "representatives") -0.10 else 0.10
    lo <- vapply(values, min, numeric(1))
    hi <- vapply(values, max, numeric(1))
    segments(xs, lo, xs, hi, col = colours[scheme_name], lwd = 2)
    segments(xs - 0.035, lo, xs + 0.035, lo, col = colours[scheme_name])
    segments(xs - 0.035, hi, xs + 0.035, hi, col = colours[scheme_name])
    points(xs, vapply(values, median, numeric(1)), pch = 19, cex = 1.4, col = colours[scheme_name])
  }
}
mtext("Points show medians across 20 draws; bars span the observed range.", side = 1,
      outer = TRUE, line = 2.2, cex = 0.85)
par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom", inset = c(0, 0.047), legend = c("One per group of similar VHH sequences", "Random rows matched for size and positives"),
       pch = 19, col = colours, horiz = TRUE, bty = "n", cex = 0.88)
invisible(dev.off())
