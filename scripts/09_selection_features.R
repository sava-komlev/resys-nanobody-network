# Describe the observed states at features linked to selection in the full MIIC graph.
source("scripts/common.R")
input <- read.csv("data/avida_network_input.csv", check.names = FALSE)
metadata <- read.csv("data/metadata_clustered.csv")
representatives <- readRDS("data/representatives.rds")$cluster90
edges <- readRDS("results/network_edges.rds")$full_MIIC
target <- "selection_label"
features <- union(edges$x[edges$y == target], edges$y[edges$x == target])
features <- union(features, getOption("resys.followup_features", character()))
stopifnot(length(features) > 0, identical(input[[target]], metadata$label),
          !anyDuplicated(metadata$cluster90[representatives]))
samples <- list(full = seq_len(nrow(input)), cluster90 = representatives)
tables <- list()
for (sample_name in names(samples)) {
  rows <- samples[[sample_name]]
  positive <- rows[input[[target]][rows] == 1]
  negative <- rows[input[[target]][rows] == 0]
  for (feature in features) {
    values <- input[[feature]]
    states <- sort(unique(values))
    count_states <- function(indices) tabulate(match(values[indices], states), length(states))
    pos <- count_states(positive)
    neg <- count_states(negative)
    stopifnot(sum(pos) == length(positive), sum(neg) == length(negative))
    tables[[length(tables) + 1L]] <- data.frame(
      dataset = sample_name, feature = feature, state = states,
      positive_count = pos, positive_total = length(positive),
      negative_count = neg, negative_total = length(negative),
      pct_of_positives = 100 * pos / length(positive),
      pct_of_negatives = 100 * neg / length(negative))
  }
  cat(sample_name, "median length in positives and negatives:",
      median(input$cdr3_length[positive]), median(input$cdr3_length[negative]), "\n")
}
save_csv(do.call(rbind, tables), "results/selection_feature_frequencies.csv")
