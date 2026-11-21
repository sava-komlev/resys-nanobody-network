# Run via sh run.sh second_donor.R after running the primary analysis.
library_root <- normalizePath("environment/library")
project <- normalizePath(".")
primary_input <- read.csv("data/avida_network_input.csv", nrows = 1)
fixed_positions <- sub("^IMGT_", "", grep("^IMGT_", names(primary_input), value = TRUE))
primary_edges <- readRDS("results/network_edges.rds")$full_MIIC
followup <- with(primary_edges, union(x[y == "selection_label"], y[x == "selection_label"]))
dir.create("second_donor", showWarnings = FALSE)
setwd("second_donor")
for (dir in c("data", "environment", "results", "logs")) dir.create(dir, showWarnings = FALSE)
dependencies <- c("scripts" = file.path(project, "scripts"),
                  "environment/library" = library_root,
                  "data/AVIDa-SARS-CoV-2.csv" = file.path(project, "data/AVIDa-SARS-CoV-2.csv"))
for (path in names(dependencies)) {
  if (!file.exists(path)) stopifnot(file.symlink(dependencies[[path]], path))
  stopifnot(normalizePath(path) == dependencies[[path]])
}
options(resys.donor = "christy", resys.fixed_positions = fixed_positions,
        resys.followup_features = followup)
source("scripts/00_inputs.R")
options(resys.domains = nrow(net))
for (script in c("01_prepare.R", "02_networks.R", "04_controls.R", "09_selection_features.R", "05_prediction_check.R"))
  source(file.path("scripts", script))
