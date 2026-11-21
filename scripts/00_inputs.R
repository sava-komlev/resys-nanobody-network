av <- read.csv("data/AVIDa-SARS-CoV-2.csv", stringsAsFactors = FALSE)
cat("AVIDa total rows:", nrow(av), "; unique supplied sequences:", length(unique(av$VHH_sequence)), "\n")
donor <- getOption("resys.donor", "puta")
sub <- subset(av, subject_name == donor & Ag_label == "WT_cell")
stopifnot(nrow(sub) == length(unique(sub$VHH_sequence)))
cat(donor, "/ WT_cell unique supplied sequences:", nrow(sub), "\n")
print(table(sub$label))
sub$variant_id <- sprintf("v%05d", seq_len(nrow(sub)) - 1L)
if (!file.exists("data/avida_imgt_H.csv")) {
  fasta <- sprintf("data/avida_%s_WT_cell.fasta", donor)
  writeLines(as.vector(rbind(paste0(">", sub$variant_id), sub$VHH_sequence)),
             fasta)
  stopifnot(system2("anarci", c("-i", fasta, "-o",
    "data/avida_imgt", "-s", "imgt", "--csv", "-r", "H", "--ncpu", "2")) == 0)
}
a <- read.csv("data/avida_imgt_H.csv", check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(a$Id), setequal(a$Id, sub$variant_id))
positions <- grep("^[0-9]+[A-Z]*$", names(a), value = TRUE)
a$domain <- gsub("-", "", apply(a[positions], 1, paste0, collapse = ""), fixed = TRUE)
a$label <- sub$label[match(a$Id, sub$variant_id)]
cat("Distinct numbered domains before QC:", length(unique(a$domain)), "\n")
conflict <- vapply(split(a$label, a$domain), function(z) length(unique(z)) > 1L, logical(1))
stopifnot(!any(conflict))
complete <- nchar(a$domain) >= 100 & nchar(a$domain) <= 160 &
  rowSums(a[c("23", "41", "104", "118")] == "-") == 0
cat("Excluded incomplete/short domains:", sum(!complete), "\n")
print(a[!complete, c("Id", "seqstart_index", "seqend_index")])
a <- a[complete, ]
cat("Duplicate domain records removed:", sum(duplicated(a$domain)), "\n")
a <- a[!duplicated(a$domain), ]
cat("Final distinct domains:", nrow(a), "\n")
print(table(a$label))
cdr_positions <- positions[as.integer(sub("[A-Z]+$", "", positions)) %in% c(27:38, 56:65, 105:117)]
qc <- do.call(rbind, lapply(cdr_positions, function(p) {
  v <- a[[p]]; f <- prop.table(table(v[v != "-"]))
  data.frame(position = p, coverage = mean(v != "-"), entropy = -sum(f * log2(f)),
             major_fraction = max(f), states_ge_1pct = sum(f >= .01))
}))
eligible <- subset(qc, coverage >= .95 & major_fraction <= .95 & states_ge_1pct >= 2)
eligible <- eligible[order(-eligible$entropy, eligible$position), ]
if (is.null(getOption("resys.fixed_positions"))) {
  stopifnot(nrow(eligible) >= 24)
  selected <- eligible$position[1:24]
} else {
  selected <- getOption("resys.fixed_positions")
  stopifnot(length(selected) == 24L, all(selected %in% cdr_positions))
}
cat("Eligible and selected CDR positions:", nrow(eligible), length(selected), "\n")
print(selected)
qc$selected <- qc$position %in% selected
write.csv(qc, "data/cdr_position_qc.csv", row.names = FALSE)
cdr3_positions <- positions[as.integer(sub("[A-Z]+$", "", positions)) %in% 105:117]
cdr3 <- gsub("-", "", apply(a[cdr3_positions], 1, paste0, collapse = ""), fixed = TRUE)
cat("Unique exact CDR3 strings (not independent-clone count):", length(unique(cdr3)), "\n")
net <- a[selected]
names(net) <- paste0("IMGT_", selected)
# A gap is an observed absence at an alignment coordinate, not an unknown residue.
# Rare categories can be pooled for inference, with a sensitivity analysis.
net$cdr3_length <- nchar(cdr3)
net$selection_label <- a$label
stopifnot(ncol(net) == 26L, !anyNA(net), nrow(net) >= 1000)
write.csv(net, "data/avida_network_input.csv", row.names = FALSE)
write.csv(data.frame(variant_id = a$Id, domain = a$domain, cdr3 = cdr3, label = a$label),
          "data/avida_variant_metadata.csv", row.names = FALSE)
cat("Network input:", nrow(net), "observations x", ncol(net), "variables\n\n")

