# Impact of sequence redundancy on network inference in nanobody screening data

Which nanobody sequence features are associated with a positive result in a spike-binding screen? Do the same associations persist when each group of similar nanobodies contributes one sequence?

This R project uses the wild-type spike cell assay from AVIDa-SARS-CoV-2. After filtering and deduplication, the analysis included 6,039 distinct VHH sequences from Puta (1,311 positive) and 4,726 from Christy (937 positive). Each network contains 26 variables: the amino acid at each of 24 positions across CDR1, CDR2 and CDR3, the length of CDR3, and the experimental selection label. Three methods were compared: pairwise mutual information, ARACNE pruning and undirected MIIC.

## Findings

- In the full datasets, MIIC linked 17 features to selection in Puta and three in Christy. At 90% CDR3 identity, none of those original links survived when each group contributed one representative. Christy gained a link to CDR3 length. [Puta network](figures/miic_network_comparison.png) and [comparison across alpacas](figures/donor_comparison.png).
- In Puta, random samples matched to the CDR3 representative sample in both size and positive count retained 5–9 selection links across 30 draws. All 30 representative draws gave zero. [Matched controls](results/matched_controls.csv).
- Grouping by similarity across the entire VHH sequence gave the same loss of original links at 90% identity in both donors. At 95%, however, Christy’s reference CDR3 network contained nine selection links, while all 20 draws based on the entire VHH sequence yielded only the length link. [Comparison across grouping thresholds](figures/domain_group_validation.png).
- Prediction was harder when entire CDR3 groups were held out. Mean ROC AUC fell from 0.890 to 0.702 in Puta and from 0.949 to 0.862 in Christy when random splits were replaced by 90% CDR3 group splits. [Puta prediction results](results/prediction_metrics.csv) and [Christy prediction results](second_donor/results/prediction_metrics.csv).

## Reproduce the analysis

The original analysis used R 4.2.1, Biostrings 2.64.1, igraph 1.6.0 and miic 2.0.4. The [package manifest](environment/package_versions.csv) records the tested versions. The input step also requires the ANARCI command-line program with IMGT numbering and HMMER on the executable search path.

From the project root, create a private R library. Install compatible versions of Biostrings (Bioconductor), igraph, Rcpp, scales, rpart, MASS and jsonlite with their dependencies in your R installation first. The setup script copies these dependencies into the project library. Then install ppcor and miic into that library. For example, with the indicated R and Bioconductor releases:

```sh
Rscript --vanilla scripts/setup_environment.R
export R_LIBS_USER="$PWD/environment/library"
Rscript --vanilla -e 'install.packages("ppcor", repos="https://cloud.r-project.org")'
Rscript --vanilla -e 'remotes::install_version("miic", version="2.0.4", lib=Sys.getenv("R_LIBS_USER"), upgrade="never", repos="https://cloud.r-project.org")'
```

The versioned installation command requires the remotes package. Download the source CSV and verify its SHA-256 digest before running the analysis:

```sh
mkdir -p data
curl -fL https://huggingface.co/datasets/COGNANO/AVIDa-SARS-CoV-2/resolve/main/AVIDa-SARS-CoV-2.csv -o data/AVIDa-SARS-CoV-2.csv
echo '02111685d5f02357021637a1232630249ff3e44b7c9f5159e471ccf9696920e6  data/AVIDa-SARS-CoV-2.csv' | shasum -a 256 -c -
```

Run the scripts in order from the project root. On the first run, step 00 numbers the sequences with ANARCI. Subsequent runs reuse the numbered alignment in data/avida_imgt_H.csv.

```sh
set -e
sh run.sh scripts/test_methods.R
for script in scripts/[0-9][0-9]_*.R; do sh run.sh "$script"; done
```

Steps 00–02 prepare the data, group sequences and reconstruct networks. Steps 03–06 run resampling, matched controls, prediction and stability checks. Step 07 builds the figures and step 08 validates the outputs. Step 09 compares amino-acid and loop-length frequencies between positive and negative sequences at features linked to selection. Scripts use fixed random seeds and overwrite generated outputs on rerun.

## Data and methods

The source data are available under CC BY-NC 4.0. The outcome records selection in a specific assay. The comparison uses the same target and assay in both donors.

- [AVIDa-SARS-CoV-2 data](https://huggingface.co/datasets/COGNANO/AVIDa-SARS-CoV-2) and [dataset paper](https://doi.org/10.52202/079017-3688)
- [ANARCI sequence numbering](https://doi.org/10.1093/bioinformatics/btv552)
- [ARACNE](https://doi.org/10.1186/1471-2105-7-S1-S7) and [MIIC](https://doi.org/10.1371/journal.pcbi.1007866)

## Second alpaca

After the primary analysis, run `sh run.sh second_donor.R` and then `sh run.sh compare_donors.R`. Christy uses the same WT cell assay and the 24 positions chosen in Puta. The script numbers Christy's sequences separately and repeats the network analysis, matched controls and within-donor prediction checks. For Christy's classifier the panel is fixed in advance. Puta's original classifier selected positions within each training fold. The comparison also checks Christy's full MIIC network after excluding exact domains found in Puta.

After both donor prediction checks, run `sh run.sh test_position_pairs.R` to compare separate-position and joint-position models on the saved held-out CDR3 groups. [Results by donor, split and position pair](results/position_pair_validation.csv).

After preparing both donors, run `sh run.sh test_domain_groups.R` to compare representatives grouped by similarity across the entire VHH sequence with random samples matched for size and positive count. [Network results](results/domain_group_validation.csv) and [group membership overlap](results/domain_group_agreement.csv).
