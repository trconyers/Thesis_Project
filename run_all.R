# Master Driver Script: Genetic Architecture of Longevity & Immunity
# Author: Taryn Conyers-Casey
# Purpose: Execute full analysis pipeline from raw data to immunity comparisons

# 0. Load global functions and setup relative paths
source("00_misc_functions.R")

# 1. Prepare and clean data
source("01_data_preparation.R")

# 2. Map SNPs to gene regions
source("02_SNPs2genes.R")

# 3. Perform Overlap Analyses
source("03a_GeneListAnalyses_Genomics.R")
source("03b_GeneListAnalyses_Transcriptomics.R")
source("03c_GeneListAnalyses_Combined.R")

# 4. Functional Enrichment
source("04a_go_enrichment_Genomics.R")
source("04b_go_enrichment_GSEA.R")

# 5. Theoretical Synthesis: Resistance vs. Tolerance
source("05_Tolerance_Resistance.R")

print("Pipeline execution complete. Plots saved to /results/Graphs and tables saved to /results/Tables.")
