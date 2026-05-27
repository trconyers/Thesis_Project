# Master Driver Script: Genetic Architecture of Longevity & Immunity
# Author: Taryn Conyers-Casey
# Purpose: Execute full analysis pipeline from raw data to immunity comparisons

# Initialize environment
source(".Rprofile")

# 0. Load global functions and setup relative paths
source("Scripts/00_misc_functions.R")

# 1. Prepare and clean data
source("Scripts/01_data_preparation.R")

# 2. Map SNPs to gene regions
source("Scripts/02_SNPs2genes.R")

# 3. Perform Overlap Analyses
source("Scripts/03a_GeneListAnalyses_Combined.R")
source("Scripts/03b_GeneListAnalyses_Genomics.R")
source("Scripts/03c_GeneListAnalyses_Transcriptomics.R")

# 4. Functional Enrichment
source("Scripts/04a_go_enrichment_Genomics.R")
source("Scripts/04b_go_enrichment_GSEA.R")

# 5. Theoretical Synthesis: Resistance vs. Tolerance
source("Scripts/05_Tolerance_Resistance.R")

print("Pipeline execution complete. Plots saved to /Results/Graphs and tables saved to /Results/Tables.")
