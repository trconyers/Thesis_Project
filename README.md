# Genetic Architecture of Longevity & Immunity in _Drosophila melanogaster_
### Reproducible Bioinformatic Pipeline for E&R and Transcriptomic Datasets

This repository provides a reproducible, relative-path R bioinformatics pipeline developed for my 2026 Master’s Thesis at California State University, Fullerton (Shahrestani Lab, Rose lineage). The system is engineered to analyze multi-generational selection data from experimentally evolved _Drosophila_ populations to map genomic and transcriptomic overlaps across divergent selection regimes.

🧬 Archival Record: My core thesis findings and experimental evolution datasets are archived at Zenodo: doi.org/10.5281/zenodo.19961092

## Core Capabilities
*   **Automated Gene-List Overlap Tests:** Quantifies and maps selection candidates across longevity, immunity, and stress-tolerance experimental regimes using custom hyper-geometric overlap scripts.
*   **Functional GO Enrichment Pipelines:** Clusters functional categories, performing over-representation analyses (ORA) and gene set enrichment analyses (GSEA).
*   **Differentiation of immune defenses:** Performs overlap analyses with "tolerance" and "resistance" genes (Howick & Lazzaro, 2017).

## Technical Philosophy & Portability
To ensure compliance with Open Science data hygiene principles, this pipeline strictly utilizes relative-path directory routing.

## Current Application
This pipeline will be optimized in a different repository to explore the predictability boundaries of polygenic adaptation. It is designed to scale directly into multi-parent mapping panels—such as the Drosophila Synthetic Population Resource (DSPR)—to dissect complex, multiallelic architectures governing host-pathogen interactions.
