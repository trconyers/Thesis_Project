# Genomic Architecture of Polygenic Adaptation in Drosophila melanogaster
### Reproducible Bioinformatic Pipeline for E&R and Transcriptomic Datasets

This repository provides a reproducible, relative-path R bioinformatics pipeline developed for my 2026 Master’s Thesis at California State University, Fullerton (Shahrestani Lab, Rose lineage). The system is engineered to analyze multi-generational selection data from experimentally evolved _Drosophila_ populations to map genomic and transcriptomic overlaps across divergent selection regimes.

🧬 Archival Record: My core thesis findings and experimental evolution datasets are archived at Zenodo: doi.org/10.5281/zenodo.19961092

## Core Capabilities
*   **Automated Gene-List Overlap Tests:** Quantifies and maps selection candidates across longevity, immunity, and stress-tolerance experimental regimes using custom hyper-geometric and matrix overlap scripts.
*   **Functional GO Enrichment Pipelines:** Utilizes automated sweeps to programmatically cluster functional categories, validating a "blended" immune mechanism model that prioritizes generalized stress tolerance over canonical immune pathways.
*   **Cross-Platform Integration:** Synchronizes genomic SNP mapping platforms with downstream transcriptomic expression datasets.

## Technical Philosophy & Portability
To ensure compliance with Open Science data hygiene principles, this pipeline strictly utilizes relative-path directory routing. This eliminates hardcoded path constraints, ensuring immediate, zero-onboarding portability across local workspaces, private servers, or High-Performance Computing (HPC) environments.

## Current Application
This pipeline is currently being optimized to explore the predictability boundaries of polygenic adaptation. It is designed to scale directly into multi-parent mapping panels—such as the Drosophila Synthetic Population Resource (DSPR)—to dissect complex, multiallelic architectures governing host-pathogen interactions.
