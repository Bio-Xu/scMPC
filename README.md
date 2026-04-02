# scMPC: Identification of Metastatic Potential Cells from Single-Cell Transcriptomic Data
**scMPC** is a computational framework designed to identify *metastatic potential cells (MPCs)* from single-cell RNA-seq data by integrating bulk transcriptomic information and ensemble machine learning strategies.

Tumor metastasis is a complex and heterogeneous process, and only a subset of tumor cells possess the ability to initiate and sustain metastatic progression. scMPC addresses this challenge by leveraging metastasis-associated gene signatures derived from bulk transcriptomic datasets and transferring this information to the single-cell level.

## Key Features
- **Transfer learning framework**
  Integrates bulk transcriptomic data with single-cell RNA-seq data to identify metastasis-associated cellular subpopulations.
- **Monte Carlo–based ensemble learning**
  Utilizes large-scale random gene set sampling and iterative screening to robustly identify metastasis-relevant gene signatures.
- **Consensus clustering strategy**
  Combines multiple clustering results to define MPCs as stable and reproducible cell populations.
- **Robustness and stability**
  Incorporates multi-metric evaluation, resampling strategies, and sensitivity analyses to ensure reliable identification of MPCs across patients.
## Input
- Gene expression count matrix (single-cell RNA-seq data)
## Output
- MPC annotations at single-cell resolution
- Metastasis-associated gene signatures (M_signatures)
- MPC activity scores for downstream analysis
## Applications
- Identification of metastasis-associated tumor cell subpopulations
- Analysis of tumor heterogeneity and evolutionary dynamics
- Integration with downstream analyses (e.g., trajectory inference, regulatory networks, tumor microenvironment interactions)

<img width="1845" height="347" alt="image" src="https://github.com/user-attachments/assets/0db67065-6cc6-46c1-92b7-e0a56682a77c" />

## Article
Integrative Multi-Modal Transcriptomic Identification of Metastatic Potential Cells Reveals Mechanistic Insights and Pro-Metastatic Ecosystems in Breast Cancer
