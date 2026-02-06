# Project #1: Utilizing multi-omic integration to interrogate molecular etiology of leukemia.

This repository contains a multi-omic integration project to interrogate the molecular etiology of leukemia.

## Project Goal
Determine what insights can be learned about the biology of leukemia from integrating various omics datasets of patients. We also want to 1) classify patients according to their molecular profiles for precision medicine purpose and 2) predict patient mutation status using machine learning. 

---

## Approach
- Integrated 4 omics modalities (somatic mutation profiles, DNA methylation, mRNA expression, and drug response measurements) and trained the MOFA2 model.
- Examined variance decomposition via the percentage of variance explained by each factor across each data modality.
- Interpreted factors using somatic mutation feature weights.
- Predicted missing IGHV mutation status using a Random Forest classifier trained on MOFA factor values.
- Performed Reactome gene set enrichment analysis to define gene sets associated with specific factors.

---

## Results
- Mean feature weights identified IGHV and trisomy12 as key features.
- Random Forest predictions successfully inferred missing IGHV mutation statuses. 
- Pathway enrichment analysis indicated that Factor 5 is associated with stress response gene sets. 

---

## Example Outputs
[Variance Decomposition Analysis](results/figures/03_plot_variance_explained_by_factor.png)

[Predicting somatic IGHV mutation using ML](results/figures/09_plot_factors_1_vs_3_IGHV_predictions.png)

---

## Tools & Packages
R, MOFA2, randomForest, ggplot2, tidyverse.

---

## Repository Structure
`analysis/MOFA_Leukemia.R` analysis script

`results/figures/` contains exported analysis figures

`environment/install_packages.R` installs required packages

---

## Data availability
Data used in this project is house in MOFAdata and downloaded automatically when running the analysis.
