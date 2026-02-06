# Install required packages for this repo

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor packages
BiocManager::install(c("MOFA2", "MOFAdata"), ask = FALSE, update = FALSE)

# CRAN packages
install.packages(c("data.table", "ggplot2", "tidyverse", "randomForest"))

message("Packages installed. Next run: source('analysis/run_mofa2_cll.R')")
