# Learning how to carry out multi-omics integration via MOFA2.
# Want to gain deeper insights from omics data in the context of cancer. 
# We are integrating 4 datasets taken from 200 patients: 
# DNA mutations, mRNA expression, DNA methylation, and drug treament response. 
# Goal: To train a MOFA model and then explore it to gain novel insights into CLL cancer. 
# Perhaps we can identify features or trends that define patients and make precision medicine possible. 

# See information on this tutorial here:
#   https://www.youtube.com/watch?v=_BfHeZ0s2i0
#   https://biofam.github.io/MOFA2/tutorials.html
#   https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html
#installed the following packages: randomForest, MOFA2, MOFAdata
###############################################################################

# Portfolio touch: create an output directory for saved figures
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Helper: reliably save a plot by opening a PNG device and drawing into it
# (avoids RStudio device timing issues when running via source())
save_plot_png <- function(filename, expr, width = 2100, height = 1500, res = 300) {
  png(filename, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  
  # Run the plotting code
  p <- eval.parent(substitute(expr))
  
  # If a ggplot-like object is returned, explicitly print it
  if (inherits(p, c("gg", "ggplot"))) print(p)
  
  invisible(p)
}

# First, we load in the libraries and datasets. 
library(MOFA2)        # Core MOFA2 package for multi-omics factor analysis
library(MOFAdata)     # Example datasets used in MOFA tutorials
library(data.table)   # Fast data import and manipulation
library(ggplot2)
library(tidyverse)

# Load in the Chronic Lymphocytic Leukemia (CLL) multi-omics dataset
utils::data("CLL_data")

# Inspect dimensions of each omics modality
# Rows = features, columns = samples
# The modalities are:
# Drugs- drug response features
# Methylation- CpG methylation features
# mRNA- gene expression
# Mutations- mutation presence/absence
lapply(CLL_data, dim)

# Next, we load in metadata of clinical information from the 200 patients
#this includes gender, age, treatment status, survival, etc. 
CLL_metadata <- fread(
  "ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt"
)
###############################################################################

# Now, we want to create an untrained MOFA object with our 4 views
MOFAobject <- create_mofa(CLL_data)
MOFAobject

# Plot the missing values and dimensionality per view
save_plot_png(
  "results/figures/01_plot_data_overview.png",
  plot_data_overview(MOFAobject),
  width = 2400, height = 1200, res = 300
)

# Next we want to define the data arguments for the MOFA object
# Do we want to scale views or groups to the same total variance? Usually, we don't.
data_opts <- get_default_data_options(MOFAobject)
data_opts

# Next we want to choose the model arguments for our object
model_opts <- get_default_model_options(MOFAobject)
# Set number of latent factors to learn
model_opts$num_factors <- 15
model_opts

# Finally, we want to choose the training arguments for our object
train_opts <- get_default_training_options(MOFAobject)
# Slower convergence = more stable solution
train_opts$convergence_mode <- "slow"
# Fix random seed for reproducibility
train_opts$seed <- 42
train_opts
###############################################################################

# Now we want to train our MOFA model! 
# Combine all the options into final MOFA object
#MOFAobject <- prepare_mofa(
#  MOFAobject,
#  data_options = data_opts,
#  model_options = model_opts,
#  training_options = train_opts
#)

# If desired, you can instead load the pretrained model below
MOFAobject <- readRDS(
  url("http://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/MOFA2_CLL.rds")
)
###############################################################################

# Now we want to begin exploring our trained MOFA model. 
# The MOFA object consists of multiple slots where relevant data/info is stored. 
# See documentation for information on all slots. 
slotNames(MOFAobject)

# The data slot contains the input omics data used to train the model
names(MOFAobject@data)

# Example matrix dimensions
dim(MOFAobject@data$Drugs$group1)

# Posterior expectations:
# Z = factor values (samples × factors)
# W = weights (features × factors)
names(MOFAobject@expectations)
dim(MOFAobject@expectations$Z$group1) # we expect this to be 200 x 15 since it's the group1 samples x factors
dim(MOFAobject@expectations$W$mRNA) # we expect this to be 5000 x 15 since it's mRNA features x factors

# Attach clinical metadata to trained model
samples_metadata(MOFAobject) <- CLL_metadata

# Sanity check- verify that factors are uncorrelated
# If there is a lot of correlation between factors, the model is a poor fit
save_plot_png(
  "results/figures/02_plot_factor_cor.png",
  plot_factor_cor(MOFAobject),
  width = 2100, height = 1800, res = 300
)

# The most important insight that MOFA generates is the variance decomposition analysis. 
# This plot shows the percentage of variance explained by each factor across each data modality (and group, if provided). 
# It summarizes the sources of variation from a complex heterogeneous data set in a single figure.
# Variance explained per factor and view
save_plot_png(
  "results/figures/03_plot_variance_explained_by_factor.png",
  plot_variance_explained(MOFAobject, max_r2 = 15),
  width = 2400, height = 1200, res = 300
)
# We can see that drugs, methylation, and mutations explain most variation in Factor 1.
# Then, drugs alone explains the variation in Factor 2. 
# This suggests that drug responses are the most important modality in our omics integration. 

# Next, we determine how much of the variance in each modality the model is able to explain. 
# Total variance explained per modality
save_plot_png(
  "results/figures/04_plot_variance_explained_total.png",
  plot_variance_explained(MOFAobject, plot_total = TRUE)[[2]],
  width = 2400, height = 1200, res = 300
)
###############################################################################

# Let's begin diving further into the model.
# There are a few strategies to characterize the molecular etiology underlying 
# the MOFA Factors and to relate them to the sample covariates:
#     1. Association analysis between the sample metadata and the Factor values.
#     2. Inspection of factor values.
#     3. Inspection of the feature weights.
#     4. Gene set enrichment analysis on the mRNA weights.

# The weights provide a score for each feature on each factor. 
# We saw that mutations capture a lot of the variation in Factor 1, so 
# let's plot the weights for the features in mutations Factor 1
# Plot mutation weights for Factor 1
save_plot_png(
  "results/figures/05_plot_weights_mutations_factor1.png",
  plot_weights(
    MOFAobject,
    view = "Mutations",
    factor = 1,
    nfeatures = 10, # top features
    scale = TRUE
  ),
  width = 2400, height = 1200, res = 300
)

# Cleaner view showing only top weights
save_plot_png(
  "results/figures/06_plot_top_weights_mutations_factor1.png",
  plot_top_weights(
    MOFAobject,
    view = "Mutations",
    factor = 1,
    nfeatures = 10,
    scale = TRUE
  ),
  width = 2400, height = 1200, res = 300
)
# This shows us that a feature called IGHV is important for/influences CLL
# Literature shows that this is actually a clinical marker for CLL!
# Additional analyses to learn more about IGHV can be done, see tutorial for more. 


# Now let's move on to Factor 3. 
# Let's do the same thing and plot the weight of the mutations features: 
save_plot_png(
  "results/figures/07_plot_weights_mutations_factor3.png",
  plot_weights(
    MOFAobject,
    view = "Mutations",
    factor = 3,
    nfeatures = 10,
    abs = FALSE
  ),
  width = 2400, height = 1200, res = 300
)
# This shows us that a feature called trisomy12 is important for/influences CLL
# Literature shows that this is actually a clinical marker for CLL!
# Additional analyses to learn more about IGHV can be done, see tutorial for more.
###############################################################################

#Now that we have characterised the etiology of the two main Factors, let’s explore them together:
# Visualize samples in Factor1 vs Factor3 space using IGHV and trisomy12
p <- plot_factors(
  MOFAobject,
  factors = c(1, 3),
  color_by = "IGHV",
  shape_by = "trisomy12",
  dot_size = 2.5,
  show_missing = TRUE
)
p <- p +
  geom_hline(yintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = "dashed") # Add quadrant guides

save_plot_png(
  "results/figures/08_plot_factors_1_vs_3_IGHV_trisomy12.png",
  print(p),
  width = 2100, height = 1500, res = 300
)
# This plot is extremely important. It classifies the patients into four different subgroups 
# depending on their (multi-omic) molecular profile. As shown in the analysis above, 
# both factors are associated with differences in the drug response assay they are are 
# strongly linked to somatic mutations (IGHV and trisomy12) that are easy to profile in clinical practice. 
# This is fantastic for the aim of personalised medicine.
###############################################################################

# The scatterplot of Factor 1 vs Factor 3 reveals that a few samples are missing the somatic mutation status. 
# In this case, the doctors were not able to classify patients into their clinical subgroups. 
# But we can now use MOFA to exploit the molecular profiles and attempt to impute the IGHV and trisomy12 status.
library(randomForest)

# Note- may need to attach the metadata to the MOFA if you didn't previously: 
CLL_metadata <- data.table::fread(
  "ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt"
)
samples_metadata(MOFAobject) <- CLL_metadata

# Convert MOFA object to dataframe
df <- as.data.frame(
  get_factors(MOFAobject, factors = c(1,2))[[1]]
)

# Train the model for IGHV
df$IGHV <- as.factor(MOFAobject@samples_metadata$IGHV)
model.ighv <- randomForest(
  IGHV ~ ., data = df[!is.na(df$IGHV),], ntree = 10
)
df$IGHV <- NULL

# Predict missing status
MOFAobject@samples_metadata$IGHV.pred <- stats::predict(model.ighv, newdata = df)

# Plot the predictions
MOFAobject@samples_metadata$IGHV.pred_logical <-
  c("True","Predicted")[
    as.numeric(is.na(MOFAobject@samples_metadata$IGHV)) + 1
  ]
p <- plot_factors(
  MOFAobject,
  factors = c(1,3),
  color_by = "IGHV.pred",
  shape_by = "IGHV.pred_logical",
  dot_size = 2.5,
  show_missing = TRUE
)

save_plot_png(
  "results/figures/09_plot_factors_1_vs_3_IGHV_predictions.png",
  print(p),
  width = 2100, height = 1500, res = 300
)
###############################################################################

#In addition to exploring the individual weights for each factor, 
# we can use enrichment analysis to look for signiificant associations of factors to genesets. 
# Here, we use the Reactome genesets.
# Load Reactome pathway gene sets
utils::data(reactomeGS)

# GSEA for positive weights
res.positive <- run_enrichment(
  MOFAobject,
  feature.sets = reactomeGS,
  view = "mRNA",
  sign = "positive"
)

# GSEA for negative weights
res.negative <- run_enrichment(
  MOFAobject,
  feature.sets = reactomeGS,
  view = "mRNA",
  sign = "negative"
)

names(res.positive)

# Plot results
save_plot_png(
  "results/figures/10_plot_enrichment_heatmap_positive.png",
  plot_enrichment_heatmap(res.positive),
  width = 2100, height = 1500, res = 300
)

save_plot_png(
  "results/figures/11_plot_enrichment_heatmap_negative.png",
  plot_enrichment_heatmap(res.negative),
  width = 2100, height = 1500, res = 300
)
# This shows us that Factor 5 has a strong gene set signature. 
# Let's see what gene set is enriched here.

save_plot_png(
  "results/figures/12_plot_enrichment_factor5.png",
  plot_enrichment(res.positive, factor = 5, max.pathways = 15),
  width = 2100, height = 1500, res = 300
)
# It seems that Factor 5 is capturing differences in the stress response of the blood cells.

# Portfolio touch: final completion message
message("Done! All plots have been saved as PNGs in: results/figures/")
