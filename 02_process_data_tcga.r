# This script loads raw TCGA BRCA data, merges subtype labels,
# processes expression data, and saves a processed data frame.

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

# Set project directory (update to local path)
project_dir <- "E:/genomics/tcga brca data"
data_dir <- file.path(project_dir, "data")

# Load the raw data object saved in previous step
load(file.path(data_dir, "tcga_brca_data.RData"))           # loads 'data'

# Get PAM50 subtype labels for BRCA patients
subtypes <- TCGAquery_subtype("BRCA") %>%
  select(patient, BRCA_Subtype_PAM50) %>%
  drop_na()

# Rename subtype column for easier reference
colnames(subtypes)[2] <- "Subtype"

# Extract expression matrix and sample patient IDs
expr <- assay(data)
samples <- colnames(expr)
patients <- substr(samples, 1, 12)                          # first 12 chars = patient barcode

# Create metadata table with patient and sample IDs
meta <- data.frame(patient = patients, sample = samples)

# Merge subtype info with sample metadata
meta <- merge(meta, subtypes, by = "patient")

# Filter expression matrix to keep samples with subtype info
expr_filtered <- expr[, meta$sample]

# Transpose expression matrix and convert to data frame
expr_df <- t(expr_filtered) %>% as.data.frame()

# Add subtype labels to the expression data frame
expr_df$Subtype <- meta$Subtype

# Print subtype distribution and preview first 6 genes
print(table(expr_df$Subtype))
print(head(expr_df)[, 1:6])

# Save processed expression data frame for downstream analysis
save(expr_df, file = file.path(data_dir, "tcga_brca_expr_df.RData"))

cat("Expression data processing complete.\n")
