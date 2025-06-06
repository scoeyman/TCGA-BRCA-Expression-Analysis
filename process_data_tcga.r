library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

# Load saved data object
load("E:/genomics/tcga brca data/tcga_brca_data.RData")  # loads 'data'

# Get subtype labels
subtypes <- TCGAquery_subtype("BRCA") %>%
  select(patient, BRCA_Subtype_PAM50) %>%
  drop_na()

# Rename column for consistency
colnames(subtypes)[2] <- "Subtype"

# Extract expression matrix and patient IDs
expr <- assay(data)
samples <- colnames(expr)
patients <- substr(samples, 1, 12)

# Merge subtype info with sample metadata
meta <- data.frame(patient = patients, sample = samples)
meta <- merge(meta, subtypes, by = "patient")

# Filter expression matrix and create data frame with subtype labels
expr_filtered <- expr[, meta$sample]
expr_df <- t(expr_filtered) %>% as.data.frame()
expr_df$Subtype <- meta$Subtype

# Quick sanity checks
print(table(expr_df$Subtype))
print(head(expr_df)[, 1:6])  # Show first few gene columns + subtype

# Save processed dataframe for analysis
save(expr_df, file = "E:/genomics/tcga brca data/tcga_brca_expr_df.RData")
