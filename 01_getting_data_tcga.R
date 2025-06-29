# This script downloads TCGA BRCA gene expression data using TCGAbiolinks and saves the raw data object for later use.

# ---- Install required packages (run once, then comment out) ----
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks")
# install.packages(c("tidyverse", "caret", "randomForest"))

library(TCGAbiolinks)
library(tidyverse)
library(caret)
library(randomForest)

# Manually set project directory here:
project_dir <- "E:/genomics/tcga brca data"                     # change this to local path

# Set data directory inside project folder
data_dir <- file.path(project_dir, "data")

# Create data directory if it doesn't exist
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# Query TCGA BRCA data
query <- GDCquery(
    project = "TCGA-BRCA", 
    data.category = "Transcriptome Profiling", 
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts")

# Download data (files will be saved in data_dir)
GDCdownload(query, files.per.chunk = 10, directory = data_dir)

# Prepare the data into a SummarizedExperiment object
data <- GDCprepare(query, directory = data_dir)

# Save the raw data object for future steps
save(data, file = file.path(data_dir, "tcga_brca_data.RData"))

cat("Data download and save complete.\n")