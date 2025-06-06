# # Step 1: Install required packages (run once, comment out after)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("TCGAbiolinks")
# install.packages(c("tidyverse", "caret", "randomForest"))

library(TCGAbiolinks)
library(tidyverse)
library(caret)
library(randomForest)

download_dir <- "E:/genomics/tcga brca data"

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, files.per.chunk = 10, directory = download_dir)
data <- GDCprepare(query, directory = download_dir)

# Save the raw SummarizedExperiment object for later
save(data, file = file.path(download_dir, "tcga_brca_data.RData"))

