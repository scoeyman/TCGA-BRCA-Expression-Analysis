# requirements.R

# Install Bioconductor manager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages
bioc_packages <- c("TCGAbiolinks", "DESeq2")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# Install CRAN packages
cran_packages <- c(
  "tidyverse",
  "caret",
  "randomForest",
  "e1071",
  "Rtsne",
  "factoextra",
  "pheatmap",
  "RColorBrewer",
  "xgboost"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}
