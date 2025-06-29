TCGA BRCA Gene Expression Analysis
----------------------------------
This repository contains a complete R-based pipeline for downloading, preprocessing, analyzing, and modeling TCGA BRCA (Breast Invasive Carcinoma) gene expression data. The pipeline uses public resources via the TCGAbiolinks package and includes dimensionality reduction, visualization, and machine learning classification.


Project Structure
-----------------
Organized into three main R scripts:
1. `01_getting_data_tcga.R`
   - Downloads TCGA BRCA RNA-Seq gene expression data using TCGAbiolinks.
   - Saves the raw data in an RData file for downstream analysis.

2. `02_process_data_tcga.R`
   - Loads the raw TCGA expression data.
   - Merges PAM50 subtype labels.
   - Filters and reshapes the data into a tidy data frame format.
   - Saves the processed expression data.

3. `03_analyze_data_tcga.R`
   - Loads processed expression data.
   - Normalizes with DESeq2 (VST).
   - Runs PCA and t-SNE for dimensionality reduction.
   - Applies K-means clustering.
   - Trains and evaluates SVM, Random Forest, and XGBoost models.
   - Saves visualizations and model metrics.


Prerequisites
-------------
Before running the scripts, install the required R packages by running: source("requirements.R")


Usage
---------------
Follow these steps in order. Set your working directory to the local path where this project is saved.
1. Run 01_getting_data_tcga to download TCGA BRCA gene expression data. This will save the raw data as tcga_brca_data.RData.
2. Run 2_preprocess_expression_data.R to extract the gene expression matrix and attach subtype labels. This will save tcga_brca_expr_df.RData.
3. Run 3_analysis_modeling.R to:
- Normalize expression data using DESeq2
- Run PCA and t-SNE for visualization
- Perform K-means clustering and compare to PAM50 subtypes
- Train classification models (SVM, Random Forest, XGBoost)
- Output plots and accuracy comparisons

Outputs include:
- PCA and t-SNE plots
- Clustering heatmap
- Model accuracy CSV
- Variable importance plots


Directory Structure
--------------------
project_folder/
│
├── data/
│   ├── tcga_brca_data.RData
│   └── tcga_brca_expr_df.RData
│
├── 01_getting_data_tcga.R
├── 2_preprocess_expression_data.R
├── 3_analysis_modeling.R
├── README.txt
│
├── pca_plot.png
├── tsne_plot.png
├── tsne_kmeans_clusters.png
├── tsne_cluster_vs_subtype_heatmap.png
├── svm_tuning_plot.png
├── svm_variable_importance.png
├── rf_variable_importance.png
├── xgb_variable_importance.png
├── model_accuracy_comparison.csv


Notes
---------------
Adjust the project_dir paths at the top of each script to match your local setup.
The scripts assume raw count data and are intended for exploratory subtype analysis and supervised classification.
Each script prints time stamps to monitor progress for long-running operations.
