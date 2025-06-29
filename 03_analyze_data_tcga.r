# This script loads processed TCGA BRCA expression data,
# performs normalization, dimensionality reduction, clustering,
# and trains classification models (SVM, Random Forest, XGBoost),
# saving plots and results for downstream analysis.

# Load required libraries 
cat("Loading libraries - started at", Sys.time(), "\n")
library(tidyverse)
library(caret)
library(randomForest)
library(DESeq2)
library(Rtsne)
library(e1071)                                                                    # for SVM
library(factoextra)                                                               # for clustering visualization
library(pheatmap)
library(RColorBrewer)
cat("Loading libraries - finished at", Sys.time(), "\n\n")

# Set project directory (update to local path)
project_dir <- "E:/genomics/tcga brca data"
data_dir <- file.path(project_dir, "data")

# Load processed expression data
cat("Loading processed expression data - started at", Sys.time(), "\n")
load(file.path(data_dir, "tcga_brca_expr_df.RData"))  # loads expr_df
cat("Loading processed expression data - finished at", Sys.time(), "\n\n")

# Prepare count matrix (genes x samples)
cat("Preparing count matrix - started at", Sys.time(), "\n")
count_matrix <- as.matrix(t(expr_df %>% select(-Subtype)))
colnames(count_matrix) <- paste0("Sample", seq_len(ncol(count_matrix)))  # generic sample names
rownames(count_matrix) <- colnames(expr_df)[-ncol(expr_df)]             # gene IDs
cat("Preparing count matrix - finished at", Sys.time(), "\n\n")

# Select top 2000 most variable genes
cat("Selecting top 2000 variable genes - started at", Sys.time(), "\n")
gene_vars <- apply(count_matrix, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:2000]
count_matrix <- count_matrix[top_genes, ]
cat("Selecting top 2000 variable genes - finished at", Sys.time(), "\n\n")

# Create sample metadata
cat("Creating sample metadata - started at", Sys.time(), "\n")
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  Subtype = factor(expr_df$Subtype))
cat("Creating sample metadata - finished at", Sys.time(), "\n\n")

# Summary info
cat("Count matrix dimensions:", dim(count_matrix), "\n")
cat("Subtype distribution:\n")
print(table(expr_df$Subtype))
cat("\n")

# Create DESeq2 dataset for normalization 
cat("Creating DESeqDataSet - started at", Sys.time(), "\n")
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Subtype)
cat("Creating DESeqDataSet - finished at", Sys.time(), "\n\n")

# Variance Stabilizing Transformation (VST)
cat("Performing VST - started at", Sys.time(), "\n")
vsd <- vst(dds, blind = TRUE)
cat("Performing VST - finished at", Sys.time(), "\n\n")

# Prepare VST dataframe 
vsd_matrix <- t(assay(vsd))
vsd_df <- as.data.frame(vsd_matrix)
vsd_df$Subtype <- sample_info$Subtype

# PCA
cat("Performing PCA - started at", Sys.time(), "\n")
pca_res <- prcomp(vsd_matrix, scale. = TRUE)
cat("Performing PCA - finished at", Sys.time(), "\n\n")

pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$Subtype <- vsd_df$Subtype

pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Subtype)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA of TCGA BRCA expression data", x = "PC1", y = "PC2") +
  theme_minimal()

ggsave(filename = file.path(project_dir, "pca_plot.png"), plot = pca_plot, width = 7, height = 5)
cat("PCA plot saved.\n\n")

# t-SNE on top 50 PCs
cat("Performing t-SNE - started at", Sys.time(), "\n")
set.seed(123)
tsne_out <- Rtsne(pca_res$x[, 1:50], perplexity = 30, verbose = TRUE, max_iter = 1000)
cat("Performing t-SNE - finished at", Sys.time(), "\n\n")

tsne_df <- as.data.frame(tsne_out$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df$Subtype <- vsd_df$Subtype

tsne_plot <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = Subtype)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "t-SNE of TCGA BRCA expression data", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

ggsave(filename = file.path(project_dir, "tsne_plot.png"), plot = tsne_plot, width = 7, height = 5)
cat("t-SNE plot saved.\n\n")

# Clustering on t-SNE coordinates
cat("Performing k-means clustering on t-SNE - started at", Sys.time(), "\n")
set.seed(123)
k_clusters <- 5
kmeans_res <- kmeans(tsne_df[, c("tSNE1", "tSNE2")], centers = k_clusters, nstart = 25)
tsne_df$Cluster <- factor(kmeans_res$cluster)

cluster_plot <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = Cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = paste0("K-means Clusters (k=", k_clusters, ") on t-SNE"), x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

ggsave(filename = file.path(project_dir, "tsne_kmeans_clusters.png"), plot = cluster_plot, width = 7, height = 5)
cat("t-SNE clustering plot saved.\n\n")

# Compare clusters to subtypes
cat("Comparing clusters to PAM50 subtypes - started at", Sys.time(), "\n")
cluster_table <- table(Cluster = tsne_df$Cluster, Subtype = tsne_df$Subtype)
print(cluster_table)

heatmap_file <- file.path(project_dir, "tsne_cluster_vs_subtype_heatmap.png")
png(filename = heatmap_file, width = 600, height = 500)
pheatmap(cluster_table, cluster_rows = FALSE, cluster_cols = FALSE, main = "t-SNE Cluster vs Subtype")
dev.off()
cat("Cluster vs subtype heatmap saved.\n\n")

# Train/test split
cat("Splitting data into training/testing - started at", Sys.time(), "\n")
set.seed(123)
trainIndex <- createDataPartition(vsd_df$Subtype, p = 0.7, list = FALSE)
train_data <- vsd_df[trainIndex, ]
test_data <- vsd_df[-trainIndex, ]

# Ensure factor levels are consistent
train_data$Subtype <- factor(train_data$Subtype)
test_data$Subtype <- factor(test_data$Subtype, levels = levels(train_data$Subtype))
cat("Data split complete.\n\n")

# Train SVM model
cat("Training SVM model - started at", Sys.time(), "\n")
set.seed(123)
svm_tune <- train(Subtype ~ ., data = train_data,
                  method = "svmRadial",
                  preProcess = c("center", "scale"),
                  tuneLength = 5,
                  trControl = trainControl(method = "cv", number = 5, verboseIter = TRUE))
cat("SVM training finished.\n")

svm_tune_plot <- ggplot(svm_tune)
ggsave(filename = file.path(project_dir, "svm_tuning_plot.png"), plot = svm_tune_plot, width = 7, height = 5)
cat("SVM tuning plot saved.\n\n")

# Predict and evaluate SVM
cat("Predicting with SVM and evaluating - started at", Sys.time(), "\n")
test_pred <- predict(svm_tune, newdata = test_data)
conf_mat <- confusionMatrix(test_pred, test_data$Subtype)
print(conf_mat)
cat("SVM evaluation done.\n\n")

# SVM variable importance
cat("Calculating SVM variable importance - started at", Sys.time(), "\n")
var_imp <- varImp(svm_tune)
print(var_imp)

var_imp_plot <- ggplot(var_imp, top = 20)
ggsave(filename = file.path(project_dir, "svm_variable_importance.png"), plot = var_imp_plot, width = 7, height = 5)
cat("SVM variable importance plot saved.\n\n")

# Train Random Forest model
cat("Training Random Forest model - started at", Sys.time(), "\n")
set.seed(123)
rf_model <- train(Subtype ~ ., data = train_data,
                  method = "rf",
                  importance = TRUE,
                  trControl = trainControl(method = "cv", number = 5))
cat("Random Forest training finished.\n")

rf_pred <- predict(rf_model, newdata = test_data)
rf_cm <- confusionMatrix(rf_pred, test_data$Subtype)
print(rf_cm)
cat("Random Forest Accuracy:", round(rf_cm$overall['Accuracy'], 4), "\n\n")

rf_imp <- varImp(rf_model)
rf_plot <- ggplot(rf_imp, top = 20)
ggsave(filename = file.path(project_dir, "rf_variable_importance.png"), plot = rf_plot, width = 7, height = 5)

# Train XGBoost model 
cat("Training XGBoost model - started at", Sys.time(), "\n")
set.seed(123)
xgb_model <- train(Subtype ~ ., data = train_data,
                   method = "xgbTree",
                   trControl = trainControl(method = "cv", number = 5),
                   verbosity = 0)
cat("XGBoost training finished.\n")

xgb_pred <- predict(xgb_model, newdata = test_data)
xgb_cm <- confusionMatrix(xgb_pred, test_data$Subtype)
print(xgb_cm)
cat("XGBoost Accuracy:", round(xgb_cm$overall['Accuracy'], 4), "\n\n")

xgb_imp <- varImp(xgb_model)
xgb_plot <- ggplot(xgb_imp, top = 20)
ggsave(filename = file.path(project_dir, "xgb_variable_importance.png"), plot = xgb_plot, width = 7, height = 5)

# Save model accuracy comparison
cat("Saving model accuracy comparison...\n")
acc_df <- data.frame(
  Model = c("SVM", "Random Forest", "XGBoost"),
  Accuracy = c(conf_mat$overall['Accuracy'],
               rf_cm$overall['Accuracy'],
               xgb_cm$overall['Accuracy'])
)
write.csv(acc_df, file.path(project_dir, "model_accuracy_comparison.csv"), row.names = FALSE)
print(acc_df)
cat("Model accuracy comparison saved.\n\n")