# Load libraries with timestamps
cat("Loading libraries - started at", Sys.time(), "\n")
library(tidyverse)
library(caret)
library(randomForest)
library(DESeq2)
library(Rtsne)
library(e1071)        # for svm
library(factoextra)   # for clustering visualization
library(pheatmap)
library(RColorBrewer)
cat("Loading libraries - finished at", Sys.time(), "\n\n")

# Load processed expression data
cat("Loading data - started at", Sys.time(), "\n")
load("E:/genomics/tcga brca data/tcga_brca_expr_df.RData")  # loads expr_df
cat("Loading data - finished at", Sys.time(), "\n\n")

# Prepare count matrix: genes x samples
cat("Preparing count matrix - started at", Sys.time(), "\n")
count_matrix <- as.matrix(t(expr_df %>% select(-Subtype)))
colnames(count_matrix) <- paste0("Sample", 1:ncol(count_matrix))  # sample names
rownames(count_matrix) <- colnames(expr_df)[-ncol(expr_df)]       # gene IDs
cat("Preparing count matrix - finished at", Sys.time(), "\n\n")

# Select top 2000 most variable genes
cat("Selecting top 2000 variable genes - started at", Sys.time(), "\n")
gene_vars <- apply(count_matrix, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:2000]
count_matrix <- count_matrix[top_genes, ]
cat("Selecting top 2000 variable genes - finished at", Sys.time(), "\n\n")

# Create sample metadata dataframe
cat("Creating sample metadata - started at", Sys.time(), "\n")
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  Subtype = factor(expr_df$Subtype)  # make sure factor
)
cat("Creating sample metadata - finished at", Sys.time(), "\n\n")

# Print dimension info and subtype counts
cat("Dimensions of count_matrix: ", dim(count_matrix), "\n")
cat("Table of Subtypes:\n")
print(table(expr_df$Subtype))
cat("\n")

# Create DESeq2 dataset for normalization and variance stabilization
cat("Creating DESeqDataSet - started at", Sys.time(), "\n")
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Subtype)
cat("Creating DESeqDataSet - finished at", Sys.time(), "\n\n")

# Perform Variance Stabilizing Transformation (VST)
cat("Performing VST - started at", Sys.time(), "\n")
vsd <- vst(dds, blind = TRUE)
cat("Performing VST - finished at", Sys.time(), "\n\n")

# Transpose VST matrix back to samples x genes
vsd_matrix <- t(assay(vsd))

# Prepare final dataframe with subtype labels
cat("Preparing vsd_df - started at", Sys.time(), "\n")
vsd_df <- as.data.frame(vsd_matrix)
vsd_df$Subtype <- sample_info$Subtype
cat("Preparing vsd_df - finished at", Sys.time(), "\n\n")

# Perform PCA on VST data
cat("Performing PCA - started at", Sys.time(), "\n")
pca_res <- prcomp(vsd_matrix, scale. = TRUE)
cat("Performing PCA - finished at", Sys.time(), "\n\n")

# Prepare PCA plot data
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$Subtype <- vsd_df$Subtype

# Plot PCA and save as PNG
cat("Plotting PCA - started at", Sys.time(), "\n")
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Subtype)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA of TCGA BRCA expression data",
       x = "PC1", y = "PC2") +
  theme_minimal()
ggsave(filename = "E:/genomics/tcga brca data/pca_plot.png", plot = pca_plot, width = 7, height = 5)
cat("Plotting PCA - finished at", Sys.time(), "\n\n")

# Perform t-SNE on top 50 PCs from PCA for better visualization
cat("Performing t-SNE - started at", Sys.time(), "\n")
pca_for_tsne <- pca_res$x[, 1:50]

set.seed(123)
tsne_out <- Rtsne(pca_for_tsne, perplexity = 30, verbose = TRUE, max_iter = 1000)
cat("Performing t-SNE - finished at", Sys.time(), "\n\n")

# Prepare t-SNE plot data
tsne_df <- as.data.frame(tsne_out$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df$Subtype <- vsd_df$Subtype

# Plot t-SNE and save as PNG
cat("Plotting t-SNE - started at", Sys.time(), "\n")
tsne_plot <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = Subtype)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "t-SNE of TCGA BRCA expression data",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()
ggsave(filename = "E:/genomics/tcga brca data/tsne_plot.png", plot = tsne_plot, width = 7, height = 5)
cat("Plotting t-SNE - finished at", Sys.time(), "\n\n")

# --- Clustering on t-SNE coordinates ---
cat("Performing clustering on t-SNE - started at", Sys.time(), "\n")
# For example: k-means with 5 clusters
set.seed(123)
k_clusters <- 5
kmeans_res <- kmeans(tsne_df[, c("tSNE1", "tSNE2")], centers = k_clusters, nstart = 25)
tsne_df$Cluster <- factor(kmeans_res$cluster)

# Plot clusters on t-SNE and save
cluster_plot <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = Cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = paste0("K-means Clusters (k=", k_clusters, ") on t-SNE"),
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()
ggsave(filename = "E:/genomics/tcga brca data/tsne_kmeans_clusters.png", plot = cluster_plot, width = 7, height = 5)
cat("Clustering on t-SNE - finished at", Sys.time(), "\n\n")

# --- Prepare data for classification ---
cat("Splitting data into training/testing - started at", Sys.time(), "\n")
set.seed(123)
trainIndex <- createDataPartition(vsd_df$Subtype, p = 0.7, list = FALSE)
train_data <- vsd_df[trainIndex, ]
test_data <- vsd_df[-trainIndex, ]

# Ensure factors have same levels
train_data$Subtype <- factor(train_data$Subtype)
test_data$Subtype <- factor(test_data$Subtype, levels = levels(train_data$Subtype))
cat("Splitting data into training/testing - finished at", Sys.time(), "\n\n")

# --- Compare Clusters to PAM50 Subtypes ---
cat("Comparing t-SNE clusters with subtypes - started at", Sys.time(), "\n")
tsne_df$Subtype <- vsd_df$Subtype
cluster_table <- table(Cluster = tsne_df$Cluster, Subtype = tsne_df$Subtype)
print(cluster_table)

# Heatmap of cluster vs subtype
heatmap_file <- "E:/genomics/tcga brca data/tsne_cluster_vs_subtype_heatmap.png"
png(filename = heatmap_file, width = 600, height = 500)
pheatmap(cluster_table, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "t-SNE Cluster vs Subtype")
dev.off()
cat("Comparing t-SNE clusters with subtypes - finished at", Sys.time(), "\n\n")

# --- Train SVM with radial kernel and cross-validation ---
cat("Training SVM model - started at", Sys.time(), "\n")
set.seed(123)
svm_tune <- train(Subtype ~ ., data = train_data,
                  method = "svmRadial",
                  preProcess = c("center", "scale"),
                  tuneLength = 5,
                  trControl = trainControl(method = "cv", number = 5, verboseIter = TRUE))
cat("Training SVM model - finished at", Sys.time(), "\n\n")

# Save SVM tuning results plot
svm_tune_plot <- ggplot(svm_tune)
ggsave(filename = "E:/genomics/tcga brca data/svm_tuning_plot.png", plot = svm_tune_plot, width = 7, height = 5)

# --- Predict on test data and evaluate ---
cat("Predicting and evaluating SVM model - started at", Sys.time(), "\n")
test_pred <- predict(svm_tune, newdata = test_data)
conf_mat <- confusionMatrix(test_pred, test_data$Subtype)
print(conf_mat)
cat("Predicting and evaluating SVM model - finished at", Sys.time(), "\n\n")

# --- Variable importance from SVM ---
cat("Calculating variable importance - started at", Sys.time(), "\n")
var_imp <- varImp(svm_tune)
print(var_imp)
var_imp_plot <- ggplot(var_imp, top = 20)
ggsave(filename = "E:/genomics/tcga brca data/svm_variable_importance.png", plot = var_imp_plot, width = 7, height = 5)
cat("Calculating variable importance - finished at", Sys.time(), "\n\n")

# --- Random Forest Model ---
cat("Training Random Forest model - started at", Sys.time(), "\n")
set.seed(123)
rf_model <- train(Subtype ~ ., data = train_data,
                  method = "rf",
                  importance = TRUE,
                  trControl = trainControl(method = "cv", number = 5))
rf_pred <- predict(rf_model, newdata = test_data)
rf_cm <- confusionMatrix(rf_pred, test_data$Subtype)
print(rf_cm)
cat("Random Forest Accuracy:", round(rf_cm$overall['Accuracy'], 4), "\n")
cat("Training Random Forest model - finished at", Sys.time(), "\n\n")

# --- Random Forest Variable Importance ---
rf_imp <- varImp(rf_model)
rf_plot <- ggplot(rf_imp, top = 20)
ggsave("E:/genomics/tcga brca data/rf_variable_importance.png", plot = rf_plot, width = 7, height = 5)

# --- XGBoost Model ---
cat("Training XGBoost model - started at", Sys.time(), "\n")
set.seed(123)
xgb_model <- train(Subtype ~ ., data = train_data,
                   method = "xgbTree",
                   trControl = trainControl(method = "cv", number = 5),
                   verbosity = 0)
xgb_pred <- predict(xgb_model, newdata = test_data)
xgb_cm <- confusionMatrix(xgb_pred, test_data$Subtype)
print(xgb_cm)
cat("XGBoost Accuracy:", round(xgb_cm$overall['Accuracy'], 4), "\n")
cat("Training XGBoost model - finished at", Sys.time(), "\n\n")

# --- XGBoost Variable Importance ---
xgb_imp <- varImp(xgb_model)
xgb_plot <- ggplot(xgb_imp, top = 20)
ggsave("E:/genomics/tcga brca data/xgb_variable_importance.png", plot = xgb_plot, width = 7, height = 5)

# --- Compare Model Accuracies ---
cat("Saving model accuracy comparison...\n")
acc_df <- data.frame(
  Model = c("SVM", "Random Forest", "XGBoost"),
  Accuracy = c(conf_mat$overall['Accuracy'],
               rf_cm$overall['Accuracy'],
               xgb_cm$overall['Accuracy'])
)
write.csv(acc_df, "E:/genomics/tcga brca data/model_accuracy_comparison.csv", row.names = FALSE)
print(acc_df)
cat("Model accuracy comparison saved.\n\n")



