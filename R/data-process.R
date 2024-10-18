options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

# Install zellkonverter, which includes dependencies for reading h5ad files
BiocManager::install(c("zellkonverter", "SingleCellExperiment"))

library(reticulate)
conda_list()
reticulate::use_condaenv(condaenv = 'r_studio',required = TRUE)

BiocManager::install("scater")

if (!requireNamespace("uwot", quietly = TRUE)) {
  install.packages("uwot")
}
    

install.packages("ggplot2")       

library(zellkonverter)  # For reading and writing h5ad files
library(SingleCellExperiment)  # For handling single-cell data

library(scater)
library(uwot)

library(ggplot2)


# Provide the path to your h5ad file
h5ad_file <- "/home/javeed/breast-cancer/HTAPP-783-SMP-4081_scRNAseq_processed.h5ad"

# Read the h5ad file into an R object (SingleCellExperiment)
sce <- readH5AD(h5ad_file)


colnames(colData(sce))

head(colData(sce)$cell_type)


table(colData(sce)$cell_type)  # Adjust 'cell_type' to the actual name



# Perform PCA (usually required before UMAP)
sce <- runPCA(sce)

# Perform UMAP using the reduced PCA dimensions
sce <- runUMAP(sce, dimred = "PCA")

# Plot the UMAP, coloring by cell type (or any other column in colData)
plotUMAP(sce, colour_by = "cell_type")  # Replace 'cell_type' with your actual metadata column name



# Plot UMAP with cells colored by cell type
plotUMAP(sce, colour_by = "cell_type")  # Adjust 'cell_type' to the actual column name

reducedDimNames(sce)


reducedDim(sce, "X_umap")




umap_data <- reducedDim(sce, "X_umap")
umap_df <- as.data.frame(umap_data)
colnames(umap_df) <- c("UMAP1", "UMAP2")

ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = sce$cell_type)) +  # Replace 'your_cell_annotation' with a column in your SCE object
  theme_minimal() +
  ggtitle("UMAP Visualization")

# View the structure of the SingleCellExperiment object
sce


counts(sce)


library(ggplot2)

# Extract UMAP coordinates from the SingleCellExperiment object
umap_data <- data.frame(
  UMAP1 = reducedDim(sce, "UMAP")[,1],
  UMAP2 = reducedDim(sce, "UMAP")[,2],
  cell_type = colData(sce)$cell_type  # Replace 'cell_type' with your actual column name
)

# Plot UMAP with ggplot2
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = 1.5) +
  labs(title = "UMAP of Single-Cell Data", x = "UMAP1", y = "UMAP2") +
  theme_minimal()


