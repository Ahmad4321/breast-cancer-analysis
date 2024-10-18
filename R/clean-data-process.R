# to Setup the download mirror
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install zellkonverter, which includes dependencies for reading h5ad files
BiocManager::install(c("zellkonverter", "SingleCellExperiment","reticulate","ggplot2"))

library(reticulate)
conda_list()
reticulate::use_condaenv(condaenv = 'r_studio',required = TRUE)


library(zellkonverter)  # For reading and writing h5ad files
library(SingleCellExperiment)  # For handling single-cell data


library(ggplot2)


# Provide the path to your h5ad file
h5ad_file <- "/home/javeed/breast-cancer/HTAPP-783-SMP-4081_scRNAseq_processed.h5ad"

# Read the h5ad file into an R object (SingleCellExperiment)
sce <- readH5AD(h5ad_file)

# Display the columns of dataset
colnames(colData(sce))


#display cell type list
head(colData(sce)$cell_type)



#display count of each cell type in dataset
table(colData(sce)$cell_type)  # Adjust 'cell_type' to the actual name


# Reduced Dimensions
umap_data <- reducedDim(sce, "X_umap")

#Make a matrix so ggplot can understand the data
umap_df <- as.data.frame(umap_data)
colnames(umap_df) <- c("UMAP1", "UMAP2")

ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = sce$cell_type)) +  # Replace 'your_cell_annotation' with a column in your SCE object
  theme_minimal() +
  ggtitle("Breast Cancer: Cell Type  Visualization")

