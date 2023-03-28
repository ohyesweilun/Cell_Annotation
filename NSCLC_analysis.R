#install.packages("tidyverse")
#install.packages("Matrix")
#install.packages("RCurl")
#install.packages("scales")
#install.packages("devtools")
#install.packages("BiocManager")
#install.packages("Seurat")
#install.packages("readr")
#install.packages("installr")
#remotes::install_github("pmbio/MuDataSeurat")

# Load library
library(Seurat)
library(dplyr)
library(Matrix)
library(readr)
library(installr)
library(MuDataSeurat)

# Set default folder director 
setwd("D:/case study/pilot_nuclei/pilot_nuclei")

# Read 10x data
pbmc.data <- Read10X(data.dir = "./")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "NSCLC", min.cells = 3, min.features = 200)
pbmc.data
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells that have unique feature counts over 2,500 or less than 200
# and cells that have >5% mitochondrial counts
pbmc.filter <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Apply global-scaling normalization method “LogNormalize” 
# that normalizes the feature expression measurements for each cell by the total expression
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
pbmc.filter <- NormalizeData(pbmc.filter, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.filter <- NormalizeData(pbmc.filter)

# Feature selection
pbmc.filter.feature <- FindVariableFeatures(pbmc.filter, selection.method = "vst", nfeatures = 2000)
pbmc.filter.feature

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc.filter.feature), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc.filter.feature)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data: 
# mean expression across cells is 0 & the variance across cells is 1
all.genes <- rownames(pbmc.filter.feature)
pbmc.filter.feature.scale <- ScaleData(pbmc.filter.feature) # features = all.genes)
pbmc.filter.feature.scale

# Run PCA to reduce dimensions
pbmc.filter.feature.scale.reduce <- RunPCA(pbmc.filter.feature.scale, features = VariableFeatures(object = pbmc.filter.feature.scale))
# Examine and visualize PCA results a few different ways
print(pbmc.filter.feature.scale.reduce[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc.filter.feature.scale.reduce, dims = 1:2, reduction = "pca")

DimPlot(pbmc.filter.feature.scale.reduce, reduction = "pca")

# Explore the heterogeneity in the dataset
DimHeatmap(pbmc.filter.feature.scale.reduce, dims = 1, cells = 500, balanced = TRUE)
# Print all 9 PCA clusters
DimHeatmap(pbmc.filter.feature.scale.reduce, dims = 1:9, cells = 500, balanced = TRUE)

# Determine the dimentionality of the dataset
pbmc.filter.feature.scale.reduce <- JackStraw(pbmc.filter.feature.scale.reduce, num.replicate = 100)
pbmc.filter.feature.scale.reduce <- ScoreJackStraw(pbmc.filter.feature.scale.reduce, dims = 1:20)
JackStrawPlot(pbmc.filter.feature.scale.reduce, dims = 1:20)

# Alternative plot
ElbowPlot(pbmc.filter.feature.scale.reduce)

# Cluster the cells (nearest neighbor graph)
pbmc.filter.feature.scale.reduce <- FindNeighbors(pbmc.filter.feature.scale.reduce, dims = 1:10)
pbmc.filter.feature.scale.reduce <- FindClusters(pbmc.filter.feature.scale.reduce, resolution = 0.5)
head(Idents(pbmc.filter.feature.scale.reduce), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
pbmc.filter.feature.scale.reduce <- RunUMAP(pbmc.filter.feature.scale.reduce, dims = 1:10)
DimPlot(pbmc.filter.feature.scale.reduce, reduction = "umap")

# Save the object for future use
saveRDS(pbmc, file = "./pbmc_tutorial.rds")

# Find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc.filter.feature.scale.reduce, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# Find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc.filter.feature.scale.reduce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
View(pbmc.markers) 

# ROC for each cluster
cluster0.markers <- FindMarkers(pbmc.filter.feature.scale.reduce, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster0.markers 

# Visualize marker expression
VlnPlot(pbmc.filter.feature.scale.reduce, features = c("AGBL1", "NCKAP5"))
# Plot raw counts as well
VlnPlot(pbmc.filter.feature.scale.reduce, features = c("AGBL1", "NCKAP5"), slot = "counts", log = TRUE)

# Feature Plot
FeaturePlot(pbmc.filter.feature.scale.reduce, features = c("AGBL1", "NCKAP5", "AOAH", "NTM", "ROBO2"))

# Heatmap
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top3
DoHeatmap(pbmc.filter.feature.scale.reduce, features = top3$gene) + NoLegend()


# Prepare the count matrix for CellTypist


"
# Read in `matrix.mtx`
counts <- readMM("./matrix.mtx.gz")

# Read in `genes.tsv`
genes <- read_tsv("./features.tsv.gz", col_names = FALSE)
gene_ids <- genes$X1

# Read in `barcodes.tsv`
cell_ids <- read_tsv("./barcodes.tsv.gz", col_names = FALSE)$X1

# Make the column names as the cell IDs and the row names as the gene IDs
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
#write.csv (counts, './count.matrix.csv')
"

# Convert Seurat to Scanpy h5ad
# step 1: Slim down a Seurat object. So you get raw counts, lognorm counts
seu = DietSeurat(
  pbmc.filter.feature.scale.reduce,
  counts = TRUE, # so, raw counts save to adata.layers['counts']
  data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
  scale.data = FALSE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
  features = rownames(pbmc.filter.feature.scale.reduce), # export all genes, not just top highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)

# single modality
MuDataSeurat::WriteH5AD(pbmc.filter, "pbmc.filter.h5ad", assay="RNA")
# multi modality, ATAC+RNA, CITE-seq 
MuDataSeurat::WriteH5MU(pbmc.filter, "pbmc.filter.h5mu")


#　Assign cell type identity to clusters
#  Cluster ID	Markers	Cell Type
# 0	AGBL1, LHFPL3	Naive CD4+ T
# 1	NCKAP5, NCKAP5-AS2	CD14+ Mono
# 2	AOAH, CD247	Memory CD4+
# 3	NTM, AC022325.1	B
# 4	ROBO2, LAMA2	CD8+ T
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
