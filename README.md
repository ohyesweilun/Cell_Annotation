# Cell_Annotation

NSCLC_analysis.R
1. Read 10x data into a Seurat object
2. Visualize QC metrics
3. Apply filtering 
4. Normalize the count and log-transform
5. Identify top10 most highly varible genes
6. Scale the data
7. Run PCA to reduce features
8. DimPlot, DimHeatmap, JackStrawPlot, ElbowPlot
9. Cluster the cells (nearest neighbor graph)
10. Run tSNE and UMAP
11. Find markers for every cluster copmaring to all remaining cells
12. Visualized marker expression- VlnPlot, FeaturePlot and Heatmap
13. Convert Seurat object to Scanpy h5ad format (Save .h5ad and .h5mu files)

CellTypist.ipynb
1. Annotate cells using CellTypist (Models: Human_Lung_Atlas & Immune_All_High)
2. Explore cell types by samples and top 3 gene markers
