## Library loading
library(Seurat)
library(tidyverse)
library(scCustomize)
set.seed(333)

# Set working directory
setwd("C:/Users/Qamar/Documents/ConM_project/Analysis/05_supervised_Analysis")
## Load the seurat Object
seurat_conM <- readRDS("~/ConM_project/Analysis/02_Data_Normalization/seurat_conM_Clustered.RDS")

## Visualization of Genes and dotplots
feat_trans <- c("TCL1A", "CD24", "CD38", "IL4R", "VPREB3", "IGLL1", "IRF8","CD72")
feat_naive <- c("IGHM", "IGHD", "CD19", "MS4A1", "TCL1A", "IL4R", "CR2")
feat_mem <- c("CD27", "CD24", "MS4A1", "BANK1", "TNFRSF13B", "AICDA")
feat_pb <- c("PRDM1", "XBP1", "IRF4", "MZB1", "CD38", "CD27")
feat_plasma <- c("PRDM1", "XBP1", "IRF4", "MZB1", "SDC1", "JCHAIN", "TNFRSF17")
feat_aty <- c("CD72", "SOX5", "FCRL5")
feat_canon_b <- unique(c(feat_trans, feat_naive, feat_mem, feat_pb, feat_plasma,feat_aty))

## Plot production
cluster_plot <- DimPlot(seurat_conM, group.by = "seurat_clusters"); cluster_plot
trans_plot <- FeaturePlot_scCustom(seurat_conM, features = feat_trans); trans_plot
naive_plot <- FeaturePlot_scCustom(seurat_conM, features = feat_naive); naive_plot
mem_plot <- FeaturePlot_scCustom(seurat_conM, features = feat_mem); mem_plot
pb_plot <- FeaturePlot_scCustom(seurat_conM, features = feat_pb); pb_plot
plasma_plot <- FeaturePlot_scCustom(seurat_conM, features = feat_plasma); plasma_plot
aty_plot <- FeaturePlot_scCustom(seurat_conM, features = feat_aty); aty_plot
clust_dot_plot <- Clustered_DotPlot(seurat_conM, features = feat_canon_b, plot_km_elbow = F); clust_dot_plot

# Create PDF
pdf(file = "Canonical_UMAP_B_Cell_Supervised_Analysis.pdf")
print(cluster_plot)
print(trans_plot)
print(naive_plot)
print(mem_plot)
print(pb_plot)
print(plasma_plot)
print(clust_dot_plot)
dev.off()

pdf(file = "Canonical_DotPlot_B_Cell_Supervised_Analysis.pdf")
print(clust_dot_plot)
dev.off()

