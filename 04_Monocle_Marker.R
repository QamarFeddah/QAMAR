##Author: Qamar Feddah
#script for performing graph based marker identification

# Load libraries
set.seed(333)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(tidyverse)
library(ggplot2)
library(scCustomize)

# read in the data
setwd("~/ConM_project/Analysis/04_Graph_based_marker_definition")
seurat_conM <- readRDS("~/ConM_project/Analysis/02_Data_Normalization/seurat_conM_Clustered.RDS")
DefaultAssay(seurat_conM) <- "SCT"

# monocle3 requires cell_data_set object
cds_conM <- as.cell_data_set(seurat_conM)
cds_conM <- preprocess_cds(cds_conM, num_dim = 30) #runs PCA and saves the first 30 PCs.
cds_conM <- reduce_dimension(cds_conM, reduction_method = "UMAP")  #computes a 2D UMAP from those PCs, and get reducedDims(cds)$UMAP
cds_conM <- cluster_cells(cds_conM, reduction_method = "UMAP", cluster_method = "leiden")  # buids a kNN (k-nearest neighbors) graph on the UMAP embedding and finds communities with the Leiden method.

## Add the gene_short_name column to the cds object
# To get cell metadata
colData(cds_conM)
# To get feature (gene) metadata
fData(cds_conM)
rownames(fData(cds_conM))
#create a new column and assign it with the row names of feature data
fData(cds_conM)$gene_short_name <- rownames(fData(cds_conM))
# to build a trajectory , we use the clustering performed previously is seurat, so we need
#to retrieve clustering info and store it inside the cds object
#monocle3 has cluster functions, which determines not only the clusters but also partitions
#partitions are super clusters
##we need to do three things:
#1-assign all cells to one partition:One partition ensures a single trajectory across all cells
recreate.partition <- c(rep(1,length(cds_conM@colData@rownames)))
names(recreate.partition) <- cds_conM@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_conM@clusters$UMAP$partitions <- recreate.partition

##2-assign the cluster info: Seurat clusters keep our grouping consistent with our upstream analysis
cluster_info <- as.numeric(seurat_conM@meta.data[["seurat_clusters"]])
names(cluster_info) <- cds_conM@colData@rownames
cds_conM@clusters$UMAP$clusters <- as.factor(cluster_info)

##3-assign umap coordinates_cell embeddings from seurat into cds object
umap_info <- seurat_conM@reductions$umap@cell.embeddings
row.names(umap_info) <- cds_conM@colData@rownames
cds_conM@int_colData@listData$reducedDims$UMAP <- umap_info

## learn the trajectory graph
cds_conM <- learn_graph(cds_conM, use_partition = FALSE)

## Differentially expressed genes graph
differential_genes <- graph_test(cds_conM, neighbor_graph = 'knn', expression_family = "negbinomial", cores = 19)

top_genes <-
  differential_genes %>%
  filter(status == 'OK') %>%
  filter(q_value <= 0.05) %>%
  arrange(desc(morans_I)) %>%
  head(n = 300L); top_genes

top_select <- row.names(subset(top_genes))
write_csv(top_genes, file = "Top_Gene_Graph_based_Markers_Monocle3.csv")
gene_module_df <- find_gene_modules(cds_conM[top_select,], resolution=1e-2)


## Plot differential genes
plot_cells(cds_conM, genes=row.names(top_genes)[1:50],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) -> deg_feature_plot_1; deg_feature_plot_1

plot_cells(cds_conM, genes=row.names(top_genes)[51:100],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) -> deg_feature_plot_2; deg_feature_plot_2

plot_cells(cds_conM, genes=row.names(top_genes)[101:150],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) -> deg_feature_plot_3; deg_feature_plot_3

plot_cells(cds_conM, genes=row.names(top_genes)[151:200],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) -> deg_feature_plot_4; deg_feature_plot_4

plot_cells(cds_conM, genes=row.names(top_genes)[201:250],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) -> deg_feature_plot_5; deg_feature_plot_5

plot_cells(cds_conM, genes=row.names(top_genes)[251:300],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) -> deg_feature_plot_6; deg_feature_plot_6

top_select <- row.names(subset(top_genes))
gene_module_df <- find_gene_modules(cds_conM[top_select,], resolution=1e-2)

write_csv(gene_module_df, file = "Gene_Modules_conM_Monocle3.csv")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_conM)), 
                                cell_group=clusters(cds_conM)[colnames(cds_conM)])
agg_mat <- aggregate_gene_expression(cds_conM, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("", colnames(agg_mat))

htm <- pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                          scale="column", clustering_method="ward.D2",
                          fontsize=10)

umap_plot <-
  plot_cells(cds_conM, color_cells_by = "cluster", show_trajectory_graph = FALSE,
             graph_label_size = 8, group_label_size = 8, cell_size = 0.75)

DefaultAssay(seurat_conM) <- "SCT"
# dot_plot <- DotPlot(seurat_conM, features = top_genes$gene_short_name[1:75]) +
#   RotatedAxis()
dotplot_cl <- Clustered_DotPlot(seurat_conM, features = top_genes$gene_short_name[1:300], 
                                assay = "SCT", plot_km_elbow = FALSE)

## Write the data plot on pdf
pdf(file = "conM_Graph_based_RNA_markers.pdf", height = 11, width = 18)
print(htm)
print(umap_plot)
print(deg_feature_plot_1)
print(deg_feature_plot_2)
print(deg_feature_plot_3)
print(deg_feature_plot_4)
print(deg_feature_plot_5)
print(deg_feature_plot_6)
dev.off()

pdf(file = "conM_Graph_based_RNA_markers_dotplot.pdf", height = 50, width = 8)
print(dotplot_cl)
dev.off()