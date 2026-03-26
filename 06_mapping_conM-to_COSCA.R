### Author Qamar Feddah
### Script map ConM data to reference data. 

## Library loading and set the working directory
library(Seurat)
library(tidyverse)
set.seed(333)

project_dir <- "../.."
setwd(project_dir)

## Load the Seurat objects
reference_seurat <- readRDS("Analysis/06_Reference mapping/reference_data.RDS")
ConM_seurat <- readRDS("Analysis/02_Data_Normalization/seurat_conM_Clustered.RDS")

## Mapping Pipeline
DefaultAssay(reference_seurat) <- "SCT"
DefaultAssay(ConM_seurat) <- "SCT"

mapping.anchors <- FindTransferAnchors(
  reference = reference_seurat,
  query = ConM_seurat,
  dims = 1:30,
  reference.reduction = "pca",
  normalization.method = "SCT",
  reference.assay = "SCT",
  query.assay = "SCT"
)

prediction_cluster <- TransferData(
  anchorset = mapping.anchors,
  refdata = reference_seurat$seurat_clusters,
  dims = 1:30
)

ConM_seurat <- AddMetaData(
  ConM_seurat,
  metadata = prediction_cluster$predicted.id,
  col.name = "reference_analogy"
)

## Visualise the ConM UMAP with the reference information
umap_conm <- DimPlot(ConM_seurat)
umap_with_reference <- DimPlot(ConM_seurat, group.by = "reference_analogy", label = TRUE)
umap_reference <- DimPlot(reference_seurat, reduction = "wnn.umap", label = TRUE)
umap_conm | umap_with_reference | umap_reference

## Annotate the cells based on all information
final_clusters <- rep(NA, times = ncol(ConM_seurat))
names(final_clusters) <- colnames(ConM_seurat)

naive_cells <- rownames(ConM_seurat@meta.data[ConM_seurat@meta.data$seurat_clusters == "2", ])
final_clusters[naive_cells] <- "Naive cells"

atypical_cells <- rownames(ConM_seurat@meta.data[
  ConM_seurat@meta.data$reference_analogy == "17" &
    ConM_seurat@meta.data$seurat_clusters == "6", ])
final_clusters[atypical_cells] <- "Atypical cells"

plasmablasts <- rownames(ConM_seurat@meta.data[ConM_seurat@meta.data$seurat_clusters == "7", ])
final_clusters[plasmablasts] <- "Plasmablasts"

act_b_cells <- rownames(ConM_seurat@meta.data[
  ConM_seurat@meta.data$reference_analogy == "10" &
    ConM_seurat@meta.data$seurat_clusters == "1", ])
final_clusters[act_b_cells] <- "Activated B cells"

act_mem_b_cells <- rownames(ConM_seurat@meta.data[ConM_seurat@meta.data$reference_analogy == "3", ])
final_clusters[act_mem_b_cells] <- "Activated B cells"

cd24_mem <- rownames(ConM_seurat@meta.data[ConM_seurat@meta.data$reference_analogy == "6", ])
final_clusters[cd24_mem] <- "Resting Memory B cells"

cd73_mem <- rownames(ConM_seurat@meta.data[ConM_seurat@meta.data$reference_analogy == "7", ])
final_clusters[cd73_mem] <- "Resting Memory B cells"

igm_mem1 <- rownames(ConM_seurat@meta.data[ConM_seurat@meta.data$reference_analogy == "8", ])
final_clusters[igm_mem1] <- "IgM Memory B cells"

igm_mem2 <- rownames(ConM_seurat@meta.data[ConM_seurat@meta.data$reference_analogy == "9", ])
final_clusters[igm_mem2] <- "IgM Memory B cells"

ConM_seurat <- AddMetaData(ConM_seurat, metadata = final_clusters, col.name = "final_clusters")

## Save the object with NAs
saveRDS(ConM_seurat, file = "Analysis/06_Reference mapping/ConM_seurat_finetune_clustering_NAs.RDS")

## Remove NA cells
valid_cells <- rownames(ConM_seurat@meta.data)[!is.na(ConM_seurat@meta.data$final_clusters)]
ConM_seurat <- subset(ConM_seurat, cells = valid_cells)

## Save the filtered object
saveRDS(ConM_seurat, file = "Analysis/06_Reference mapping/ConM_seurat_finetune_cluster.RDS")

## Produce the interpretation plots
fine_tuned_clustering_umap_plot <- DimPlot(ConM_seurat, group.by = "final_clusters") +
  ggtitle("Clustering after Fine Tuning"); fine_tuned_clustering_umap_plot

## Produce the PDF
pdf(file = "Analysis/06_Reference mapping/ConM_Mapping_to_Reference.pdf", width = 10, height = 6)
print(umap_conm | umap_with_reference | umap_reference)
print(fine_tuned_clustering_umap_plot)
dev.off()
