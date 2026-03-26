### Author Qamar Feddah
### Script map ConM data to reference data. 

## Library loading and set the working directory
library(Seurat)
library(tidyverse)
set.seed(333)

project_dir <- "../.."
setwd(project_dir)

## Load the Seurat objects
reference_seurat <- readRDS("Analysis/06_Reference mapping/reference_data")
query_seurat <- readRDS("Analysis/02_Data_Normalization/seurat_conM_Clustered.RDS")

## Mapping Pipeline
DefaultAssay(seurat_reference) <- "SCT"
DefaultAssay(seurat_ConM) <- "SCT"
mapping.anchors <- FindTransferAnchors(reference = seurat_reference, query = seurat_ConM, dims = 1:30,
                                       reference.reduction = "pca", normalization.method = "SCT",
                                       reference.assay = "SCT", query.assay = "SCT")
prediction_cluster <- TransferData(anchorset = mapping.anchors, refdata = seurat_reference$seurat_clusters, dims = 1:30)
seurat_ConM <- AddMetaData(seurat_ConM, metadata = prediction_cluster$predicted.id, col.name = "reference_analogy")

## Vizualise the ConM UMAP wuth the reference information
umap_conm <- DimPlot(seurat_ConM)
umap_with_reference <- DimPlot(seurat_ConM, group.by = "reference_analogy", label = TRUE)
umap_reference <-DimPlot(seurat_reference, reduction = "wnn.umap", label = TRUE)
umap_conm | umap_with_reference | umap_reference

## Annotatew the cells based on all information
final_clusters <- rep(NA, times = ncol(seurat_ConM))
names(final_clusters) <- colnames(seurat_ConM)
naive_cells <- rownames(seurat_ConM@meta.data[seurat_ConM@meta.data$seurat_clusters == "2",])
final_clusters[naive_cells] <- "Naive cells"
atypical_cells <- rownames(seurat_ConM@meta.data[seurat_ConM@meta.data$reference_analogy == "17" &
                                                   seurat_ConM@meta.data$seurat_clusters == "6",])
final_clusters[atypical_cells] <- "Atypical cells"
plasmablasts <- rownames(seurat_ConM@meta.data[seurat_ConM@meta.data$seurat_clusters == "7",])
final_clusters[plasmablasts] <- "Plasmablasts"
act_b_cells <- rownames(seurat_ConM@meta.data[seurat_ConM@meta.data$reference_analogy == "10" &
                                                seurat_ConM@meta.data$seurat_clusters == "1",])
final_clusters[act_b_cells] <- "Activated B cells"
act_mem_b_cells <- rownames(seurat_ConM@meta.data[seurat_ConM@meta.data$reference_analogy == "3",])
final_clusters[act_mem_b_cells] <- "Activated B cells"
cd24_mem <- rownames(seurat_ConM@meta.data[seurat_ConM@meta.data$reference_analogy == "6",])
final_clusters[cd24_mem] <- "Resting Memory B cells"
cd73_mem <- rownames(seurat_ConM@meta.data[seurat_ConM@meta.data$reference_analogy == "7",])
final_clusters[cd73_mem] <- "Resting Memory B cells"
igm_mem1 <- rownames(seurat_ConM@meta.data[seurat_ConM@meta.data$reference_analogy == "8",])
final_clusters[igm_mem1] <- "IgM Memory B cells"
igm_mem2 <- rownames(seurat_ConM@meta.data[seurat_ConM@meta.data$reference_analogy == "9",])
final_clusters[igm_mem2] <- "IgM Memory B cells"
seurat_ConM <- AddMetaData(seurat_ConM, metadata = final_clusters, col.name = "final_clusters")

## Save the object
saveRDS(seurat_ConM, file = "seurat_ConM_finetune_clustering_NAs.RDS")saveRDS(query_seurat, file = "Analysis/06_Reference mapping/query_seurat_finetune_clustering_NAs.RDS")

## Remove NA cells
valid_cells <- rownames(seurat_ConM@meta.data)[!is.na(seurat_ConM@meta.data$final_clusters)]
seurat_ConM <- subset(seurat_ConM, cells = valid_cells)

## Save the object
saveRDS(query_seurat, file = "Analysis/06_Reference mapping/ConM_seurat_finetune_cluster.RDS")

## Produce the interpretation plots
fine_tuned_clustering_umap_plot <- DimPlot(seurat_ConM, group.by = "final_clusters") +
  ggtitle("Clustering after Fine Tuning"); fine_tuned_clustering_umap_plot

## Produce the PDF
pdf(file = "Analysis/06_Reference mapping/ConM_Mapping_to_Reference.pdf", width = 10, height = 6)
print(umap_query | umap_with_reference | umap_reference)
print(fine_tuned_clustering_umap_plot)
dev.off()
