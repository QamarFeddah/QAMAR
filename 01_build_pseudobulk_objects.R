### Author Qamar Feddah
### Build pseudobulk objects

source("00_setup.R")

seurat_conM <- readRDS("C:/Users/Qamar/Documents/ConM_project/Analysis/06_Reference mapping/seurat_ConM_finetune_cluster.RDS")
conM_meta <- seurat_conM@meta.data %>% as.data.frame()

sample_meta_data <- read.csv("C:/Users/Qamar/Documents/ConM_project/Raw_data/Metadata.csv" , sep = ";", dec = ",")
#("Metadata.csv", sep = ";", dec = ",")

names(conM_meta)[names(conM_meta) == "orig.ident"] <- "Lane"
names(conM_meta)[names(conM_meta) == "Sample_ID"] <- "Hashtag_ID"

conM_meta <- conM_meta %>%
  mutate(Identifier = paste(conM_meta$Lane, conM_meta$Hashtag_ID, sep = "_"))

sample_meta_data$Lane <- paste("Lane", sample_meta_data$Lane, sep = "")
sample_meta_data$Hashtag <- gsub(pattern = "_", replacement = "-", sample_meta_data$Hashtag)

sample_meta_data <- sample_meta_data %>%
  mutate(Identifier = paste(sample_meta_data$Lane, sample_meta_data$Hashtag, sep = "_"))

conM_meta <- conM_meta %>% rownames_to_column(var = "Barcode")

unified_meta_data <- left_join(conM_meta, sample_meta_data, by = "Identifier") %>%
  select(-`Lane.y`, -Hashtag)

names(unified_meta_data)[names(unified_meta_data) == "Lane.x"] <- "Lane"

unified_meta_data <- unified_meta_data %>%
  column_to_rownames(var = "Barcode")

cat_cols <- c("Sex","Randomisation","Age_group","BMI_group","Visit")
for (cn in cat_cols) {
  if (cn %in% colnames(unified_meta_data)) {
    unified_meta_data[[cn]] <- as.character(unified_meta_data[[cn]])
    unified_meta_data[[cn]][is.na(unified_meta_data[[cn]])] <- "NA"
  }
}
## Pass th metadata
seurat_conM@meta.data <- unified_meta_data

## Create the Pseudo-bulk expression matrix
count_matrix <-
  AggregateExpression(
    seurat_conM,
    group.by = c("Donor","Visit"),
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )$RNA %>%
  as.matrix()

count_matrix <- count_matrix[,-47]

## Create the column data frame
colData <- data.frame(Sample = colnames(count_matrix))
row.names(colData) <- NULL

sample_meta_data$Donor_Visit <- paste(sample_meta_data$Donor, sample_meta_data$Visit, sep = "_")

colData <-
  left_join(colData, sample_meta_data, by = c("Sample" = "Donor_Visit")) %>%
  as.data.frame() %>%
  distinct(Sample, .keep_all = TRUE) %>%
  select(Sample,Visit, Sex, Randomisation, Age_group, BMI_group,IgG, IgG1, IgG2, IgG3, IgG4, IgA, ID50) %>%
  column_to_rownames(var = "Sample")


## Build the DESeq object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ 1)
dds <- estimateSizeFactors(dds) ##Normalize for sequencing depth and Stabilize variance
vsd <- vst(dds, blind = TRUE)  # blind = TRUE for exploratory PCA
pca_data <- prcomp(t(assay(vsd)))  # transpose so PCA is on samples

## Run PCA
pca_df <- as.data.frame(pca_data$x)
pca_df$Visit <- as.factor(colData$Visit)
pca_df$Sex <- as.factor(colData$Sex)
pca_df$Randomisation <- as.factor(colData$Randomisation)
pca_df$Age_group <- as.factor(colData$Age_group)
pca_df$BMI_group <- as.factor(colData$BMI_group)

pca_df$IgG  <- colData$IgG
pca_df$IgG1 <- colData$IgG1
pca_df$IgG2 <- colData$IgG2
pca_df$IgG3 <- colData$IgG3
pca_df$IgG4 <- colData$IgG4
pca_df$IgA  <- colData$IgA

saveRDS(seurat_conM, "outputs/seurat_conM.rds")
saveRDS(sample_meta_data, "outputs/sample_meta_data.rds")
saveRDS(unified_meta_data, "outputs/unified_meta_data.rds")
saveRDS(count_matrix, "outputs/count_matrix.rds")
saveRDS(colData, "outputs/colData.rds")
saveRDS(dds, "outputs/dds_pca_base.rds")
saveRDS(vsd, "outputs/vsd.rds")
saveRDS(pca_data, "outputs/pca_data.rds")
saveRDS(pca_df, "outputs/pca_df.rds")
