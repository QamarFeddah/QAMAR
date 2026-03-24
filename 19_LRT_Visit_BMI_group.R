### LRT tests

source("00_setup.R")

count_matrix <- readRDS("outputs/count_matrix.rds")
colData <- readRDS("outputs/colData.rds")
dir.create("outputs/LRT_Visit_BMI")

## Compare the reduced vs the full model
dds_full <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ Visit + BMI_group)
dds_full <- DESeq(dds_full)

dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ Visit)
res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)

write.csv(res_df, "outputs/LRT_Visit_BMI/LRT_Visit_plus_BMI_group_reduced_Visit.csv")

########
## Compare the reduced vs the full model
dds_full <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ Visit + BMI_group)
dds_full <- DESeq(dds_full)

dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~BMI_group)
res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)

write.csv(res_df, "outputs/LRT_Visit_BMI/LRT_Visit_plus_BMI_group_reduced_BMI_group.csv")