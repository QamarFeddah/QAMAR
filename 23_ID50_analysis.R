### Author Qamar Feddah
### ID50 related analysis in COnM dataset

## Library loading and set the working directory
library(Seurat)
library(tidyverse)
library(DESeq2)
library(viridis)
library(EnhancedVolcano)
library(ggrepel)
setwd("~/ConM_project/Analysis/07_confounder_assessment+DESeq2_ConM")
set.seed(333)


dir.create("outputs/ID50_results_in_ConM")

#load seurat object
seurat_conM <- readRDS("C:/Users/Qamar/Documents/ConM_project/Analysis/07_confounder_assessment+DESeq2_ConM/outputs/seurat_conM.rds")

sample_meta_data <- read.csv("C:/Users/Qamar/Documents/ConM_project/Raw_data/Metadata.csv" , sep = ";", dec = ",")
#("Metadata.csv", sep = ";", dec = ",")

### Clean ID50
sample_meta_data$ID50 <- as.character(sample_meta_data$ID50)

sample_meta_data$ID50[sample_meta_data$ID50 == ""] <- NA

sample_meta_data$ID50_num <- ifelse(
  sample_meta_data$ID50 == "<20",
  20,
  as.numeric(sample_meta_data$ID50)
)

sample_meta_data$log10_ID50 <- log10(sample_meta_data$ID50_num)

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

## Create the column data frame
colData <- data.frame(Sample = colnames(count_matrix))
row.names(colData) <- NULL

sample_meta_data$Donor_Visit <- paste(sample_meta_data$Donor, sample_meta_data$Visit, sep = "_")

colData <-
  left_join(colData, sample_meta_data, by = c("Sample" = "Donor_Visit")) %>%
  as.data.frame() %>%
  distinct(Sample, .keep_all = TRUE) %>%
  select(Sample,Visit, Sex, Randomisation, Age_group, BMI_group,IgG, IgG1, IgG2, IgG3, IgG4, IgA, ID50, ID50_num, log10_ID50) %>%
  column_to_rownames(var = "Sample")
colData <- colData[colnames(count_matrix), , drop = FALSE]

## Build the DESeq object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ 1)
dds <- estimateSizeFactors(dds) ##Normalize for sequencing depth and Stabilize variance
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)  # blind = TRUE for exploratory PCA
pca_data <- prcomp(t(assay(vsd)))  # transpose so PCA is on samples

## make sure variable types are usable
colData$log10_ID50 <- as.numeric(as.character(colData$log10_ID50))
colData$Sex <- as.factor(colData$Sex)
colData$Visit <- as.factor(colData$Visit)
colData$Randomisation <- as.factor(colData$Randomisation)
colData$Age_group <- as.factor(colData$Age_group)
colData$BMI_group <- as.factor(colData$BMI_group)
## Convert immunoglobulins to numeric and create normalized versions
colData$IgG  <- as.numeric(gsub(",", ".", as.character(colData$IgG)))
colData$IgG1 <- as.numeric(gsub(",", ".", as.character(colData$IgG1)))
colData$IgG2 <- as.numeric(gsub(",", ".", as.character(colData$IgG2)))
colData$IgG3 <- as.numeric(gsub(",", ".", as.character(colData$IgG3)))
colData$IgG4 <- as.numeric(gsub(",", ".", as.character(colData$IgG4)))
colData$IgA  <- as.numeric(gsub(",", ".", as.character(colData$IgA)))
colData$NormIgG  <- as.numeric(scale(log10(colData$IgG  + 1)))
colData$NormIgG1 <- as.numeric(scale(log10(colData$IgG1 + 1)))
colData$NormIgG2 <- as.numeric(scale(log10(colData$IgG2 + 1)))
colData$NormIgG3 <- as.numeric(scale(log10(colData$IgG3 + 1)))
colData$NormIgG4 <- as.numeric(scale(log10(colData$IgG4 + 1)))
colData$NormIgA  <- as.numeric(scale(log10(colData$IgA  + 1)))



############################################################
## Sex given log10_ID50
############################################################

keep_samples <- complete.cases(colData[, c("log10_ID50", "Sex")])
colData_sex <- colData[keep_samples, , drop = FALSE]
count_matrix_sex <- count_matrix[, rownames(colData_sex), drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(countData = count_matrix_sex, colData = colData_sex,
                                   design = ~ Sex + log10_ID50)

dds_full <- DESeq(dds_full)
dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ log10_ID50)

res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)

write.csv(
  res_df,
  file.path("outputs/ID50_results_in_ConM", "LRT_Sex_given_log10_ID50.csv"),
  row.names = FALSE
)
mod_comp_plot_log10_ID50 <-
  ggplot(res_df, aes(x = log10(baseMean + 1),
                     y = -log10(padj))) +
  geom_point(aes(color = padj < 0.1),
             alpha = 0.6, size = 1.8, show.legend = FALSE) +
  scale_color_manual(values = c("grey70", viridis(1, option = "A"))) +
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed", linewidth = 0.5, color = "red") +
  labs(
    title = "LRT for Effect of Sex given log10_ID50 in ACT",
    subtitle = paste0("Genes with FDR < 0.1: ", sum(!is.na(res_df$padj) & res_df$padj < 0.1)),
    x = "log10(mean expression + 1)",
    y = "-log10 adjusted p-value"
  ) +
  geom_text_repel(
    data = res_df %>% filter(!is.na(padj), padj < 0.1),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2),
    panel.grid.major.y = element_line(linewidth = 0.2)
  ) ; mod_comp_plot_log10_ID50

############################################################
## Randomisation given log10_ID50
############################################################

keep_samples <- complete.cases(colData[, c("log10_ID50", "Randomisation")])
colData_rand <- colData[keep_samples, , drop = FALSE]
count_matrix_rand <- count_matrix[, rownames(colData_rand), drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(countData = count_matrix_rand, colData = colData_rand,
                                   design = ~ Randomisation + log10_ID50)

dds_full <- DESeq(dds_full)
dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ log10_ID50)

res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)

write.csv(
  res_df,
  file.path("outputs/ID50_results_in_ConM", "LRT_Randomisation_given_log10_ID50_all.csv"),
  row.names = FALSE
)


############################################################
## Age_group given log10_ID50
############################################################

keep_samples <- complete.cases(colData[, c("log10_ID50", "Age_group")])
colData_age <- colData[keep_samples, , drop = FALSE]
count_matrix_age <- count_matrix[, rownames(colData_age), drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(countData = count_matrix_age, colData  = colData_age,
                                   design    = ~ Age_group + log10_ID50)

dds_full <- DESeq(dds_full)
dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ log10_ID50)

res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)
write.csv(
  res_df,
  file.path("outputs/ID50_results_in_ConM", "LRT_Age_group_given_log10_ID50.csv"),
  row.names = FALSE
)

############################################################
## BMI_group given log10_ID50
############################################################

keep_samples <- complete.cases(colData[, c("log10_ID50", "BMI_group")])
colData_bmi <- colData[keep_samples, , drop = FALSE]
count_matrix_bmi <- count_matrix[, rownames(colData_bmi), drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(countData = count_matrix_bmi, colData = colData_bmi,
                                   design  = ~ BMI_group + log10_ID50)

dds_full <- DESeq(dds_full)
dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ log10_ID50)

res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)
write.csv(
  res_df,
  file.path("outputs/ID50_results_in_ConM/", "LRT_BMI_group_given_log10_ID50.csv"),
  row.names = FALSE
)

############################################################
## log10_ID50 given Sex
############################################################
keep_samples <- complete.cases(colData[, c("log10_ID50", "Sex")])
colData_sex <- colData[keep_samples, , drop = FALSE]
count_matrix_sex <- count_matrix[, rownames(colData_sex), drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(countData = count_matrix_sex, colData = colData_sex,
                                   design = ~ Sex + log10_ID50)

dds_full <- DESeq(dds_full)
dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ Sex)

res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)

write.csv(
  res_df,
  file.path("outputs/ID50_results_in_ConM", "LRT_log10_ID50_given_Sex.csv"),
  row.names = FALSE
)

############################################################
## log10_ID50 given Age_group
############################################################
keep_samples <- complete.cases(colData[, c("log10_ID50", "Age_group")])
colData_Age_group <- colData[keep_samples, , drop = FALSE]
count_matrix_Age_group <- count_matrix[, rownames(colData_Age_group), drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(countData = count_matrix_Age_group, colData = colData_Age_group,
                                   design = ~ Age_group + log10_ID50)

dds_full <- DESeq(dds_full)
dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ Age_group)

res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)

write.csv(
  res_df,
  file.path("outputs/ID50_results_in_ConM", "LRT_log10_ID50_given_Age_group.csv"),
  row.names = FALSE
)
############################################################
## log10_ID50 given BMI_group
############################################################
keep_samples <- complete.cases(colData[, c("log10_ID50", "BMI_group")])
colData_BMI_group <- colData[keep_samples, , drop = FALSE]
count_matrix_BMI_group <- count_matrix[, rownames(colData_BMI_group), drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(countData = count_matrix_BMI_group, colData = colData_BMI_group,
                                   design = ~ BMI_group + log10_ID50)

dds_full <- DESeq(dds_full)
dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ BMI_group)

res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)

write.csv(
  res_df,
  file.path("outputs/ID50_results_in_ConM", "LRT_log10_ID50_given_BMI_group.csv"),
  row.names = FALSE
)
############################################################
## log10_ID50 given Randomisation
############################################################
keep_samples <- complete.cases(colData[, c("log10_ID50", "Randomisation")])
colData_Randomisation <- colData[keep_samples, , drop = FALSE]
count_matrix_Randomisation <- count_matrix[, rownames(colData_Randomisation), drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(countData = count_matrix_Randomisation, colData = colData_Randomisation,
                                   design = ~ Randomisation + log10_ID50)

dds_full <- DESeq(dds_full)
dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ Randomisation)

res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)

write.csv(
  res_df,
  file.path("outputs/ID50_results_in_ConM", "LRT_log10_ID50_given_Randomisation.csv"),
  row.names = FALSE
)

############################################################
## DESeq2 test: gene expression associated with log10(ID50)
## Model: ~ log10_ID50
############################################################
keep_samples <- complete.cases(colData[, "log10_ID50"])
colData <- colData[keep_samples, , drop = FALSE]
count_matrix <- count_matrix[, rownames(colData), drop = FALSE]

## Create a DESeq Data Set
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ log10_ID50)

## Remove genes that have counts less than 10 in total 
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes, ]

## Run the DESeq2 workflow
dds <- DESeq(dds)

res <- results(dds, alpha = 0.1)
summary(res)

subset(res, padj < 0.1 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) %>%
  as.data.frame() %>%
  write.csv(file.path("outputs/ID50_results_in_ConM/UP_Down_regulated_Genes.csv"))
res_dataframe <- as.data.frame(res)

## Create the ranked vector of genes
res_dataframe <- res_dataframe[order(res_dataframe$stat, decreasing = TRUE),]
genes <- row.names(res_dataframe)
ranked_genes <- res_dataframe$stat
names(ranked_genes) <- genes


vlcn_plot <-
  EnhancedVolcano(res_dataframe,
                  lab = rownames(res_dataframe),
                  x = "log2FoldChange",
                  y = "padj",
                  title = "Gene expression associated with log10_ID50 in The ConM dataset",
                  pCutoff = 0.1,
                  FCcutoff = 0.5,
                  pointSize = 2.0,
                  labSize = 5.0)

pdf(file = "outputs/ID50_results_in_ConM/ID50_DEG_Volcano_ACT.pdf", width = 12, height = 9)
print(vlcn_plot)
dev.off()
