### Author Qamar Feddah
### DESeq for IgG1

source("00_setup.R")

library(ggrepel)
library(EnhancedVolcano)
set.seed(333)

count_matrix <- readRDS("outputs/count_matrix.rds")
colData <- readRDS("outputs/colData.rds")

dir.create("outputs/DESeq_for_IgG1", showWarnings = FALSE, recursive = TRUE)

## Make sure variables exist and are in correct type
colData$IgG1 <- as.numeric(gsub(",", ".", as.character(colData$IgG1)))


############################################################
## 2) DESeq: association with IgG1 (treat IgG1 as continuous)
############################################################

## Normalize IgG1
colData <-  colData %>%
  mutate(NormIgG1 = scale(log10(IgG1+1)))

## Create a DESeq Data Set
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ NormIgG1)

## Remove genes that have counts less than 10 in total 
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes,]

## Run the DESeq2 workflow
dds <- DESeq(dds)
res <- results(dds, alpha = 0.1)
summary(res)
subset(res, padj < 0.1 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) %>%
  as.data.frame() %>%
  write.csv(file.path("outputs/DESeq_for_IgG1/UP_Down_regulated_Genes.csv"))
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
                  title = "Gene expression associated with IgG1",
                  pCutoff = 0.1,
                  FCcutoff = 0.5,
                  pointSize = 2.0,
                  labSize = 5.0)

pdf(file = "outputs/DESeq_for_IgG1/ConM_IgG1_DEG_Volcano.pdf", width = 12, height = 9)
print(vlcn_plot)
dev.off()
