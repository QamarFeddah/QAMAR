### Author Qamar Feddah
### DESeq Randomisation

source("00_setup.R")

count_matrix <- readRDS("outputs/count_matrix.rds")
colData <- readRDS("outputs/colData.rds")

dir.create("outputs/DESeq_Randomisation/")
## Create a DESeq Data Set
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ Randomisation)

## Remove genes that have counts less than 10 in total and set the base condition
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes,]
dds$Randomisation <- factor(dds$Randomisation)
dds$Randomisation <- relevel(dds$Randomisation, ref = "1")

## Run the DESeq2 workflow
dds <- DESeq(dds)
res <- results(dds, alpha = 0.1)
summary(res)
subset(res, padj < 0.1 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)) %>%
  as.data.frame() %>%
  write.csv(file.path("outputs/DESeq_Randomisation/UP_Down_regulated_Genes.csv"))
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
                  title = "Gene expression associated with Randomisation",
                  subtitle = "baseline group 1",
                  pCutoff = 0.1,
                  FCcutoff = 0.5,
                  pointSize = 2.0,
                  labSize = 5.0)

pdf(file = "outputs/DESeq_Randomisation/ConM_Randomisation_DEG_Volcano.pdf", width = 12, height = 9)
print(vlcn_plot)
dev.off()
