### Author Qamar Feddah
### LRT tests: Ig effects adjusted for Visit

library(DESeq2)
library(tidyverse)
library(viridis)
library(ggrepel)
set.seed(333)

count_matrix <- readRDS("outputs/count_matrix.rds")
colData <- readRDS("outputs/colData.rds")

dir.create("outputs/LRT_Ig_Visit", showWarnings = FALSE, recursive = TRUE)

## Make sure Visit is a factor
colData$Visit <- as.factor(colData$Visit)

## List of Ig variables
ig_vars <- c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA")

for (ig in ig_vars) {
  
  message("Running LRT for ", ig)
  
  ## Convert Ig to numeric (comma decimal safe)
  colData[[ig]] <- as.numeric(gsub(",", ".", as.character(colData[[ig]])))
  
  ## Remove samples with missing Ig
  keep_samples <- !is.na(colData[[ig]])
  cm <- count_matrix[, keep_samples]
  cd <- colData[keep_samples, ]
  
  ## LRT
  dds_full <- DESeqDataSetFromMatrix(
    countData = cm,
    colData = cd,
    design = as.formula(paste("~ Visit +", ig))
  )
  
  dds_full <- DESeq(dds_full)
  
  dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ Visit)
  res_lrt <- results(dds_lrt, alpha = 0.1)
  
  res_df <- as.data.frame(res_lrt)
  res_df$gene <- rownames(res_df)
  
  write.csv(
    res_df,
    file = paste0("outputs/LRT_Ig_Visit/LRT_", ig, "_adj_Visit.csv")
  )
  
  ## Plot
  p <-
    ggplot(res_df, aes(log10(baseMean + 1), -log10(padj))) +
    geom_point(aes(color = padj < 0.1),
               size = 1.8, alpha = 0.6, show.legend = FALSE) +
    scale_color_manual(values = c("grey70", viridis(1))) +
    geom_hline(yintercept = -log10(0.1),
               linetype = "dashed", linewidth = 0.5) +
    labs(
      title = paste("LRT for", ig, "adjusted for Visit"),
      subtitle = paste(
        "Genes with FDR < 0.1:",
        sum(!is.na(res_df$padj) & res_df$padj < 0.1)
      ),
      x = "log10(mean expression + 1)",
      y = "-log10 adjusted p-value"
    ) +
    geom_text_repel(
      data = res_df %>% filter(!is.na(padj), padj < 0.1),
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    theme_minimal(base_size = 14)
  
  pdf(
    file = paste0("outputs/LRT_Ig_Visit/LRT_", ig, "_adj_Visit.pdf"),
    width = 10,
    height = 7
  )
  print(p)
  dev.off()
}
