### Author Qamar Feddah
### IgG1 LRT (Sex effect adjusted for IgG1) + DESeq on IgG1 + plots

source("00_setup.R")

library(ggrepel)
library(EnhancedVolcano)
set.seed(333)

count_matrix <- readRDS("outputs/count_matrix.rds")
colData <- readRDS("outputs/colData.rds")

dir.create("outputs/LRT_for_Sex_IgG1", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/LRT_for_Sex_IgG1/plots", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/LRT_for_Sex_IgG1/tables", showWarnings = FALSE, recursive = TRUE)

## Make sure variables exist and are in correct type
colData$Sex <- as.factor(colData$Sex)
colData$IgG1 <- as.numeric(gsub(",", ".", as.character(colData$IgG1)))

############################################################
## 1) 
## Full model: ~ IgG1 + Sex
## Reduced model: ~ IgG1
############################################################

dds_full <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ IgG1 + Sex)
dds_full <- DESeq(dds_full)

dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ Sex)
res_lrt <- results(dds_lrt, alpha = 0.1)
summary(res_lrt)

res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)

write.csv(res_df, "outputs/LRT_for_Sex_IgG1/tables/LRT_Sex_adj_IgG1.csv")

mod_comp_plot_IgG1 <-
  ggplot(res_df, aes(x = log10(baseMean + 1),
                     y = -log10(padj))) +
  geom_point(aes(color = padj < 0.1),
             alpha = 0.6, size = 1.8, show.legend = FALSE) +
  scale_color_manual(values = c("grey70", viridis(1, option = "A"))) +
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed", linewidth = 0.5, color = "red") +
  labs(
    title = "LRT for Effect of Sex adjusted for IgG1",
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
  )

pdf("outputs/LRT_for_Sex_IgG1/plots/LRT_Sex_adj_IgG1.pdf", width = 10, height = 7)
print(mod_comp_plot_IgG1)
dev.off()

