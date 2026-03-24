### Author Qamar Feddah
### LRT tests

source("00_setup.R")

count_matrix <- readRDS("outputs/count_matrix.rds")
colData <- readRDS("outputs/colData.rds")

## Compare the reduced vs the full model
dds_full <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~ Visit + Sex)
dds_full <- DESeq(dds_full)

dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ Visit)
res_lrt <- results(dds_lrt, alpha = 0.1)
res_df <- as.data.frame(res_lrt)
res_df$gene <- rownames(res_df)
summary(res_lrt)

write.csv(res_df, "outputs/tables/LRT_Visit_plus_Sex_reduced_Visit.csv")

## Visualize with plot
mod_comp_plot <-
  ggplot(res_df, aes(x = log10(baseMean + 1), 
                     y = -log10(padj))) +
  geom_point(aes(color = padj < 0.1),
             alpha = 0.6, size = 1.8, show.legend = FALSE) +
  scale_color_manual(values = c("grey70", viridis(1, option = "A"))) +
  geom_hline(yintercept = -log10(0.1), 
             linetype = "dashed", linewidth = 0.5, color = "red") +
  geom_text_repel(
    data = res_df %>% filter(!is.na(padj), padj < 0.1),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "LRT for Effect of Sex: No Global Contribution Detected",
    subtitle = "Only 7 genes show padj < 0.1 across 17033 tested",
    x = "log10(mean expression + 1)",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(size = 0.2),
    panel.grid.major.y = element_line(size = 0.2)); mod_comp_plot

pdf("outputs/plots/LRT_Sex_Visit_plot.pdf", width = 10, height = 7)
print(mod_comp_plot)
dev.off()
