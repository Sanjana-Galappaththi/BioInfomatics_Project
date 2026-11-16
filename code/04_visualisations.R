# 04 - Volcano plot, MA plot and heatmap of top genes
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

results_df <- read.csv("results/differential_expression.csv")
vsd <- readRDS("results/vsd.rds")
dds <- readRDS("results/dds_deseq.rds")
res <- readRDS("results/res.rds")

# Volcano plot
p <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Downregulated" = "blue",
                                "Not Significant" = "grey70",
                                "Upregulated" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "Log2 Fold Change (treated vs untreated)",
       y = "-Log10 Adjusted P-value") +
  theme_minimal()
ggsave("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/results/figures/04_volcano_plot.png", p, width = 8, height = 6, dpi = 150)

# MA plot
png("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/results/figures/05_ma_plot.png", width = 900, height = 700, res = 120)
plotMA(res, ylim = c(-5, 5), main = "MA plot: dexamethasone vs untreated")
dev.off()

# Heatmap of the 30 genes with smallest adjusted p-values
top_idx <- head(order(res$padj, na.last = NA), 30)
mat <- assay(vsd)[top_idx, ]
labels <- results_df$symbol[match(rownames(mat), results_df$gene_id)]
labels[is.na(labels)] <- rownames(mat)[is.na(labels)]
rownames(mat) <- make.unique(labels)

png("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/results/figures/06_heatmap_top30.png", width = 1000, height = 900, res = 120)
pheatmap(mat, scale = "row",
         annotation_col = as.data.frame(colData(dds)[, c("dex", "cell")]),
         main = "Top 30 Differentially Expressed Genes")
dev.off()
