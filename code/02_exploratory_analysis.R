#install.packages("pheatmap")

# 02 - Exploratory data analysis
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

dds <- readRDS("results/dds_filtered.rds")
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, "results/vsd.rds")

# Box plot of normalised expression per sample
png("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/results/figures/01_boxplot_expression.png", width = 1000, height = 700, res = 120)
boxplot(as.data.frame(assay(vsd)),
        main = "Gene Expression Distribution",
        xlab = "Samples", ylab = "Normalised Expression (VST)",
        las = 2, col = "lightblue")
dev.off()

# PCA
pca_data <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = dex, shape = cell)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Airway Gene Expression") +
  theme_minimal()
ggsave("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/results/figures/02_pca.png", p, width = 7, height = 5, dpi = 150)
write.csv(pca_data, "results/pca_coordinates.csv", row.names = TRUE)

# Sample-to-sample distance heatmap
library(pheatmap)
d <- as.matrix(dist(t(assay(vsd))))
rownames(d) <- paste(vsd$dex, vsd$cell, sep = "-")
colnames(d) <- NULL
png("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/results/figures/03_sample_distance_heatmap.png", width = 900, height = 800, res = 120)
pheatmap(d, main = "Sample-to-sample distances")
dev.off()
