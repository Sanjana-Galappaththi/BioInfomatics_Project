install.packages("dplyr")

# 03 - Differential gene expression with DESeq2
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

dds <- readRDS("results/dds_filtered.rds")
dds <- DESeq(dds)
res <- results(dds, contrast = c("dex", "trt", "untrt"))

cat("\nDESeq2 summary (alpha = 0.1 default):\n"); summary(res)

results_df <- as.data.frame(res)
results_df$gene_id <- rownames(results_df)

# add gene symbols (Ensembl -> symbol)
results_df$symbol <- mapIds(org.Hs.eg.db,
                            keys = results_df$gene_id,
                            keytype = "ENSEMBL", column = "SYMBOL",
                            multiVals = "first")

results_df <- results_df[!is.na(results_df$padj), ]

results_df$Significance <- "Not Significant"
results_df$Significance[results_df$padj < 0.05 & results_df$log2FoldChange > 1]  <- "Upregulated"
results_df$Significance[results_df$padj < 0.05 & results_df$log2FoldChange < -1] <- "Downregulated"

significant_genes <- results_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)

cat("\nSignificant genes (padj < 0.05, |log2FC| > 1):", nrow(significant_genes), "\n")
print(table(results_df$Significance))
print(head(significant_genes, 10))

write.csv(results_df, "results/differential_expression.csv", row.names = FALSE)
write.csv(significant_genes, "results/significant_genes.csv", row.names = FALSE)
saveRDS(dds, "results/dds_deseq.rds")
saveRDS(res, "results/res.rds")
