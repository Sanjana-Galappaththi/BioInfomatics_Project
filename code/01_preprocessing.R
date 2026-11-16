# 01 - Data import and preprocessing
# Run from the project root folder.
suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
})

dir.create("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/results/figures", recursive = TRUE, showWarnings = FALSE)

data("airway")
# make untreated the reference level
airway$dex <- relevel(airway$dex, ref = "untrt")

dds <- DESeqDataSet(airway, design = ~ cell + dex)

cat("Sample information:\n"); print(colData(dds)[, c("SampleName", "cell", "dex")])
cat("\nFirst rows of count matrix:\n"); print(head(counts(dds)))
cat("\nGenes before filtering:", nrow(dds), "\n")

# remove genes with very low counts (>= 10 reads in at least 2 samples)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")
cat("Missing values in counts:", sum(is.na(counts(dds))), "\n")
print(table(dds$dex))

# save a summary table and the filtered object for the next scripts
write.csv(as.data.frame(colData(dds)[, c("SampleName", "cell", "dex")]),
          "results/sample_information.csv")
saveRDS(dds, "results/dds_filtered.rds")
