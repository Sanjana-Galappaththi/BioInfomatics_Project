# Run once: installs all required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("airway", "DESeq2", "org.Hs.eg.db", "AnnotationDbi"))
install.packages(c("ggplot2", "dplyr", "pheatmap", "rmarkdown", "knitr"))
