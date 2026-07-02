############################################################
# Annotation, Merging and Common DEG Identification
# Project: Active vs Latent Tuberculosis Transcriptomics
# Author: Mrunalline Atmika
############################################################

library(GEOquery)
library(dplyr)
library(readr)
library(VennDiagram)

dir.create("results/annotated", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# Load DEG tables
deg1 <- read.csv("results/GSE19491/DEG_table.csv")
deg3 <- read.csv("results/GSE28623/DEG_table.csv")
deg4 <- read.csv("results/GSE25534/DEG_table.csv")

# Download platform annotation files
gpl1 <- Table(getGEO("GPL6947"))
gpl3 <- Table(getGEO("GPL4133"))
gpl4 <- Table(getGEO("GPL1708"))

# Extract probe-to-gene annotations
ann1 <- gpl1 %>%
  select(ID, GeneSymbol = Symbol)

ann3 <- gpl3 %>%
  select(ID, GeneSymbol = GENE_SYMBOL)

ann4 <- gpl4 %>%
  select(ID, GeneSymbol = GENE_SYMBOL)

# Annotate DEG tables
deg1_annotated <- deg1 %>%
  left_join(ann1, by = c("X" = "ID")) %>%
  mutate(GSE = "GSE19491", GPL = "GPL6947")

deg3_annotated <- deg3 %>%
  left_join(ann3, by = c("X" = "ID")) %>%
  mutate(GSE = "GSE28623", GPL = "GPL4133")

deg4_annotated <- deg4 %>%
  left_join(ann4, by = c("X" = "ID")) %>%
  mutate(GSE = "GSE25534", GPL = "GPL1708")

# Save annotated DEG tables
write.csv(deg1_annotated, "results/annotated/Annotated_DEG_GSE19491.csv", row.names = FALSE)
write.csv(deg3_annotated, "results/annotated/Annotated_DEG_GSE28623.csv", row.names = FALSE)
write.csv(deg4_annotated, "results/annotated/Annotated_DEG_GSE25534.csv", row.names = FALSE)

# Merge annotated DEG tables
merged_deg <- bind_rows(
  deg1_annotated,
  deg3_annotated,
  deg4_annotated
)

# Clean missing gene symbols
merged_deg_clean <- merged_deg %>%
  filter(!is.na(GeneSymbol), trimws(GeneSymbol) != "")

write.csv(merged_deg_clean, "results/annotated/Merged_Annotated_DEGs_Cleaned.csv", row.names = FALSE)

# Identify significant DEGs
sig_merged <- merged_deg_clean %>%
  mutate(
    logFC = as.numeric(logFC),
    adj.P.Val = as.numeric(adj.P.Val)
  ) %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)

write.csv(sig_merged, "results/annotated/Significant_Annotated_Merged_DEGs.csv", row.names = FALSE)

# Separate upregulated and downregulated genes
upregulated_genes <- sig_merged %>%
  filter(logFC > 1)

downregulated_genes <- sig_merged %>%
  filter(logFC < -1)

write.csv(upregulated_genes, "results/annotated/Upregulated_DEGs.csv", row.names = FALSE)
write.csv(downregulated_genes, "results/annotated/Downregulated_DEGs.csv", row.names = FALSE)

# Extract significant genes from each dataset
genes_GSE19491 <- sig_merged %>%
  filter(GSE == "GSE19491") %>%
  pull(GeneSymbol) %>%
  unique()

genes_GSE28623 <- sig_merged %>%
  filter(GSE == "GSE28623") %>%
  pull(GeneSymbol) %>%
  unique()

genes_GSE25534 <- sig_merged %>%
  filter(GSE == "GSE25534") %>%
  pull(GeneSymbol) %>%
  unique()

# Identify common DEGs across all datasets
common_genes <- Reduce(intersect, list(
  genes_GSE19491,
  genes_GSE28623,
  genes_GSE25534
))

write.csv(
  data.frame(GeneSymbol = common_genes),
  "results/annotated/Common_DEGs_All_Datasets.csv",
  row.names = FALSE
)

# Generate Venn diagram
venn.diagram(
  x = list(
    GSE19491 = genes_GSE19491,
    GSE28623 = genes_GSE28623,
    GSE25534 = genes_GSE25534
  ),
  filename = "figures/DEG_Venn_Diagram.png",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

cat("Annotation, merging and common DEG analysis complete.\n")
cat("Number of common genes:", length(common_genes), "\n")
