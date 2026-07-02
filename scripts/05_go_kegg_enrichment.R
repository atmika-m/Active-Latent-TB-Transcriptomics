############################################################
# GO and KEGG Enrichment Analysis
# Project: Active vs Latent Tuberculosis Transcriptomics
# Author: Mrunalline Atmika
############################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# Load significant annotated DEGs
sig_merged <- read.csv("results/annotated/Significant_Annotated_Merged_DEGs.csv")

# Extract unique gene symbols
gene_symbols <- unique(sig_merged$GeneSymbol)

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

write.csv(gene_entrez, "results/enrichment/GeneSymbol_to_EntrezID.csv", row.names = FALSE)

# GO Biological Process enrichment
go_bp <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(as.data.frame(go_bp), "results/enrichment/GO_BP_Enrichment.csv", row.names = FALSE)

go_plot <- dotplot(go_bp, showCategory = 20) +
  ggtitle("GO Biological Process Enrichment")

ggsave("figures/GO_BP_Dotplot.png", go_plot, width = 10, height = 7)

# KEGG pathway enrichment
kegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

write.csv(as.data.frame(kegg), "results/enrichment/KEGG_Enrichment.csv", row.names = FALSE)

kegg_plot <- dotplot(kegg, showCategory = 20) +
  ggtitle("KEGG Pathway Enrichment")

ggsave("figures/KEGG_Dotplot.png", kegg_plot, width = 10, height = 7)

cat("GO and KEGG enrichment analysis complete.\n")
