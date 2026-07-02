############################################################
# Differential Expression Analysis: GSE25534
# Project: Active vs Latent Tuberculosis Transcriptomics
# Author: Mrunalline Atmika
############################################################

# Load required libraries
library(GEOquery)
library(limma)
library(tidyverse)

# Create output directory
dir.create("results/GSE25534", recursive = TRUE, showWarnings = FALSE)

# Download GEO dataset
gse <- getGEO("GSE25534", GSEMatrix = TRUE)

# Extract expression matrix and metadata
expr <- exprs(gse[[1]])
metadata <- pData(gse[[1]])

# Assign disease groups
metadata$group <- NA
metadata$group[grepl("TB", metadata$source_name_ch1, ignore.case = TRUE)] <- "ActiveTB"
metadata$group[grepl("LTBI", metadata$source_name_ch1, ignore.case = TRUE)] <- "LTBI"

# Keep only Active TB and LTBI samples
metadata_clean <- metadata %>%
  filter(group %in% c("ActiveTB", "LTBI"))

expr_clean <- expr[, rownames(metadata_clean)]

# Create design matrix
group <- factor(metadata_clean$group, levels = c("LTBI", "ActiveTB"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Differential expression analysis using limma
fit <- lmFit(expr_clean, design)
contrast <- makeContrasts(ActiveTB - LTBI, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Extract DEG table
deg <- topTable(fit2, number = Inf, adjust.method = "BH")

# Save outputs
write.csv(deg, "results/GSE25534/DEG_table.csv")
write.csv(expr_clean, "results/GSE25534/expression_matrix.csv")
write.csv(metadata_clean, "results/GSE25534/metadata_clean.csv")

# Summary
cat("GSE25534 analysis complete.\n")
cat("Group counts:\n")
print(table(metadata_clean$group))
