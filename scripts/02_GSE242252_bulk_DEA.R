#!/usr/bin/env Rscript
###############################################################################
# 02_GSE242252_bulk_DEA.R
# HZ急性期 vs 恢复期 差异表达分析 (Bulk RNA-seq, 全血)
# 包括: PCA质控, 性染色体基因检查, DESeq2 DEA
###############################################################################
suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(tidyverse)
})

set.seed(42)

# ── 0. Paths ──────────────────────────────────────────────────────────────────
DATA_DIR  <- "/media/cmj/MechanicalDisk/yjs/VZV-geo/data"
RES_DIR  <- "/media/cmj/MechanicalDisk/yjs/VZV-geo/results/GSE242252"
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load data ─────────────────────────────────────────────────────────────
message("Loading GSE242252...")
gse <- getGEO("GSE242252", destdir = DATA_DIR, GSEMatrix = TRUE,
              getGPL = FALSE, AnnotGPL = FALSE)
eset <- gse[[1]]
pd   <- pData(eset)
expr <- exprs(eset)

message("Samples: ", ncol(expr))
message("Genes: ", nrow(expr))
cat("\nColumn names in pData:\n")
cat(paste0("  ", colnames(pd)), sep = "\n")

# ── 2. Build metadata ────────────────────────────────────────────────────────
# Inspect relevant columns
cat("\nSample titles:\n")
print(pd$title)

# ── 3. PCA Quality Control ───────────────────────────────────────────────────
message("\n===== PCA QC =====")

# Full gene set PCA
expr_log <- log2(expr + 1)
pca_full <- prcomp(t(expr_log), center = TRUE, scale. = TRUE)
pca_var  <- round(summary(pca_full)$importance[2, ] * 100, 1)

# Annotate with available metadata
pca_df <- data.frame(
  PC1 = pca_full$x[, 1],
  PC2 = pca_full$x[, 2],
  PC3 = pca_full$x[, 3],
  sample = colnames(expr),
  title = pd$title,
  stringsAsFactors = FALSE
)

# Try to extract sex from metadata if available
sex_col <- grep("sex|gender", colnames(pd), ignore.case = TRUE, value = TRUE)
if (length(sex_col) > 0) {
  pca_df$sex <- pd[[sex_col[1]]]
} else {
  # Infer sex from expression of sex chromosome genes
  # XIST (female-specific), RPS4Y1 / DDX3Y (male-specific)
  xist_row <- which(rownames(expr) == "XIST")
  rps4y1_row <- which(rownames(expr) == "RPS4Y1")
  ddx3y_row <- which(rownames(expr) == "DDX3Y")
  if (length(xist_row) > 0) {
    xist_expr <- as.numeric(expr[xist_row, ])
    pca_df$XIST <- xist_expr
    message("XIST expression range: ", min(xist_expr), " - ", max(xist_expr))
  }
  if (length(rps4y1_row) > 0) {
    rps4y1_expr <- as.numeric(expr[rps4y1_row, ])
    pca_df$RPS4Y1 <- rps4y1_expr
    message("RPS4Y1 expression range: ", min(rps4y1_expr), " - ", max(rps4y1_expr))
  }
  if (length(ddx3y_row) > 0) {
    ddx3y_expr <- as.numeric(expr[ddx3y_row, ])
    pca_df$DDX3Y <- ddx3y_expr
    message("DDX3Y expression range: ", min(ddx3y_expr), " - ", max(ddx3y_expr))
  }
}

# PC1 vs PC2
pdf(file.path(RES_DIR, "PCA_full_PC1vsPC2.pdf"), width = 8, height = 6)
p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)"),
       title = "GSE242252 PCA (All Genes)") +
  theme_minimal(base_size = 14)
print(p)
dev.off()
message("PCA(all genes) plot saved.")

# PC2 vs PC3
pdf(file.path(RES_DIR, "PCA_full_PC2vsPC3.pdf"), width = 8, height = 6)
p2 <- ggplot(pca_df, aes(x = PC2, y = PC3)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = paste0("PC2 (", pca_var[2], "%)"),
       y = paste0("PC3 (", pca_var[3], "%)"),
       title = "GSE242252 PCA (All Genes) - PC2 vs PC3") +
  theme_minimal(base_size = 14)
print(p2)
dev.off()

# Print top loadings for PC1
pc1_loadings <- sort(pca_full$rotation[, 1], decreasing = TRUE)
cat("\nTop 20 PC1 positive loadings:\n")
print(head(pc1_loadings, 20))
cat("\nTop 20 PC1 negative loadings:\n")
print(tail(pc1_loadings, 20))

# ── 4. Identify sex chromosome genes ─────────────────────────────────────────
message("\n===== Sex Chromosome Gene Check =====")
library(org.Hs.eg.db)

# Get gene symbols on X and Y chromosomes
gene_ids <- keys(org.Hs.eg.db, keytype = "SYMBOL")
chr_info <- select(org.Hs.eg.db, keys = gene_ids,
                   columns = c("CHR"), keytype = "SYMBOL")
chr_x <- chr_info$SYMBOL[chr_info$CHR == "X"]
chr_y <- chr_info$SYMBOL[chr_info$CHR == "Y"]

genes_x_in_data <- intersect(rownames(expr), chr_x)
genes_y_in_data <- intersect(rownames(expr), chr_y)

message("X chromosome genes in data: ", length(genes_x_in_data))
message("Y chromosome genes in data: ", length(genes_y_in_data))

# ── 5. PCA without sex chromosome genes ──────────────────────────────────────
message("\n===== PCA (no sex chr) =====")
keep_genes <- setdiff(rownames(expr), c(genes_x_in_data, genes_y_in_data))
expr_nosex <- expr[keep_genes, ]
expr_nosex_log <- log2(expr_nosex + 1)

pca_nosex <- prcomp(t(expr_nosex_log), center = TRUE, scale. = TRUE)
pca_nosex_var <- round(summary(pca_nosex)$importance[2, ] * 100, 1)

pca_nosex_df <- data.frame(
  PC1 = pca_nosex$x[, 1],
  PC2 = pca_nosex$x[, 2],
  PC3 = pca_nosex$x[, 3],
  sample = colnames(expr),
  stringsAsFactors = FALSE
)

pdf(file.path(RES_DIR, "PCA_nosexchr_PC1vsPC2.pdf"), width = 8, height = 6)
p3 <- ggplot(pca_nosex_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = paste0("PC1 (", pca_nosex_var[1], "%)"),
       y = paste0("PC2 (", pca_nosex_var[2], "%)"),
       title = "GSE242252 PCA (No Sex Chr Genes)") +
  theme_minimal(base_size = 14)
print(p3)
dev.off()

pdf(file.path(RES_DIR, "PCA_nosexchr_PC2vsPC3.pdf"), width = 8, height = 6)
p4 <- ggplot(pca_nosex_df, aes(x = PC2, y = PC3)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = paste0("PC2 (", pca_nosex_var[2], "%)"),
       y = paste0("PC3 (", pca_nosex_var[3], "%)"),
       title = "GSE242252 PCA (No Sex Chr Genes) - PC2 vs PC3") +
  theme_minimal(base_size = 14)
print(p4)
dev.off()
message("PCA (no sex chr) plots saved.")

# ── 6. Build proper metadata for DESeq2 ──────────────────────────────────────
# From the paper: paired design, acute HZ vs 1 year post, n=26 pairs
# Need to identify which samples are acute vs convalescent
cat("\nSample titles for manual annotation:\n")
for (i in seq_along(pd$title)) {
  cat(sprintf("  [%d] %s\n", i, pd$title[i]))
}

# Write sample info for manual inspection
write.csv(pd, file.path(RES_DIR, "sample_metadata.csv"), row.names = FALSE)
message("\nSample metadata saved to sample_metadata.csv for review.")

# ── 7. Save processed data ───────────────────────────────────────────────────
saveRDS(list(expr = expr, pd = pd,
             pca_full = pca_full, pca_nosex = pca_nosex,
             genes_x = genes_x_in_data, genes_y = genes_y_in_data),
        file.path(RES_DIR, "GSE242252_processed.rds"))

message("\n===== GSE242252 preprocessing done =====")
message("Output in: ", RES_DIR)
