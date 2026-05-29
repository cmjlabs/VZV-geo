#!/usr/bin/env Rscript
###############################################################################
# GSE242252 DESeq2: HZ急性期 vs 恢复期 + PCA质控 (含性染色体分析)
###############################################################################
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(org.Hs.eg.db)
  library(RColorBrewer)
})

set.seed(42)

RES_DIR <- "/media/cmj/MechanicalDisk/yjs/VZV-geo/results/GSE242252"
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load count matrix and metadata ────────────────────────────────────────
message("Loading data...")
counts <- as.matrix(read.csv(file.path(RES_DIR, "filtered_counts.csv"), row.names = 1))
meta  <- read.csv(file.path(RES_DIR, "deseq2_metadata.csv"), row.names = 1)

message(sprintf("Count matrix: %d genes x %d samples", nrow(counts), ncol(counts)))
print(table(meta$group))

# ── 2. Focus on HZ paired samples (acute vs convalescent) ────────────────────
hz_meta <- subset(meta, condition_label == "Herpes_Zoster")
hz_samples <- rownames(hz_meta)
counts_hz <- counts[, hz_samples, drop = FALSE]

# Ensure proper factor levels
hz_meta$timepoint <- factor(hz_meta$timepoint, levels = c("acute", "convalescent"))
hz_meta$patient_id <- factor(hz_meta$patient_id)

message(sprintf("HZ subset: %d genes x %d samples", nrow(counts_hz), ncol(counts_hz)))
print(table(hz_meta$timepoint))

# ── 3. DESeq2 with paired design ─────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_hz,
  colData   = hz_meta,
  design    = ~ patient_id + timepoint  # paired by patient
)

# Pre-filter: keep genes with >= 10 counts in at least 26 samples (half the cohort)
keep <- rowSums(counts(dds) >= 10) >= 26
dds <- dds[keep, ]
message(sprintf("After DESeq2 pre-filter: %d genes", nrow(dds)))

dds <- DESeq2::DESeq(dds, parallel = FALSE)

# ── 4. Results: acute vs convalescent ────────────────────────────────────────
res <- results(dds, contrast = c("timepoint", "acute", "convalescent"),
               alpha = 0.05)
res <- res[order(res$pvalue), ]
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

summary(res)
message(sprintf("DEGs (FDR < 0.05): %d up, %d down",
                sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE),
                sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)))

# ── 5. Annotate gene symbols ─────────────────────────────────────────────────
gene_ids <- rownames(dds)
gene_symbols <- tryCatch({
  mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL",
         keytype = "ENSEMBL", multiVals = "first")
}, error = function(e) {
  message("org.Hs.eg.db not available, using Ensembl IDs only")
  setNames(gene_ids, gene_ids)
})

res_df$symbol <- gene_symbols[res_df$gene_id]

# Save full results
write.csv(res_df, file.path(RES_DIR, "DE_HZ_acute_vs_convalescent.csv"), row.names = FALSE)

# Top DEGs
top_up <- head(res_df[res_df$log2FoldChange > 0 & !is.na(res_df$padj), ], 50)
top_down <- head(res_df[res_df$log2FoldChange < 0 & !is.na(res_df$padj), ], 50)
message("\nTop 10 upregulated (acute vs convalescent):")
print(head(top_up[, c("symbol", "log2FoldChange", "padj")], 10))
message("\nTop 10 downregulated (acute vs convalescent):")
print(head(top_down[, c("symbol", "log2FoldChange", "padj")], 10))

# ── 6. PCA Analysis ──────────────────────────────────────────────────────────
# 6a. PCA with all genes
vsd <- vst(dds, blind = FALSE)
pca_all <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)
pca_var <- round(summary(pca_all)$importance[2, 1:10] * 100, 1)

pca_df <- as.data.frame(pca_all$x)
pca_df$patient_id <- hz_meta$patient_id
pca_df$timepoint <- hz_meta$timepoint

# 6b. Identify sex chromosome genes in the VST matrix
vst_mat <- assay(vsd)
vst_genes <- rownames(vst_mat)

# Get sex chr genes
suppressPackageStartupMessages(library(org.Hs.eg.db))
gene_chr <- tryCatch({
  select(org.Hs.eg.db, keys = vst_genes, columns = c("CHR"), keytype = "ENSEMBL")
}, error = function(e) NULL)

if (!is.null(gene_chr)) {
  chr_x_genes <- gene_chr$ENSEMBL[gene_chr$CHR == "X" & !is.na(gene_chr$CHR)]
  chr_y_genes <- gene_chr$ENSEMBL[gene_chr$CHR == "Y" & !is.na(gene_chr$CHR)]
  sex_genes <- c(chr_x_genes, chr_y_genes)
  message(sprintf("Sex chr genes in VST: X=%d, Y=%d",
                  length(intersect(vst_genes, chr_x_genes)),
                  length(intersect(vst_genes, chr_y_genes))))

  # 6c. PCA without sex chromosome genes
  vst_nosex <- vst_mat[setdiff(vst_genes, sex_genes), ]
  pca_nosex <- prcomp(t(vst_nosex), center = TRUE, scale. = TRUE)
  pca_nosex_var <- round(summary(pca_nosex)$importance[2, 1:10] * 100, 1)

  pca_nosex_df <- as.data.frame(pca_nosex$x)
  pca_nosex_df$patient_id <- hz_meta$patient_id
  pca_nosex_df$timepoint <- hz_meta$timepoint

  # 6d. Plot PC1 vs PC2 (no sex chr)
  pdf(file.path(RES_DIR, "PCA_HZ_nosexchr.pdf"), width = 10, height = 7)
  p <- ggplot(pca_nosex_df, aes(x = PC1, y = PC2, color = timepoint)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95, alpha = 0.3) +
    labs(x = paste0("PC1 (", pca_nosex_var[1], "%)"),
         y = paste0("PC2 (", pca_nosex_var[2], "%)"),
         title = "GSE242252: HZ Acute vs Convalescent (No Sex Chr Genes)") +
    scale_color_manual(values = c("acute" = "#E41A1C", "convalescent" = "#377EB8")) +
    theme_minimal(base_size = 14)
  print(p)
  dev.off()

  # 6e. Plot PC2 vs PC3
  pdf(file.path(RES_DIR, "PCA_HZ_nosexchr_PC2vsPC3.pdf"), width = 10, height = 7)
  p2 <- ggplot(pca_nosex_df, aes(x = PC2, y = PC3, color = timepoint)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95, alpha = 0.3) +
    labs(x = paste0("PC2 (", pca_nosex_var[2], "%)"),
         y = paste0("PC3 (", pca_nosex_var[3], "%)"),
         title = "GSE242252: PC2 vs PC3 (No Sex Chr Genes)") +
    scale_color_manual(values = c("acute" = "#E41A1C", "convalescent" = "#377EB8")) +
    theme_minimal(base_size = 14)
  print(p2)
  dev.off()

  # Check if PC1 is dominated by sex
  pc1_top <- head(sort(abs(pca_all$rotation[, 1]), decreasing = TRUE), 30)
  pc1_genes <- names(pc1_top)
  pc1_sex <- intersect(pc1_genes, sex_genes)
  message(sprintf("Sex chr genes in PC1 top 30 loadings: %d / 30", length(pc1_sex)))
  if (length(pc1_sex) > 0) {
    message("  -> ", paste(pc1_sex, collapse = ", "))
  }
}

# ── 7. Volcano plot ──────────────────────────────────────────────────────────
res_plot <- res_df[!is.na(res_df$padj), ]
res_plot$sig <- "NS"
res_plot$sig[res_plot$padj < 0.05 & res_plot$log2FoldChange > 1] <- "Up"
res_plot$sig[res_plot$padj < 0.05 & res_plot$log2FoldChange < -1] <- "Down"

# Label top genes
top_genes <- res_plot[order(res_plot$padj), ]
top_genes <- top_genes[abs(top_genes$log2FoldChange) > 2, ]
label_genes <- head(top_genes, 15)

pdf(file.path(RES_DIR, "Volcano_HZ_acute_vs_convalescent.pdf"), width = 10, height = 8)
p3 <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = c("Down" = "#377EB8", "NS" = "grey70", "Up" = "#E41A1C")) +
  labs(x = "log2 Fold Change (acute vs convalescent)",
       y = "-log10(adjusted p-value)",
       title = "HZ Acute vs Convalescent (GSE242252)") +
  theme_minimal(base_size = 14)
if (nrow(label_genes) > 0) {
  p3 <- p3 + ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = symbol), size = 3, max.overlaps = 20
  )
}
print(p3)
dev.off()

# ── 8. Save all ──────────────────────────────────────────────────────────────
saveRDS(list(dds = dds, res = res, vsd = vsd, pca_all = pca_all,
             pca_nosex = if (exists("pca_nosex")) pca_nosex else NULL),
        file.path(RES_DIR, "GSE242252_DESeq2_full.rds"))

message("\n===== GSE242252 DESeq2 analysis complete =====")
message("Output in: ", RES_DIR)
