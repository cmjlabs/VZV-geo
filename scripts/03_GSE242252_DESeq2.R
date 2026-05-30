#!/usr/bin/env Rscript
###############################################################################
# GSE242252 DESeq2: HZ acute vs convalescent + PCA QC
# (Gene symbol mapping done separately in Python)
###############################################################################
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(42)

# Detect project root from script location (portable: no hardcoded paths)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0 || nchar(script_path) == 0) {
  stop("Please run via: Rscript scripts/03_GSE242252_DESeq2.R (not source() in RStudio)")
}
PROJ_ROOT <- dirname(dirname(normalizePath(script_path)))
RES_DIR <- file.path(PROJ_ROOT, "results/GSE242252")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load data ─────────────────────────────────────────────────────────────
message("Loading data...")
counts <- as.matrix(read.csv(file.path(RES_DIR, "filtered_counts.csv"), row.names = 1))
meta  <- read.csv(file.path(RES_DIR, "deseq2_metadata.csv"), row.names = 1)

message(sprintf("Count matrix: %d genes x %d samples", nrow(counts), ncol(counts)))
print(table(meta$group))

# ── 2. Focus on HZ paired samples ────────────────────────────────────────────
# Metadata contains 3 sample types:
#   - Herpes_Zoster (HZ): 26 acute + 23 convalescent → used for paired DEA
#   - Control (CO): 28 healthy donor samples → excluded (different biology)
#   - External_Control (ERC): 3 cross-study reference samples (CO053/CO109/CO271)
#     → excluded (external batch, no matched HZ timepoints)
hz_meta <- subset(meta, condition_label == "Herpes_Zoster")

# Identify patients with BOTH acute and convalescent samples (paired design requirement)
patient_tps <- table(hz_meta$patient_id, hz_meta$timepoint)
paired_patients <- rownames(patient_tps)[patient_tps[, "acute"] > 0 & patient_tps[, "convalescent"] > 0]
hz_meta <- hz_meta[hz_meta$patient_id %in% paired_patients, ]
excluded_patients <- setdiff(rownames(patient_tps), paired_patients)
message(sprintf("Excluded %d patients with only one timepoint: %s",
                length(excluded_patients),
                paste(excluded_patients, collapse = ", ")))

hz_samples <- rownames(hz_meta)
counts_hz <- counts[, hz_samples, drop = FALSE]

hz_meta$timepoint <- factor(hz_meta$timepoint, levels = c("acute", "convalescent"))
hz_meta$patient_id <- factor(hz_meta$patient_id)

message(sprintf("HZ paired subset: %d genes x %d samples (%d patients)",
                nrow(counts_hz), ncol(counts_hz), length(paired_patients)))
print(table(hz_meta$timepoint))

# ── 3. DESeq2 paired design ──────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_hz,
  colData   = hz_meta,
  design    = ~ patient_id + timepoint
)

# Pre-filter: >= 10 counts in >= half the samples per group (standard practice)
min_per_group <- ceiling(min(table(hz_meta$timepoint)) / 2)
keep <- rowSums(counts(dds) >= 10) >= min_per_group
message(sprintf("Pre-filter: >= 10 counts in >= %d samples (half of smallest group)", min_per_group))
dds <- dds[keep, ]
message(sprintf("After DESeq2 pre-filter: %d genes", nrow(dds)))

dds <- DESeq2::DESeq(dds, parallel = FALSE)

# ── 4. Results: acute vs convalescent ────────────────────────────────────────
res <- results(dds, contrast = c("timepoint", "acute", "convalescent"),
               alpha = 0.05)
res <- res[order(res$pvalue), ]

message("\n=== DESeq2 Results Summary ===")
print(summary(res))

n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
message(sprintf("DEGs (FDR<0.05): %d up, %d down", n_up, n_down))

# Save full results
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
write.csv(res_df, file.path(RES_DIR, "DE_HZ_acute_vs_convalescent.csv"), row.names = FALSE)

# Top DEG table
top_deg <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
top_deg <- top_deg[order(-abs(top_deg$log2FoldChange)), ]
message("\n=== Top 20 DEGs (by |log2FC|) ===")
print(head(top_deg[, c("gene_id", "log2FoldChange", "padj")], 20))

# ── 5. VST transformation and PCA ────────────────────────────────────────────
vsd <- vst(dds, blind = FALSE)
pca_all <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)
pca_var <- round(summary(pca_all)$importance[2, 1:10] * 100, 1)

pca_df <- as.data.frame(pca_all$x)
pca_df$patient_id <- hz_meta$patient_id
pca_df$timepoint <- hz_meta$timepoint

# ── 6. Identify sex chr genes using chromosome info from Ensembl IDs ─────────
# Load gene-to-chromosome mapping from the Ensembl annotation
# The Ensembl IDs may have version info (ENSGxxxxx.version)
vst_genes <- rownames(assay(vsd))

# Try to load org.Hs.eg.db if available; otherwise skip sex chr analysis
has_orgdb <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
if (has_orgdb) {
  library(org.Hs.eg.db)
  gene_chr <- tryCatch({
    AnnotationDbi::select(org.Hs.eg.db, keys = vst_genes,
                          columns = c("CHR"), keytype = "ENSEMBL")
  }, error = function(e) NULL)

  if (!is.null(gene_chr)) {
    chr_x_genes <- gene_chr$ENSEMBL[gene_chr$CHR == "X" & !is.na(gene_chr$CHR)]
    chr_y_genes <- gene_chr$ENSEMBL[gene_chr$CHR == "Y" & !is.na(gene_chr$CHR)]
    sex_genes <- unique(c(chr_x_genes, chr_y_genes))
    sex_genes_in_vst <- intersect(vst_genes, sex_genes)
    message(sprintf("Sex chr genes in VST: %d (X=%d, Y=%d)",
                    length(sex_genes_in_vst),
                    length(intersect(vst_genes, chr_x_genes)),
                    length(intersect(vst_genes, chr_y_genes))))

    # PCA without sex chr
    vst_nosex <- assay(vsd)[setdiff(vst_genes, sex_genes), ]
    pca_nosex <- prcomp(t(vst_nosex), center = TRUE, scale. = TRUE)
    pca_nosex_var <- round(summary(pca_nosex)$importance[2, 1:10] * 100, 1)

    # Check PC1 sex gene dominance
    pc1_loadings <- pca_all$rotation[, 1]
    pc1_top50 <- names(sort(abs(pc1_loadings), decreasing = TRUE)[1:50])
    pc1_sex_overlap <- intersect(pc1_top50, sex_genes)
    message(sprintf("Sex chr genes in PC1 top 50 loadings: %d", length(pc1_sex_overlap)))
    if (length(pc1_sex_overlap) > 0) {
      message("  -> First 10: ", paste(head(pc1_sex_overlap, 10), collapse = ", "))
    }
  }
} else {
  message("org.Hs.eg.db not available, skipping sex chromosome analysis")
  message("Run gene-level PCA-only below (no sex chr stratification possible)")
}

# ── 7. PCA Plots ─────────────────────────────────────────────────────────────
# 7a. All genes: PC1 vs PC2
pdf(file.path(RES_DIR, "PCA_HZ_allgenes_PC1vsPC2.pdf"), width = 10, height = 7)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = timepoint)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, alpha = 0.2) +
  labs(x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)"),
       title = "GSE242252 PCA: HZ Acute vs Convalescent (All Genes)") +
  scale_color_manual(values = c("acute" = "#E41A1C", "convalescent" = "#377EB8"),
                     name = "Timepoint") +
  theme_minimal(base_size = 14)
print(p)
dev.off()

# 7b. No sex chr (if available)
if (exists("pca_nosex")) {
  pca_nosex_df <- as.data.frame(pca_nosex$x)
  pca_nosex_df$patient_id <- hz_meta$patient_id
  pca_nosex_df$timepoint <- hz_meta$timepoint

  pdf(file.path(RES_DIR, "PCA_HZ_nosexchr_PC1vsPC2.pdf"), width = 10, height = 7)
  p2 <- ggplot(pca_nosex_df, aes(x = PC1, y = PC2, color = timepoint)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95, alpha = 0.2) +
    labs(x = paste0("PC1 (", pca_nosex_var[1], "%)"),
         y = paste0("PC2 (", pca_nosex_var[2], "%)"),
         title = "GSE242252 PCA: No Sex Chr Genes") +
    scale_color_manual(values = c("acute" = "#E41A1C", "convalescent" = "#377EB8"),
                       name = "Timepoint") +
    theme_minimal(base_size = 14)
  print(p2)
  dev.off()

  # PC2 vs PC3 (no sex chr)
  pdf(file.path(RES_DIR, "PCA_HZ_nosexchr_PC2vsPC3.pdf"), width = 10, height = 7)
  p3 <- ggplot(pca_nosex_df, aes(x = PC2, y = PC3, color = timepoint)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95, alpha = 0.2) +
    labs(x = paste0("PC2 (", pca_nosex_var[2], "%)"),
         y = paste0("PC3 (", pca_nosex_var[3], "%)"),
         title = "GSE242252 PCA: PC2 vs PC3 (No Sex Chr)") +
    scale_color_manual(values = c("acute" = "#E41A1C", "convalescent" = "#377EB8"),
                       name = "Timepoint") +
    theme_minimal(base_size = 14)
  print(p3)
  dev.off()
}

# ── 8. Volcano plot ──────────────────────────────────────────────────────────
# Load gene symbol annotations (from Python mygene mapping)
annot_file <- file.path(RES_DIR, "DE_HZ_annotated.csv")
if (file.exists(annot_file)) {
  annot <- read.csv(annot_file)
  res_df$symbol <- annot$symbol[match(res_df$gene_id, annot$gene_id)]
} else {
  res_df$symbol <- res_df$gene_id  # fallback to Ensembl ID
  message("DE_HZ_annotated.csv not found, using Ensembl IDs as labels")
}

res_plot <- res_df[!is.na(res_df$padj), ]
res_plot$sig <- "NS"
res_plot$sig[res_plot$padj < 0.05 & res_plot$log2FoldChange > 1] <- "Up (FDR<0.05, LFC>1)"
res_plot$sig[res_plot$padj < 0.05 & res_plot$log2FoldChange < -1] <- "Down (FDR<0.05, LFC<-1)"

n_up_lfc1 <- sum(res_plot$padj < 0.05 & res_plot$log2FoldChange > 1, na.rm = TRUE)
n_down_lfc1 <- sum(res_plot$padj < 0.05 & res_plot$log2FoldChange < -1, na.rm = TRUE)

# Label top DEGs by |log2FC|
top_pos <- res_plot[res_plot$padj < 0.05 & res_plot$log2FoldChange > 0.5, ]
top_pos <- top_pos[order(-top_pos$log2FoldChange), ]
top_neg <- res_plot[res_plot$padj < 0.05 & res_plot$log2FoldChange < -0.5, ]
top_neg <- top_neg[order(top_neg$log2FoldChange), ]
label_genes <- rbind(head(top_pos, 10), head(top_neg, 10))

x_max <- max(abs(res_plot$log2FoldChange), na.rm = TRUE) * 1.05

pdf(file.path(RES_DIR, "Volcano_HZ_acute_vs_convalescent.pdf"), width = 11, height = 8)
p4 <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.4, size = 0.6) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
  scale_x_continuous(limits = c(-x_max, x_max)) +
  scale_color_manual(
    values = c("Down (FDR<0.05, LFC<-1)" = "#377EB8",
               "NS" = "grey80",
               "Up (FDR<0.05, LFC>1)" = "#E41A1C"),
    name = "") +
  labs(x = "log2 Fold Change (acute vs convalescent)",
       y = expression(-log[10](adjusted~p-value)),
       title = "GSE242252: HZ Acute vs Convalescent",
       subtitle = paste0("DESeq2 paired, ", length(paired_patients), " patients | ",
                         n_up, " up, ", n_down, " down (FDR<0.05); ",
                         n_up_lfc1, " up, ", n_down_lfc1, " down (|LFC|>1)")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
if (nrow(label_genes) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
  p4 <- p4 + ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = symbol), size = 3.2, max.overlaps = 25, box.padding = 0.5
  )
}
print(p4)
dev.off()

# ── 9. Save ──────────────────────────────────────────────────────────────────
saveRDS(list(dds = dds, res = res, vsd = vsd, pca_all = pca_all,
             pca_nosex = if (exists("pca_nosex")) pca_nosex else NULL),
        file.path(RES_DIR, "GSE242252_DESeq2_full.rds"))

message("\n===== DESeq2 analysis complete =====")
message("Output in: ", RES_DIR)
