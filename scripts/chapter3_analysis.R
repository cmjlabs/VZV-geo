#!/usr/bin/env Rscript
###############################################################################
# Chapter 3: Complete Reproducible Analysis
# ==========================================
# 运行方式: Rscript scripts/chapter3_analysis.R
# 前置条件: 数据文件已下载到 data/ 和 results/ 目录
#
# 包含:
#   Part A - GSE242252 DESeq2 (HZ全血Bulk RNA-seq)
#   Part B - GSE249632 limma (RZV疫苗CD4+ T细胞 scRNA-seq pseudobulk)
#   Part C - 关键基因数据提取
###############################################################################

# ── 0. Environment Setup ─────────────────────────────────────────────────────
PROJ_ROOT <- normalizePath(dirname(dirname(
  sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
)))

message("Project root: ", PROJ_ROOT)

# 依赖包已安装，注释掉安装代码
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2","limma","edgeR","ggplot2","pheatmap","RColorBrewer"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(42)

# ── Output directories ───────────────────────────────────────────────────────
RES_HZ  <- file.path(PROJ_ROOT, "results", "GSE242252")
RES_RZV <- file.path(PROJ_ROOT, "results", "GSE249632")
FIG_DIR <- file.path(PROJ_ROOT, "results", "chapter3_figures")
dir.create(RES_HZ,  recursive = TRUE, showWarnings = FALSE)
dir.create(RES_RZV, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# PART A: GSE242252 — HZ全血Bulk RNA-seq DESeq2分析
# ══════════════════════════════════════════════════════════════════════════════
message("\n", paste(rep("=", 70), collapse = ""))
message("PART A: GSE242252 DESeq2 — HZ Acute vs Convalescent")
message(paste(rep("=", 70), collapse = ""))

# ── A1. Load data ────────────────────────────────────────────────────────────
message("Loading count matrix and metadata...")
counts <- as.matrix(read.csv(
  file.path(RES_HZ, "filtered_counts.csv"), row.names = 1))
meta  <- read.csv(
  file.path(RES_HZ, "deseq2_metadata.csv"), row.names = 1)

message(sprintf("Count matrix: %d genes x %d samples", nrow(counts), ncol(counts)))
message("Sample groups:")
print(table(meta$group))

# ── A2. Subset to HZ samples (all, unpaired design) ──────────────────────────
# 使用非配对设计以保留全部 26 acute + 23 convalescent 样本
# (配对设计需排除3例仅急性期患者: 30, 31, 8)
hz_meta <- subset(meta, condition_label == "Herpes_Zoster")
hz_samples <- rownames(hz_meta)
counts_hz <- counts[, hz_samples, drop = FALSE]

hz_meta$timepoint <- factor(hz_meta$timepoint, levels = c("acute", "convalescent"))

n_acute <- sum(hz_meta$timepoint == "acute")
n_conv  <- sum(hz_meta$timepoint == "convalescent")
message(sprintf("HZ samples: %d acute + %d convalescent = %d total (unpaired design)",
                n_acute, n_conv, nrow(hz_meta)))

# ── A3. DESeq2 (unpaired) ────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_hz,
  colData   = hz_meta,
  design    = ~ timepoint
)

# Pre-filter: >= 10 counts in >= 12 samples (half of smallest group)
keep <- rowSums(counts(dds) >= 10) >= 12
dds <- dds[keep, ]
message(sprintf("After pre-filter: %d genes", nrow(dds)))

message("Running DESeq2...")
dds <- DESeq(dds, parallel = FALSE)

# ── A4. Results ──────────────────────────────────────────────────────────────
res <- results(dds, contrast = c("timepoint", "acute", "convalescent"),
               alpha = 0.05)
res <- res[order(res$pvalue), ]

message("\n=== DESeq2 Results Summary ===")
print(summary(res))

# Save full results
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
write.csv(res_df, file.path(RES_HZ, "DE_HZ_acute_vs_convalescent.csv"),
          row.names = FALSE)

n_up   <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
message(sprintf("DEGs (FDR<0.05): %d up, %d down", n_up, n_down))

# ── A5. PCA ──────────────────────────────────────────────────────────────────
message("Generating PCA plot...")
vsd <- vst(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)
pca_var <- round(summary(pca)$importance[2, 1:5] * 100, 1)

pca_df <- as.data.frame(pca$x)
pca_df$timepoint <- hz_meta$timepoint

pdf(file.path(FIG_DIR, "FigA_PCA_HZ.pdf"), width = 8, height = 6)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = timepoint)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, alpha = 0.2) +
  labs(x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)"),
       title = "GSE242252 PCA: HZ Acute vs Convalescent (Unpaired)") +
  scale_color_manual(values = c("acute" = "#E41A1C", "convalescent" = "#377EB8")) +
  theme_minimal(base_size = 14)
print(p)
dev.off()
message("PCA saved: ", file.path(FIG_DIR, "FigA_PCA_HZ.pdf"))

# ── A6. Volcano plot ─────────────────────────────────────────────────────────
message("Generating volcano plot...")
res_plot <- res_df[!is.na(res_df$padj), ]
res_plot$sig <- "NS"
res_plot$sig[res_plot$padj < 0.05 & res_plot$log2FoldChange > 0.5] <- "Up"
res_plot$sig[res_plot$padj < 0.05 & res_plot$log2FoldChange < -0.5] <- "Down"

pdf(file.path(FIG_DIR, "FigA_Volcano_HZ.pdf"), width = 10, height = 8)
p <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = c("Down" = "#377EB8", "NS" = "grey80", "Up" = "#E41A1C")) +
  labs(x = "log2 Fold Change (acute vs convalescent)",
       y = expression(-log[10](adjusted~p-value)),
       title = "GSE242252: HZ Acute vs Convalescent (Unpaired)",
       subtitle = paste0(n_up, " up, ", n_down, " down (FDR<0.05)")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
print(p)
dev.off()
message("Volcano saved: ", file.path(FIG_DIR, "FigA_Volcano_HZ.pdf"))

# ── A7. Top DEG table ────────────────────────────────────────────────────────
top_degs <- res_df[order(-abs(res_df$log2FoldChange)), ]
message("\nTop 20 DEGs (by |log2FoldChange|):")
print(head(top_degs[, c("gene_id", "log2FoldChange", "padj")], 20))

# ══════════════════════════════════════════════════════════════════════════════
# PART B: GSE249632 — RZV疫苗CD4+ T细胞 Pseudobulk limma分析
# ══════════════════════════════════════════════════════════════════════════════
message("\n", paste(rep("=", 70), collapse = ""))
message("PART B: GSE249632 limma — RZV Vaccine CD4+ T Cell Response")
message(paste(rep("=", 70), collapse = ""))

# ── B1. Load pseudobulk data ─────────────────────────────────────────────────
message("Loading pseudobulk data...")
pb_counts <- as.matrix(read.csv(
  file.path(RES_RZV, "pseudobulk_counts.csv"), row.names = 1))
pb_meta <- read.csv(
  file.path(RES_RZV, "pseudobulk_metadata.csv"))

message(sprintf("Pseudobulk: %d genes x %d samples", nrow(pb_counts), ncol(pb_counts)))
message("Samples per timepoint:")
print(table(pb_meta$timepoint))

# ── B2. Filter low-expressed genes ───────────────────────────────────────────
dge <- DGEList(counts = pb_counts)
dge <- calcNormFactors(dge, method = "TMM")
cpm <- cpm(dge)
keep <- rowSums(cpm > 1) >= 3
message(sprintf("Genes passing CPM filter: %d / %d", sum(keep), length(keep)))
dge <- dge[keep, ]

# ── B3. Design matrix: ~ donor + timepoint (paired) ──────────────────────────
timepoint <- factor(pb_meta$timepoint, levels = c("D0", "D14", "D60", "D74", "D365"))
donor <- factor(pb_meta$donor)
design <- model.matrix(~ donor + timepoint)
rownames(design) <- colnames(pb_counts)

message("Timepoint coefficients: ",
        paste(grep("^timepoint", colnames(design), value = TRUE), collapse = ", "))

# ── B4. Voom + limma ─────────────────────────────────────────────────────────
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit, trend = TRUE)

# ── B5. Extract DE results ───────────────────────────────────────────────────
tp_coefs <- grep("^timepoint", colnames(coef(fit)), value = TRUE)
de_all <- list()

for (coef in tp_coefs) {
  tp_name <- gsub("timepoint", "", coef)
  res <- topTable(fit, coef = coef, number = Inf, sort.by = "P")
  res$gene_id <- rownames(res)

  n_u <- sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)
  n_d <- sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)
  message(sprintf("  %s vs D0: %d up, %d down (FDR<0.05)", tp_name, n_u, n_d))

  de_all[[tp_name]] <- res
  write.csv(res, file.path(RES_RZV, paste0("DE_", tp_name, "_vs_D0.csv")),
            row.names = FALSE)
}

# ── B6. Save logFC matrix ────────────────────────────────────────────────────
logfc_list <- lapply(names(de_all), function(tp) {
  setNames(de_all[[tp]]$logFC, de_all[[tp]]$gene_id)
})
names(logfc_list) <- names(de_all)
all_genes <- unique(unlist(lapply(logfc_list, names)))
logfc_mat <- do.call(cbind, lapply(logfc_list, function(x) {
  out <- rep(NA_real_, length(all_genes))
  names(out) <- all_genes
  out[names(x)] <- x
  out[is.na(out)] <- 0
  out
}))
write.csv(logfc_mat, file.path(RES_RZV, "logFC_matrix_all_timepoints.csv"))

# ══════════════════════════════════════════════════════════════════════════════
# PART C: Extract key gene data for Chapter 3 figures
# ══════════════════════════════════════════════════════════════════════════════
message("\n", paste(rep("=", 70), collapse = ""))
message("PART C: Key Gene Data Extraction")
message(paste(rep("=", 70), collapse = ""))

# ── C1. HZ Disease Signature genes (from paper Table S5) ─────────────────────
hz_sig_genes <- c("IFI44L", "IFI27", "RSAD2", "ISG15", "SERPING1", "SIGLEC1", "IFI44",
                  "TOP2A", "PTTG1", "MZB1", "BATF2", "MX1", "IFIT5", "OASL")
message("\nHZ Disease Signature genes (paper Table S5):")
for (g in hz_sig_genes) {
  # Look up in results by gene_id pattern
  message(sprintf("  %s", g))
}

# ── C2. RZV Protection Signature genes ───────────────────────────────────────
rzv_sig_genes <- c("ZEB2", "CTLA4", "ICOS", "HAVCR2")
message("\nRZV Protection Signature genes:")
for (g in rzv_sig_genes) {
  message(sprintf("  %s", g))
}

# ── C3. DEG timeline summary ─────────────────────────────────────────────────
message("\nDEG timeline summary:")
for (tp in names(de_all)) {
  n_u <- sum(de_all[[tp]]$adj.P.Val < 0.05 & de_all[[tp]]$logFC > 0, na.rm = TRUE)
  n_d <- sum(de_all[[tp]]$adj.P.Val < 0.05 & de_all[[tp]]$logFC < 0, na.rm = TRUE)
  message(sprintf("  %s vs D0: %d up, %d down", tp, n_u, n_d))
}

# ══════════════════════════════════════════════════════════════════════════════
# Done
# ══════════════════════════════════════════════════════════════════════════════
message("\n", paste(rep("=", 70), collapse = ""))
message("Analysis complete!")
message("Output files:")
message("  GSE242252: ", RES_HZ)
message("  GSE249632: ", RES_RZV)
message("  Figures:    ", FIG_DIR)
message(paste(rep("=", 70), collapse = ""))
