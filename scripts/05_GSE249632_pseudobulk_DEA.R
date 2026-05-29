#!/usr/bin/env Rscript
###############################################################################
# GSE249632: Pseudobulk DE analysis - RZV vaccine CD4+ T cell response
# Timepoints: D14, D60, D74, D365 vs D0
# Method: limma-trend with donor as blocking factor (paired)
###############################################################################
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
})

set.seed(42)

RES_DIR <- "/media/cmj/MechanicalDisk/yjs/VZV-geo/results/GSE249632"
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load pseudobulk data ──────────────────────────────────────────────────
message("Loading pseudobulk counts...")
pb_counts <- as.matrix(read.csv(file.path(RES_DIR, "pseudobulk_counts.csv"),
                                row.names = 1))
pb_meta <- read.csv(file.path(RES_DIR, "pseudobulk_metadata.csv"))

message(sprintf("Pseudobulk matrix: %d genes x %d samples", nrow(pb_counts), ncol(pb_counts)))
print(table(pb_meta$timepoint))

# ── 2. Filter low-expressed genes ────────────────────────────────────────────
# CPM > 1 in at least 3 samples
dge <- DGEList(counts = pb_counts)
dge <- calcNormFactors(dge, method = "TMM")

cpm <- cpm(dge)
keep <- rowSums(cpm > 1) >= 3
message(sprintf("Genes passing CPM filter: %d / %d", sum(keep), length(keep)))
dge <- dge[keep, ]

# ── 3. Design matrix: ~ donor + timepoint ────────────────────────────────────
timepoint <- factor(pb_meta$timepoint, levels = c("D0", "D14", "D60", "D74", "D365"))
donor <- factor(pb_meta$donor)

design <- model.matrix(~ donor + timepoint)
rownames(design) <- colnames(pb_counts)

message("Design matrix dimensions: ", nrow(design), " x ", ncol(design))
cat("Timepoint coefficients:\n")
print(colnames(design)[grep("^timepoint", colnames(design))])

# ── 4. Voom + limma ──────────────────────────────────────────────────────────
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit, trend = TRUE)

# ── 5. Extract DE results for each timepoint vs D0 ──────────────────────────
tp_coefs <- grep("^timepoint", colnames(coef(fit)), value = TRUE)

de_all <- list()
for (coef in tp_coefs) {
  tp_name <- gsub("timepoint", "", coef)
  res <- topTable(fit, coef = coef, number = Inf, sort.by = "P")
  res$gene_id <- rownames(res)

  n_up <- sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)
  n_down <- sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)
  message(sprintf("  %s vs D0: %d up, %d down (FDR<0.05)", tp_name, n_up, n_down))

  de_all[[tp_name]] <- res
  write.csv(res, file.path(RES_DIR, paste0("DE_", tp_name, "_vs_D0.csv")),
            row.names = FALSE)
}

# ── 6. Save combined log2FC matrix for cross-dataset comparison ──────────────
# Extract logFC for all timepoints
logfc_list <- list()
for (tp in names(de_all)) {
  df <- de_all[[tp]]
  lfc <- setNames(df$logFC, df$gene_id)
  logfc_list[[tp]] <- lfc
}

# Combine into a single matrix
all_genes <- unique(unlist(lapply(logfc_list, names)))
logfc_mat <- do.call(cbind, lapply(logfc_list, function(x) x[all_genes]))
rownames(logfc_mat) <- all_genes
logfc_mat[is.na(logfc_mat)] <- 0

write.csv(logfc_mat, file.path(RES_DIR, "logFC_matrix_all_timepoints.csv"))

# Also extract p-values
pval_list <- list()
for (tp in names(de_all)) {
  df <- de_all[[tp]]
  pv <- setNames(df$adj.P.Val, df$gene_id)
  pval_list[[tp]] <- pv
}
pval_mat <- do.call(cbind, lapply(pval_list, function(x) x[all_genes]))
rownames(pval_mat) <- all_genes
write.csv(pval_mat, file.path(RES_DIR, "padj_matrix_all_timepoints.csv"))

# ── 7. Heatmap of top variable genes across timepoints ───────────────────────
# Select genes significant in at least 1 timepoint
sig_genes <- unique(unlist(lapply(de_all, function(x) {
  rownames(x)[x$adj.P.Val < 0.05 & abs(x$logFC) > 0.5 & !is.na(x$adj.P.Val)]
})))
message(sprintf("Genes significant (FDR<0.05 & |logFC|>0.5) in >=1 tp: %d",
                length(sig_genes)))

if (length(sig_genes) > 200) {
  # Take top by variance across timepoints
  lfc_sig <- logfc_mat[sig_genes, , drop = FALSE]
  var_rank <- order(-apply(lfc_sig, 1, var))
  sig_genes <- sig_genes[var_rank[1:min(200, length(var_rank))]]
}

if (length(sig_genes) > 2) {
  # Get mean expression per timepoint for heatmap
  expr_means <- matrix(NA, nrow = length(sig_genes), ncol = length(levels(timepoint)))
  rownames(expr_means) <- sig_genes
  colnames(expr_means) <- levels(timepoint)

  for (tp in levels(timepoint)) {
    samples <- colnames(pb_counts)[pb_meta$timepoint == tp]
    expr_means[, tp] <- rowMeans(cpm[sig_genes, samples, drop = FALSE])
  }

  # Z-score normalize
  expr_z <- t(scale(t(expr_means)))

  pdf(file.path(RES_DIR, "Heatmap_DEGs_across_timepoints.pdf"), width = 8, height = 14)
  pheatmap(expr_z,
           cluster_rows = TRUE, cluster_cols = FALSE,
           show_rownames = (length(sig_genes) <= 50),
           main = paste0("RZV CD4+ T Cell Response (", length(sig_genes), " DEGs)"),
           color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
           fontsize = 8)
  dev.off()
}

# ── 8. Line plot for top DEGs ────────────────────────────────────────────────
top_genes <- head(sig_genes, 12)
if (length(top_genes) > 0) {
  pdf(file.path(RES_DIR, "TopDEGs_lineplot.pdf"), width = 12, height = 8)

  plot_data <- data.frame()
  for (g in top_genes) {
    for (tp in levels(timepoint)) {
      samples <- colnames(pb_counts)[pb_meta$timepoint == tp]
      vals <- cpm[g, samples]
      plot_data <- rbind(plot_data, data.frame(
        gene = g, timepoint = tp,
        mean = mean(vals), se = sd(vals) / sqrt(length(vals))
      ))
    }
  }

  p <- ggplot(plot_data, aes(x = timepoint, y = mean, group = gene, color = gene)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1, alpha = 0.4) +
    labs(x = "Timepoint", y = "Mean log2 CPM",
         title = "Top DEGs Across Vaccination Timepoints (GSE249632)") +
    theme_minimal(base_size = 13)
  print(p)
  dev.off()
}

# ── 9. Save ──────────────────────────────────────────────────────────────────
saveRDS(list(dge = dge, v = v, fit = fit, de_all = de_all,
             logfc_mat = logfc_mat, pval_mat = pval_mat, pb_meta = pb_meta),
        file.path(RES_DIR, "GSE249632_pseudobulk_DEA.rds"))

message("\n===== GSE249632 pseudobulk DEA complete =====")
message("Output in: ", RES_DIR)
