#!/usr/bin/env Rscript
###############################################################################
# 03_GSE249632_scRNA_DEA.R
# RZV疫苗: gE特异性CD4+ T细胞 scRNA-seq 差异表达分析
# 时间点: D0, D14, D60, D74, D365
# 方法: Pseudobulk per donor/timepoint, 然后 limma-trend
###############################################################################
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(tidyverse)
  library(Matrix)
})

set.seed(42)

RES_DIR <- "/media/cmj/MechanicalDisk/yjs/VZV-geo/results/GSE249632"
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load data ─────────────────────────────────────────────────────────────
message("Loading GSE249632...")
gse <- getGEO("GSE249632", GSEMatrix = TRUE, getGPL = FALSE, AnnotGPL = FALSE)
eset <- gse[[1]]
pd   <- pData(eset)
expr <- exprs(eset)

message("Samples (cells): ", ncol(expr))
message("Genes: ", nrow(expr))

# Save metadata for inspection
cat("\nColumn names in pData:\n")
cat(paste0("  ", colnames(pd)), sep = "\n")

# ── 2. Understand sample structure ───────────────────────────────────────────
cat("\nSample titles (first 20):\n")
for (i in 1:min(20, length(pd$title))) {
  cat(sprintf("  [%d] %s\n", i, pd$title[i]))
}

# Try to extract timepoint/donor info from metadata
# SMART-Seq data may have specific columns
write.csv(pd, file.path(RES_DIR, "sample_metadata.csv"), row.names = FALSE)
message("\nFull metadata saved to sample_metadata.csv")

# ── 3. Parse time point and donor information ────────────────────────────────
# From the paper: time points are D0, D14, D60, D74, D365
# Need to extract from sample titles or characteristics

# Try extracting from characteristics columns
char_cols <- grep("characteristics", colnames(pd), value = TRUE)
for (cc in char_cols) {
  cat("\n", cc, " (first 5):\n", sep = "")
  print(head(pd[[cc]], 5))
}

# ── 4. Build annotation from sample titles ───────────────────────────────────
# Look for patterns like D0, D14, D60, D74, D365 in sample titles
sample_info <- data.frame(
  sample_id = colnames(expr),
  title = pd$title,
  stringsAsFactors = FALSE
)

# Extract time point
time_pattern <- "D(\\d+)|day\\s*(\\d+)|Day\\s*(\\d+)"
time_matches <- regmatches(pd$title, gregexpr(time_pattern, pd$title,
                                               ignore.case = TRUE, perl = TRUE))
sample_info$timepoint <- sapply(time_matches, function(x) {
  if (length(x) > 0) x[1] else NA
})
sample_info$timepoint <- toupper(sample_info$timepoint)
sample_info$timepoint <- gsub("DAY", "D", sample_info$timepoint)

cat("\nParsed time points:\n")
print(table(sample_info$timepoint))

# Extract donor ID
donor_pattern <- "Donor\\s*(\\d+)|donor\\s*(\\d+)|D(\\d+)_"
donor_matches <- regmatches(pd$title, gregexpr(donor_pattern, pd$title,
                                                ignore.case = TRUE, perl = TRUE))
sample_info$donor <- sapply(donor_matches, function(x) {
  if (length(x) > 0) x[1] else NA
})

cat("\nParsed donors:\n")
print(table(sample_info$donor))

# ── 5. Pseudobulk aggregation ────────────────────────────────────────────────
# Aggregate by donor × timepoint
sample_info$group <- paste(sample_info$donor, sample_info$timepoint, sep = "_")

groups_unique <- unique(sample_info$group)
message("\nUnique donor×timepoint groups: ", length(groups_unique))

pseudobulk <- matrix(NA, nrow = nrow(expr), ncol = length(groups_unique))
colnames(pseudobulk) <- groups_unique
rownames(pseudobulk) <- rownames(expr)

pb_meta <- data.frame(
  group = groups_unique,
  n_cells = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(groups_unique)) {
  cells <- which(sample_info$group == groups_unique[i])
  pb_meta$n_cells[i] <- length(cells)
  if (length(cells) == 1) {
    pseudobulk[, i] <- expr[, cells]
  } else {
    pseudobulk[, i] <- rowSums(expr[, cells, drop = FALSE])
  }
}

# Parse donor and timepoint from group name
pb_meta$donor <- gsub("_D\\d+$", "", pb_meta$group)
pb_meta$timepoint <- gsub(".*_(D\\d+)$", "\\1", pb_meta$group)

message("Pseudobulk matrix: ", nrow(pseudobulk), " genes × ", ncol(pseudobulk), " samples")
cat("\nCells per pseudobulk sample:\n")
print(summary(pb_meta$n_cells))

# ── 6. Filter low-expressed genes ────────────────────────────────────────────
# Keep genes with reasonable expression in at least N samples
# Use CPM > 1 in at least 3 samples
cpm <- t(t(pseudobulk) / colSums(pseudobulk)) * 1e6
keep_genes <- rowSums(cpm > 1) >= 3
message("Genes passing CPM filter: ", sum(keep_genes), " / ", length(keep_genes))

pseudobulk_filt <- pseudobulk[keep_genes, ]

# ── 7. limma-trend differential expression ───────────────────────────────────
# For each post-vaccination timepoint vs D0

dge <- DGEList(counts = pseudobulk_filt)
dge <- calcNormFactors(dge, method = "TMM")

# Design: ~ donor + timepoint (paired by donor)
timepoint <- factor(pb_meta$timepoint, levels = c("D0", "D14", "D60", "D74", "D365"))
donor <- factor(pb_meta$donor)

design <- model.matrix(~ donor + timepoint)
rownames(design) <- pb_meta$group

# Voom
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit, trend = TRUE)

# Extract results for each timepoint vs D0
time_coefs <- grep("^timepoint", colnames(coef(fit)), value = TRUE)

de_results <- list()
for (coef_name in time_coefs) {
  tp <- gsub("timepoint", "", coef_name)
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  res$gene <- rownames(res)
  de_results[[tp]] <- res

  n_up <- sum(res$adj.P.Val < 0.05 & res$logFC > 0)
  n_down <- sum(res$adj.P.Val < 0.05 & res$logFC < 0)
  message(sprintf("  %s vs D0: %d up, %d down (FDR<0.05)", tp, n_up, n_down))
}

# ── 8. Save results ──────────────────────────────────────────────────────────
saveRDS(list(
  pseudobulk = pseudobulk,
  pseudobulk_filt = pseudobulk_filt,
  pb_meta = pb_meta,
  de_results = de_results,
  v = v,
  fit = fit
), file.path(RES_DIR, "GSE249632_pseudobulk_DEA.rds"))

# Write DE tables
for (tp in names(de_results)) {
  write.csv(de_results[[tp]],
            file.path(RES_DIR, paste0("DE_", tp, "_vs_D0.csv")),
            row.names = FALSE)
}

message("\n===== GSE249632 pseudobulk DEA done =====")
message("Output in: ", RES_DIR)
