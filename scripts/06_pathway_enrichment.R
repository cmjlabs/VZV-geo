#!/usr/bin/env Rscript
###############################################################################
# Pathway enrichment: GO + KEGG on HZ DEGs
# Use clusterProfiler to define functional gene modules for stratified scoring
###############################################################################
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(enrichplot)
})

set.seed(42)

RES_DIR <- "/media/cmj/MechanicalDisk/yjs/VZV-geo/results/GSE242252"
OUT_DIR <- "/media/cmj/MechanicalDisk/yjs/VZV-geo/results/pathway_enrichment"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load DEG results ──────────────────────────────────────────────────────
degs <- read.csv(file.path(RES_DIR, "DE_HZ_annotated.csv"))
message(sprintf("Loaded %d genes from DESeq2 results", nrow(degs)))

# Filter to significant DEGs
sig_degs <- degs[!is.na(degs$padj) & degs$padj < 0.05, ]
message(sprintf("Significant DEGs (FDR<0.05): %d", nrow(sig_degs)))

# ── 2. GO enrichment (upregulated genes) ─────────────────────────────────────
up_genes <- sig_degs$ensembl_id_clean[sig_degs$log2FoldChange > 0.5 & !is.na(sig_degs$ensembl_id_clean)]
down_genes <- sig_degs$ensembl_id_clean[sig_degs$log2FoldChange < -0.5 & !is.na(sig_degs$ensembl_id_clean)]

message(sprintf("Up-regulated (LFC>0.5, FDR<0.05): %d genes", length(up_genes)))
message(sprintf("Down-regulated (LFC<-0.5, FDR<0.05): %d genes", length(down_genes)))

# GO BP enrichment
go_up <- enrichGO(
  gene          = up_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

go_down <- enrichGO(
  gene          = down_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# Simplify GO results
if (!is.null(go_up) && nrow(go_up) > 0) {
  go_up_sim <- simplify(go_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
  go_up_df <- as.data.frame(go_up_sim)
  write.csv(go_up_df, file.path(OUT_DIR, "GO_BP_upregulated.csv"), row.names = FALSE)

  message(sprintf("GO BP terms (up): %d (after simplify: %d)", nrow(go_up), nrow(go_up_sim)))
  cat("\nTop GO BP terms (up-regulated):\n")
  print(head(go_up_df[, c("Description", "p.adjust", "Count")], 15))

  # Dotplot
  pdf(file.path(OUT_DIR, "GO_BP_upregulated_dotplot.pdf"), width = 14, height = 8)
  print(dotplot(go_up_sim, showCategory = 20, font.size = 10) +
        ggtitle("GO BP: Upregulated in HZ Acute vs Convalescent"))
  dev.off()
}

if (!is.null(go_down) && nrow(go_down) > 0) {
  go_down_sim <- simplify(go_down, cutoff = 0.7, by = "p.adjust", select_fun = min)
  go_down_df <- as.data.frame(go_down_sim)
  write.csv(go_down_df, file.path(OUT_DIR, "GO_BP_downregulated.csv"), row.names = FALSE)

  message(sprintf("\nGO BP terms (down): %d (after simplify: %d)", nrow(go_down), nrow(go_down_sim)))
  cat("\nTop GO BP terms (down-regulated):\n")
  print(head(go_down_df[, c("Description", "p.adjust", "Count")], 15))

  pdf(file.path(OUT_DIR, "GO_BP_downregulated_dotplot.pdf"), width = 14, height = 8)
  print(dotplot(go_down_sim, showCategory = 20, font.size = 10) +
        ggtitle("GO BP: Downregulated in HZ Acute vs Convalescent"))
  dev.off()
}

# ── 3. KEGG enrichment ───────────────────────────────────────────────────────
entrez_up <- tryCatch(
  bitr(up_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
  error = function(e) NULL
)
entrez_down <- tryCatch(
  bitr(down_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
  error = function(e) NULL
)

if (!is.null(entrez_up)) {
  kegg_up <- enrichKEGG(
    gene         = entrez_up$ENTREZID,
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    kegg_up_df <- as.data.frame(kegg_up)
    write.csv(kegg_up_df, file.path(OUT_DIR, "KEGG_upregulated.csv"), row.names = FALSE)
    message(sprintf("KEGG terms (up): %d", nrow(kegg_up)))
    print(head(kegg_up_df[, c("Description", "p.adjust", "Count")], 10))
  }
}

if (!is.null(entrez_down)) {
  kegg_down <- enrichKEGG(
    gene         = entrez_down$ENTREZID,
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    kegg_down_df <- as.data.frame(kegg_down)
    write.csv(kegg_down_df, file.path(OUT_DIR, "KEGG_downregulated.csv"), row.names = FALSE)
    message(sprintf("KEGG terms (down): %d", nrow(kegg_down)))
    print(head(kegg_down_df[, c("Description", "p.adjust", "Count")], 10))
  }
}

# ── 4. Define functional gene modules for stratified scoring ─────────────────
# Based on top GO terms and biological knowledge, define modules

define_module <- function(go_result, keywords, name) {
  if (is.null(go_result) || nrow(go_result) == 0) return(NULL)
  matched_terms <- go_result[grep(keywords, go_result$Description, ignore.case = TRUE), ]
  if (nrow(matched_terms) == 0) return(NULL)
  genes <- unique(unlist(strsplit(matched_terms$geneID, "/")))
  message(sprintf("  Module '%s': %d genes from %d GO terms", name, length(genes), nrow(matched_terms)))
  list(name = name, genes = genes, terms = matched_terms$Description)
}

modules_up <- list()
if (!is.null(go_up_sim)) {
  modules_up[["Type_I_IFN"]]     <- define_module(go_up_sim, "interferon.*type I|type I.*interferon|interferon-alpha|interferon-beta|response to interferon", "Type I Interferon Response")
  modules_up[["ISG15_Conjugation"]] <- define_module(go_up_sim, "ISG15", "ISG15 Conjugation")
  modules_up[["Viral_Transcription"]] <- define_module(go_up_sim, "viral.*transcri|translation.*viral|viral.*gene", "Viral Transcription/Translation")
  modules_up[["T_Cell_Activation"]] <- define_module(go_up_sim, "T.*cell.*activ|T.*cell.*receptor|TCR", "T Cell Activation")
  modules_up[["B_Cell_Activation"]] <- define_module(go_up_sim, "B.*cell.*activ|B.*cell.*receptor|BCR", "B Cell Activation")
  modules_up[["Complement"]]     <- define_module(go_up_sim, "complement|humoral", "Complement/Humoral Immunity")
  modules_up[["Cell_Cycle"]]     <- define_module(go_up_sim, "cell.*cycle|mitotic|division|proliferation|mitosis", "Cell Cycle & Proliferation")
  modules_up[["Antigen_Presentation"]] <- define_module(go_up_sim, "antigen.*process|MHC|peptide.*antigen", "Antigen Processing/Presentation")
}

# Remove NULL modules
modules_up <- modules_up[!sapply(modules_up, is.null)]

# Also define a combined "HZ_Disease_Signature" from all up-regulated DEGs
all_up_genes <- sig_degs$symbol[sig_degs$log2FoldChange > 0.5 & !is.na(sig_degs$padj) & sig_degs$padj < 0.05]
modules_up[["HZ_Disease_All"]] <- list(
  name = "HZ Disease Signature (All Up)",
  genes = all_up_genes,
  terms = "All up-regulated DEGs (LFC>0.5, FDR<0.05)"
)

modules_down <- list()
if (!is.null(go_down_sim)) {
  modules_down[["Immune_Suppression"]] <- define_module(go_down_sim, "suppress|negative.*regulat|inhibit", "Immune Suppression")
  modules_down[["Cell_Adhesion"]]      <- define_module(go_down_sim, "adhesion|junction|integrin", "Cell Adhesion")
}

modules_down <- modules_down[!sapply(modules_down, is.null)]

# ── 5. Save gene modules ─────────────────────────────────────────────────────
saveRDS(list(modules_up = modules_up, modules_down = modules_down,
             go_up = go_up_sim, go_down = go_down_sim,
             all_degs = sig_degs),
        file.path(OUT_DIR, "pathway_modules.rds"))

# Write gene lists
for (mod in c(modules_up, modules_down)) {
  writeLines(mod$genes, file.path(OUT_DIR, paste0("module_", gsub(" ", "_", mod$name), ".txt")))
}

message(sprintf("\nTotal modules defined: %d up, %d down", length(modules_up), length(modules_down)))
cat("\nModule summary:\n")
for (m in c(modules_up, modules_down)) {
  cat(sprintf("  %-30s : %3d genes\n", m$name, length(m$genes)))
}

message("\n===== Pathway enrichment complete =====")
message("Output in: ", OUT_DIR)
