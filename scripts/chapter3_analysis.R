#!/usr/bin/env Rscript
###############################################################################
# 第三章: 公共转录组数据分析 — 完整可复现代码
# ============================================
# 运行方式: Rscript scripts/chapter3_analysis.R
#          或在RStudio中全选(Ctrl+A) → Run(Ctrl+Enter)
# 前置条件: results/GSE242252/filtered_counts.csv, deseq2_metadata.csv
#           results/GSE249632/pseudobulk_counts.csv, pseudobulk_metadata.csv
#
# 三个分析模块:
#   Part A — GSE242252: HZ患者全血Bulk RNA-seq, DESeq2差异分析
#   Part B — GSE249632: RZV疫苗CD4+ T细胞scRNA-seq pseudobulk, limma差异分析
#   Part C — 提取关键基因数据供图表制作
###############################################################################

# ===========================================================================
# 0. 环境设置: 定位项目根目录、加载依赖包、创建输出文件夹
# ===========================================================================

# 自动检测项目根目录(脚本所在目录的上一级)
# 在RStudio中运行时, 请手动设置下面这行为你的实际项目路径
# 例如: PROJ_ROOT <- "/media/cmj/MechanicalDisk/yjs/VZV-geo"
script_path <- tryCatch({
  sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
}, error = function(e) "")

if (length(script_path) == 0 || nchar(script_path) == 0) {
  # RStudio环境下无法自动检测, 使用当前工作目录或手动设置
  # 请取消下面这行注释并填入你的实际项目路径:
  PROJ_ROOT <- "/media/cmj/MechanicalDisk/yjs/VZV-geo"
  # 如果当前工作目录就是项目根目录, 也可以直接: PROJ_ROOT <- getwd()
} else {
  PROJ_ROOT <- normalizePath(dirname(dirname(script_path)))
}
message("项目根目录: ", PROJ_ROOT)

# 依赖包已安装，注释掉安装代码。首次运行需取消注释
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2","limma","edgeR","ggplot2","pheatmap","RColorBrewer"))

# 加载分析所需的R包
suppressPackageStartupMessages({
  library(DESeq2)        # 差异表达分析(Bulk RNA-seq)
  library(limma)         # 线性模型差异分析(scRNA-seq pseudobulk)
  library(edgeR)         # 表达数据标准化(TMM)
  library(ggplot2)       # 绑图(PCA、火山图)
  library(pheatmap)      # 热图
  library(RColorBrewer)  # 调色板
})

set.seed(42)  # 设置随机种子，确保结果可复现

# 定义输出目录
RES_HZ  <- file.path(PROJ_ROOT, "results", "GSE242252")   # HZ Bulk数据结果
RES_RZV <- file.path(PROJ_ROOT, "results", "GSE249632")   # RZV scRNA数据结果
FIG_DIR <- file.path(PROJ_ROOT, "results", "rnaseq")  # 图片输出
dir.create(RES_HZ,  recursive = TRUE, showWarnings = FALSE)
dir.create(RES_RZV, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)


# #############################################################################
#                                                                             #
#   PART A: GSE242252 — HZ患者全血Bulk RNA-seq差异表达分析                    #
#   数据来源: Vandoren et al. 2024, J Infect Dis                              #
#   设计: 26例HZ患者，急性发作期 vs 恢复期(约1年后)                            #
#   方法: DESeq2非配对设计 (保留全部样本)                                      #
#                                                                             #
# #############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("PART A: GSE242252 — HZ急性期 vs 恢复期 差异表达分析")
message(paste(rep("=", 70), collapse = ""))

# --- A1. 读入count矩阵和样本元数据 ---
# filtered_counts.csv: 基因×样本的原始count矩阵(已过滤低表达基因)
# deseq2_metadata.csv: 每列的样本信息(condition, patient_id, timepoint)
message("读入count矩阵和样本元数据...")
counts <- as.matrix(read.csv(
  file.path(RES_HZ, "filtered_counts.csv"), row.names = 1))
meta  <- read.csv(
  file.path(RES_HZ, "deseq2_metadata.csv"), row.names = 1)

message(sprintf("Count矩阵: %d 基因 × %d 样本", nrow(counts), ncol(counts)))
message("样本分组情况:")
print(table(meta$group))

# --- A2. 提取HZ样本(非配对设计) ---
# 使用非配对设计保留全部 26 acute + 23 convalescent = 49 样本
# 配对设计需排除3例仅急性期患者(patient 30, 31, 8), 会丢失统计效力
# 原文(Vandoren 2024)同样使用非配对比较
hz_meta <- subset(meta, condition_label == "Herpes_Zoster")
hz_samples <- rownames(hz_meta)
counts_hz <- counts[, hz_samples, drop = FALSE]

# 将timepoint设为因子, acute作为分子, convalescent作为分母
hz_meta$timepoint <- factor(hz_meta$timepoint, levels = c("acute", "convalescent"))

n_acute <- sum(hz_meta$timepoint == "acute")
n_conv  <- sum(hz_meta$timepoint == "convalescent")
message(sprintf("HZ样本: %d 急性期 + %d 恢复期 = %d 总计 (非配对设计)",
                n_acute, n_conv, nrow(hz_meta)))

# --- A3. 构建DESeq2对象并运行差异分析 ---
# 设计公式: ~ timepoint (非配对，比较急性vs恢复)
dds <- DESeqDataSetFromMatrix(
  countData = counts_hz,
  colData   = hz_meta,
  design    = ~ timepoint
)

# 预过滤: 至少在12个样本中count≥10 (12 = 较小样本组23的一半)
keep <- rowSums(counts(dds) >= 10) >= 12
dds <- dds[keep, ]
message(sprintf("预过滤后: %d 个基因", nrow(dds)))

message("运行DESeq2差异分析(估计size factor → 离散度 → 模型拟合)...")
dds <- DESeq(dds, parallel = FALSE)

# --- A4. 提取差异表达结果 ---
# contrast = c("timepoint", "acute", "convalescent")
# 含义: log2(acute/convalescent), 正值=急性期上调
res <- results(dds, contrast = c("timepoint", "acute", "convalescent"),
               alpha = 0.05)
res <- res[order(res$pvalue), ]  # 按p值排序

message("\n=== DESeq2 结果汇总 ===")
print(summary(res))

# 保存完整结果到CSV
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# --- 添加基因Symbol列 (通过org.Hs.eg.db映射Ensembl ID → Gene Symbol) ---
library(org.Hs.eg.db)
gene_ids <- gsub("\\..*", "", res_df$gene_id)  # 去掉Ensembl ID版本号
res_df$symbol <- mapIds(org.Hs.eg.db, keys = gene_ids,
                        column = "SYMBOL", keytype = "ENSEMBL",
                        multiVals = "first")
message(sprintf("已添加基因Symbol列 (匹配 %d / %d 个基因)",
                sum(!is.na(res_df$symbol)), nrow(res_df)))

write.csv(res_df, file.path(RES_HZ, "DE_HZ_acute_vs_convalescent.csv"),
          row.names = FALSE)

# 统计显著差异基因数量 (FDR < 0.05)
n_up   <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
message(sprintf("显著差异基因 (FDR<0.05): %d 上调, %d 下调", n_up, n_down))

# --- A5. PCA主成分分析 ---
# 用vst(variance stabilizing transformation)标准化后做PCA
# 观察急性期和恢复期样本是否在主成分空间分离
message("绘制PCA图...")
vsd <- vst(dds, blind = FALSE)               # VST标准化
pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)  # PCA
pca_var <- round(summary(pca)$importance[2, 1:5] * 100, 1)  # 各PC解释方差百分比

pca_df <- as.data.frame(pca$x)
pca_df$timepoint <- hz_meta$timepoint

pdf(file.path(FIG_DIR, "FigA_PCA_HZ.pdf"), width = 8, height = 6)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = timepoint)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, alpha = 0.2) +
  labs(x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)"),
       title = "GSE242252 PCA: HZ急性期 vs 恢复期 (非配对)") +
  scale_color_manual(values = c("acute" = "#E41A1C", "convalescent" = "#377EB8")) +
  theme_minimal(base_size = 14)
print(p)
dev.off()
message("PCA图已保存: ", file.path(FIG_DIR, "FigA_PCA_HZ.pdf"))

# --- A6. 火山图 ---
# X轴=log2FoldChange, Y轴=-log10(adjusted p-value)
# 红=上调(FDR<0.05, LFC>0.5), 蓝=下调, 灰=不显著
message("绘制火山图...")
res_plot <- res_df[!is.na(res_df$padj), ]
res_plot$sig <- "不显著"
res_plot$sig[res_plot$padj < 0.05 & res_plot$log2FoldChange > 0.58] <- "上调"
res_plot$sig[res_plot$padj < 0.05 & res_plot$log2FoldChange < -0.58] <- "下调"

pdf(file.path(FIG_DIR, "FigA_Volcano_HZ.pdf"), width = 10, height = 8)
p <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.5) +           # LFC=0线
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", alpha = 0.3) +     # |LFC|=1线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) + # FDR=0.05线
  scale_color_manual(values = c("下调" = "#377EB8", "不显著" = "grey80", "上调" = "#E41A1C")) +
  labs(x = "log2 倍数变化 (急性期 vs 恢复期)",
       y = expression(-log[10](调整后p值)),
       title = "GSE242252: HZ急性期 vs 恢复期 (非配对)",
       subtitle = paste0(n_up, " 上调, ", n_down, " 下调 (FDR<0.05)")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
print(p)
dev.off()
message("火山图已保存: ", file.path(FIG_DIR, "FigA_Volcano_HZ.pdf"))

# --- A7. 打印Top差异基因 ---
top_degs <- res_df[order(-abs(res_df$log2FoldChange)), ]
message("\n按|log2FC|排序的Top 20差异基因:")
print(head(top_degs[, c("gene_id", "log2FoldChange", "padj")], 20))


# #############################################################################
#                                                                             #
#   PART B: GSE249632 — RZV疫苗 CD4+ T细胞 Pseudobulk limma分析               #
#   数据来源: GEO GSE249632                                                    #
#   设计: 7例健康接种者, gE四聚体+CD4+ T细胞, 纵向追踪                        #
#   时间点: D0(基线), D14, D60, D74, D365                                     #
#   方法: pseudobulk聚合 → TMM标准化 → voom → limma配对设计                    #
#                                                                             #
# #############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("PART B: GSE249632 — RZV疫苗CD4+ T细胞应答 差异表达分析")
message(paste(rep("=", 70), collapse = ""))

# --- B1. 读入pseudobulk数据 ---
# pseudobulk_counts.csv: 基因 × donor_timepoint组的count矩阵
# pseudobulk_metadata.csv: 每列信息(donor, timepoint, n_cells)
message("读入pseudobulk数据...")
pb_counts <- as.matrix(read.csv(
  file.path(RES_RZV, "pseudobulk_counts.csv"), row.names = 1))
pb_meta <- read.csv(
  file.path(RES_RZV, "pseudobulk_metadata.csv"))

message(sprintf("Pseudobulk矩阵: %d 基因 × %d 样本组", nrow(pb_counts), ncol(pb_counts)))
message("各时间点样本数:")
print(table(pb_meta$timepoint))

# --- B2. 低表达基因过滤 ---
# 用edgeR的TMM标准化, 保留CPM>1的基因(在至少3个样本中)
dge <- DGEList(counts = pb_counts)
dge <- calcNormFactors(dge, method = "TMM")    # TMM标准化因子
cpm <- cpm(dge)                                 # 计算CPM
keep <- rowSums(cpm > 1) >= 3                   # CPM>1在≥3个样本
message(sprintf("通过CPM过滤的基因: %d / %d", sum(keep), length(keep)))
dge <- dge[keep, ]

# --- B3. 构建设计矩阵 ---
# ~ donor + timepoint: 配对设计(donor为blocking factor)
# 比较每个疫苗接种后时间点 vs D0基线
timepoint <- factor(pb_meta$timepoint, levels = c("D0", "D14", "D60", "D74", "D365"))
donor <- factor(pb_meta$donor)
design <- model.matrix(~ donor + timepoint)
rownames(design) <- colnames(pb_counts)

message("时间点系数: ",
        paste(grep("^timepoint", colnames(design), value = TRUE), collapse = ", "))

# --- B4. voom标准化 + limma线性模型 ---
# voom: 将count数据转换为log2-CPM并估计权重
# lmFit + eBayes: 经验贝叶斯调节的线性模型
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit, trend = TRUE)

# --- B5. 提取每个时间点的差异表达结果 ---
# 循环4个时间点系数(D14, D60, D74, D365), 每个都vs D0
tp_coefs <- grep("^timepoint", colnames(coef(fit)), value = TRUE)
de_all <- list()

for (coef in tp_coefs) {
  tp_name <- gsub("timepoint", "", coef)  # 提取时间点名称(如"D14")
  res <- topTable(fit, coef = coef, number = Inf, sort.by = "P")
  res$gene_id <- rownames(res)

  n_u <- sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)
  n_d <- sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)
  message(sprintf("  %s vs D0: %d 上调, %d 下调 (FDR<0.05)", tp_name, n_u, n_d))

  de_all[[tp_name]] <- res
  write.csv(res, file.path(RES_RZV, paste0("DE_", tp_name, "_vs_D0.csv")),
            row.names = FALSE)
}

# --- B6. 保存合并的logFC矩阵(用于后续象限图等跨数据集比较) ---
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
message("logFC矩阵已保存: ", file.path(RES_RZV, "logFC_matrix_all_timepoints.csv"))


# #############################################################################
#                                                                             #
#   PART C: 提取第三章关键基因数据                                             #
#   用于后续制作火山图、基因轨迹图、疾病/保护Signature表                       #
#                                                                             #
# #############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("PART C: 提取关键基因数据")
message(paste(rep("=", 70), collapse = ""))

# --- C1. HZ疾病Signature基因 (来自原文 Supplementary Table S5) ---
# 7个最显著上调的ISGs + 增殖/浆细胞标志基因
hz_sig_genes <- c("IFI44L", "IFI27", "RSAD2", "ISG15", "SERPING1", "SIGLEC1", "IFI44",
                  "TOP2A", "PTTG1", "MZB1", "BATF2", "MX1", "IFIT5", "OASL")
message("\nHZ疾病Signature基因 (Vandoren et al. 2024 Table S5):")
message(paste(hz_sig_genes, collapse = ", "))

# --- C2. RZV保护Signature基因 ---
# 在疫苗后持续上调的T细胞分化/调控/记忆基因
rzv_sig_genes <- c("ZEB2", "CTLA4", "ICOS", "HAVCR2")
message("\nRZV保护Signature基因:")
message(paste(rzv_sig_genes, collapse = ", "))

# --- C3. 各时间点DEG数量汇总(用于DEG时间线柱状图) ---
message("\n各时间点差异基因数量汇总 (用于图3.4):")
for (tp in names(de_all)) {
  n_u <- sum(de_all[[tp]]$adj.P.Val < 0.05 & de_all[[tp]]$logFC > 0, na.rm = TRUE)
  n_d <- sum(de_all[[tp]]$adj.P.Val < 0.05 & de_all[[tp]]$logFC < 0, na.rm = TRUE)
  message(sprintf("  %s vs D0: %d ↑ / %d ↓", tp, n_u, n_d))
}


# #############################################################################
#                                                                             #
#   分析完成                                                                   #
#                                                                             #
# #############################################################################

message("\n", paste(rep("=", 70), collapse = ""))
message("全部分析完成!")
message("输出文件:")
message("  GSE242252 DESeq2结果: ", file.path(RES_HZ, "DE_HZ_acute_vs_convalescent.csv"))
message("  GSE249632 DE表:       ", file.path(RES_RZV, "DE_D[14,60,74,365]_vs_D0.csv"))
message("  logFC合并矩阵:        ", file.path(RES_RZV, "logFC_matrix_all_timepoints.csv"))
message("  PCA图:                ", file.path(FIG_DIR, "FigA_PCA_HZ.pdf"))
message("  火山图:               ", file.path(FIG_DIR, "FigA_Volcano_HZ.pdf"))
message(paste(rep("=", 70), collapse = ""))
