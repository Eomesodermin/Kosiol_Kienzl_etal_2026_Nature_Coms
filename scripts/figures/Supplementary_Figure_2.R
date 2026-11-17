# scripts/figures/Supplementary_Figure_2.R
#
# Supplementary Figure 2a-h: type I interferon-stimulated gene (ISG) and
# cGAS/STING signatures, day 4, Vehicle vs PDS
#   a  ISG pseudobulk heatmap, myeloid cells (DC, Mac, Mon, Neu)
#   b  ISG pseudobulk heatmap, CD45- cells (Tumor, Fibroblasts)
#   c  ISG module-score violins, myeloid cells
#   d  ISG module-score violins, CD45- cells
#   e  cGAS/STING pseudobulk heatmap, myeloid cells
#   f  cGAS/STING pseudobulk heatmap, CD45- cells
#   g  cGAS/STING module-score violins, myeloid cells
#   h  cGAS/STING module-score violins, CD45- cells
# (Panel i, phospho-gamma-H2AX flow cytometry, is not produced here.)
#
# Statistics (c, d, g, h): two-sided Wilcoxon rank-sum test on per-cell module
# scores, Holm-corrected across cell types within each panel. Exact values for
# all comparisons are written to Supplementary_Figure_2_module_score_statistics.csv.
#
# Input : CD45pos_seurat.rds, CD45neg_seurat.rds (see _shared.R)
# Output: results/Supplementary_Figure_2/
# Run   : Rscript scripts/figures/Supplementary_Figure_2.R

this.dir <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else file.path(getwd(), "scripts", "figures") })
source(file.path(this.dir, "_shared.R"))
out <- fig_dir("Supplementary_Figure_2")

## ---- Gene signatures (curated; intersected with detected genes below) ----
isg.genes <- c(
  "Isg15", "Ifit1", "Ifit2", "Ifit3", "Irf7", "Stat1", "Stat2", "Oas1a",
  "Mx2", "Rsad2", "Usp18", "Ifi44", "Ifi44l", "Oasl1", "Gbp2", "Cxcl10",
  "Ifih1", "Ifi27", "Ifi27l2b", "Ifi30", "Ifi35", "Ifitm1", "Ifitm2", "Ifitm3",
  "Isg20", "Ddx60", "Eif2ak2", "Ly6e", "Oas1g", "Parp12", "Parp14", "Rtp4",
  "Oas1b", "Oas2", "Oas3", "Oasl2", "Ddx58")
cgas.genes <- c(
  "Cgas", "Ifi204", "Ddx41", "Tmem173", "Tbk1", "Ikbke", "Irf3", "Irf7",
  "Nfkb1", "Rela", "Trex1", "Enpp1", "Ifnb1", "Il6", "Tnf", "Cxcl10",
  "Ccl2", "Ccl5", "Aim2", "Pycard", "Casp1")

## ---- CD45+ myeloid cells, day 4 ----
log_msg("Supp Fig 2: loading objects")
cd45pos.seurat <- load_cd45pos()
myeloid.keep <- c("DC", "Macrophages", "Monocytes", "Neutrophils")
Idents(cd45pos.seurat) <- cd45pos.seurat@meta.data$Time
cd45pos.d4 <- subset(cd45pos.seurat, idents = "Day4")
Idents(cd45pos.d4) <- cd45pos.d4@meta.data$cell_type
cd45pos.d4 <- subset(cd45pos.d4, idents = myeloid.keep)
DefaultAssay(cd45pos.d4) <- "SCT"
cd45pos.d4$Treatment <- factor(cd45pos.d4$Treatment, levels = c("Vehicle", "PDS"))
cd45pos.d4$cell_type <- factor(cd45pos.d4$cell_type, levels = myeloid.keep)
myeloid.col.order <- c("DC_Vehicle", "Macrophages_Vehicle", "Monocytes_Vehicle", "Neutrophils_Vehicle",
                       "DC_PDS", "Macrophages_PDS", "Monocytes_PDS", "Neutrophils_PDS")
rm(cd45pos.seurat)

## ---- CD45- tumor cells (all Tumor_* clusters) and fibroblasts (day 4 only) ----
cd45neg.seurat <- load_cd45neg()
cd45neg.seurat$cell_type <- ifelse(grepl("^Tumor", cd45neg.seurat$cluster_anno), "Tumor",
                                   as.character(cd45neg.seurat$cluster_anno))
cd45neg.d4 <- subset(cd45neg.seurat, subset = cell_type %in% c("Tumor", "Fibroblasts"))
cd45neg.d4$Treatment <- factor(cd45neg.d4$Treatment, levels = c("Vehicle", "PDS"))
cd45neg.d4$cell_type <- factor(cd45neg.d4$cell_type, levels = c("Tumor", "Fibroblasts"))
cd45neg.col.order <- c("Tumor_Vehicle", "Fibroblasts_Vehicle", "Tumor_PDS", "Fibroblasts_PDS")
rm(cd45neg.seurat)

stats.all <- list()

## ---- a, b: ISG heatmaps ----
res <- plot_pseudobulk_heatmap(cd45pos.d4, isg.genes, c("cell_type", "Treatment"), myeloid.col.order, "ISG Signature - Myeloid Cells")
save_pdf(res$plot, file.path(out, "Supplementary_Figure_2a_ISG_heatmap_myeloid.pdf"), width = 8, height = 9)
res <- plot_pseudobulk_heatmap(cd45neg.d4, isg.genes, c("cell_type", "Treatment"), cd45neg.col.order, "ISG Signature - Tumor & Fibroblasts")
save_pdf(res$plot, file.path(out, "Supplementary_Figure_2b_ISG_heatmap_CD45neg.pdf"), width = 6, height = 9)

## ---- c, d: ISG module scores ----
cd45pos.d4 <- AddModuleScore(cd45pos.d4, features = list(intersect(isg.genes, rownames(cd45pos.d4))), name = "ISG", assay = "SCT")
cd45pos.d4$ISG_score <- cd45pos.d4$ISG1; cd45pos.d4$ISG1 <- NULL
res <- plot_violin_holm(cd45pos.d4, "ISG_score", "cell_type", "Treatment", "ISG Score", "ISG Module Score - Myeloid Cells", panel = "2c")
ggsave(file.path(out, "Supplementary_Figure_2c_ISG_violin_myeloid.pdf"), res$plot, width = 10, height = 6, units = "in")
stats.all$c <- res$stats

cd45neg.d4 <- AddModuleScore(cd45neg.d4, features = list(intersect(isg.genes, rownames(cd45neg.d4))), name = "ISG", assay = "SCT")
cd45neg.d4$ISG_score <- cd45neg.d4$ISG1; cd45neg.d4$ISG1 <- NULL
res <- plot_violin_holm(cd45neg.d4, "ISG_score", "cell_type", "Treatment", "ISG Score", "ISG Module Score - Tumor & Fibroblasts", panel = "2d")
ggsave(file.path(out, "Supplementary_Figure_2d_ISG_violin_CD45neg.pdf"), res$plot, width = 8, height = 6, units = "in")
stats.all$d <- res$stats

## ---- e, f: cGAS/STING heatmaps ----
res <- plot_pseudobulk_heatmap(cd45pos.d4, cgas.genes, c("cell_type", "Treatment"), myeloid.col.order, "cGAS/STING Signature - Myeloid Cells")
save_pdf(res$plot, file.path(out, "Supplementary_Figure_2e_cGAS_heatmap_myeloid.pdf"), width = 8, height = 7)
res <- plot_pseudobulk_heatmap(cd45neg.d4, cgas.genes, c("cell_type", "Treatment"), cd45neg.col.order, "cGAS/STING Signature - Tumor & Fibroblasts")
save_pdf(res$plot, file.path(out, "Supplementary_Figure_2f_cGAS_heatmap_CD45neg.pdf"), width = 6, height = 7)

## ---- g, h: cGAS/STING module scores ----
cd45pos.d4 <- AddModuleScore(cd45pos.d4, features = list(intersect(cgas.genes, rownames(cd45pos.d4))), name = "cGAS", assay = "SCT")
cd45pos.d4$cGAS_score <- cd45pos.d4$cGAS1; cd45pos.d4$cGAS1 <- NULL
res <- plot_violin_holm(cd45pos.d4, "cGAS_score", "cell_type", "Treatment", "cGAS/STING Score", "cGAS/STING Module Score - Myeloid Cells", panel = "2g")
ggsave(file.path(out, "Supplementary_Figure_2g_cGAS_violin_myeloid.pdf"), res$plot, width = 10, height = 6, units = "in")
stats.all$g <- res$stats

cd45neg.d4 <- AddModuleScore(cd45neg.d4, features = list(intersect(cgas.genes, rownames(cd45neg.d4))), name = "cGAS", assay = "SCT")
cd45neg.d4$cGAS_score <- cd45neg.d4$cGAS1; cd45neg.d4$cGAS1 <- NULL
res <- plot_violin_holm(cd45neg.d4, "cGAS_score", "cell_type", "Treatment", "cGAS/STING Score", "cGAS/STING Module Score - Tumor & Fibroblasts", panel = "2h")
ggsave(file.path(out, "Supplementary_Figure_2h_cGAS_violin_CD45neg.pdf"), res$plot, width = 8, height = 6, units = "in")
stats.all$h <- res$stats

## ---- exact statistics for all comparisons ----
write.csv(dplyr::bind_rows(stats.all), file.path(out, "Supplementary_Figure_2_module_score_statistics.csv"), row.names = FALSE)

log_msg("Supplementary Figure 2 done -> ", out)
