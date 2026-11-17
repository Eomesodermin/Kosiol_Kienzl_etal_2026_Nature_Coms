# scripts/figures/Figure_5.R
#
# Figure 5 (single-cell RNA-seq panels, CD45+ cells, day 4 after treatment)
#   5c  UMAP of CD45+ cells, Vehicle
#   5d  UMAP of CD45+ cells, PDS
#   5e  pseudobulk heatmap of pro-inflammatory cytokines/chemokines per cell type x treatment
#   5f  GO/KEGG enrichment of PDS-vs-Vehicle DEGs in monocytes
#   5g  GO/KEGG enrichment of PDS-vs-Vehicle DEGs in neutrophils
# (Fig 5a-b and 5h-i are multiplex/ELISA data and are not produced here.)
#
# Input : CD45pos_seurat.rds (see _shared.R)
# Output: results/Figure_5/
# Run   : Rscript scripts/figures/Figure_5.R   (5f-g need internet access for Enrichr)

this.dir <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else file.path(getwd(), "scripts", "figures") })
source(file.path(this.dir, "_shared.R"))
out <- fig_dir("Figure_5")

log_msg("Figure 5: loading CD45+ object")
day4.seurat <- make_day4_cd45pos(load_cd45pos())

## ---- 5c / 5d: UMAPs (downsampled to 7,000 cells per treatment for a fair visual comparison) ----
Idents(day4.seurat) <- day4.seurat@meta.data$Treatment
small.seurat <- subset(day4.seurat, downsample = 7000)
for (arm in c("Vehicle", "PDS")) {
  Idents(small.seurat) <- small.seurat@meta.data$Treatment
  p <- scCustomize::DimPlot_scCustom(subset(small.seurat, idents = arm),
    reduction = "umap", label = TRUE, figure_plot = TRUE, group.by = "cell_type", pt.size = 1)
  panel <- if (arm == "Vehicle") "5c" else "5d"
  save_pdf(p, file.path(out, paste0("Figure_", panel, "_UMAP_", arm, ".pdf")))
}
rm(small.seurat)

## ---- 5e: pseudobulk heatmap, cytokine/chemokine panel ----
DefaultAssay(day4.seurat) <- "SCT"
pseudo <- AggregateExpression(day4.seurat, assays = "SCT", return.seurat = TRUE, group.by = c("cell_type", "Treatment"))
custom.genes <- c("Il1b", "Il6", "Il12a", "Il12b", "Ifng", "Tnf", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
                  "Cxcl1", "Cxcl9", "Cxcl10")   # Cxcl13 is not detected in the CD45+ data
p <- DoHeatmap(pseudo, features = custom.genes, draw.lines = FALSE, angle = 90, raster = FALSE) +
  ggplot2::scale_fill_gradient2(low = "blue", high = "red2")
save_pdf(p, file.path(out, "Figure_5e_heatmap_cytokines_chemokines.pdf"))
rm(pseudo)

## ---- 5f / 5g: DEGs -> Enrichr GO/KEGG dot plots ----
degs <- compute_day4_degs(day4.seurat, c("Monocytes", "Neutrophils"))
for (ct in names(degs)) write.csv(degs[[ct]], file.path(out, paste0("Figure_5fg_DEGs_", ct, "_PDS_vs_Vehicle.csv")))

run_celltype_go_dotplot(degs[["Monocytes"]], "Monocytes",
  GO.terms = c("cytokine-mediated signaling pathway", "inflammatory response", "positive regulation of inflammatory response",
    "cellular response to type I interferon", "positive regulation of cytokine production", "translation",
    "gene expression", "ribosome biogenesis"),
  KEGG.terms = c("Lysosome", "Phagosome", "NOD-like receptor signaling pathway", "Apoptosis",
    "Toll-like receptor signaling pathway", "TNF signaling pathway", "Ribosome", "Oxidative phosphorylation"),
  save.path = file.path(out, "Figure_5f_Monocytes_GO_KEGG_dotplot.pdf"))

run_celltype_go_dotplot(degs[["Neutrophils"]], "Neutrophils",
  GO.terms = c("cytokine-mediated signaling pathway", "neutrophil degranulation", "neutrophil mediated immunity",
    "cellular response to type I interferon", "phagosome acidification", "phagosome maturation",
    "translation", "gene expression", "ribosome biogenesis"),
  KEGG.terms = c("Lysosome", "Phagosome", "NOD-like receptor signaling pathway", "Apoptosis",
    "NF-kappa B signaling pathway", "TNF signaling pathway", "Antigen processing and presentation",
    "RIG-I-like receptor signaling pathway", "Ribosome", "Oxidative phosphorylation"),
  save.path = file.path(out, "Figure_5g_Neutrophils_GO_KEGG_dotplot.pdf"))

log_msg("Figure 5 done -> ", out)
