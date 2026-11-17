# scripts/figures/Supplementary_Figure_5.R
#
# Supplementary Figure 5 (single-cell RNA-seq panels, day 4 after treatment)
#   5b  % of each immune cell subtype among CD45+ cells per treatment (CSV; plotted in GraphPad Prism)
#   5c  pseudobulk heatmap of cytokines/chemokines across CD45- clusters
#   5d  GO/KEGG enrichment of PDS-vs-Vehicle DEGs in macrophages
#   5e  GO/KEGG enrichment of PDS-vs-Vehicle DEGs in dendritic cells
# (Supp Fig 5a and 5f-i are ELISA/survival/in vitro data and are not produced here.)
#
# Input : CD45pos_seurat.rds, CD45neg_seurat.rds (see _shared.R)
# Output: results/Supplementary_Figure_5/
# Run   : Rscript scripts/figures/Supplementary_Figure_5.R   (5d-e need internet access for Enrichr)

this.dir <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else file.path(getwd(), "scripts", "figures") })
source(file.path(this.dir, "_shared.R"))
out <- fig_dir("Supplementary_Figure_5")

log_msg("Supp Fig 5: loading CD45+ object")
day4.seurat <- make_day4_cd45pos(load_cd45pos())

## ---- 5b: cell-type composition (% of cells per treatment) ----
dist.table <- table(day4.seurat@meta.data$cell_type, day4.seurat@meta.data$Treatment)
freq.table <- as.data.frame(sweep(dist.table, 2, colSums(dist.table), FUN = "/") * 100)
colnames(freq.table) <- c("cell_type", "Treatment", "percent_of_treatment")
freq.table$Treatment <- factor(freq.table$Treatment, levels = c("Vehicle", "PDS"))
write.csv(freq.table, file.path(out, "Supplementary_Figure_5b_celltype_frequency_day4.csv"), row.names = FALSE)

## ---- 5d / 5e: DEGs -> Enrichr GO/KEGG dot plots ----
degs <- compute_day4_degs(day4.seurat, c("Macrophages", "DC"))
for (ct in names(degs)) write.csv(degs[[ct]], file.path(out, paste0("Supplementary_Figure_5de_DEGs_", ct, "_PDS_vs_Vehicle.csv")))
rm(day4.seurat)

run_celltype_go_dotplot(degs[["Macrophages"]], "Macrophages",
  GO.terms = c("cytokine-mediated signaling pathway", "regulation of inflammatory response", "inflammatory response",
    "toll-like receptor signaling pathway", "pattern recognition receptor signaling pathway",
    "cellular response to Interleukin-1", "translation", "gene expression", "ribosome biogenesis"),
  KEGG.terms = c("Lysosome", "Phagosome", "NOD-like receptor signaling pathway", "Apoptosis",
    "Toll-like receptor signaling pathway", "TNF signaling pathway", "Antigen processing and presentation",
    "Ribosome", "Oxidative phosphorylation"),
  save.path = file.path(out, "Supplementary_Figure_5d_Macrophages_GO_KEGG_dotplot.pdf"))

run_celltype_go_dotplot(degs[["DC"]], "DC",
  GO.terms = c("cytokine-mediated signaling pathway", "inflammatory response", "cellular response to Interleukin-1",
    "positive regulation of inflammatory response", "positive regulation of cytokine production",
    "ribosome biogenesis", "translation", "gene expression"),
  KEGG.terms = c("Phagosome", "NOD-like receptor signaling pathway", "Apoptosis", "Toll-like receptor signaling pathway",
    "TNF signaling pathway", "NF-kappa B signaling pathway", "RIG-I-like receptor signaling pathway",
    "Ribosome", "Oxidative phosphorylation"),
  save.path = file.path(out, "Supplementary_Figure_5e_DC_GO_KEGG_dotplot.pdf"))

## ---- 5c: CD45- cytokine/chemokine pseudobulk heatmap (all CD45- clusters) ----
log_msg("Supp Fig 5: loading CD45- object")
CD45neg.seurat <- load_cd45neg()
pseudo <- AggregateExpression(CD45neg.seurat, assays = "SCT", return.seurat = TRUE, group.by = c("cluster_anno"))
custom.genes.neg <- c("Il1b", "Il6", "Il12a", "Ifng", "Tnf", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
                      "Cxcl1", "Cxcl9", "Cxcl10", "Cxcl13")   # Il12b is not detected in the CD45- data
p <- DoHeatmap(pseudo, features = custom.genes.neg, draw.lines = FALSE, angle = 90, raster = FALSE) +
  scico::scale_fill_scico(palette = "batlow", direction = 1, na.value = "white")
save_pdf(p, file.path(out, "Supplementary_Figure_5c_heatmap_cytokines_chemokines_CD45neg.pdf"))

log_msg("Supplementary Figure 5 done -> ", out)
