# scripts/figures/_shared.R
#
# Shared set-up for the per-figure scripts in scripts/figures/.
# Sourced at the top of Figure_5.R, Supplementary_Figure_2.R and
# Supplementary_Figure_5.R; not meant to be run on its own.
#
# Provides
#   repo.root           repository root (working directory is set to it)
#   obj.dir             folder holding the annotated Seurat objects
#   load_cd45pos(), load_cd45neg()
#   fig_dir(name)       creates and returns results/<name>/
#   save_pdf()          headless print-to-PDF
#   compute_day4_degs() Day-4 PDS-vs-Vehicle DEGs per myeloid cell type
#   run_celltype_go_dotplot()   Enrichr GO/KEGG dot plot (Fig 5f-g, Supp Fig 5d-e)
#   plot_pseudobulk_heatmap()   ordered pseudobulk heatmap (Supp Fig 2a,b,e,f)
#   plot_violin_holm()          module-score violin + Holm-adjusted Wilcoxon (Supp Fig 2c,d,g,h)
#
# INPUT OBJECTS
#   CD45pos_seurat.rds and CD45neg_seurat.rds, produced by
#   scripts/02_CD45Pos_data_annotation.Rmd and scripts/03_CD45Neg_data_annotation.Rmd
#   and archived on Zenodo (doi:10.5281/zenodo.21391578). By default they are
#   read from saves/; point elsewhere with the environment variable
#   KOSIOL_SEURAT_DIR, e.g.  KOSIOL_SEURAT_DIR=/path/to/objects Rscript scripts/figures/run_all.R

## ---- working directory: repository root ----
.get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", args[hit[1]]))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  file.path(getwd(), "scripts", "figures")  # interactive fallback: cwd = repo root
}
if (!exists("repo.root")) {
  repo.root <- normalizePath(file.path(.get_script_dir(), "..", ".."))
}
setwd(repo.root)
set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)        # 5.5.0
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(scCustomize)   # 3.0.1 (DimPlot_scCustom)
  library(scico)         # 1.5.0 (batlow palette)
  library(enrichR)       # 3.4
  library(scales)
  library(SeuratExtend)  # 1.2.4 (VlnPlot2)
  library(ggpubr)        # 0.6.3 (compare_means, stat_pvalue_manual)
})

obj.dir <- Sys.getenv("KOSIOL_SEURAT_DIR", unset = "saves")

load_cd45pos <- function() {
  obj <- LoadSeuratRds(file.path(obj.dir, "CD45pos_seurat.rds"))
  DefaultAssay(obj) <- "SCT"
  obj
}
load_cd45neg <- function() {
  obj <- LoadSeuratRds(file.path(obj.dir, "CD45neg_seurat.rds"))
  DefaultAssay(obj) <- "SCT"
  obj
}

fig_dir <- function(name) {
  d <- file.path("results", name)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

save_pdf <- function(plot_obj, path, width = 10, height = 8) {
  pdf(path, width = width, height = height)
  print(plot_obj)
  dev.off()
}

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))


## ---- Day-4 CD45+ object (basophils removed; only 2 basophils at Day 4) ----
make_day4_cd45pos <- function(CD45pos.seurat) {
  CD45pos.seurat@meta.data$Treatment <- factor(CD45pos.seurat@meta.data$Treatment, levels = c("Vehicle", "PDS"))
  Idents(CD45pos.seurat) <- CD45pos.seurat@meta.data$Time
  day4 <- subset(CD45pos.seurat, idents = "Day4")
  Idents(day4) <- day4@meta.data$cell_type
  subset(day4, idents = "Basophils", invert = TRUE)
}


## ---- Differential expression: PDS vs Vehicle within each Day-4 cell type ----
## Wilcoxon rank-sum (Seurat FindMarkers default) on the SCT assay.
compute_day4_degs <- function(day4.seurat, celltypes) {
  day4.seurat <- PrepSCTFindMarkers(day4.seurat)
  DefaultAssay(day4.seurat) <- "SCT"
  day4.seurat@meta.data$Treatment_cluster <- paste0(day4.seurat@meta.data$Treatment, "_", day4.seurat@meta.data$cell_type)
  Idents(day4.seurat) <- "Treatment_cluster"
  degs <- list()
  for (ct in celltypes) {
    log_msg("FindMarkers Day4 ", ct, " (PDS vs Vehicle)")
    degs[[ct]] <- FindMarkers(day4.seurat, assay = "SCT",
      ident.1 = paste0("PDS_", ct), ident.2 = paste0("Vehicle_", ct),
      only.pos = FALSE, recorrect_umi = FALSE, verbose = FALSE)
  }
  degs
}


## ---- Enrichr over-representation analysis -> curated GO/KEGG dot plot ----
## Up- (PDS) and down-regulated (Vehicle) DEGs (p_val_adj < 0.05) are queried
## against GO_Biological_Process_2021 and KEGG_2019_Mouse; terms are filtered
## to a curated list and only terms significant in the PDS direction are plotted.
## Requires internet access to maayanlab.cloud/Enrichr.
dbs_use <- c("GO_Biological_Process_2021", "KEGG_2019_Mouse")

filter_terms <- function(data, go_terms, group_name, prefix = "") {
  data %>%
    dplyr::select(c(Term, Overlap, Adjusted.P.value, Odds.Ratio, Combined.Score)) %>%
    filter(grepl(paste(go_terms, collapse = "|"), Term, ignore.case = TRUE)) %>%
    mutate(Term = paste0(prefix, Term), group = group_name)
}

process_data <- function(input.1, input.2, input.3, input.4, pval_threshold = 0.05) {
  tab <- bind_rows(input.1, input.2, input.3, input.4) %>%
    mutate(Overlap_Ratio = as.numeric(sub("/.*", "", Overlap)) / as.numeric(sub(".*/", "", Overlap))) %>%
    mutate(color_value = ifelse(Adjusted.P.value > pval_threshold, NA, Combined.Score)) %>%
    mutate(color_value_log = log1p(color_value))
  tab$group <- factor(tab$group, levels = c("PDS", "Vehicle"))
  tab %>% arrange(group, desc(Combined.Score)) %>% mutate(Term_order = factor(Term, levels = unique(Term)))
}

plot_pds_sig <- function(term.table, title.var = "", pval_threshold = 0.05, save.path, width.var = 12, height.var = 8) {
  sig.terms <- term.table %>% filter(group == "PDS" & Adjusted.P.value < pval_threshold) %>% pull(Term)
  tab <- term.table %>% filter(Term %in% sig.terms)
  p1 <- ggplot(tab, aes(x = group, y = Term_order, size = Overlap_Ratio)) +
    geom_point(aes(color = color_value_log), na.rm = TRUE) +
    scale_color_gradient(low = "blue", high = "red", na.value = "gray",
      name = "Combined Score", labels = scales::label_log()) +
    scale_size_continuous(range = c(0, 10)) +
    scale_y_discrete(limits = rev) +
    theme_minimal() +
    labs(title = title.var, x = "Treatment group (Day 4)", y = NULL)
  ggsave(save.path, p1, width = width.var, height = height.var)
  invisible(tab)
}

run_celltype_go_dotplot <- function(deg.table, celltype, GO.terms, KEGG.terms, save.path) {
  log_msg("Enrichr GO/KEGG dot plot: ", celltype)
  enrichR::setEnrichrSite(site = "Enrichr")
  up.genes <- deg.table %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0) %>% rownames()
  dn.genes <- deg.table %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC < 0) %>% rownames()
  up <- enrichR::enrichr(up.genes, dbs_use)
  dn <- enrichR::enrichr(dn.genes, dbs_use)
  GO.rx <- paste0("^", GO.terms, " ")
  KEGG.rx <- paste0("^", KEGG.terms, "$")
  tab <- process_data(
    filter_terms(up$GO_Biological_Process_2021, GO.rx, "PDS"),
    filter_terms(up$KEGG_2019_Mouse, KEGG.rx, "PDS", prefix = "KEGG: "),
    filter_terms(dn$GO_Biological_Process_2021, GO.rx, "Vehicle"),
    filter_terms(dn$KEGG_2019_Mouse, KEGG.rx, "Vehicle", prefix = "KEGG: "))
  plot_pds_sig(tab, title.var = paste0(celltype, " PDS vs. Vehicle"), save.path = save.path)
  write.csv(tab %>% dplyr::select(-Term_order), sub("\\.pdf$", "_terms.csv", save.path), row.names = FALSE)
}


## ---- Ordered pseudobulk heatmap (blue-white-red2) ----
plot_pseudobulk_heatmap <- function(seu, features, group.by, col.order, title) {
  pseudo <- AggregateExpression(seu, assays = "SCT", return.seurat = TRUE, group.by = group.by)
  features.present <- intersect(features, rownames(pseudo))
  col.order <- col.order[col.order %in% colnames(pseudo)]
  Idents(pseudo) <- factor(Idents(pseudo), levels = col.order)
  pseudo <- pseudo[, col.order]
  p <- DoHeatmap(pseudo, features = features.present, draw.lines = FALSE, angle = 90, raster = FALSE) +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red2") +
    ggplot2::labs(title = title)
  list(plot = p, n.features = length(features.present), n.total = length(features))
}


## ---- Module-score violin with Holm-adjusted Wilcoxon annotation ----
## Two-sided Wilcoxon rank-sum test (ggpubr::compare_means) of per-cell module
## scores, Vehicle vs PDS within each cell type; Holm correction across the
## cell types of one panel. Each cell is treated as an independent observation
## (one pooled 10x library per treatment arm). Non-significant comparisons
## (Holm p >= 0.05) are not drawn. Returns list(plot, stats) where stats holds
## the exact values for every comparison, including non-significant ones.
format_p_label <- function(p) {
  ifelse(p < 1e-4, "p < 0.0001", paste0("p = ", formatC(signif(p, 2), format = "fg", digits = 2)))
}

plot_violin_holm <- function(seu, feature, group.by, split.by, ylab, title, panel,
                             step.increase = 0.12, tip.length = 0.03) {
  p <- VlnPlot2(seu, features = feature, group.by = group.by, split.by = split.by,
    stat.method = "none", pt = FALSE, hide.outlier = TRUE) +
    labs(title = title, y = ylab, x = "Cell Type")

  scores <- p$data
  group_by_arg <- if ("f2" %in% colnames(scores)) c("feature", "f") else "feature"
  formula <- if ("f2" %in% colnames(scores)) value ~ f2 else value ~ f
  stat.test <- ggpubr::compare_means(formula, data = scores, method = "wilcox.test",
    group.by = group_by_arg, p.adjust.method = "holm")

  # exact statistics for every comparison (written out by the caller)
  md <- seu@meta.data
  stats <- bind_rows(lapply(levels(md[[group.by]]), function(ct) {
    v <- md[md[[group.by]] == ct & md[[split.by]] == "Vehicle", feature]
    d <- md[md[[group.by]] == ct & md[[split.by]] == "PDS", feature]
    w <- wilcox.test(d, v)
    data.frame(panel = panel, score = feature, cell_type = ct,
      n_cells_vehicle = length(v), n_cells_PDS = length(d),
      median_vehicle = median(v), median_PDS = median(d),
      W_PDS_vs_vehicle = unname(w$statistic),
      AUC_PDS_gt_vehicle = unname(w$statistic) / (length(v) * length(d)),
      p_raw = w$p.value)
  }))
  stats$p_holm <- p.adjust(stats$p_raw, method = "holm")
  stats$label <- ifelse(stats$p_holm < 0.05, format_p_label(stats$p_holm), "ns (not shown)")

  stat.test <- dplyr::filter(stat.test, p.adj < 0.05)
  if (nrow(stat.test) == 0) return(list(plot = p, stats = stats))

  stat.test$groups <- apply(stat.test, 1, function(x) c(x[["group1"]], x[["group2"]]), simplify = FALSE)
  stat.test <- SeuratExtend:::vlnplot2_Stat_add_y(stat.test, scores = scores, step.increase = step.increase)
  stat.test$annotation <- format_p_label(stat.test$p.adj)
  list(plot = p + ggpubr::stat_pvalue_manual(stat.test, label = "annotation", tip.length = tip.length),
       stats = stats)
}
