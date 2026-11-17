[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24311376.svg)](https://doi.org/10.5281/zenodo.24311376)

# G-quadruplex stabilization by pyridostatin augments the response to anti-PD-1 therapy by modulating the tumor immune microenvironment in mice

### Accepted in principle at Nature Communications — publication link to follow

This repository holds the **single-cell RNA-seq analysis code** behind the scRNA-seq panels of the paper (Figure 5c–g, Supplementary Figure 2a–h, Supplementary Figure 5b–e).

## Authors & Affiliations

**Nils Kosiol**<sup>1,$</sup>, **Melanie Kienzl**<sup>2,$</sup>, **Hendrik Gerritzen**<sup>1,2</sup>, **Dillon Corvino**<sup>3</sup>, **Saul Rodriguez Roque**<sup>4</sup>, **Rebecca Linke**<sup>1,4</sup>, **Chrystel Flores**<sup>1</sup>, **Ann-Kathrin Baumgart**<sup>5</sup>, **Franziska Schneppenheim**<sup>3</sup>, **Sandra Köll**<sup>2</sup>, **Oliver Schanz**<sup>1</sup>, **Peter Brossart**<sup>1</sup>, **Christian Kurts**<sup>5</sup>, **Julia Kargl**<sup>6,7</sup>, **Tobias Bald**<sup>3</sup>, **Katrin Paeschke**<sup>4,\*</sup>, **Annkristin Heine**<sup>1,2,\*,#</sup>

<sup>1</sup> Medical Clinic III for Oncology, Hematology, Rheumatology and Immune-Oncology, CIO ABCD, University Hospital Bonn, Bonn, Germany  
<sup>2</sup> Clinical Division of Hematology, Department of Internal Medicine, Medical University of Graz, Graz, Austria  
<sup>3</sup> Institute of Experimental Oncology, University Hospital Bonn, Bonn, Germany  
<sup>4</sup> Institute of Clinical Chemistry and Clinical Pharmacology, University Hospital Bonn, Bonn, Germany  
<sup>5</sup> Institute of Molecular Medicine and Experimental Immunology, University Hospital Bonn, Bonn, Germany  
<sup>6</sup> Division of Pharmacology, Otto Loewi Research Center, Medical University of Graz, Graz, Austria  
<sup>7</sup> BioTechMed-Graz, Graz, Austria  

<sup>$</sup> NK and MK contributed equally · <sup>\*</sup> KP and AH contributed equally · <sup>#</sup> Corresponding author

## Contact Information

**Correspondence:** Annkristin Heine — [annkristin.heine@medunigraz.at](mailto:annkristin.heine@medunigraz.at) / [annkristin.heine@ukbonn.de](mailto:annkristin.heine@ukbonn.de)  
**Code:** Dillon Corvino — [dcorvino@uni-bonn.de](mailto:dcorvino@uni-bonn.de)

---

## Summary

Immune checkpoint inhibitors fail in many patients because effector T cells do not infiltrate the tumour. This study tests whether stabilising G-quadruplex (G4) DNA structures with the G4 ligand **pyridostatin (PDS)** can turn checkpoint-resistant tumours into responsive ones. In mouse models, intratumoral PDS slows tumour growth and prolongs survival in a CD8⁺ T cell-dependent manner. It increases tumour-cell antigen presentation and immunogenic cell death, drives cytokine and chemokine production by intratumoral myeloid cells, and boosts the cytotoxicity of tumour-infiltrating T and NK cells. PDS also sensitises tumours to anti-PD-1 therapy.

## Scope of this repository

Single-cell RNA-seq (10x Genomics 5′ gene expression) of **CD45⁺ and CD45⁻ cells sorted from B16OVA tumours 4 days after intratumoral Vehicle or PDS treatment**. Each library contains cells pooled from five tumours, one library per treatment and compartment. Bulk RNA-seq, flow cytometry, cytokine and survival data belong to the co-authors and are not part of this repository.

### Manuscript panels produced here

| Figure panel | Content | Script | Output (`results/…`) |
|---|---|---|---|
| Fig 5c, 5d | UMAP of CD45⁺ cells, Vehicle and PDS (7,000 cells per arm shown) | `scripts/figures/Figure_5.R` | `Figure_5/Figure_5c_UMAP_Vehicle.pdf`, `Figure_5d_UMAP_PDS.pdf` |
| Fig 5e | Pseudobulk heatmap of pro-inflammatory cytokines/chemokines per cell type × treatment | `scripts/figures/Figure_5.R` | `Figure_5/Figure_5e_heatmap_cytokines_chemokines.pdf` |
| Fig 5f, 5g | GO/KEGG enrichment of PDS-vs-Vehicle DEGs in monocytes (f) and neutrophils (g) | `scripts/figures/Figure_5.R` | `Figure_5/Figure_5f_…`, `Figure_5g_…` (+ DEG and term tables) |
| Supp Fig 2a, 2b | ISG pseudobulk heatmaps, myeloid cells and CD45⁻ cells | `scripts/figures/Supplementary_Figure_2.R` | `Supplementary_Figure_2/…2a…`, `…2b…` |
| Supp Fig 2c, 2d | ISG module-score violins with Holm-adjusted Wilcoxon tests | `scripts/figures/Supplementary_Figure_2.R` | `…2c…`, `…2d…` |
| Supp Fig 2e, 2f | cGAS/STING pseudobulk heatmaps | `scripts/figures/Supplementary_Figure_2.R` | `…2e…`, `…2f…` |
| Supp Fig 2g, 2h | cGAS/STING module-score violins | `scripts/figures/Supplementary_Figure_2.R` | `…2g…`, `…2h…` |
| Supp Fig 2 (stats) | Exact n, medians, W, raw and Holm-adjusted p for every comparison in 2c/d/g/h | `scripts/figures/Supplementary_Figure_2.R` | `Supplementary_Figure_2/Supplementary_Figure_2_module_score_statistics.csv` |
| Supp Fig 5b | Cell-type composition of CD45⁺ cells per treatment (table; plotted in Prism) | `scripts/figures/Supplementary_Figure_5.R` | `Supplementary_Figure_5/Supplementary_Figure_5b_celltype_frequency_day4.csv` |
| Supp Fig 5c | Cytokine/chemokine pseudobulk heatmap across CD45⁻ clusters | `scripts/figures/Supplementary_Figure_5.R` | `Supplementary_Figure_5/Supplementary_Figure_5c_…` |
| Supp Fig 5d, 5e | GO/KEGG enrichment of PDS-vs-Vehicle DEGs in macrophages (d) and DCs (e) | `scripts/figures/Supplementary_Figure_5.R` | `Supplementary_Figure_5/Supplementary_Figure_5d_…`, `…5e…` |

## Repository structure

```
scripts/
├── 01_preprocessing_QC_integration.Rmd   # Cell Ranger output -> SoupX, QC, SCTransform, CCA integration, clustering
├── 02_CD45Pos_data_annotation.Rmd        # CD45+ annotation (SingleR + markers) -> CD45pos_seurat.rds
├── 03_CD45Neg_data_annotation.Rmd        # CD45- clean-up, annotation, SCEVAN tumour calling -> CD45neg_seurat.rds
├── Setup.R                               # shared set-up sourced by 01-03
├── variables/Colour_scheme_variable.R
├── figures/
│   ├── _shared.R                         # paths, package loading, plotting/statistics helpers
│   ├── Figure_5.R                        # Fig 5c-g
│   ├── Supplementary_Figure_2.R          # Supp Fig 2a-h (+ exact statistics table)
│   ├── Supplementary_Figure_5.R          # Supp Fig 5b-e
│   └── run_all.R                         # runs the three figure scripts
└── archive/                              # exploratory and superseded code, not used for any panel
```

`data/`, `saves/` and `results/` are not version-controlled (see Data availability).

## Reproducing the figures

1. Download the annotated Seurat objects `CD45pos_seurat.rds` and `CD45neg_seurat.rds` from Zenodo (doi:[10.5281/zenodo.21391578](https://doi.org/10.5281/zenodo.21391578)) into `saves/`. Alternatively, point the environment variable `KOSIOL_SEURAT_DIR` at the folder holding them.
2. From the repository root run:

```bash
Rscript scripts/figures/run_all.R            # all panels
Rscript scripts/figures/Supplementary_Figure_2.R   # or a single figure
```

Panels 5f–g and Supp 5d–e query the Enrichr web service and need internet access. To rebuild the objects from Cell Ranger output, run scripts 01 → 02 → 03 in order.

**Software** (versions used for the manuscript panels): R 4.5, Seurat 5.5.0, SeuratObject 5.4.0, SeuratExtend 1.2.4, ggpubr 0.6.3, scCustomize 3.x, enrichR 3.4, scico 1.5.0, ggplot2; upstream: Cell Ranger 7.1.0, SoupX 1.6.2, clustree 0.5.1, SingleR/celldex, SCEVAN.

## Methods (single-cell analysis)

Reads were aligned with Cell Ranger 7.1.0 to mm10 (refdata-gex-mm10-2020-A), and ambient RNA was removed with SoupX 1.6.2. Cells were retained with 200–5,500 detected features, 1,000–45,000 UMIs, <10% mitochondrial reads and log10(genes)/UMI > 0.8 (CD45⁺) or > 0.75 (CD45⁻). Data were normalised with SCTransform (v2, regressing out mitochondrial content), integrated by canonical correlation analysis and clustered at resolution 0.4, chosen with clustree. Clusters were annotated by canonical markers. Two day-4 basophils were excluded from day-4 analyses.

**Differential expression.** PDS vs Vehicle within each cell type used `FindMarkers()` (Wilcoxon rank-sum). Genes with adjusted p < 0.05 were split by sign of log2FC and tested for over-representation with enrichR against GO Biological Process 2021 and KEGG 2019 Mouse. A curated list of terms was plotted.

**Heatmaps.** Expression was pseudobulked with `AggregateExpression()` (SCT assay) and plotted with `DoHeatmap()`.

**Module scores.** ISG and cGAS/STING scores were computed per cell with `AddModuleScore()` on curated gene panels intersected with detected genes. Vehicle and PDS were compared within each cell type with a two-sided Wilcoxon rank-sum test (`ggpubr::compare_means()`), with Holm correction across cell types within each panel. Each cell is treated as an independent observation (one pooled library per treatment arm).

## Data availability

| Resource | Accession / DOI | Content |
|---|---|---|
| ENA | [PRJEB97458](https://www.ebi.ac.uk/ena/browser/view/PRJEB97458) | Raw scRNA-seq reads (analysed here; Figure 5) |
| ENA | [PRJEB97406](https://www.ebi.ac.uk/ena/browser/view/PRJEB97406) | Bulk RNA-seq, in vitro CD8⁺ T cells (Figure 4; not analysed here) |
| ENA | [PRJEB97321](https://www.ebi.ac.uk/ena/browser/view/PRJEB97321) | Bulk RNA-seq, tumour cells day 9 (Figure 2; not analysed here) |
| Zenodo | [10.5281/zenodo.21391578](https://doi.org/10.5281/zenodo.21391578) | Cell Ranger output and processed Seurat objects |

## Citation

If you use this code, please cite the paper (reference to follow on publication) and the archived code:

> Corvino D. *Code for: G-quadruplex stabilization by pyridostatin augments the response to anti-PD-1 therapy by modulating the tumor immune microenvironment in mice* (single-cell RNA-seq analysis). GitHub repository Eomesodermin/Kosiol_Kienzl_etal_2026_Nature_Coms, v1.0.0. Zenodo, 2026. https://doi.org/10.5281/zenodo.24311376

The DOI above refers to v1.0.0, the version used for the manuscript; [10.5281/zenodo.24311375](https://doi.org/10.5281/zenodo.24311375) always resolves to the latest release.

---
Analysis code by **Dillon Corvino** · [GitHub](https://github.com/Eomesodermin) · [dilloncorvino.com](https://dilloncorvino.com)
