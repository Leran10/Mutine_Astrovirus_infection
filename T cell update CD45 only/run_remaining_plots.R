library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(patchwork)

output.dir <- "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Harshad/annotation/T cell update CD45 only"

# Load the saved CD45-only T cell object
T.cell <- LoadH5Seurat(file.path(output.dir, "combined.all.final.T.cell.CD45only.h5Seurat"),
                        assays = c("RNA", "SCT"),
                        reductions = c("pca", "tsne", "umap"),
                        graphs = FALSE, images = FALSE, verbose = TRUE)

# Create directories
dir.create(file.path(output.dir, "Violin plots"), showWarnings = FALSE)
dir.create(file.path(output.dir, "DotPlots"), showWarnings = FALSE)
dir.create(file.path(output.dir, "UMAP"), showWarnings = FALSE)
dir.create(file.path(output.dir, "tSNE"), showWarnings = FALSE)
dir.create(file.path(output.dir, "Cell counts"), showWarnings = FALSE)

subtype.colors <- c("CD4+ T cell"  = "#E41A1C",
                     "CD8aa T cell" = "#377EB8",
                     "CD8ab T cell" = "#4DAF4A",
                     "Mixed CD4/CD8" = "#FF7F00",
                     "DN T cell"    = "#984EA3")

genes.of.interest <- c("Cd4", "Cd8a", "Cd8b1", "Ifng", "Gzmb", "Sell", "Cd44", "Cxcr3", "Klrg1", "Pdcd1")
marker.genes <- c("Cd3e", "Cd4", "Cd8a", "Cd8b1")

# ============================================================
# 1. Cd4/Cd8 DotPlot and VlnPlot at res 0.5
# ============================================================
cat("1. Cd4/Cd8 dotplot and vlnplot at res 0.5\n")

Idents(T.cell) <- "SCT_snn_res.0.5"

p.dot <- DotPlot(T.cell, features = marker.genes, assay = "RNA") +
         ggtitle("Cd4/Cd8 markers per cluster (res 0.5, CD45 only)") +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.dot, file = file.path(output.dir, "DotPlots", "Cd4_Cd8_dotplot_0.5.pdf"), width = 7, height = 5, dpi = 600)
ggsave(p.dot, file = file.path(output.dir, "SCT_snn_res.0.5", "Cd4_Cd8_dotplot_0.5.pdf"), width = 7, height = 5, dpi = 600)

p.vln <- VlnPlot(T.cell, features = marker.genes, assay = "RNA", ncol = 4)
ggsave(p.vln, file = file.path(output.dir, "Violin plots", "Cd4_Cd8_vlnplot_0.5.pdf"), width = 16, height = 5, dpi = 600)
ggsave(p.vln, file = file.path(output.dir, "SCT_snn_res.0.5", "Cd4_Cd8_vlnplot_0.5.pdf"), width = 16, height = 5, dpi = 600)

# ============================================================
# 2. DotPlot genes of interest by subtype + sample
# ============================================================
cat("2. DotPlot genes of interest by subtype + sample\n")

Idents(T.cell) <- "Tcell_subtype"
T.cell@meta.data$subtype_sample <- paste0(T.cell@meta.data$Tcell_subtype, " - ", T.cell@meta.data$modified.ident)
Idents(T.cell) <- "subtype_sample"

p.dot2 <- DotPlot(T.cell, features = genes.of.interest, assay = "RNA") +
           coord_flip() +
           ggtitle("Gene expression: WT naive vs WT+FFT by T cell subtype (CD45 only)") +
           theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(p.dot2, file = file.path(output.dir, "DotPlots", "DotPlot_genes_of_interest_by_subtype_sample.pdf"),
       width = 12, height = 7, dpi = 600)

# Subtype markers dotplot (copy to DotPlots folder too)
Idents(T.cell) <- "Tcell_subtype"
p.dot3 <- DotPlot(T.cell, features = marker.genes, assay = "RNA") +
           ggtitle("Key markers by T cell subtype (CD45 only)") +
           theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.dot3, file = file.path(output.dir, "DotPlots", "Tcell_subtype_markers_dotplot.pdf"), width = 7, height = 5, dpi = 600)

# ============================================================
# 3. UMAP plots — subtype combined + perSample
# ============================================================
cat("3. UMAP subtype plots\n")

Idents(T.cell) <- "Tcell_subtype"

p.umap <- DimPlot(T.cell, reduction = "umap", cols = subtype.colors, label = TRUE, repel = TRUE) +
           ggtitle("T cell subtypes - UMAP (CD45 only)")
ggsave(p.umap, file = file.path(output.dir, "UMAP", "Tcell_subtype.umap.pdf"), width = 7, height = 5, dpi = 600)

p.umap.split <- DimPlot(T.cell, reduction = "umap", split.by = "modified.ident",
                          cols = subtype.colors, label = FALSE) +
                ggtitle("T cell subtypes - per sample (CD45 only)")
ggsave(p.umap.split, file = file.path(output.dir, "UMAP", "Tcell_subtype.umap.perSample.pdf"),
       width = 10, height = 5, dpi = 600)

# ============================================================
# 4. FeaturePlots — UMAP versions for Cd3e, Cd4, Cd8a, Cd8b1
# ============================================================
cat("4. UMAP FeaturePlots\n")

for (gene in marker.genes) {
  p <- FeaturePlot(T.cell, features = gene, cols = c("#F9F4F4", "#0F04F7"), reduction = "umap", order = TRUE) +
       ggtitle(paste0(gene, " - UMAP"))
  ggsave(p, file = file.path(output.dir, "FeaturePlots", paste0(gene, ".featureplot.umap.pdf")), width = 6, height = 5, dpi = 600)

  p.split <- FeaturePlot(T.cell, features = gene, cols = c("#F9F4F4", "#0F04F7"),
                          split.by = "modified.ident", reduction = "umap", order = TRUE)
  ggsave(p.split, file = file.path(output.dir, "FeaturePlots", paste0(gene, ".featureplot.umap.perSample.pdf")),
         width = 10, height = 4, dpi = 600)
}

# ============================================================
# 5. VlnPlots — 10 genes by CLUSTER (individual + combined)
# ============================================================
cat("5. VlnPlots by cluster\n")

Idents(T.cell) <- "SCT_snn_res.0.5"

vln.list <- list()
for (gene in genes.of.interest) {
  p <- VlnPlot(T.cell, features = gene, assay = "RNA", split.by = "modified.ident",
               cols = c("#2787cf", "#d24339")) +
       ggtitle(paste0(gene, " - by cluster"))
  ggsave(p, file = file.path(output.dir, "Violin plots", paste0("VlnPlot_cluster_", gene, ".pdf")),
         width = 10, height = 5, dpi = 600)
  vln.list[[gene]] <- p
}

p.combined <- wrap_plots(vln.list, ncol = 2)
ggsave(p.combined, file = file.path(output.dir, "Violin plots", "VlnPlot_cluster_all_genes_combined.pdf"),
       width = 16, height = 25, dpi = 600)

# ============================================================
# 6. VlnPlots — 10 genes by SUBTYPE (individual + combined)
# ============================================================
cat("6. VlnPlots by subtype\n")

Idents(T.cell) <- "Tcell_subtype"

vln.list2 <- list()
for (gene in genes.of.interest) {
  p <- VlnPlot(T.cell, features = gene, assay = "RNA", split.by = "modified.ident",
               cols = c("#2787cf", "#d24339")) +
       ggtitle(paste0(gene, " - by subtype"))
  ggsave(p, file = file.path(output.dir, "Violin plots", paste0("VlnPlot_subtype_", gene, ".pdf")),
         width = 10, height = 5, dpi = 600)
  vln.list2[[gene]] <- p
}

p.combined2 <- wrap_plots(vln.list2, ncol = 2)
ggsave(p.combined2, file = file.path(output.dir, "Violin plots", "VlnPlot_subtype_all_genes_combined.pdf"),
       width = 14, height = 25, dpi = 600)

# ============================================================
# 7. VlnPlots — per subtype group (10 genes each)
# ============================================================
cat("7. VlnPlots per subtype group\n")

subtypes <- c("CD4+ T cell", "CD8aa T cell", "CD8ab T cell", "DN T cell")

for (sub in subtypes) {
  temp <- subset(T.cell, subset = Tcell_subtype == sub)
  if (ncol(temp) < 10) next
  Idents(temp) <- "modified.ident"

  p.vln <- VlnPlot(temp, features = genes.of.interest, assay = "RNA",
                    split.by = "modified.ident", ncol = 5) +
           patchwork::plot_annotation(title = paste0(sub, " (CD45 only)"))

  sub.clean <- gsub("[/+ ]", "_", sub)
  ggsave(p.vln, file = file.path(output.dir, "Violin plots", paste0("VlnPlot_genes_", sub.clean, ".pdf")),
         width = 20, height = 8, dpi = 600)
}

# ============================================================
# 8. Copy cell counts and proportions to Cell counts folder
# ============================================================
cat("8. Copying cell counts\n")

file.copy(file.path(output.dir, "Tcell_subtype_cell_counts.csv"),
          file.path(output.dir, "Cell counts", "Tcell_subtype_cell_counts.csv"), overwrite = TRUE)
file.copy(file.path(output.dir, "Tcell_subtype_proportions_barplot.pdf"),
          file.path(output.dir, "Cell counts", "Tcell_subtype_proportions_barplot.pdf"), overwrite = TRUE)

cat("\nDone! All remaining plots generated.\n")
