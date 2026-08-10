library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)

output.dir <- "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Harshad/annotation/T cell update CD45 only"

T.cell <- LoadH5Seurat(file.path(output.dir, "combined.all.final.T.cell.CD45only.h5Seurat"),
                        assays = c("RNA", "SCT"),
                        reductions = c("pca", "tsne", "umap"),
                        graphs = FALSE, images = FALSE, verbose = TRUE)

marker.genes <- c("Cd3e", "Cd4", "Cd8a", "Cd8b1")

# ============================================================
# 1. Re-cluster at res 0.3 and 0.8
# ============================================================
cat("1. Re-clustering at res 0.3 and 0.8\n")

DefaultAssay(T.cell) <- "SCT"
T.cell <- FindNeighbors(T.cell, reduction = "pca", dims = 1:50)
T.cell <- FindClusters(T.cell, resolution = c(0.3, 0.5, 0.8))

for (res in c(0.3, 0.8)) {
  res.col <- paste0("SCT_snn_res.", res)
  res.dir <- file.path(output.dir, res.col)
  dir.create(res.dir, showWarnings = FALSE, recursive = TRUE)

  Idents(T.cell) <- res.col
  cat("\nClusters at resolution", res, ":\n")
  print(table(Idents(T.cell)))
}

# ============================================================
# 2. Per-resolution plots: tSNE, UMAP, dotplot, vlnplot, bubble
# ============================================================

for (res in c(0.3, 0.5, 0.8)) {
  res.col <- paste0("SCT_snn_res.", res)
  res.dir <- file.path(output.dir, res.col)
  dir.create(res.dir, showWarnings = FALSE, recursive = TRUE)

  Idents(T.cell) <- res.col
  cat("\n--- Resolution", res, "---\n")

  # tSNE + UMAP cluster plots
  p.tsne <- DimPlot(T.cell, label = TRUE, reduction = "tsne") + ggtitle(paste0("Clusters - Resolution ", res))
  p.umap <- DimPlot(T.cell, label = TRUE, reduction = "umap") + ggtitle(paste0("Clusters - Resolution ", res))
  ggsave(p.tsne, file = file.path(res.dir, "p.tsne.pdf"), width = 6, height = 5, dpi = 600)
  ggsave(p.umap, file = file.path(res.dir, "p.umap.pdf"), width = 6, height = 5, dpi = 600)

  # Copies in tSNE/ and UMAP/ folders
  ggsave(p.tsne, file = file.path(output.dir, "tSNE", paste0("p.tsne.", res, ".pdf")), width = 6, height = 5, dpi = 600)
  ggsave(p.umap, file = file.path(output.dir, "UMAP", paste0("p.umap.", res, ".pdf")), width = 6, height = 5, dpi = 600)

  # Cd4/Cd8 dotplot
  p.dot <- DotPlot(T.cell, features = marker.genes, assay = "RNA") +
           ggtitle(paste0("Cd4/Cd8 markers - Resolution ", res)) +
           theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(p.dot, file = file.path(res.dir, paste0("Cd4_Cd8_dotplot_", res, ".pdf")), width = 7, height = 5, dpi = 600)
  ggsave(p.dot, file = file.path(output.dir, "DotPlots", paste0("Cd4_Cd8_dotplot_", res, ".pdf")), width = 7, height = 5, dpi = 600)

  # Cd4/Cd8 vlnplot
  p.vln <- VlnPlot(T.cell, features = marker.genes, assay = "RNA", ncol = 4)
  ggsave(p.vln, file = file.path(res.dir, paste0("Cd4_Cd8_vlnplot_", res, ".pdf")), width = 16, height = 5, dpi = 600)
  ggsave(p.vln, file = file.path(output.dir, "Violin plots", paste0("Cd4_Cd8_vlnplot_", res, ".pdf")), width = 16, height = 5, dpi = 600)

  # Bubble plot (FindAllMarkers + top 20)
  cat("  FindAllMarkers...\n")
  DE.markers <- FindAllMarkers(T.cell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(DE.markers, file = file.path(res.dir, paste0("DE.", res, ".csv")), row.names = FALSE)

  top20 <- DE.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC) %>% ungroup()
  T.cell <- ScaleData(T.cell, features = unique(top20$gene), assay = "RNA")

  p.bubble <- DotPlot(T.cell, features = rev(unique(top20$gene)), assay = "RNA") +
              coord_flip() +
              ggtitle(paste0("Top 20 DE genes per cluster - Resolution ", res, " (CD45 only)")) +
              theme(axis.text.y = element_text(size = 6))

  n.genes <- length(unique(top20$gene))
  plot.height <- max(15, n.genes * 0.18)

  ggsave(p.bubble, file = file.path(res.dir, paste0("top20_gene_expression_", res, ".pdf")),
         width = 10, height = plot.height, dpi = 600)
  ggsave(p.bubble, file = file.path(output.dir, "Bubble plots", paste0("top20_gene_expression_", res, ".pdf")),
         width = 10, height = plot.height, dpi = 600)

  cat("  Done.\n")
}

# ============================================================
# 3. Bubble plot - 11 clusters at res 0.5 (dedicated version)
# ============================================================
cat("\n3. Bubble plot for 11 clusters at res 0.5\n")

Idents(T.cell) <- "SCT_snn_res.0.5"
# Already have DE.0.5 from the main analysis, but re-do for the specific "11clusters" label
DE.markers.0.5 <- FindAllMarkers(T.cell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DE.markers.0.5, file = file.path(output.dir, "SCT_snn_res.0.5", "DE.0.5.all_clusters.csv"), row.names = FALSE)

top20.0.5 <- DE.markers.0.5 %>% group_by(cluster) %>% top_n(20, avg_log2FC) %>% ungroup()
T.cell <- ScaleData(T.cell, features = unique(top20.0.5$gene), assay = "RNA")

p.bubble.11 <- DotPlot(T.cell, features = rev(unique(top20.0.5$gene)), assay = "RNA") +
               coord_flip() +
               ggtitle("Top 20 DE genes per cluster - 11 clusters res 0.5 (CD45 only)") +
               theme(axis.text.y = element_text(size = 6))

n.genes.11 <- length(unique(top20.0.5$gene))
plot.height.11 <- max(15, n.genes.11 * 0.18)

ggsave(p.bubble.11, file = file.path(output.dir, "Bubble plots", "top20_gene_expression_11clusters_0.5.pdf"),
       width = 10, height = plot.height.11, dpi = 600)

# Also copy bubble plot to SCT_snn_res.0.5
ggsave(p.bubble.11, file = file.path(output.dir, "SCT_snn_res.0.5", "top20_gene_expression_0.5.pdf"),
       width = 10, height = plot.height.11, dpi = 600)

# ============================================================
# 4. Root-level featureplot copies
# ============================================================
cat("\n4. Root-level featureplot copies\n")

for (gene in marker.genes) {
  # tSNE
  p <- FeaturePlot(T.cell, features = gene, cols = c("#F9F4F4", "#0F04F7"), reduction = "tsne", order = TRUE) +
       ggtitle(paste0(gene, " - tSNE"))
  ggsave(p, file = file.path(output.dir, paste0(gene, ".featureplot.tsne.pdf")), width = 6, height = 5, dpi = 600)

  p.split.tsne <- FeaturePlot(T.cell, features = gene, cols = c("#F9F4F4", "#0F04F7"),
                               split.by = "modified.ident", reduction = "tsne", order = TRUE)
  ggsave(p.split.tsne, file = file.path(output.dir, paste0(gene, ".featureplot.tsne.perSample.pdf")),
         width = 10, height = 4, dpi = 600)

  # UMAP
  p.umap <- FeaturePlot(T.cell, features = gene, cols = c("#F9F4F4", "#0F04F7"), reduction = "umap", order = TRUE) +
            ggtitle(paste0(gene, " - UMAP"))
  ggsave(p.umap, file = file.path(output.dir, paste0(gene, ".featureplot.umap.pdf")), width = 6, height = 5, dpi = 600)

  # Cd3e also gets UMAP perSample at root
  if (gene == "Cd3e") {
    p.split.umap <- FeaturePlot(T.cell, features = gene, cols = c("#F9F4F4", "#0F04F7"),
                                 split.by = "modified.ident", reduction = "umap", order = TRUE)
    ggsave(p.split.umap, file = file.path(output.dir, paste0(gene, ".featureplot.umap.perSample.pdf")),
           width = 10, height = 4, dpi = 600)
  }
}

# ============================================================
# 5. Subtype markers vlnplot copy to Violin plots/
# ============================================================
cat("5. Subtype markers vlnplot copy\n")

Idents(T.cell) <- "Tcell_subtype"
p.vln.sub <- VlnPlot(T.cell, features = marker.genes, assay = "RNA", ncol = 4)
ggsave(p.vln.sub, file = file.path(output.dir, "Violin plots", "Tcell_subtype_markers_vlnplot.pdf"),
       width = 16, height = 5, dpi = 600)

# ============================================================
# 6. Save updated object with all resolutions
# ============================================================
cat("6. Saving updated object\n")
SaveH5Seurat(T.cell, filename = file.path(output.dir, "combined.all.final.T.cell.CD45only.h5Seurat"), overwrite = TRUE)

cat("\nDone! All remaining plots generated.\n")
