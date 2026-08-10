library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(ggrepel)

output.dir <- "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Harshad/annotation/T cell update CD45 only"

T.cell <- LoadH5Seurat(file.path(output.dir, "combined.all.final.T.cell.CD45only.h5Seurat"),
                        assays = c("RNA", "SCT"),
                        reductions = c("pca", "tsne", "umap"),
                        graphs = FALSE, images = FALSE, verbose = TRUE)

# ============================================================
# 1. FIX: Remove Cd3e from subtype marker dotplot and vlnplot
#    Only keep Cd4, Cd8a, Cd8b1
# ============================================================
cat("1. Fixing subtype marker dotplot and vlnplot (removing Cd3e)\n")

marker.genes <- c("Cd4", "Cd8a", "Cd8b1")
Idents(T.cell) <- "Tcell_subtype"

p.dot <- DotPlot(T.cell, features = marker.genes, assay = "RNA") +
         ggtitle("Key markers by T cell subtype (CD45 only)") +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.dot, file = file.path(output.dir, "Tcell_subtype_markers_dotplot.pdf"), width = 7, height = 5, dpi = 600)
ggsave(p.dot, file = file.path(output.dir, "DotPlots", "Tcell_subtype_markers_dotplot.pdf"), width = 7, height = 5, dpi = 600)

p.vln <- VlnPlot(T.cell, features = marker.genes, assay = "RNA", ncol = 3)
ggsave(p.vln, file = file.path(output.dir, "Tcell_subtype_markers_vlnplot.pdf"), width = 12, height = 5, dpi = 600)
ggsave(p.vln, file = file.path(output.dir, "Violin plots", "Tcell_subtype_markers_vlnplot.pdf"), width = 12, height = 5, dpi = 600)

# ============================================================
# 2. FIX: Flip volcano plots so WT naive on left, WT+FFT on right
#    Negate avg_log2FC so positive = up in WT+FFT
# ============================================================
cat("2. Flipping volcano plots\n")

subtypes <- c("CD8aa T cell", "CD8ab T cell", "CD4+ T cell", "DN T cell")

for (subtype in subtypes) {
  subtype.clean <- gsub("[+ /]", "_", subtype)
  de.file <- file.path(output.dir, "DE results",
                       paste0("DE_genomewide_naive_vs_fft_", subtype.clean, ".csv"))

  if (!file.exists(de.file)) {
    cat("  Skipping", subtype, "- DE file not found\n")
    next
  }

  de.results <- read.csv(de.file)

  # Flip: negate log2FC so positive = up in WT+FFT
  de.results$flipped_log2FC <- -de.results$avg_log2FC

  de.results$significance <- "Not significant"
  de.results$significance[de.results$p_val_adj < 0.05 & de.results$flipped_log2FC > 0] <- "Up in WT + FFT"
  de.results$significance[de.results$p_val_adj < 0.05 & de.results$flipped_log2FC < 0] <- "Up in WT naive"
  de.results$significance <- factor(de.results$significance,
                                     levels = c("Up in WT naive", "Up in WT + FFT", "Not significant"))

  sig.genes <- de.results %>% filter(p_val_adj < 0.05) %>% arrange(p_val_adj) %>% head(20)
  n.sig <- sum(de.results$p_val_adj < 0.05, na.rm = TRUE)
  n.up.naive <- sum(de.results$p_val_adj < 0.05 & de.results$flipped_log2FC < 0, na.rm = TRUE)
  n.up.fft <- sum(de.results$p_val_adj < 0.05 & de.results$flipped_log2FC > 0, na.rm = TRUE)

  p <- ggplot(de.results, aes(x = flipped_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up in WT naive" = "#E41A1C",
                                   "Up in WT + FFT" = "#377EB8",
                                   "Not significant" = "grey70")) +
    geom_text_repel(data = sig.genes, aes(x = flipped_log2FC, label = gene),
                    size = 3, max.overlaps = 20, color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggtitle(paste0(subtype, " - WT naive vs WT + FFT (CD45 only)"),
            subtitle = paste0(n.sig, " significant genes (", n.up.naive, " up in naive, ", n.up.fft, " up in FFT)")) +
    xlab("log2 Fold Change\n<-- WT naive | WT + FFT -->") +
    ylab("-log10(adjusted p-value)") +
    theme_bw() +
    theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10))

  ggsave(p, file = file.path(output.dir, "Volcano plots", paste0("volcano_", subtype.clean, ".pdf")),
         width = 8, height = 7, dpi = 600)

  cat("  Saved flipped volcano for", subtype, "\n")
}

cat("\nDone!\n")
