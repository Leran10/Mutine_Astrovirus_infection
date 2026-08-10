library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)

output.dir <- "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Harshad/annotation/T cell update CD45 only"

# Load the saved CD45-only T cell object
T.cell <- LoadH5Seurat(file.path(output.dir, "combined.all.final.T.cell.CD45only.h5Seurat"),
                        assays = c("RNA", "SCT"),
                        reductions = c("pca", "tsne", "umap"),
                        graphs = FALSE, images = FALSE, verbose = TRUE)

dir.create(file.path(output.dir, "Bubble plots"), showWarnings = FALSE)

# Use resolution 0.5 clusters
Idents(T.cell) <- "SCT_snn_res.0.5"

cat("Finding markers for resolution 0.5...\n")
DE.markers <- FindAllMarkers(T.cell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(DE.markers, file = file.path(output.dir, "SCT_snn_res.0.5", "DE.0.5.csv"), row.names = FALSE)

top20 <- DE.markers %>%
  group_by(cluster) %>%
  top_n(20, avg_log2FC) %>%
  ungroup()

T.cell <- ScaleData(T.cell, features = unique(top20$gene), assay = "RNA")

p.bubble <- DotPlot(T.cell, features = rev(unique(top20$gene)), assay = "RNA") +
            coord_flip() +
            ggtitle("Top 20 DE genes per cluster - Resolution 0.5 (CD45 only)") +
            theme(axis.text.y = element_text(size = 6))

n.genes <- length(unique(top20$gene))
plot.height <- max(15, n.genes * 0.18)

ggsave(p.bubble, file = file.path(output.dir, "Bubble plots", "top20_gene_expression_0.5.pdf"),
       width = 10, height = plot.height, dpi = 600)

cat("Saved bubble plot:", n.genes, "unique genes across", length(unique(top20$cluster)), "clusters\n")
cat("Done!\n")
