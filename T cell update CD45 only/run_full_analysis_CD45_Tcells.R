###############################################################################
# Full T cell analysis — CD45-only subset
#
# This script replicates the T cell update pipeline but uses ONLY T cells
# from the CD45+ clusters (clusters 14,7,3,13,25,0,16,6,10,23), excluding
# 136 T cells that were in Epcam cluster 20.
#
# CD45-only T cells: 6,429 (2,047 WT naive + 4,382 WT+FFT)
# vs original:       6,565 (2,106 WT naive + 4,459 WT+FFT)
###############################################################################

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(data.table)
library(tidyverse)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(fgsea)
library(msigdbr)

# ============================================================
# 1. LOAD DATA & SUBSET TO CD45 T CELLS
# ============================================================

cat("============================================================\n")
cat("1. Loading data and subsetting to CD45 T cells\n")
cat("============================================================\n\n")

combined.all.final <- LoadH5Seurat(
  "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Harshad/annotation/combined.all.annotated.final.h5seurat",
  assays = "RNA",
  reductions = c("pca", "tsne", "umap"),
  graphs = FALSE, images = FALSE, verbose = TRUE)

# Rename sample groups
combined.all.final$modified.ident <- ifelse(combined.all.final$modified.ident == "Wtnaive", "WT naïve",
                                     ifelse(combined.all.final$modified.ident == "WTFFT", "WT + FFT",
                                     ifelse(combined.all.final$modified.ident == "Rag2Il2rg", "Rag2-/-Il2rg-/-",
                                            "Rag2-/-Il2rg-/-Ifnlr1-/-")))
combined.all.final$modified.ident <- factor(combined.all.final$modified.ident, levels = c("WT naïve", "WT + FFT"))

# Subset to CD45 clusters first, THEN filter to T cells
cd45.clusters <- c("14", "7", "3", "13", "25", "0", "16", "6", "10", "23")
CD45.subset <- subset(combined.all.final, seurat_clusters_updated %in% cd45.clusters)
T.cell <- subset(CD45.subset, custom.cell.annotation %like% "T cell")

# Remove Rag2 knockout samples
T.cell <- subset(T.cell, subset = modified.ident %in% c("WT naïve", "WT + FFT"))

cat("CD45-only T cells (WT only):", ncol(T.cell), "\n")
cat("Per sample:\n")
print(table(T.cell$modified.ident))

# Set output directory
output.dir <- "~/Handley Lab Dropbox/leran wang/Baldridge/scRNA/Harshad/annotation/T cell update CD45 only"
dir.create(file.path(output.dir, "SCT_snn_res.0.5"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output.dir, "tSNE"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output.dir, "FeaturePlots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output.dir, "Volcano plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output.dir, "DE results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output.dir, "Pathway analysis"), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. SCTransform + DIMENSIONALITY REDUCTION + CLUSTERING
# ============================================================

cat("\n============================================================\n")
cat("2. SCTransform + dim reduction + clustering\n")
cat("============================================================\n\n")

T.cell <- SCTransform(T.cell, vars.to.regress = "percent.mito", verbose = FALSE)
T.cell <- RunPCA(T.cell, verbose = FALSE, npcs = 100)
T.cell <- RunUMAP(T.cell, dims = 1:50, verbose = FALSE, reduction = "pca", n.neighbors = 5)
T.cell <- RunTSNE(object = T.cell, dims = 1:50, reduction = "pca", perplexity = 50)

DefaultAssay(T.cell) <- "SCT"
T.cell <- FindNeighbors(T.cell, reduction = "pca", dims = 1:50)
T.cell <- FindClusters(T.cell, resolution = 0.5)

cat("Clusters at resolution 0.5:\n")
print(table(T.cell$SCT_snn_res.0.5))

# ============================================================
# 3. MARKER QUANTIFICATION PER CLUSTER
# ============================================================

cat("\n============================================================\n")
cat("3. Marker quantification per cluster\n")
cat("============================================================\n\n")

Idents(T.cell) <- "SCT_snn_res.0.5"
clusters <- levels(Idents(T.cell))
subtype.genes <- c("Cd3e", "Cd4", "Cd8a", "Cd8b1")
marker.expr <- GetAssayData(T.cell, assay = "RNA")

stats <- data.frame(cluster = clusters, n_cells = NA,
                    pct_Cd3e = NA, pct_Cd4 = NA, pct_Cd8a = NA, pct_Cd8b1 = NA,
                    avg_Cd3e = NA, avg_Cd4 = NA, avg_Cd8a = NA, avg_Cd8b1 = NA)

for (j in seq_along(clusters)) {
  cells <- WhichCells(T.cell, idents = clusters[j])
  stats$n_cells[j] <- length(cells)
  for (gene in subtype.genes) {
    expr <- marker.expr[gene, cells]
    stats[[paste0("pct_", gene)]][j] <- round(sum(expr > 0) / length(expr) * 100, 2)
    stats[[paste0("avg_", gene)]][j] <- round(mean(expr), 4)
  }
}

cat("Marker expression per cluster (res 0.5):\n")
print(stats)
write.csv(stats, file = file.path(output.dir, "Cd4_Cd8_expression_per_cluster_0.5.csv"), row.names = FALSE)

# ============================================================
# 4. SUBTYPE LABELING
# ============================================================

cat("\n============================================================\n")
cat("4. Assigning T cell subtype labels\n")
cat("============================================================\n\n")

T.cell@meta.data$Tcell_subtype <- sapply(
  as.character(T.cell@meta.data$SCT_snn_res.0.5),
  function(cl) {
    s <- stats[stats$cluster == cl, ]
    if (s$pct_Cd4 > 30) return("CD4+ T cell")
    if (s$pct_Cd8b1 > 30) return("CD8ab T cell")
    if (s$pct_Cd8a < 30 & s$pct_Cd4 < 10) return("DN T cell")
    if (s$pct_Cd4 >= 10 & s$pct_Cd4 <= 30 & s$pct_Cd8a >= 20) return("Mixed CD4/CD8")
    return("CD8aa T cell")
  }
)

T.cell@meta.data$Tcell_subtype <- factor(
  T.cell@meta.data$Tcell_subtype,
  levels = c("CD4+ T cell", "CD8aa T cell", "CD8ab T cell", "Mixed CD4/CD8", "DN T cell"))

cat("T cell subtype counts:\n")
print(table(T.cell$Tcell_subtype))
cat("\nSubtype counts per sample:\n")
print(table(T.cell$Tcell_subtype, T.cell$modified.ident))

# ============================================================
# 5. PLOTS — CLUSTERS, SUBTYPES, FEATURE PLOTS
# ============================================================

cat("\n============================================================\n")
cat("5. Generating plots\n")
cat("============================================================\n\n")

# --- Cluster tSNE/UMAP ---
p.tsne.cluster <- DimPlot(T.cell, label = TRUE, reduction = "tsne") + ggtitle("Clusters - Resolution 0.5")
p.umap.cluster <- DimPlot(T.cell, label = TRUE, reduction = "umap") + ggtitle("Clusters - Resolution 0.5")
ggsave(p.tsne.cluster, file = file.path(output.dir, "SCT_snn_res.0.5", "p.tsne.pdf"), width = 6, height = 5, dpi = 600)
ggsave(p.umap.cluster, file = file.path(output.dir, "SCT_snn_res.0.5", "p.umap.pdf"), width = 6, height = 5, dpi = 600)

# --- FeaturePlots for markers ---
for (gene in subtype.genes) {
  p <- FeaturePlot(T.cell, features = gene, cols = c("#F9F4F4", "#0F04F7"), reduction = "tsne", order = TRUE) +
       ggtitle(paste0(gene, " - tSNE"))
  ggsave(p, file = file.path(output.dir, "FeaturePlots", paste0(gene, ".featureplot.tsne.pdf")), width = 6, height = 5, dpi = 600)

  p.split <- FeaturePlot(T.cell, features = gene, cols = c("#F9F4F4", "#0F04F7"),
                          split.by = "modified.ident", reduction = "tsne", order = TRUE)
  ggsave(p.split, file = file.path(output.dir, "FeaturePlots", paste0(gene, ".featureplot.tsne.perSample.pdf")),
         width = 10, height = 5, dpi = 600)
}

# --- Subtype tSNE/UMAP ---
Idents(T.cell) <- "Tcell_subtype"

subtype.colors <- c("CD4+ T cell"  = "#E41A1C",
                     "CD8aa T cell" = "#377EB8",
                     "CD8ab T cell" = "#4DAF4A",
                     "Mixed CD4/CD8" = "#FF7F00",
                     "DN T cell"    = "#984EA3")

p.tsne <- DimPlot(T.cell, reduction = "tsne", cols = subtype.colors, label = TRUE, repel = TRUE) +
          ggtitle("T cell subtypes - tSNE (CD45 only)")
ggsave(p.tsne, file = file.path(output.dir, "tSNE", "Tcell_subtype.tsne.pdf"), width = 7, height = 5, dpi = 600)

p.tsne.split <- DimPlot(T.cell, reduction = "tsne", split.by = "modified.ident",
                         cols = subtype.colors, label = FALSE) +
                ggtitle("T cell subtypes - per sample (CD45 only)") +
                coord_fixed()
ggsave(p.tsne.split, file = file.path(output.dir, "tSNE", "Tcell_subtype.tsne.perSample.pdf"),
       width = 10, height = 5, dpi = 600)
ggsave(p.tsne.split, file = file.path(output.dir, "tSNE", "Tcell_subtype.tsne.perSample.svg"),
       width = 10, height = 5)

# --- DotPlot of key markers ---
marker.genes <- c("Cd3e", "Cd4", "Cd8a", "Cd8b1")
p.dot <- DotPlot(T.cell, features = marker.genes, assay = "RNA") +
         ggtitle("Key markers by T cell subtype (CD45 only)") +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.dot, file = file.path(output.dir, "Tcell_subtype_markers_dotplot.pdf"), width = 7, height = 5, dpi = 600)

# --- VlnPlot of key markers ---
p.vln <- VlnPlot(T.cell, features = marker.genes, assay = "RNA", ncol = 4)
ggsave(p.vln, file = file.path(output.dir, "Tcell_subtype_markers_vlnplot.pdf"), width = 16, height = 5, dpi = 600)

# --- Cell counts & proportions barplot ---
count.table <- table(T.cell$Tcell_subtype, T.cell$modified.ident)
write.csv(as.data.frame.matrix(count.table), file = file.path(output.dir, "Tcell_subtype_cell_counts.csv"))

count.df <- as.data.frame(count.table)
colnames(count.df) <- c("Tcell_subtype", "Sample", "Count")

p.bar <- ggplot(count.df, aes(x = Sample, y = Count, fill = Tcell_subtype)) +
         geom_bar(stat = "identity", position = "fill") +
         scale_fill_manual(values = subtype.colors) +
         ylab("Proportion") +
         ggtitle("T cell subtype proportions per sample (CD45 only)") +
         theme_bw() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p.bar, file = file.path(output.dir, "Tcell_subtype_proportions_barplot.pdf"), width = 7, height = 5, dpi = 600)

cat("Plots saved.\n")

# ============================================================
# 6. FeaturePlots for 10 genes of interest
# ============================================================

cat("\n============================================================\n")
cat("6. FeaturePlots for 10 genes of interest\n")
cat("============================================================\n\n")

genes.of.interest <- c("Cd4", "Cd8a", "Cd8b1", "Ifng", "Gzmb", "Sell", "Cd44", "Cxcr3", "Klrg1", "Pdcd1")

for (gene in genes.of.interest) {
  if (!gene %in% rownames(T.cell)) {
    cat("Gene", gene, "not found, skipping\n")
    next
  }
  p.tsne <- FeaturePlot(T.cell, features = gene, cols = c("#F9F4F4", "#0F04F7"),
                         split.by = "modified.ident", reduction = "tsne", order = TRUE)
  ggsave(p.tsne, file = file.path(output.dir, "FeaturePlots", paste0(gene, ".featureplot.tsne.perSample.pdf")),
         width = 10, height = 5, dpi = 600)
}

cat("FeaturePlots saved.\n")

# ============================================================
# 7. DE ANALYSIS — GENOME-WIDE PER SUBTYPE + VOLCANO PLOTS
# ============================================================

cat("\n============================================================\n")
cat("7. DE analysis + volcano plots\n")
cat("============================================================\n\n")

subtypes <- c("CD8aa T cell", "CD8ab T cell", "CD4+ T cell", "DN T cell")

for (subtype in subtypes) {
  cat("Processing", subtype, "...\n")

  cells <- WhichCells(T.cell, expression = Tcell_subtype == subtype)
  sub.obj <- subset(T.cell, cells = cells)
  Idents(sub.obj) <- "modified.ident"

  if (length(unique(Idents(sub.obj))) < 2) {
    cat("  Skipping", subtype, "- only one group present\n")
    next
  }

  de.results <- tryCatch({
    FindMarkers(sub.obj, ident.1 = "WT naïve", ident.2 = "WT + FFT",
                min.pct = 0.1, logfc.threshold = 0, assay = "RNA", test.use = "wilcox")
  }, error = function(e) {
    cat("  Error in DE for", subtype, ":", e$message, "\n")
    return(NULL)
  })

  if (is.null(de.results)) next

  de.results$gene <- rownames(de.results)
  subtype.clean <- gsub("[+ /]", "_", subtype)

  write.csv(de.results, file = file.path(output.dir, "DE results",
            paste0("DE_genomewide_naive_vs_fft_", subtype.clean, ".csv")), row.names = FALSE)

  # Volcano plot
  de.results$significance <- "Not significant"
  de.results$significance[de.results$p_val_adj < 0.05 & de.results$avg_log2FC > 0] <- "Up in WT naive"
  de.results$significance[de.results$p_val_adj < 0.05 & de.results$avg_log2FC < 0] <- "Up in WT + FFT"
  de.results$significance <- factor(de.results$significance,
                                     levels = c("Up in WT naive", "Up in WT + FFT", "Not significant"))

  sig.genes <- de.results %>% filter(p_val_adj < 0.05) %>% arrange(p_val_adj) %>% head(20)
  n.sig <- sum(de.results$p_val_adj < 0.05, na.rm = TRUE)
  n.up.naive <- sum(de.results$p_val_adj < 0.05 & de.results$avg_log2FC > 0, na.rm = TRUE)
  n.up.fft <- sum(de.results$p_val_adj < 0.05 & de.results$avg_log2FC < 0, na.rm = TRUE)

  p <- ggplot(de.results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up in WT naive" = "#E41A1C", "Up in WT + FFT" = "#377EB8",
                                   "Not significant" = "grey70")) +
    geom_text_repel(data = sig.genes, aes(label = gene), size = 3, max.overlaps = 20, color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggtitle(paste0(subtype, " — WT naive vs WT + FFT (CD45 only)"),
            subtitle = paste0(n.sig, " significant genes (", n.up.naive, " up in naive, ", n.up.fft, " up in FFT)")) +
    xlab("log2 Fold Change (positive = higher in WT naive)") +
    ylab("-log10(adjusted p-value)") +
    theme_bw() +
    theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10))

  ggsave(p, file = file.path(output.dir, "Volcano plots", paste0("volcano_", subtype.clean, ".pdf")),
         width = 8, height = 7, dpi = 600)

  cat("  ", n.sig, "significant genes\n")
}

# ============================================================
# 8. 10-GENE DE ANALYSIS + HEATMAP
# ============================================================

cat("\n============================================================\n")
cat("8. 10-gene DE analysis + heatmap\n")
cat("============================================================\n\n")

de.results.10genes <- list()

for (sub in subtypes) {
  temp <- subset(T.cell, subset = Tcell_subtype == sub)
  if (length(unique(temp$modified.ident)) < 2) next
  Idents(temp) <- "modified.ident"

  tryCatch({
    de <- FindMarkers(temp, ident.1 = "WT naïve", ident.2 = "WT + FFT",
                      features = genes.of.interest, only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
    de$gene <- rownames(de)
    de$subtype <- sub
    de.results.10genes[[sub]] <- de
  }, error = function(e) {
    cat("Error in 10-gene DE for", sub, ":", e$message, "\n")
  })
}

de.10genes <- bind_rows(de.results.10genes)
write.csv(de.10genes, file = file.path(output.dir, "DE results",
          "DE_genes_of_interest_naive_vs_fft_by_subtype.csv"), row.names = FALSE)

# Heatmap
de.10genes$subtype <- gsub("CD4\\+ T cell", "CD4+", de.10genes$subtype)
de.10genes$subtype <- gsub("CD8aa T cell", "CD8aa", de.10genes$subtype)
de.10genes$subtype <- gsub("CD8ab T cell", "CD8ab", de.10genes$subtype)
de.10genes$subtype <- gsub("DN T cell", "DN", de.10genes$subtype)

de.10genes$subtype <- factor(de.10genes$subtype, levels = c("CD8aa", "CD8ab", "CD4+", "DN"))
de.10genes$gene <- factor(de.10genes$gene, levels = c("Cd4", "Cd8a", "Cd8b1", "Gzmb", "Ifng",
                                                       "Pdcd1", "Cd44", "Sell", "Cxcr3", "Klrg1"))

de.10genes$sig_label <- ifelse(de.10genes$p_val_adj < 0.001, "***",
                        ifelse(de.10genes$p_val_adj < 0.01, "**",
                        ifelse(de.10genes$p_val_adj < 0.05, "*", "")))

de.10genes$log2FC_flipped <- -de.10genes$avg_log2FC

p1 <- ggplot(de.10genes, aes(x = subtype, y = gene, fill = log2FC_flipped)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_label), size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0,
                       name = "log2FC\n(+ = up in WT+FFT\n- = up in WT naive)") +
  labs(title = "DE of 10 Genes of Interest (CD45 only)",
       subtitle = "WT naive vs WT+FFT per T cell subtype\n*** p_adj < 0.001  |  ** p_adj < 0.01  |  * p_adj < 0.05",
       x = "T cell subtype", y = "Gene") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        panel.grid = element_blank())

ggsave(file.path(output.dir, "DE results", "10_genes_of_interest_heatmap.pdf"), p1, width = 8, height = 8)

cat("10-gene heatmap saved.\n")

# ============================================================
# 9. PATHWAY ANALYSIS — ORA + GSEA PER SUBTYPE
# ============================================================

cat("\n============================================================\n")
cat("9. Pathway analysis\n")
cat("============================================================\n\n")

subtype_labels <- c("CD8aa T cell", "CD8ab T cell", "CD4+ T cell", "DN T cell")
subtype_dirs <- c("CD8aa", "CD8ab", "CD4", "DN")
subtype_files <- c("CD8aa_T_cell", "CD8ab_T_cell", "CD4__T_cell", "DN_T_cell")

for (d in subtype_dirs) {
  dir.create(file.path(output.dir, "Pathway analysis", d), showWarnings = FALSE, recursive = TRUE)
}

# Helper function for GSEA barplot
make_gsea_barplot <- function(fgsea_sig, title_prefix, outfile) {
  top <- head(fgsea_sig[order(fgsea_sig$padj), ], 20)
  top$pathway_clean <- top$pathway
  top$pathway_clean <- gsub("^GOBP_", "", top$pathway_clean)
  top$pathway_clean <- gsub("^KEGG_", "", top$pathway_clean)
  top$pathway_clean <- gsub("^HALLMARK_", "", top$pathway_clean)
  top$pathway_clean <- gsub("_", " ", top$pathway_clean)
  top$pathway_clean <- tolower(top$pathway_clean)
  top$pathway_clean <- paste0(toupper(substring(top$pathway_clean, 1, 1)), substring(top$pathway_clean, 2))
  top$pathway_clean <- ifelse(nchar(top$pathway_clean) > 60,
                               paste0(substring(top$pathway_clean, 1, 57), "..."), top$pathway_clean)
  top$direction <- ifelse(top$NES > 0, "Up in WT + FFT", "Up in WT naive")

  p <- ggplot(top, aes(x = NES, y = reorder(pathway_clean, NES), fill = direction)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c("Up in WT naive" = "#E41A1C", "Up in WT + FFT" = "#377EB8")) +
    labs(title = title_prefix, x = "Normalized Enrichment Score (NES)\n<-- WT naive | WT+FFT -->",
         y = "", fill = "Direction") +
    theme_classic() +
    theme(plot.title = element_text(size = 13, face = "bold"),
          axis.text.y = element_text(size = 13), axis.text.x = element_text(size = 13),
          axis.title = element_text(size = 14), legend.text = element_text(size = 12),
          legend.title = element_text(size = 13))
  plot_height <- max(4, nrow(top) * 0.5 + 2)
  ggsave(outfile, p, width = 8, height = plot_height)
}

for (i in seq_along(subtype_files)) {
  label <- subtype_labels[i]
  outdir <- file.path(output.dir, "Pathway analysis", subtype_dirs[i])
  de.file <- file.path(output.dir, "DE results", paste0("DE_genomewide_naive_vs_fft_", subtype_files[i], ".csv"))

  if (!file.exists(de.file)) {
    cat("Skipping", label, "- DE file not found\n")
    next
  }

  de.results <- read.csv(de.file, row.names = NULL)
  cat("\n", label, ":", nrow(de.results), "total genes\n")

  # Convert gene symbols to Entrez IDs
  gene_map <- bitr(de.results$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  de.merged <- merge(de.results, gene_map, by.x = "gene", by.y = "SYMBOL")
  sig.genes <- de.merged[de.merged$p_val_adj < 0.05, ]
  cat("Significant genes:", nrow(sig.genes), "\n")

  # --- ORA ---
  if (nrow(sig.genes) >= 5) {
    ego <- enrichGO(gene = sig.genes$ENTREZID, universe = de.merged$ENTREZID,
                    OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
    if (!is.null(ego) && nrow(ego) > 0) {
      write.csv(as.data.frame(ego), file.path(outdir, paste0("ORA_GO_all_significant_", subtype_dirs[i], ".csv")), row.names = FALSE)
      n_show <- min(20, nrow(ego))
      p <- dotplot(ego, showCategory = n_show) +
        ggtitle(paste0("ORA GO BP — ", label, " (CD45 only)")) +
        theme(plot.title = element_text(size = 11, face = "bold"))
      ggsave(file.path(outdir, paste0("ORA_GO_top20_", subtype_dirs[i], ".pdf")), p, width = 12, height = 10)
      cat("ORA GO BP:", nrow(ego), "significant terms\n")
    }

    ekegg <- enrichKEGG(gene = sig.genes$ENTREZID, universe = de.merged$ENTREZID,
                        organism = "mmu", pAdjustMethod = "BH", pvalueCutoff = 0.05)
    if (!is.null(ekegg) && nrow(ekegg) > 0) {
      write.csv(as.data.frame(ekegg), file.path(outdir, paste0("ORA_KEGG_all_significant_", subtype_dirs[i], ".csv")), row.names = FALSE)
      n_show <- min(20, nrow(ekegg))
      p <- dotplot(ekegg, showCategory = n_show) +
        ggtitle(paste0("ORA KEGG — ", label, " (CD45 only)")) +
        theme(plot.title = element_text(size = 11, face = "bold"))
      ggsave(file.path(outdir, paste0("ORA_KEGG_top20_", subtype_dirs[i], ".pdf")), p, width = 12, height = 10)
      cat("ORA KEGG:", nrow(ekegg), "significant pathways\n")
    }
  }

  # --- GSEA ---
  de.ranked <- de.merged[!is.na(de.merged$avg_log2FC), ]
  de.ranked$rank_score <- -de.ranked$avg_log2FC
  gene_ranks <- de.ranked$rank_score
  names(gene_ranks) <- de.ranked$ENTREZID
  gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)

  # GSEA: GO BP
  go_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
  go_pathways <- split(go_sets$ncbi_gene, go_sets$gs_name)
  go_pathways <- lapply(go_pathways, as.character)
  fgsea_go <- fgsea(pathways = go_pathways, stats = gene_ranks, minSize = 15, maxSize = 500)
  fgsea_go_sig <- fgsea_go[fgsea_go$padj < 0.05, ]
  if (nrow(fgsea_go_sig) > 0) {
    fgsea_go_save <- fgsea_go_sig
    fgsea_go_save$leadingEdge <- sapply(fgsea_go_save$leadingEdge, paste, collapse = ", ")
    write.csv(as.data.frame(fgsea_go_save), file.path(outdir, paste0("GSEA_GO_all_significant_", subtype_dirs[i], ".csv")), row.names = FALSE)
    make_gsea_barplot(fgsea_go_sig, paste0("GSEA GO BP — ", label, " (CD45 only)"),
                      file.path(outdir, paste0("GSEA_GO_top20_", subtype_dirs[i], ".pdf")))
    cat("GSEA GO BP:", nrow(fgsea_go_sig), "significant pathways\n")
  }

  # GSEA: KEGG
  kegg_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  kegg_pathways <- split(kegg_sets$ncbi_gene, kegg_sets$gs_name)
  kegg_pathways <- lapply(kegg_pathways, as.character)
  fgsea_kegg <- fgsea(pathways = kegg_pathways, stats = gene_ranks, minSize = 15, maxSize = 500)
  fgsea_kegg_sig <- fgsea_kegg[fgsea_kegg$padj < 0.05, ]
  if (nrow(fgsea_kegg_sig) > 0) {
    fgsea_kegg_save <- fgsea_kegg_sig
    fgsea_kegg_save$leadingEdge <- sapply(fgsea_kegg_save$leadingEdge, paste, collapse = ", ")
    write.csv(as.data.frame(fgsea_kegg_save), file.path(outdir, paste0("GSEA_KEGG_all_significant_", subtype_dirs[i], ".csv")), row.names = FALSE)
    make_gsea_barplot(fgsea_kegg_sig, paste0("GSEA KEGG — ", label, " (CD45 only)"),
                      file.path(outdir, paste0("GSEA_KEGG_top20_", subtype_dirs[i], ".pdf")))
    cat("GSEA KEGG:", nrow(fgsea_kegg_sig), "significant pathways\n")
  }

  # GSEA: Hallmark
  hallmark_sets <- msigdbr(species = "Mus musculus", collection = "H")
  hallmark_pathways <- split(hallmark_sets$ncbi_gene, hallmark_sets$gs_name)
  hallmark_pathways <- lapply(hallmark_pathways, as.character)
  fgsea_hallmark <- fgsea(pathways = hallmark_pathways, stats = gene_ranks, minSize = 15, maxSize = 500)
  fgsea_hallmark_sig <- fgsea_hallmark[fgsea_hallmark$padj < 0.05, ]
  if (nrow(fgsea_hallmark_sig) > 0) {
    fgsea_hallmark_save <- fgsea_hallmark_sig
    fgsea_hallmark_save$leadingEdge <- sapply(fgsea_hallmark_save$leadingEdge, paste, collapse = ", ")
    write.csv(as.data.frame(fgsea_hallmark_save), file.path(outdir, paste0("GSEA_Hallmark_all_significant_", subtype_dirs[i], ".csv")), row.names = FALSE)
    make_gsea_barplot(fgsea_hallmark_sig, paste0("GSEA Hallmark — ", label, " (CD45 only)"),
                      file.path(outdir, paste0("GSEA_Hallmark_top20_", subtype_dirs[i], ".pdf")))
    cat("GSEA Hallmark:", nrow(fgsea_hallmark_sig), "significant pathways\n")
  }
}

# ============================================================
# 10. SAVE OBJECT
# ============================================================

cat("\n============================================================\n")
cat("10. Saving annotated Seurat object\n")
cat("============================================================\n\n")

SaveH5Seurat(T.cell, filename = file.path(output.dir, "combined.all.final.T.cell.CD45only.h5Seurat"), overwrite = TRUE)

# ============================================================
# 11. COMPARISON SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("COMPARISON: CD45-only vs Original T cell update\n")
cat("============================================================\n\n")

cat("Cell counts:\n")
cat("  Original:  6,565 total (2,106 WT naive + 4,459 WT+FFT)\n")
cat("  CD45 only:", ncol(T.cell), "total\n")
cat("  Difference: 136 fewer cells (from Epcam cluster 20)\n\n")

cat("Subtype distribution (CD45 only):\n")
print(table(T.cell$Tcell_subtype, T.cell$modified.ident))

cat("\n\nDone! All outputs saved to:\n")
cat(output.dir, "\n")
