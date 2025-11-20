################################################################################
# Script: 03_integration_clustering.R
# Description: Integration, normalization, clustering, and cell type annotation
# Input: Filtered Seurat objects from step 02
# Output: Integrated Seurat object with clusters and cell type annotations
# Author: Bioinformatics Pipeline
# Date: 2025-11-14
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)  # For integration
  library(clustree)  # For cluster tree visualization
})

# Set working directory
work_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
setwd(work_dir)

# Create output directories - SEPARATE for harmony to avoid conflicts
OUTPUT_DIR <- "results/clustering/harmony"
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)

cat("\n================================================================\n")
cat("scRNA-seq Integration and Clustering (WITH Harmony)\n")
cat("================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

################################################################################
# 1. Load filtered data
################################################################################
cat("Step 1: Loading filtered Seurat objects...\n")
cat("----------------------------------------------------------------\n")

seurat_obj <- readRDS("results/qc/combined_filtered.rds")

cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  Total features:", nrow(seurat_obj), "\n")
cat("  Conditions:", table(seurat_obj$condition), "\n\n")

################################################################################
# 2. Normalization and feature selection
################################################################################
cat("Step 2: Normalizing data and selecting variable features...\n")
cat("----------------------------------------------------------------\n")

# Standard log-normalization
seurat_obj <- NormalizeData(seurat_obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj,
                                   selection.method = "vst",
                                   nfeatures = 3000)

# Identify top 20 most variable genes
top20 <- head(VariableFeatures(seurat_obj), 20)
cat("  Top 20 variable genes:\n")
cat("   ", paste(top20, collapse = ", "), "\n\n")

# Plot variable features
p1 <- VariableFeaturePlot(seurat_obj)
p2 <- LabelPoints(plot = p1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)

ggsave(file.path(OUTPUT_DIR, "plots", "01_variable_features.png"),
       plot = p2, width = 12, height = 6, dpi = 300)

cat("  Variable features identified and plotted\n\n")

################################################################################
# 3. Scaling and PCA
################################################################################
cat("Step 3: Scaling data and running PCA...\n")
cat("----------------------------------------------------------------\n")

# Scale data (regress out mitochondrial content)
seurat_obj <- ScaleData(seurat_obj,
                        vars.to.regress = c("percent.mt", "nCount_RNA"))

# Run PCA
seurat_obj <- RunPCA(seurat_obj,
                     features = VariableFeatures(object = seurat_obj),
                     npcs = 50,
                     verbose = FALSE)

# Visualize PCA results
p3 <- DimPlot(seurat_obj, reduction = "pca", group.by = "condition")
p4 <- ElbowPlot(seurat_obj, ndims = 50)

p_pca <- p3 + p4
ggsave(file.path(OUTPUT_DIR, "plots", "02_pca_overview.png"),
       plot = p_pca, width = 12, height = 5, dpi = 300)

# PCA heatmaps for top PCs
png(file.path(OUTPUT_DIR, "plots", "03_pca_heatmaps.png"), width = 12, height = 15, units = "in", res = 300)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE, ncol = 3)
dev.off()

cat("  PCA complete. Top 50 PCs computed\n\n")

################################################################################
# 4. Batch correction with Harmony
################################################################################
cat("Step 4: Running Harmony integration...\n")
cat("----------------------------------------------------------------\n")

# Run Harmony to integrate samples
# Note: Minimal parameters to avoid conflicts - reduction defaults to "pca"
# This is the ONLY way that works with latest Harmony version
seurat_obj <- RunHarmony(seurat_obj, "condition")

# Visualize Harmony correction
p5 <- DimPlot(seurat_obj, reduction = "pca", group.by = "condition") +
      ggtitle("Before Harmony")
p6 <- DimPlot(seurat_obj, reduction = "harmony", group.by = "condition") +
      ggtitle("After Harmony")

p_harmony <- p5 + p6
ggsave(file.path(OUTPUT_DIR, "plots", "04_harmony_integration.png"),
       plot = p_harmony, width = 14, height = 5, dpi = 300)

cat("  Harmony integration complete\n\n")

################################################################################
# 5. UMAP and clustering
################################################################################
cat("Step 5: Running UMAP and clustering...\n")
cat("----------------------------------------------------------------\n")

# Run UMAP on Harmony embeddings
seurat_obj <- RunUMAP(seurat_obj,
                      reduction = "harmony",
                      dims = 1:30,
                      verbose = FALSE)

# Find neighbors and clusters at multiple resolutions
seurat_obj <- FindNeighbors(seurat_obj,
                            reduction = "harmony",
                            dims = 1:30,
                            verbose = FALSE)

# Test multiple resolutions for clustering
resolutions <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2)
seurat_obj <- FindClusters(seurat_obj,
                           resolution = resolutions,
                           verbose = FALSE)

# Visualize clustering at different resolutions using clustree
png(file.path(OUTPUT_DIR, "plots", "05_clustree.png"), width = 10, height = 12, units = "in", res = 300)
clustree(seurat_obj, prefix = "RNA_snn_res.")
dev.off()

cat("  Clustering complete at multiple resolutions\n")
cat("  Recommended resolution: 0.5 (can be adjusted)\n\n")

# Set default identity to resolution 0.5
Idents(seurat_obj) <- "RNA_snn_res.0.5"
seurat_obj$seurat_clusters <- seurat_obj$RNA_snn_res.0.5

cat("  Number of clusters at res=0.5:", length(unique(seurat_obj$seurat_clusters)), "\n\n")

################################################################################
# 6. Visualization
################################################################################
cat("Step 6: Generating visualization plots...\n")
cat("----------------------------------------------------------------\n")

# UMAP by cluster
p7 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
      ggtitle("Clusters (res=0.5)")

# UMAP by condition
p8 <- DimPlot(seurat_obj, reduction = "umap", group.by = "condition") +
      ggtitle("Condition")

# UMAP split by condition
p9 <- DimPlot(seurat_obj, reduction = "umap", split.by = "condition", label = TRUE) +
      ggtitle("Clusters by Condition")

# Combined plots
p_umap1 <- p7 + p8
ggsave(file.path(OUTPUT_DIR, "plots", "06_umap_clusters_condition.png"),
       plot = p_umap1, width = 14, height = 6, dpi = 300)

ggsave(file.path(OUTPUT_DIR, "plots", "07_umap_split_condition.png"),
       plot = p9, width = 14, height = 6, dpi = 300)

# QC metrics on UMAP
p10 <- FeaturePlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3)
ggsave(file.path(OUTPUT_DIR, "plots", "08_umap_qc_metrics.png"),
       plot = p10, width = 18, height = 5, dpi = 300)

cat("  UMAP plots generated\n\n")

################################################################################
# 7. Find cluster markers
################################################################################
cat("Step 7: Finding marker genes for each cluster...\n")
cat("----------------------------------------------------------------\n")

# Join layers before finding markers (required for Seurat v5)
seurat_obj <- JoinLayers(seurat_obj)

# Find markers for all clusters
all_markers <- FindAllMarkers(seurat_obj,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              verbose = FALSE)

# Add gene column from rownames (CRITICAL: FindAllMarkers returns genes as rownames!)
all_markers$gene <- rownames(all_markers)

# Save all markers
write.csv(all_markers, file.path(OUTPUT_DIR, "cluster_markers_all.csv"), row.names = FALSE)

# Get top markers per cluster
# Note: Column name might be 'cluster' or another name depending on Seurat version
if ("cluster" %in% colnames(all_markers)) {
  top10_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

  top5_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
} else {
  # For newer Seurat versions, markers are returned with rownames as cluster info
  all_markers$cluster <- sub("\\..*", "", rownames(all_markers))

  top10_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

  top5_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
}

write.csv(top10_markers, file.path(OUTPUT_DIR, "cluster_markers_top10.csv"), row.names = FALSE)

cat("  Marker genes identified for all clusters\n\n")

# Heatmap of top markers (using gene column, not rownames)
png(file.path(OUTPUT_DIR, "plots", "09_top_markers_heatmap.png"), width = 14, height = 10, units = "in", res = 300)
DoHeatmap(seurat_obj, features = top5_markers$gene, size = 3) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

# Dot plot of top markers (using gene column, not rownames)
png(file.path(OUTPUT_DIR, "plots", "10_top_markers_dotplot.png"), width = 16, height = 10, units = "in", res = 300)
DotPlot(seurat_obj, features = unique(top5_markers$gene)) +
  RotatedAxis()
dev.off()

cat("  Marker visualization plots generated\n\n")

################################################################################
# 8. Automatic cell type annotation
################################################################################
cat("Step 8: Annotating cell types based on markers...\n")
cat("----------------------------------------------------------------\n")

# Define brain organoid marker genes
# These are common markers - adjust based on your specific organoid protocol
brain_markers <- list(
  "Radial Glia/NPCs" = c("PAX6", "SOX2", "NES", "HES1", "VIM"),
  "Intermediate Progenitors" = c("EOMES", "TBR2"),
  "Excitatory Neurons" = c("NEUROD6", "TBR1", "SATB2", "BCL11B", "SLC17A7"),
  "Inhibitory Neurons" = c("GAD1", "GAD2", "DLX2", "DLX5", "SLC32A1"),
  "Cajal-Retzius Cells" = c("RELN", "CALB2", "TRPC3"),
  "Astrocytes" = c("GFAP", "AQP4", "S100B", "ALDH1L1"),
  "Oligodendrocyte Precursors" = c("PDGFRA", "OLIG2", "CSPG4"),
  "Oligodendrocytes" = c("MBP", "MOG", "PLP1"),
  "Microglia" = c("CX3CR1", "P2RY12", "TMEM119", "AIF1"),
  "Choroid Plexus" = c("TTR", "CLIC6"),
  "Cycling Cells" = c("MKI67", "TOP2A", "CDK1")
)

# Visualize canonical markers
for (cell_type in names(brain_markers)) {
  genes <- brain_markers[[cell_type]]
  genes_present <- genes[genes %in% rownames(seurat_obj)]

  if (length(genes_present) > 0) {
    cat("  Visualizing", cell_type, "markers:", paste(genes_present, collapse = ", "), "\n")

    p <- FeaturePlot(seurat_obj, features = genes_present, ncol = 3)

    filename <- file.path(OUTPUT_DIR, "plots",
                         paste0("markers_", gsub("[/ ]", "_", tolower(cell_type)), ".png"))
    ggsave(filename, plot = p, width = 15, height = ceiling(length(genes_present)/3) * 4, dpi = 300)
  }
}

# Create a dot plot with all markers
all_marker_genes <- unique(unlist(brain_markers))
all_marker_genes <- all_marker_genes[all_marker_genes %in% rownames(seurat_obj)]

png(file.path(OUTPUT_DIR, "plots", "11_canonical_markers_dotplot.png"), width = 16, height = 10, units = "in", res = 300)
DotPlot(seurat_obj, features = all_marker_genes) +
  RotatedAxis() +
  ggtitle("Canonical Brain Cell Type Markers")
dev.off()

cat("\n  Cell type marker plots generated\n")
cat("  NOTE: Manual annotation recommended based on marker expression\n\n")

################################################################################
# 9. Cell cycle scoring
################################################################################
cat("Step 9: Cell cycle scoring...\n")
cat("----------------------------------------------------------------\n")

# Cell cycle genes
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

# Score cell cycle
seurat_obj <- CellCycleScoring(seurat_obj,
                               s.features = s_genes,
                               g2m.features = g2m_genes,
                               set.ident = FALSE)

# Visualize cell cycle scores
p11 <- DimPlot(seurat_obj, reduction = "umap", group.by = "Phase") +
       ggtitle("Cell Cycle Phase")

p12 <- FeaturePlot(seurat_obj, features = c("S.Score", "G2M.Score"), ncol = 2)

p_cellcycle <- p11 / p12
ggsave(file.path(OUTPUT_DIR, "plots", "12_cell_cycle.png"),
       plot = p_cellcycle, width = 12, height = 10, dpi = 300)

cat("  Cell cycle scoring complete\n\n")

################################################################################
# 10. Cluster composition analysis
################################################################################
cat("Step 10: Analyzing cluster composition...\n")
cat("----------------------------------------------------------------\n")

# Calculate cluster proportions by condition
cluster_composition <- as.data.frame.matrix(table(seurat_obj$seurat_clusters, seurat_obj$condition))
cluster_composition$cluster <- rownames(cluster_composition)
cluster_composition$total <- cluster_composition$KO + cluster_composition$WT
cluster_composition$KO_pct <- 100 * cluster_composition$KO / cluster_composition$total
cluster_composition$WT_pct <- 100 * cluster_composition$WT / cluster_composition$total

write.csv(cluster_composition, file.path(OUTPUT_DIR, "cluster_composition.csv"), row.names = FALSE)

cat("\n  Cluster Composition:\n")
print(cluster_composition[, c("cluster", "KO", "WT", "KO_pct", "WT_pct")])

# Plot cluster composition
library(reshape2)
composition_melt <- melt(cluster_composition[, c("cluster", "KO", "WT")],
                         id.vars = "cluster")
colnames(composition_melt) <- c("Cluster", "Condition", "Cells")

p13 <- ggplot(composition_melt, aes(x = Cluster, y = Cells, fill = Condition)) +
       geom_bar(stat = "identity", position = "dodge") +
       theme_minimal() +
       ggtitle("Cluster Composition by Condition") +
       theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "plots", "13_cluster_composition.png"),
       plot = p13, width = 10, height = 6, dpi = 300)

cat("\n  Cluster composition analysis complete\n\n")

################################################################################
# 11. Save integrated object
################################################################################
cat("Step 11: Saving integrated Seurat object...\n")
cat("----------------------------------------------------------------\n")

saveRDS(seurat_obj, file = "results/clustering/seurat_integrated_harmony.rds")

cat("  Integrated object saved to results/clustering/seurat_integrated_harmony.rds\n\n")

################################################################################
# 12. Generate summary report
################################################################################
cat("Step 12: Generating summary report...\n")
cat("----------------------------------------------------------------\n")

summary_stats <- data.frame(
  Metric = c("Total Cells", "Total Genes", "Number of Clusters",
             "KO Cells", "WT Cells", "Variable Features"),
  Value = c(ncol(seurat_obj), nrow(seurat_obj),
            length(unique(seurat_obj$seurat_clusters)),
            sum(seurat_obj$condition == "KO"),
            sum(seurat_obj$condition == "WT"),
            length(VariableFeatures(seurat_obj)))
)

write.csv(summary_stats, file.path(OUTPUT_DIR, "integration_summary.csv"), row.names = FALSE)

cat("\n  Summary Statistics:\n")
print(summary_stats)

cat("\n================================================================\n")
cat("Integration and Clustering Complete\n")
cat("================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nOutputs:\n")
cat("  - Integrated object: results/clustering/seurat_integrated.rds\n")
cat("  - Cluster markers: results/clustering/cluster_markers_*.csv\n")
cat("  - Plots: results/clustering/plots/\n")
cat("  - Summary statistics: results/clustering/integration_summary.csv\n")
cat("\nNext steps:\n")
cat("  1. Review UMAP plots and cluster markers\n")
cat("  2. Manually annotate cell types if needed\n")
cat("  3. Proceed to differential expression analysis (step 04)\n")
cat("================================================================\n\n")
