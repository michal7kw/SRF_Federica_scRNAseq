################################################################################
# Script: regenerate_heatmaps.R
# Description: Fix and regenerate marker gene heatmaps
# Bug Fix: Use rownames(markers) instead of markers$gene
# Author: Bioinformatics Pipeline
# Date: 2025-11-20
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# Set working directory
work_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
setwd(work_dir)

cat("\n================================================================\n")
cat("Regenerating Marker Gene Heatmaps (Bug Fix)\n")
cat("================================================================\n\n")

################################################################################
# Process both no_harmony and harmony versions
################################################################################

process_heatmaps <- function(method_name, seurat_path, output_dir) {

  cat(paste0("Processing ", method_name, " version...\n"))
  cat("----------------------------------------------------------------\n")

  # Check if Seurat object exists
  if (!file.exists(seurat_path)) {
    cat("  ✗ Seurat object not found:", seurat_path, "\n")
    cat("  Skipping this version.\n\n")
    return(FALSE)
  }

  # Load Seurat object
  cat("  Loading Seurat object...\n")
  seurat_obj <- readRDS(seurat_path)

  # Load or recalculate markers
  markers_file <- file.path(output_dir, "cluster_markers_all.csv")

  if (file.exists(markers_file)) {
    cat("  Loading existing markers from CSV...\n")
    all_markers <- read.csv(markers_file)

    # If gene column doesn't exist, it might be in rownames
    # Try to reconstruct it
    if (!"gene" %in% colnames(all_markers)) {
      # Check if first column looks like gene names
      if (all(grepl("^[A-Z0-9-]+$", head(all_markers[,1], 100)))) {
        all_markers$gene <- all_markers[,1]
      } else {
        cat("  ⚠ Warning: Cannot identify gene column in CSV\n")
        cat("  Recalculating markers...\n")
        all_markers <- NULL
      }
    }
  } else {
    all_markers <- NULL
  }

  # Recalculate markers if needed
  if (is.null(all_markers)) {
    cat("  Calculating markers (this may take a few minutes)...\n")

    # Join layers first (Seurat v5 requirement)
    seurat_obj <- JoinLayers(seurat_obj)

    all_markers <- FindAllMarkers(seurat_obj,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25,
                                  verbose = FALSE)

    # Add gene column from rownames
    all_markers$gene <- rownames(all_markers)

    # Save the updated markers
    write.csv(all_markers, markers_file, row.names = FALSE)
  }

  # Get top markers per cluster
  if ("cluster" %in% colnames(all_markers)) {
    top10_markers <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC)

    top5_markers <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC)
  } else {
    # Extract cluster from rownames if needed
    if (is.null(all_markers$cluster)) {
      all_markers$cluster <- sub("\\..*", "", rownames(all_markers))
    }

    top10_markers <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC)

    top5_markers <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC)
  }

  cat("  Top markers identified:\n")
  cat("    Clusters:", length(unique(top5_markers$cluster)), "\n")
  cat("    Top 5 markers per cluster:", nrow(top5_markers), "genes\n")
  cat("    Top 10 markers per cluster:", nrow(top10_markers), "genes\n\n")

  # Create plots directory if needed
  dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

  ################################################################################
  # FIX 1: Heatmap using correct gene names from 'gene' column
  ################################################################################
  cat("  Creating heatmap (TOP 5 markers)...\n")

  # Get unique genes in correct order
  genes_for_heatmap <- unique(top5_markers$gene)

  # Remove any NA or empty values
  genes_for_heatmap <- genes_for_heatmap[!is.na(genes_for_heatmap) & genes_for_heatmap != ""]

  # Keep only genes that exist in the Seurat object
  genes_for_heatmap <- genes_for_heatmap[genes_for_heatmap %in% rownames(seurat_obj)]

  cat("    Genes to plot:", length(genes_for_heatmap), "\n")

  if (length(genes_for_heatmap) > 0) {
    # Set cluster identity for ordering
    Idents(seurat_obj) <- "seurat_clusters"

    # Create heatmap
    p_heatmap <- DoHeatmap(seurat_obj,
                           features = genes_for_heatmap,
                           size = 3,
                           angle = 90) +
      scale_fill_gradientn(colors = c("blue", "white", "red")) +
      theme(text = element_text(size = 8))

    # Save as PNG
    ggsave(file.path(output_dir, "plots", "09_top_markers_heatmap.png"),
           plot = p_heatmap,
           width = 14,
           height = max(10, length(genes_for_heatmap) * 0.2),  # Dynamic height
           dpi = 300)

    cat("    ✓ Heatmap saved: 09_top_markers_heatmap.png\n")
  } else {
    cat("    ✗ ERROR: No valid genes found for heatmap!\n")
  }

  ################################################################################
  # FIX 2: Dot plot
  ################################################################################
  cat("  Creating dot plot (TOP 5 markers)...\n")

  if (length(genes_for_heatmap) > 0) {
    p_dotplot <- DotPlot(seurat_obj,
                         features = genes_for_heatmap) +
      RotatedAxis() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    ggsave(file.path(output_dir, "plots", "10_top_markers_dotplot.png"),
           plot = p_dotplot,
           width = max(12, length(genes_for_heatmap) * 0.3),  # Dynamic width
           height = 8,
           dpi = 300)

    cat("    ✓ Dot plot saved: 10_top_markers_dotplot.png\n")
  } else {
    cat("    ✗ ERROR: No valid genes found for dot plot!\n")
  }

  ################################################################################
  # BONUS: Create heatmap with TOP 10 markers as well
  ################################################################################
  cat("  Creating heatmap (TOP 10 markers)...\n")

  genes_for_heatmap_10 <- unique(top10_markers$gene)
  genes_for_heatmap_10 <- genes_for_heatmap_10[!is.na(genes_for_heatmap_10) & genes_for_heatmap_10 != ""]
  genes_for_heatmap_10 <- genes_for_heatmap_10[genes_for_heatmap_10 %in% rownames(seurat_obj)]

  cat("    Genes to plot:", length(genes_for_heatmap_10), "\n")

  if (length(genes_for_heatmap_10) > 0) {
    p_heatmap_10 <- DoHeatmap(seurat_obj,
                              features = genes_for_heatmap_10,
                              size = 2.5,
                              angle = 90) +
      scale_fill_gradientn(colors = c("blue", "white", "red")) +
      theme(text = element_text(size = 7))

    ggsave(file.path(output_dir, "plots", "09b_top10_markers_heatmap.png"),
           plot = p_heatmap_10,
           width = 16,
           height = max(12, length(genes_for_heatmap_10) * 0.15),
           dpi = 300)

    cat("    ✓ Top 10 heatmap saved: 09b_top10_markers_heatmap.png\n")
  }

  cat("\n  ✓", method_name, "heatmaps regenerated successfully!\n\n")
  return(TRUE)
}

################################################################################
# Process both versions
################################################################################

# Try no_harmony version
no_harm_success <- process_heatmaps(
  method_name = "No-Harmony",
  seurat_path = "results/clustering/no_harmony/seurat_integrated.rds",
  output_dir = "results/clustering/no_harmony"
)

# Try harmony version
harm_success <- process_heatmaps(
  method_name = "Harmony",
  seurat_path = "results/clustering/harmony/seurat_integrated.rds",
  output_dir = "results/clustering/harmony"
)

# Try old file structure as fallback
if (!no_harm_success && !harm_success) {
  cat("Trying old file structure...\n")
  cat("----------------------------------------------------------------\n")

  if (file.exists("results/clustering/seurat_integrated_no_harmony.rds")) {
    process_heatmaps(
      method_name = "No-Harmony (old path)",
      seurat_path = "results/clustering/seurat_integrated_no_harmony.rds",
      output_dir = "results/clustering"
    )
  }

  if (file.exists("results/clustering/seurat_integrated_harmony.rds")) {
    process_heatmaps(
      method_name = "Harmony (old path)",
      seurat_path = "results/clustering/seurat_integrated_harmony.rds",
      output_dir = "results/clustering"
    )
  }
}

cat("\n================================================================\n")
cat("Heatmap Regeneration Complete\n")
cat("================================================================\n\n")

cat("Summary:\n")
cat("  Bug fixed: Now using 'gene' column instead of non-existent rownames\n")
cat("  Output files:\n")
cat("    - 09_top_markers_heatmap.png (top 5 markers per cluster)\n")
cat("    - 09b_top10_markers_heatmap.png (top 10 markers per cluster)\n")
cat("    - 10_top_markers_dotplot.png (top 5 markers dotplot)\n\n")
