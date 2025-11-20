################################################################################
# Script: compare_integration_methods.R
# Description: Compare Harmony vs No-Harmony integration results
# Author: Bioinformatics Pipeline
# Date: 2025-11-14
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
})

work_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
setwd(work_dir)

# Create comparison output directory
dir.create("results/clustering/comparison", recursive = TRUE, showWarnings = FALSE)

cat("\n================================================================\n")
cat("Comparing Integration Methods: Harmony vs No Harmony\n")
cat("================================================================\n\n")

################################################################################
# 1. Load both Seurat objects
################################################################################
cat("Step 1: Loading Seurat objects...\n")
cat("----------------------------------------------------------------\n")

# Check which files exist - try new directory structure first, then old
# New structure (separate directories)
harmony_file_new <- "results/clustering/harmony/seurat_integrated.rds"
no_harmony_file_new <- "results/clustering/no_harmony/seurat_integrated.rds"

# Old structure (same directory, different names)
harmony_file_old <- "results/clustering/seurat_integrated_harmony.rds"
no_harmony_file_old <- "results/clustering/seurat_integrated_no_harmony.rds"
default_file_old <- "results/clustering/seurat_integrated.rds"

# Determine which files exist
harmony_file <- if (file.exists(harmony_file_new)) harmony_file_new else harmony_file_old
no_harmony_file <- if (file.exists(no_harmony_file_new)) no_harmony_file_new else no_harmony_file_old

files_available <- c(
  harmony = file.exists(harmony_file),
  no_harmony = file.exists(no_harmony_file),
  default = file.exists(default_file_old)
)

cat("  Files found:\n")
if (files_available["harmony"]) cat("    ✓ Harmony version:", harmony_file, "\n")
if (files_available["no_harmony"]) cat("    ✓ No-Harmony version:", no_harmony_file, "\n")
if (files_available["default"] && !files_available["harmony"] && !files_available["no_harmony"]) {
  cat("    ✓ Default version:", default_file_old, "\n")
}

# Load available objects
if (files_available["harmony"] && files_available["no_harmony"]) {
  seurat_harmony <- readRDS(harmony_file)
  seurat_no_harmony <- readRDS(no_harmony_file)
  comparison_possible <- TRUE
} else if (files_available["harmony"] && files_available["default"]) {
  seurat_harmony <- readRDS(harmony_file)
  seurat_no_harmony <- readRDS(default_file_old)
  comparison_possible <- TRUE
} else if (files_available["no_harmony"] && files_available["default"]) {
  seurat_harmony <- readRDS(default_file_old)
  seurat_no_harmony <- readRDS(no_harmony_file)
  comparison_possible <- TRUE
} else {
  cat("\n  ERROR: Need at least 2 Seurat objects to compare!\n")
  cat("  Run both integration methods first.\n\n")
  cat("  Expected files:\n")
  cat("    New structure:\n")
  cat("      - results/clustering/harmony/seurat_integrated.rds\n")
  cat("      - results/clustering/no_harmony/seurat_integrated.rds\n")
  cat("    OR old structure:\n")
  cat("      - results/clustering/seurat_integrated_harmony.rds\n")
  cat("      - results/clustering/seurat_integrated_no_harmony.rds\n\n")
  quit(status = 1)
}

cat("  Both objects loaded successfully\n\n")

################################################################################
# 2. Basic statistics comparison
################################################################################
cat("Step 2: Comparing basic statistics...\n")
cat("----------------------------------------------------------------\n")

stats_comparison <- data.frame(
  Metric = c("Total Cells", "KO Cells", "WT Cells", "Total Genes",
             "Number of Clusters", "Harmony Used"),
  No_Harmony = c(
    ncol(seurat_no_harmony),
    sum(seurat_no_harmony$condition == "KO"),
    sum(seurat_no_harmony$condition == "WT"),
    nrow(seurat_no_harmony),
    length(unique(seurat_no_harmony$seurat_clusters)),
    "No"
  ),
  With_Harmony = c(
    ncol(seurat_harmony),
    sum(seurat_harmony$condition == "KO"),
    sum(seurat_harmony$condition == "WT"),
    nrow(seurat_harmony),
    length(unique(seurat_harmony$seurat_clusters)),
    "Yes"
  )
)

cat("\n")
print(stats_comparison)
cat("\n")

write.csv(stats_comparison, "results/clustering/comparison/stats_comparison.csv",
          row.names = FALSE)

################################################################################
# 3. UMAP comparison
################################################################################
cat("Step 3: Comparing UMAP visualizations...\n")
cat("----------------------------------------------------------------\n")

# UMAPs colored by condition
p1 <- DimPlot(seurat_no_harmony, reduction = "umap", group.by = "condition") +
      ggtitle("No Harmony - By Condition") +
      theme(legend.position = "bottom")

p2 <- DimPlot(seurat_harmony, reduction = "umap", group.by = "condition") +
      ggtitle("With Harmony - By Condition") +
      theme(legend.position = "bottom")

# UMAPs colored by cluster
p3 <- DimPlot(seurat_no_harmony, reduction = "umap", label = TRUE, repel = TRUE) +
      ggtitle("No Harmony - By Cluster") +
      NoLegend()

p4 <- DimPlot(seurat_harmony, reduction = "umap", label = TRUE, repel = TRUE) +
      ggtitle("With Harmony - By Cluster") +
      NoLegend()

# Combine plots
p_umap_comparison <- (p1 | p2) / (p3 | p4)

ggsave("results/clustering/comparison/01_umap_comparison.png",
       plot = p_umap_comparison, width = 14, height = 12, dpi = 300)

cat("  UMAP comparison plot saved\n\n")

################################################################################
# 4. Cluster composition comparison
################################################################################
cat("Step 4: Comparing cluster composition...\n")
cat("----------------------------------------------------------------\n")

# No Harmony composition
comp_no_harmony <- as.data.frame.matrix(
  table(seurat_no_harmony$seurat_clusters, seurat_no_harmony$condition)
)
comp_no_harmony$cluster <- rownames(comp_no_harmony)
comp_no_harmony$method <- "No Harmony"
comp_no_harmony$KO_pct <- 100 * comp_no_harmony$KO /
                                 (comp_no_harmony$KO + comp_no_harmony$WT)

# With Harmony composition
comp_harmony <- as.data.frame.matrix(
  table(seurat_harmony$seurat_clusters, seurat_harmony$condition)
)
comp_harmony$cluster <- rownames(comp_harmony)
comp_harmony$method <- "With Harmony"
comp_harmony$KO_pct <- 100 * comp_harmony$KO /
                              (comp_harmony$KO + comp_harmony$WT)

# Combine
composition_comparison <- rbind(
  comp_no_harmony[, c("cluster", "KO", "WT", "KO_pct", "method")],
  comp_harmony[, c("cluster", "KO", "WT", "KO_pct", "method")]
)

write.csv(composition_comparison,
          "results/clustering/comparison/cluster_composition_comparison.csv",
          row.names = FALSE)

# Plot composition comparison
library(reshape2)
comp_melt <- melt(composition_comparison[, c("cluster", "KO", "WT", "method")],
                  id.vars = c("cluster", "method"))

p5 <- ggplot(comp_melt, aes(x = cluster, y = value, fill = variable)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~method, ncol = 1) +
      theme_minimal() +
      labs(title = "Cluster Composition Comparison",
           x = "Cluster", y = "Number of Cells", fill = "Condition") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/clustering/comparison/02_cluster_composition_comparison.png",
       plot = p5, width = 12, height = 8, dpi = 300)

cat("  Cluster composition comparison saved\n\n")

################################################################################
# 5. Batch mixing metrics
################################################################################
cat("Step 5: Calculating batch mixing metrics...\n")
cat("----------------------------------------------------------------\n")

# Function to calculate mixing metric
calculate_mixing <- function(seurat_obj) {
  # Get UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, reduction = "umap")
  condition <- seurat_obj$condition

  # For each cell, find nearest neighbors
  nn_results <- list()
  for (i in 1:min(100, nrow(umap_coords))) {  # Sample 100 cells for speed
    # Calculate distances to all other cells
    distances <- sqrt(rowSums((sweep(umap_coords, 2, umap_coords[i, ]))^2))
    # Get 20 nearest neighbors (excluding self)
    nn_idx <- order(distances)[2:21]
    # Calculate proportion of neighbors from same condition
    same_condition <- sum(condition[nn_idx] == condition[i]) / 20
    nn_results[[i]] <- same_condition
  }

  mean_same_condition <- mean(unlist(nn_results))
  return(mean_same_condition)
}

mixing_no_harmony <- calculate_mixing(seurat_no_harmony)
mixing_harmony <- calculate_mixing(seurat_harmony)

cat("\n  Batch Mixing Scores (20 nearest neighbors):\n")
cat("    No Harmony: ", round(mixing_no_harmony, 3), "\n")
cat("    With Harmony: ", round(mixing_harmony, 3), "\n\n")
cat("  Interpretation:\n")
cat("    1.0 = Perfect separation (all neighbors same condition)\n")
cat("    0.5 = Perfect mixing (50% neighbors from each condition)\n")
cat("    0.0 = Complete integration (all neighbors different condition)\n\n")

# Optimal range: 0.4-0.6 suggests good mixing while preserving biology
if (mixing_no_harmony >= 0.4 && mixing_no_harmony <= 0.6) {
  cat("  ✓ No Harmony shows good mixing - likely no batch effect\n")
} else if (mixing_no_harmony > 0.7) {
  cat("  ⚠ No Harmony shows poor mixing - possible batch effect\n")
}

if (mixing_harmony >= 0.4 && mixing_harmony <= 0.6) {
  cat("  ✓ With Harmony shows good mixing\n")
} else if (mixing_harmony < 0.4) {
  cat("  ⚠ Harmony may be over-correcting (too much mixing)\n")
}

cat("\n")

mixing_comparison <- data.frame(
  Method = c("No Harmony", "With Harmony"),
  Mixing_Score = c(mixing_no_harmony, mixing_harmony),
  Interpretation = c(
    ifelse(mixing_no_harmony >= 0.4 & mixing_no_harmony <= 0.6,
           "Good mixing",
           ifelse(mixing_no_harmony > 0.7, "Poor mixing", "Over-mixing")),
    ifelse(mixing_harmony >= 0.4 & mixing_harmony <= 0.6,
           "Good mixing",
           ifelse(mixing_harmony > 0.7, "Poor mixing", "Over-mixing"))
  )
)

write.csv(mixing_comparison,
          "results/clustering/comparison/mixing_metrics.csv",
          row.names = FALSE)

################################################################################
# 6. Marker gene preservation
################################################################################
cat("Step 6: Comparing marker gene patterns...\n")
cat("----------------------------------------------------------------\n")

# Use some canonical brain markers
canonical_markers <- c("PAX6", "SOX2", "NEUROD6", "GAD1", "GFAP", "MBP")
canonical_markers <- canonical_markers[canonical_markers %in% rownames(seurat_no_harmony)]

if (length(canonical_markers) > 0) {
  # Feature plots for both methods
  p6 <- FeaturePlot(seurat_no_harmony, features = canonical_markers, ncol = 3) +
        plot_annotation(title = "No Harmony")

  p7 <- FeaturePlot(seurat_harmony, features = canonical_markers, ncol = 3) +
        plot_annotation(title = "With Harmony")

  ggsave("results/clustering/comparison/03_markers_no_harmony.png",
         plot = p6, width = 15, height = 10, dpi = 300)

  ggsave("results/clustering/comparison/04_markers_with_harmony.png",
         plot = p7, width = 15, height = 10, dpi = 300)

  cat("  Marker gene comparison plots saved\n\n")
}

################################################################################
# 7. Silhouette score comparison
################################################################################
cat("Step 7: Calculating silhouette scores...\n")
cat("----------------------------------------------------------------\n")

# Function to calculate mean silhouette score
calc_silhouette <- function(seurat_obj) {
  library(cluster)
  umap_coords <- Embeddings(seurat_obj, reduction = "umap")
  clusters <- as.numeric(as.character(seurat_obj$seurat_clusters))

  # Calculate silhouette (sample if too many cells)
  if (nrow(umap_coords) > 1000) {
    sample_idx <- sample(1:nrow(umap_coords), 1000)
    sil <- silhouette(clusters[sample_idx], dist(umap_coords[sample_idx, ]))
  } else {
    sil <- silhouette(clusters, dist(umap_coords))
  }

  return(mean(sil[, 3]))
}

sil_no_harmony <- calc_silhouette(seurat_no_harmony)
sil_harmony <- calc_silhouette(seurat_harmony)

cat("\n  Silhouette Scores (cluster separation):\n")
cat("    No Harmony: ", round(sil_no_harmony, 3), "\n")
cat("    With Harmony: ", round(sil_harmony, 3), "\n\n")
cat("  Interpretation:\n")
cat("    > 0.5: Strong cluster separation\n")
cat("    0.25-0.5: Reasonable cluster separation\n")
cat("    < 0.25: Weak cluster separation\n\n")

silhouette_comparison <- data.frame(
  Method = c("No Harmony", "With Harmony"),
  Silhouette_Score = c(sil_no_harmony, sil_harmony),
  Quality = c(
    ifelse(sil_no_harmony > 0.5, "Strong",
           ifelse(sil_no_harmony > 0.25, "Reasonable", "Weak")),
    ifelse(sil_harmony > 0.5, "Strong",
           ifelse(sil_harmony > 0.25, "Reasonable", "Weak"))
  )
)

write.csv(silhouette_comparison,
          "results/clustering/comparison/silhouette_scores.csv",
          row.names = FALSE)

################################################################################
# 8. Generate recommendation
################################################################################
cat("Step 8: Generating recommendation...\n")
cat("----------------------------------------------------------------\n")

# Decision criteria
criteria <- data.frame(
  Criterion = c(
    "Batch Mixing Score",
    "Silhouette Score",
    "Cluster Quality"
  ),
  No_Harmony = c(
    round(mixing_no_harmony, 3),
    round(sil_no_harmony, 3),
    length(unique(seurat_no_harmony$seurat_clusters))
  ),
  With_Harmony = c(
    round(mixing_harmony, 3),
    round(sil_harmony, 3),
    length(unique(seurat_harmony$seurat_clusters))
  )
)

cat("\n  Decision Criteria:\n")
print(criteria)
cat("\n")

# Generate recommendation
recommendation <- "No Harmony"  # Default recommendation
reasons <- c()

# Check mixing
if (mixing_no_harmony >= 0.4 && mixing_no_harmony <= 0.6) {
  reasons <- c(reasons, "✓ No Harmony shows good mixing (no batch effect detected)")
} else if (mixing_no_harmony > 0.7) {
  recommendation <- "With Harmony"
  reasons <- c(reasons, "⚠ Strong batch effect detected (poor mixing without Harmony)")
}

# Check silhouette
if (sil_no_harmony > sil_harmony) {
  reasons <- c(reasons, "✓ No Harmony has better cluster separation")
} else {
  reasons <- c(reasons, "✓ Harmony has better cluster separation")
}

# Over-correction check
if (mixing_harmony < 0.3) {
  recommendation <- "No Harmony"
  reasons <- c(reasons, "⚠ Harmony may be over-correcting (excessive mixing)")
}

cat("\n================================================================\n")
cat("RECOMMENDATION\n")
cat("================================================================\n\n")

cat("  Recommended method:", recommendation, "\n\n")
cat("  Reasons:\n")
for (reason in reasons) {
  cat("   ", reason, "\n")
}

cat("\n  Additional considerations:\n")
cat("    • No Harmony preserves more biological variation\n")
cat("    • Use Harmony only if strong batch effects exist\n")
cat("    • For KO vs WT comparison, preserving variation is crucial\n\n")

# Save recommendation
recommendation_report <- data.frame(
  Recommended_Method = recommendation,
  Mixing_Score_NoHarmony = mixing_no_harmony,
  Mixing_Score_Harmony = mixing_harmony,
  Silhouette_NoHarmony = sil_no_harmony,
  Silhouette_Harmony = sil_harmony,
  Reasoning = paste(reasons, collapse = "; ")
)

write.csv(recommendation_report,
          "results/clustering/comparison/recommendation.csv",
          row.names = FALSE)

cat("================================================================\n")
cat("Comparison Complete\n")
cat("================================================================\n")
cat("\nOutputs saved to: results/clustering/comparison/\n\n")
cat("Review these files:\n")
cat("  1. 01_umap_comparison.png - Visual comparison\n")
cat("  2. 02_cluster_composition_comparison.png - Cluster composition\n")
cat("  3. recommendation.csv - Final recommendation\n")
cat("  4. mixing_metrics.csv - Batch mixing scores\n")
cat("  5. silhouette_scores.csv - Cluster quality scores\n\n")
