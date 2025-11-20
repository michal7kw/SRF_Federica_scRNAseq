################################################################################
# Script: 04_differential_expression_harmony.R
# Description: Differential expression analysis between KO and WT conditions
# Input: Harmony integrated Seurat object from step 03
# Output: Differential expression results, volcano plots, heatmaps
# Author: Bioinformatics Pipeline
# Date: 2025-11-20
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(ggrepel)
})

# Set working directory
work_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
setwd(work_dir)

# Create output directories (separate for no-harmony)
OUTPUT_DIR <- "results/diff_expression/harmony"
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

cat("\n================================================================\n")
cat("Differential Expression Analysis: KO vs WT (Harmony Version)\n")
cat("================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

################################################################################
# 1. Load integrated data
################################################################################
cat("Step 1: Loading Harmony integrated Seurat object...\n")
cat("----------------------------------------------------------------\n")

seurat_obj <- readRDS("results/clustering/harmony/seurat_integrated.rds")

cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  KO cells:", sum(seurat_obj$condition == "KO"), "\n")
cat("  WT cells:", sum(seurat_obj$condition == "WT"), "\n")
cat("  Clusters:", length(unique(seurat_obj$seurat_clusters)), "\n\n")

################################################################################
# 2. Global differential expression: KO vs WT (all cells)
################################################################################
cat("Step 2: Global differential expression (KO vs WT, all cells)...\n")
cat("----------------------------------------------------------------\n")

Idents(seurat_obj) <- "condition"

# Find DEGs between KO and WT across all cells
global_degs <- FindMarkers(seurat_obj,
                           ident.1 = "KO",
                           ident.2 = "WT",
                           test.use = "wilcox",
                           logfc.threshold = 0.25,
                           min.pct = 0.1,
                           verbose = FALSE)

global_degs$gene <- rownames(global_degs)
global_degs <- global_degs[order(global_degs$p_val_adj), ]

# Add significance labels
global_degs$significance <- "NS"
global_degs$significance[global_degs$p_val_adj < 0.05 & global_degs$avg_log2FC > 0.5] <- "Up in KO"
global_degs$significance[global_degs$p_val_adj < 0.05 & global_degs$avg_log2FC < -0.5] <- "Down in KO"

# Save results
write.csv(global_degs, file.path(OUTPUT_DIR, "tables", "global_degs_KO_vs_WT.csv"), row.names = FALSE)

# Summary
n_up <- sum(global_degs$significance == "Up in KO")
n_down <- sum(global_degs$significance == "Down in KO")
cat("  Total DEGs (padj < 0.05, |log2FC| > 0.5):", n_up + n_down, "\n")
cat("  Up-regulated in KO:", n_up, "\n")
cat("  Down-regulated in KO:", n_down, "\n\n")

################################################################################
# IMPROVED VOLCANO PLOT
################################################################################
cat("  Creating enhanced volcano plot...\n")

# Select genes to label (top 15 up and down by log2FC)
genes_to_label <- c(
  head(global_degs %>% filter(significance == "Up in KO") %>% arrange(desc(avg_log2FC)) %>% pull(gene), 15),
  head(global_degs %>% filter(significance == "Down in KO") %>% arrange(avg_log2FC) %>% pull(gene), 15)
)

# Create improved volcano plot with custom colors
p1 <- EnhancedVolcano(global_degs,
                     lab = global_degs$gene,
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     selectLab = genes_to_label,
                     title = 'KO vs WT (All Cells - Harmony)',
                     subtitle = paste0(n_up, ' up-regulated, ', n_down, ' down-regulated'),
                     pCutoff = 0.05,
                     FCcutoff = 0.5,
                     pointSize = 2.5,
                     labSize = 4.5,
                     labCol = 'black',
                     labFace = 'bold',
                     boxedLabels = TRUE,
                     col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                     colAlpha = 0.7,
                     legendPosition = 'right',
                     legendLabSize = 12,
                     legendIconSize = 5.0,
                     drawConnectors = TRUE,
                     widthConnectors = 0.5,
                     colConnectors = 'black',
                     max.overlaps = 50,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE,
                     border = 'full',
                     borderWidth = 1.5,
                     borderColour = 'black')

ggsave(file.path(OUTPUT_DIR, "plots", "01_volcano_global_KO_vs_WT.png"),
       plot = p1, width = 14, height = 12, dpi = 300)

cat("  Global volcano plot saved\n\n")

# Top DEGs table
top_global_degs <- global_degs %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  arrange(p_val_adj) %>%
  head(50)

write.csv(top_global_degs, file.path(OUTPUT_DIR, "tables", "top50_global_degs.csv"), row.names = FALSE)

cat("  Top 10 upregulated genes in KO:\n")
up_genes <- global_degs %>%
  filter(avg_log2FC > 0.5, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(10)
print(up_genes[, c("gene", "avg_log2FC", "p_val_adj")])

cat("\n  Top 10 downregulated genes in KO:\n")
down_genes <- global_degs %>%
  filter(avg_log2FC < -0.5, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) %>%
  head(10)
print(down_genes[, c("gene", "avg_log2FC", "p_val_adj")])

cat("\n")

################################################################################
# 3. Cluster-specific differential expression
################################################################################
cat("Step 3: Cluster-specific differential expression...\n")
cat("----------------------------------------------------------------\n")

Idents(seurat_obj) <- "seurat_clusters"
clusters <- unique(seurat_obj$seurat_clusters)

# Initialize list to store results
cluster_degs_list <- list()
all_cluster_degs <- data.frame()

# Loop through each cluster
for (cluster in clusters) {
  cat("  Analyzing cluster", cluster, "...\n")

  # Subset cells in this cluster
  cluster_cells <- WhichCells(seurat_obj, idents = cluster)

  if (length(cluster_cells) < 10) {
    cat("    Skipping (too few cells)\n")
    next
  }

  # Get KO and WT cells in this cluster
  ko_cells <- intersect(cluster_cells, WhichCells(seurat_obj, expression = condition == "KO"))
  wt_cells <- intersect(cluster_cells, WhichCells(seurat_obj, expression = condition == "WT"))

  if (length(ko_cells) < 3 || length(wt_cells) < 3) {
    cat("    Skipping (insufficient cells in one condition)\n")
    next
  }

  # Perform DE test
  tryCatch({
    cluster_de <- FindMarkers(seurat_obj,
                             ident.1 = ko_cells,
                             ident.2 = wt_cells,
                             test.use = "wilcox",
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             verbose = FALSE)

    if (nrow(cluster_de) > 0) {
      cluster_de$gene <- rownames(cluster_de)
      cluster_de$cluster <- cluster
      cluster_de$significance <- "NS"
      cluster_de$significance[cluster_de$p_val_adj < 0.05 & cluster_de$avg_log2FC > 0.5] <- "Up in KO"
      cluster_de$significance[cluster_de$p_val_adj < 0.05 & cluster_de$avg_log2FC < -0.5] <- "Down in KO"

      cluster_degs_list[[as.character(cluster)]] <- cluster_de
      all_cluster_degs <- rbind(all_cluster_degs, cluster_de)

      n_sig <- sum(cluster_de$p_val_adj < 0.05 & abs(cluster_de$avg_log2FC) > 0.5)
      cat("    DEGs found:", n_sig, "\n")
    }
  }, error = function(e) {
    cat("    Error in DE analysis:", e$message, "\n")
  })
}

cat("\n")

# Save cluster-specific results
write.csv(all_cluster_degs, file.path(OUTPUT_DIR, "tables", "cluster_specific_degs.csv"), row.names = FALSE)

# Save individual cluster results
for (cluster in names(cluster_degs_list)) {
  filename <- file.path(OUTPUT_DIR, "tables", paste0("cluster_", cluster, "_degs.csv"))
  write.csv(cluster_degs_list[[cluster]], filename, row.names = FALSE)
}

cat("  Cluster-specific DEG tables saved\n\n")

################################################################################
# 4. Visualizations
################################################################################
cat("Step 4: Generating visualization plots...\n")
cat("----------------------------------------------------------------\n")

# Feature plots for top global DEGs
top_genes <- c(head(up_genes$gene, 6), head(down_genes$gene, 6))
top_genes <- top_genes[top_genes %in% rownames(seurat_obj)]

if (length(top_genes) > 0) {
  p2 <- FeaturePlot(seurat_obj, features = top_genes, ncol = 4)
  ggsave(file.path(OUTPUT_DIR, "plots", "02_feature_plots_top_degs.png"),
         plot = p2, width = 16, height = 8, dpi = 300)

  # Split by condition
  p3 <- FeaturePlot(seurat_obj, features = top_genes[1:min(6, length(top_genes))],
                    split.by = "condition", ncol = 4)
  ggsave(file.path(OUTPUT_DIR, "plots", "03_feature_plots_split_condition.png"),
         plot = p3, width = 16, height = 12, dpi = 300)
}

# Violin plots for top DEGs
p4 <- VlnPlot(seurat_obj, features = top_genes[1:min(9, length(top_genes))],
              group.by = "condition", ncol = 3, pt.size = 0)
ggsave(file.path(OUTPUT_DIR, "plots", "04_violin_plots_top_degs.png"),
       plot = p4, width = 15, height = 12, dpi = 300)

# Dot plot for top DEGs by cluster
p5 <- DotPlot(seurat_obj, features = top_genes, group.by = "seurat_clusters", split.by = "condition") +
      RotatedAxis() +
      ggtitle("Top DEGs by Cluster and Condition")
ggsave(file.path(OUTPUT_DIR, "plots", "05_dotplot_top_degs_by_cluster.png"),
       plot = p5, width = 14, height = 8, dpi = 300)

cat("  Feature plots generated\n\n")

# Heatmap of top DEGs
top_heatmap_genes <- global_degs %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(50) %>%
  pull(gene)

# Filter genes that are in scale.data
top_heatmap_genes <- top_heatmap_genes[top_heatmap_genes %in% rownames(seurat_obj)]

if (length(top_heatmap_genes) > 0) {
  png(file.path(OUTPUT_DIR, "plots", "06_heatmap_top50_degs.png"),
      width = 12, height = 14, units = "in", res = 300)
  DoHeatmap(seurat_obj,
            features = top_heatmap_genes,
            group.by = "condition",
            size = 3) +
    scale_fill_gradientn(colors = c("blue", "white", "red"))
  dev.off()
}

cat("  Heatmap generated\n\n")

################################################################################
# 5. Summary statistics
################################################################################
cat("Step 5: Generating summary statistics...\n")
cat("----------------------------------------------------------------\n")

# DEG counts per cluster
cluster_deg_summary <- all_cluster_degs %>%
  group_by(cluster) %>%
  summarise(
    Total_DEGs = sum(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5),
    Up_in_KO = sum(significance == "Up in KO"),
    Down_in_KO = sum(significance == "Down in KO")
  )

write.csv(cluster_deg_summary, file.path(OUTPUT_DIR, "tables", "cluster_deg_summary.csv"), row.names = FALSE)

cat("\n  DEGs per cluster:\n")
print(cluster_deg_summary)

# Create bar plot of DEG counts
library(reshape2)
cluster_deg_melt <- melt(cluster_deg_summary, id.vars = "cluster")
p6 <- ggplot(cluster_deg_melt, aes(x = cluster, y = value, fill = variable)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      ggtitle("DEG Counts by Cluster (Harmony)") +
      xlab("Cluster") +
      ylab("Number of DEGs") +
      scale_fill_manual(values = c("Total_DEGs" = "grey40",
                                     "Up_in_KO" = "red2",
                                     "Down_in_KO" = "royalblue")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "plots", "07_deg_counts_by_cluster.png"),
       plot = p6, width = 10, height = 6, dpi = 300)

cat("\n")

################################################################################
# 6. Generate comprehensive summary
################################################################################
cat("Step 6: Generating comprehensive summary...\n")
cat("----------------------------------------------------------------\n")

summary_report <- data.frame(
  Analysis = c("Integration Method", "Global DE (all cells)", "Total clusters analyzed",
               "Clusters with DEGs", "Total unique DEGs"),
  Result = c("No Harmony",
             paste0(n_up + n_down, " (", n_up, " up, ", n_down, " down)"),
             length(unique(seurat_obj$seurat_clusters)),
             length(unique(all_cluster_degs$cluster[all_cluster_degs$p_val_adj < 0.05])),
             length(unique(all_cluster_degs$gene[all_cluster_degs$p_val_adj < 0.05])))
)

write.csv(summary_report, file.path(OUTPUT_DIR, "summary_report.csv"), row.names = FALSE)

cat("\n  Summary Report:\n")
print(summary_report)

cat("\n================================================================\n")
cat("Differential Expression Analysis Complete (Harmony)\n")
cat("================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nOutputs:\n")
cat("  - Global DEGs:", file.path(OUTPUT_DIR, "tables", "global_degs_KO_vs_WT.csv"), "\n")
cat("  - Cluster DEGs:", file.path(OUTPUT_DIR, "tables", "cluster_specific_degs.csv"), "\n")
cat("  - Plots:", file.path(OUTPUT_DIR, "plots/"), "\n")
cat("  - Summary:", file.path(OUTPUT_DIR, "summary_report.csv"), "\n")
cat("\nNext steps:\n")
cat("  - Review top DEGs and their biological relevance\n")
cat("  - Proceed to pathway enrichment analysis (step 05)\n")
cat("================================================================\n\n")
