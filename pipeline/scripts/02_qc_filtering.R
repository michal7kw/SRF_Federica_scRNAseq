################################################################################
# Script: 02_qc_filtering.R
# Description: Quality control and filtering of scRNA-seq data
# Input: CellRanger filtered matrices
# Output: Filtered Seurat objects, QC plots
# Author: Bioinformatics Pipeline
# Date: 2025-11-14
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

# Set working directory
work_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
setwd(work_dir)

# Create output directories
dir.create("results/qc/plots", recursive = TRUE, showWarnings = FALSE)

cat("\n================================================================\n")
cat("scRNA-seq Quality Control and Filtering\n")
cat("================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

################################################################################
# 1. Load CellRanger outputs
################################################################################
cat("Step 1: Loading CellRanger outputs...\n")
cat("----------------------------------------------------------------\n")

# Define paths to CellRanger outputs
ko_data <- Read10X(data.dir = "results/cellranger/KO/outs/filtered_feature_bc_matrix")
wt_data <- Read10X(data.dir = "results/cellranger/WT/outs/filtered_feature_bc_matrix")

# Create Seurat objects
ko <- CreateSeuratObject(counts = ko_data,
                         project = "BrainOrganoid",
                         min.cells = 3,
                         min.features = 200)
ko$condition <- "KO"

wt <- CreateSeuratObject(counts = wt_data,
                         project = "BrainOrganoid",
                         min.cells = 3,
                         min.features = 200)
wt$condition <- "WT"

cat("  KO cells:", ncol(ko), "\n")
cat("  WT cells:", ncol(wt), "\n\n")

################################################################################
# 2. Calculate QC metrics
################################################################################
cat("Step 2: Calculating QC metrics...\n")
cat("----------------------------------------------------------------\n")

# Calculate mitochondrial percentage
ko[["percent.mt"]] <- PercentageFeatureSet(ko, pattern = "^MT-")
wt[["percent.mt"]] <- PercentageFeatureSet(wt, pattern = "^MT-")

# Calculate ribosomal percentage
ko[["percent.ribo"]] <- PercentageFeatureSet(ko, pattern = "^RP[SL]")
wt[["percent.ribo"]] <- PercentageFeatureSet(wt, pattern = "^RP[SL]")

# Calculate hemoglobin percentage (contamination check)
ko[["percent.hb"]] <- PercentageFeatureSet(ko, pattern = "^HB[^(P)]")
wt[["percent.hb"]] <- PercentageFeatureSet(wt, pattern = "^HB[^(P)]")

cat("  QC metrics calculated for both samples\n\n")

################################################################################
# 3. Visualize QC metrics before filtering
################################################################################
cat("Step 3: Generating QC plots before filtering...\n")
cat("----------------------------------------------------------------\n")

# Combine for easier plotting
combined_unfiltered <- merge(ko, y = wt, add.cell.ids = c("KO", "WT"))

# Violin plots
p1 <- VlnPlot(combined_unfiltered,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "condition",
              ncol = 3,
              pt.size = 0)

ggsave("results/qc/plots/01_violin_qc_unfiltered.png",
       plot = p1, width = 15, height = 5, dpi = 300)

# Scatter plots to identify outliers
p2 <- FeatureScatter(combined_unfiltered,
                     feature1 = "nCount_RNA",
                     feature2 = "percent.mt",
                     group.by = "condition")

p3 <- FeatureScatter(combined_unfiltered,
                     feature1 = "nCount_RNA",
                     feature2 = "nFeature_RNA",
                     group.by = "condition")

p_combined <- p2 + p3
ggsave("results/qc/plots/02_scatter_qc_unfiltered.png",
       plot = p_combined, width = 12, height = 5, dpi = 300)

cat("  QC plots saved to results/qc/plots/\n\n")

################################################################################
# 4. Print QC statistics
################################################################################
cat("Step 4: QC Statistics Summary\n")
cat("----------------------------------------------------------------\n")

print_qc_stats <- function(seurat_obj, sample_name) {
  cat("\n", sample_name, ":\n", sep = "")
  cat("  Cells: ", ncol(seurat_obj), "\n", sep = "")
  cat("  Features (genes): ", nrow(seurat_obj), "\n", sep = "")
  cat("  Median UMI counts: ", median(seurat_obj$nCount_RNA), "\n", sep = "")
  cat("  Median genes/cell: ", median(seurat_obj$nFeature_RNA), "\n", sep = "")
  cat("  Median % MT: ", round(median(seurat_obj$percent.mt), 2), "%\n", sep = "")
  cat("  Median % Ribo: ", round(median(seurat_obj$percent.ribo), 2), "%\n", sep = "")
}

cat("\nBEFORE FILTERING:\n")
print_qc_stats(ko, "KO")
print_qc_stats(wt, "WT")

################################################################################
# 5. Apply filtering criteria
################################################################################
cat("\n\nStep 5: Applying filtering criteria...\n")
cat("----------------------------------------------------------------\n")

# Define filtering thresholds
# These are conservative defaults - adjust based on your QC plots
filter_params <- list(
  min_features = 200,      # Minimum genes per cell
  max_features = 7000,     # Maximum genes per cell (exclude potential doublets)
  min_counts = 500,        # Minimum UMI counts per cell
  max_counts = 50000,      # Maximum UMI counts per cell
  max_mt = 15              # Maximum mitochondrial percentage
)

cat("  Filtering parameters:\n")
cat("    Min genes/cell:", filter_params$min_features, "\n")
cat("    Max genes/cell:", filter_params$max_features, "\n")
cat("    Min UMI/cell:", filter_params$min_counts, "\n")
cat("    Max UMI/cell:", filter_params$max_counts, "\n")
cat("    Max % MT:", filter_params$max_mt, "%\n\n")

# Apply filters to KO
ko_filtered <- subset(ko,
                      subset = nFeature_RNA > filter_params$min_features &
                               nFeature_RNA < filter_params$max_features &
                               nCount_RNA > filter_params$min_counts &
                               nCount_RNA < filter_params$max_counts &
                               percent.mt < filter_params$max_mt)

# Apply filters to WT
wt_filtered <- subset(wt,
                      subset = nFeature_RNA > filter_params$min_features &
                               nFeature_RNA < filter_params$max_features &
                               nCount_RNA > filter_params$min_counts &
                               nCount_RNA < filter_params$max_counts &
                               percent.mt < filter_params$max_mt)

cat("AFTER FILTERING:\n")
print_qc_stats(ko_filtered, "KO")
print_qc_stats(wt_filtered, "WT")

# Calculate retention rate
ko_retention <- round(100 * ncol(ko_filtered) / ncol(ko), 1)
wt_retention <- round(100 * ncol(wt_filtered) / ncol(wt), 1)

cat("\n  Retention rates:\n")
cat("    KO: ", ko_retention, "% (", ncol(ko_filtered), "/", ncol(ko), ")\n", sep = "")
cat("    WT: ", wt_retention, "% (", ncol(wt_filtered), "/", ncol(wt), ")\n\n", sep = "")

################################################################################
# 6. Visualize QC metrics after filtering
################################################################################
cat("Step 6: Generating QC plots after filtering...\n")
cat("----------------------------------------------------------------\n")

# Combine filtered data
combined_filtered <- merge(ko_filtered, y = wt_filtered,
                          add.cell.ids = c("KO", "WT"))

# Violin plots after filtering
p4 <- VlnPlot(combined_filtered,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "condition",
              ncol = 3,
              pt.size = 0)

ggsave("results/qc/plots/03_violin_qc_filtered.png",
       plot = p4, width = 15, height = 5, dpi = 300)

cat("  Post-filtering QC plots saved\n\n")

################################################################################
# 7. Save filtered objects
################################################################################
cat("Step 7: Saving filtered Seurat objects...\n")
cat("----------------------------------------------------------------\n")

saveRDS(ko_filtered, file = "results/qc/ko_filtered.rds")
saveRDS(wt_filtered, file = "results/qc/wt_filtered.rds")
saveRDS(combined_filtered, file = "results/qc/combined_filtered.rds")

cat("  Filtered objects saved to results/qc/\n\n")

################################################################################
# 8. Generate QC summary report
################################################################################
cat("Step 8: Generating QC summary report...\n")
cat("----------------------------------------------------------------\n")

# Create summary table
qc_summary <- data.frame(
  Sample = c("KO", "WT"),
  Cells_Raw = c(ncol(ko), ncol(wt)),
  Cells_Filtered = c(ncol(ko_filtered), ncol(wt_filtered)),
  Retention_Percent = c(ko_retention, wt_retention),
  Median_Genes = c(median(ko_filtered$nFeature_RNA),
                   median(wt_filtered$nFeature_RNA)),
  Median_UMI = c(median(ko_filtered$nCount_RNA),
                 median(wt_filtered$nCount_RNA)),
  Median_MT_Percent = c(round(median(ko_filtered$percent.mt), 2),
                        round(median(wt_filtered$percent.mt), 2))
)

write.csv(qc_summary, "results/qc/qc_summary.csv", row.names = FALSE)

cat("\n  QC Summary Table:\n")
print(qc_summary)

cat("\n================================================================\n")
cat("Quality Control Complete\n")
cat("================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nOutputs:\n")
cat("  - Filtered objects: results/qc/*.rds\n")
cat("  - QC plots: results/qc/plots/\n")
cat("  - Summary table: results/qc/qc_summary.csv\n")
cat("================================================================\n\n")
