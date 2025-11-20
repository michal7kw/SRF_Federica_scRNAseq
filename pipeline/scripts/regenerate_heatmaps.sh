#!/bin/bash
#SBATCH --job-name=regen_heatmaps
#SBATCH --output=logs/regenerate_heatmaps.out
#SBATCH --error=logs/regenerate_heatmaps.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

################################################################################
# Script: regenerate_heatmaps.sh
# Description: Fix and regenerate empty heatmap plots
# Bug: Original code used top5_markers$gene which was NULL
# Fix: Use gene column from markers dataframe
################################################################################

set -euo pipefail

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/scRNA_env

# Set working directory
WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
cd "${WORK_DIR}"

echo "================================================================"
echo "Regenerating Marker Gene Heatmaps"
echo "================================================================"
echo "Start time: $(date)"
echo ""
echo "Bug: Original heatmaps were empty because code used top5_markers\$gene"
echo "     but FindAllMarkers() returns genes as rownames, not a column"
echo ""
echo "Fix: Extract gene names correctly from markers dataframe"
echo "================================================================"
echo ""

# Run R script
Rscript scripts/regenerate_heatmaps.R

echo ""
echo "================================================================"
echo "Heatmap Regeneration Complete"
echo "================================================================"
echo "End time: $(date)"
echo ""
echo "Check output files:"
echo "  - results/clustering/*/plots/09_top_markers_heatmap.png"
echo "  - results/clustering/*/plots/09b_top10_markers_heatmap.png"
echo "  - results/clustering/*/plots/10_top_markers_dotplot.png"
echo ""
