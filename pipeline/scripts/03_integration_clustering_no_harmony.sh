#!/bin/bash
#SBATCH --job-name=clustering_no_harmony
#SBATCH --output=logs/03_integration_clustering_no_harmony.out
#SBATCH --error=logs/03_integration_clustering_no_harmony.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

################################################################################
# Script: 03_integration_clustering_no_harmony.sh
# Description: Run clustering without Harmony integration
# Dependencies: QC filtered data from step 02
################################################################################

set -euo pipefail

# Load required modules
# module purge
# module load R/4.3.0 || module load R || true

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/scRNA_env

# Set working directory
WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
cd "${WORK_DIR}"

echo "================================================================"
echo "Running Integration and Clustering (NO Harmony)"
echo "================================================================"
echo "Start time: $(date)"

# Run R script
Rscript scripts/03_integration_clustering_no_harmony.R

echo "================================================================"
echo "Clustering Complete (No Harmony)"
echo "================================================================"
echo "End time: $(date)"
