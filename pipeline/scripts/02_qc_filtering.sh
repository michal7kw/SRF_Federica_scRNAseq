#!/bin/bash
#SBATCH --job-name=qc_filtering
#SBATCH --output=logs/02_qc_filtering.out
#SBATCH --error=logs/02_qc_filtering.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

################################################################################
# Script: 02_qc_filtering.sh
# Description: Run QC and filtering on scRNA-seq data
# Dependencies: CellRanger output from step 01
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
echo "Running QC and Filtering"
echo "================================================================"
echo "Start time: $(date)"

# Run R script
Rscript scripts/02_qc_filtering.R

echo "================================================================"
echo "QC and Filtering Complete"
echo "================================================================"
echo "End time: $(date)"
