#!/bin/bash
#SBATCH --job-name=diff_expression
#SBATCH --output=logs/04_differential_expression.out
#SBATCH --error=logs/04_differential_expression.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

################################################################################
# Script: 04_differential_expression.sh
# Description: Run differential expression analysis
# Dependencies: Integrated data from step 03
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
echo "Running Differential Expression Analysis"
echo "================================================================"
echo "Start time: $(date)"

# Run R script
Rscript scripts/04_differential_expression.R

echo "================================================================"
echo "Differential Expression Analysis Complete"
echo "================================================================"
echo "End time: $(date)"
