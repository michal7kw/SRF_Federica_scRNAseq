#!/bin/bash
#SBATCH --job-name=compare_methods
#SBATCH --output=logs/compare_integration_methods.out
#SBATCH --error=logs/compare_integration_methods.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

################################################################################
# Script: compare_integration_methods.sh
# Description: Compare Harmony vs No-Harmony results
################################################################################

set -euo pipefail

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/scRNA_env

# Set working directory
WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
cd "${WORK_DIR}"

echo "================================================================"
echo "Comparing Integration Methods"
echo "================================================================"
echo "Start time: $(date)"

# Run comparison script
Rscript scripts/compare_integration_methods.R

echo "================================================================"
echo "Comparison Complete"
echo "================================================================"
echo "End time: $(date)"
