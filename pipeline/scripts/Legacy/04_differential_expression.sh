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
#
# ⚠️  DEPRECATED - USE SPECIFIC VERSIONS INSTEAD:
# - 04_differential_expression_no_harmony.sh (recommended)
# - 04_differential_expression_harmony.sh
################################################################################

set -euo pipefail

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/scRNA_env

# Set working directory
WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
cd "${WORK_DIR}"

echo "================================================================"
echo "⚠️  WARNING: This script is DEPRECATED"
echo "================================================================"
echo ""
echo "Please use one of the improved versions instead:"
echo "  - 04_differential_expression_no_harmony.sh (recommended)"
echo "  - 04_differential_expression_harmony.sh"
echo ""
echo "Redirecting to no-harmony version..."
echo "================================================================"
echo ""

# Run the no-harmony version instead
sbatch scripts/04_differential_expression_no_harmony.sh

echo "Job submitted! Check logs/04_differential_expression_no_harmony.out"
