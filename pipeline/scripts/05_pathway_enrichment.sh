#!/bin/bash
#SBATCH --job-name=pathway_enrichment
#SBATCH --output=logs/05_pathway_enrichment.out
#SBATCH --error=logs/05_pathway_enrichment.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

################################################################################
# Script: 05_pathway_enrichment.sh
# Description: Run pathway enrichment analysis
# Dependencies: Differential expression results from step 04
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
echo "Running Pathway Enrichment Analysis"
echo "================================================================"
echo "Start time: $(date)"

# Run R script
Rscript scripts/05_pathway_enrichment.R

echo "================================================================"
echo "Pathway Enrichment Analysis Complete"
echo "================================================================"
echo "End time: $(date)"
