#!/bin/bash
#SBATCH --job-name=install_packages
#SBATCH --output=logs/00_install_packages.out
#SBATCH --error=logs/00_install_packages.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

################################################################################
# Script: 00_install_packages.sh
# Description: Install all required R packages for the pipeline
# Note: Run this once before starting the analysis
################################################################################

set -euo pipefail

# Load required modules
# module purge
# module load R/4.3.0 || module load R || true

# Set working directory
WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
cd "${WORK_DIR}"

echo "================================================================"
echo "Installing R Packages"
echo "================================================================"
echo "Start time: $(date)"

# Run R script
# Rscript scripts/00_install_packages.R

echo "================================================================"
echo "Package Installation Complete"
echo "================================================================"
echo "End time: $(date)"
