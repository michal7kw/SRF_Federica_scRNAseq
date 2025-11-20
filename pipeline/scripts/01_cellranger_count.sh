#!/bin/bash
#SBATCH --job-name=cellranger_count
#SBATCH --output=logs/01_cellranger_count.out
#SBATCH --error=logs/01_cellranger_count.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

################################################################################
# Script: 01_cellranger_count.sh
# Description: Process raw scRNA-seq data using CellRanger count
# Platform: 10x Genomics 3' Gene Expression v3
# Author: Bioinformatics Pipeline
# Date: 2025-11-14
################################################################################

set -euo pipefail

# Load required modules
# module purge
# module load cellranger/7.1.0 || module load cellranger || true

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/scRNA_env

# Define paths
WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq"
RAW_DATA_DIR="${WORK_DIR}/Project_SessaA_2358_scRNAseq"
OUTPUT_DIR="${WORK_DIR}/pipeline/results/cellranger"
TRANSCRIPTOME="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/refdata-gex-GRCh38-2020-A"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

echo "================================================================"
echo "CellRanger Count Pipeline - Brain Organoid scRNA-seq"
echo "================================================================"
echo "Start time: $(date)"
echo "Platform: 10x Genomics 3' Gene Expression v3"
echo "Reference: Human GRCh38"
echo "----------------------------------------------------------------"

################################################################################
# Sample 1: KO (Knockout)
################################################################################
echo ""
echo "Processing Sample: KO"
echo "----------------------------------------------------------------"

cellranger count \
    --id=KO \
    --transcriptome="${TRANSCRIPTOME}" \
    --fastqs="${RAW_DATA_DIR}/Sample_KO" \
    --sample=KO \
    --chemistry=SC3Pv3 \
    --create-bam=true \
    --localcores=16 \
    --localmem=64 \
    --expect-cells=10000

echo "KO sample processing completed: $(date)"

################################################################################
# Sample 2: WT (Wild-type)
################################################################################
echo ""
echo "Processing Sample: WT"
echo "----------------------------------------------------------------"

cellranger count \
    --id=WT \
    --transcriptome="${TRANSCRIPTOME}" \
    --fastqs="${RAW_DATA_DIR}/Sample_WT" \
    --sample=WT \
    --chemistry=SC3Pv3 \
    --create-bam=true \
    --localcores=16 \
    --localmem=64 \
    --expect-cells=10000

echo "WT sample processing completed: $(date)"

################################################################################
# Generate summary statistics
################################################################################
echo ""
echo "================================================================"
echo "CellRanger Processing Complete"
echo "================================================================"
echo "End time: $(date)"
echo ""
echo "Output directories:"
echo "  KO: ${OUTPUT_DIR}/KO"
echo "  WT: ${OUTPUT_DIR}/WT"
echo ""
echo "Key outputs:"
echo "  - Web summary: */outs/web_summary.html"
echo "  - Filtered matrix: */outs/filtered_feature_bc_matrix/"
echo "  - BAM file: */outs/possorted_genome_bam.bam"
echo "================================================================"

# Create symbolic links for easier access
ln -sf "${OUTPUT_DIR}/KO/outs/web_summary.html" "${OUTPUT_DIR}/KO_web_summary.html"
ln -sf "${OUTPUT_DIR}/WT/outs/web_summary.html" "${OUTPUT_DIR}/WT_web_summary.html"

echo "Pipeline step 01 completed successfully!"
