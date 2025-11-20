#!/bin/bash

################################################################################
# Configuration Template for scRNA-seq Pipeline
# Copy this file and modify with your specific paths and parameters
################################################################################

# ==============================================================================
# PATHS AND DIRECTORIES
# ==============================================================================

# Reference genome path (REQUIRED - UPDATE THIS!)
export TRANSCRIPTOME="/path/to/refdata-gex-GRCh38-2020-A"

# Working directory (automatically set, but can override)
export WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq"

# Raw data directory
export RAW_DATA_DIR="${WORK_DIR}/Project_SessaA_2358_scRNAseq"

# ==============================================================================
# CELLRANGER PARAMETERS
# ==============================================================================

# Expected number of cells per sample
export EXPECTED_CELLS=10000

# Chemistry type (10x Genomics)
export CHEMISTRY="SC3Pv3"

# Number of cores for CellRanger
export CELLRANGER_CORES=16

# Memory for CellRanger (GB)
export CELLRANGER_MEM=64

# ==============================================================================
# QC FILTERING PARAMETERS
# ==============================================================================

# Minimum genes per cell
export MIN_FEATURES=200

# Maximum genes per cell (doublet threshold)
export MAX_FEATURES=7000

# Minimum UMI counts per cell
export MIN_COUNTS=500

# Maximum UMI counts per cell
export MAX_COUNTS=50000

# Maximum mitochondrial percentage
export MAX_MT_PERCENT=15

# ==============================================================================
# CLUSTERING PARAMETERS
# ==============================================================================

# Number of variable features to use
export N_VARIABLE_FEATURES=3000

# Number of PCA dimensions
export N_PCS=30

# Clustering resolution (can be comma-separated for multiple)
export CLUSTERING_RESOLUTION="0.1,0.2,0.3,0.5,0.8,1.0,1.2"

# Default resolution to use
export DEFAULT_RESOLUTION=0.5

# ==============================================================================
# DIFFERENTIAL EXPRESSION PARAMETERS
# ==============================================================================

# Log2 fold-change threshold
export LOGFC_THRESHOLD=0.25

# Minimum percentage of cells expressing gene
export MIN_PCT=0.1

# Adjusted p-value cutoff
export PADJ_CUTOFF=0.05

# ==============================================================================
# PATHWAY ENRICHMENT PARAMETERS
# ==============================================================================

# Organism (for KEGG)
export ORGANISM="hsa"  # Human

# P-value cutoff for enrichment
export ENRICHMENT_PVAL=0.05

# Q-value cutoff
export ENRICHMENT_QVAL=0.2

# ==============================================================================
# SLURM PARAMETERS
# ==============================================================================

# Account name
export SLURM_ACCOUNT="kubacki.michal"

# Partition
export SLURM_PARTITION="workq"

# Email for notifications (optional)
# export SLURM_EMAIL="your.email@example.com"
# export SLURM_EMAIL_TYPE="END,FAIL"

# ==============================================================================
# MODULE LOADING
# ==============================================================================

# Adjust these based on your HPC module system
export R_MODULE="R/4.3.0"
export CELLRANGER_MODULE="cellranger/7.1.0"

# ==============================================================================
# ANALYSIS OPTIONS
# ==============================================================================

# Perform batch correction with Harmony (true/false)
export USE_HARMONY=true

# Perform cell cycle regression (true/false)
export REGRESS_CELL_CYCLE=false

# Generate extended plots (true/false)
export GENERATE_EXTENDED_PLOTS=true

################################################################################
# END OF CONFIGURATION
################################################################################

echo "Configuration loaded successfully!"
echo "Transcriptome: ${TRANSCRIPTOME}"
echo "Expected cells: ${EXPECTED_CELLS}"
echo "Working directory: ${WORK_DIR}"
