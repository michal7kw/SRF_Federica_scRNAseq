#!/bin/bash

################################################################################
# Master Pipeline Script for Brain Organoid scRNA-seq Analysis
# Author: Bioinformatics Pipeline
# Date: 2025-11-14
################################################################################

set -euo pipefail

WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
cd "${WORK_DIR}"

# Create logs directory
mkdir -p logs

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  Brain Organoid scRNA-seq Analysis Pipeline${NC}"
echo -e "${BLUE}  Platform: 10x Genomics 3' Gene Expression v3${NC}"
echo -e "${BLUE}  Samples: KO vs WT${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help              Show this help message"
    echo "  -s, --step STEP         Run specific step (0-5, or 'all')"
    echo "  -i, --interactive       Interactive mode (prompts before each step)"
    echo ""
    echo "Steps:"
    echo "  0: Install R packages (run once)"
    echo "  1: CellRanger count (preprocessing)"
    echo "  2: Quality control and filtering"
    echo "  3: Integration and clustering"
    echo "  4: Differential expression analysis"
    echo "  5: Pathway enrichment analysis"
    echo "  all: Run all steps sequentially (1-5)"
    echo ""
    echo "Examples:"
    echo "  $0 --step 0              # Install packages"
    echo "  $0 --step all            # Run complete pipeline"
    echo "  $0 --step 3              # Run only step 3"
    echo "  $0 --interactive         # Interactive mode"
    exit 1
}

# Function to submit a job and optionally wait for it
submit_job() {
    local script=$1
    local step_name=$2
    local wait_for=$3

    echo -e "${YELLOW}Submitting step: ${step_name}${NC}" >&2
    echo "Script: ${script}" >&2

    if [ -f "scripts/${script}" ]; then
        if [ -n "${wait_for}" ]; then
            job_id=$(sbatch --dependency=afterok:${wait_for} scripts/${script} | awk '{print $4}')
        else
            job_id=$(sbatch scripts/${script} | awk '{print $4}')
        fi
        echo -e "${GREEN}Job submitted: ${job_id}${NC}" >&2
        echo "" >&2
        echo "${job_id}"
    else
        echo -e "${RED}Error: Script not found: scripts/${script}${NC}" >&2
        exit 1
    fi
}

# Function to run interactive mode
interactive_mode() {
    echo -e "${BLUE}Running in interactive mode${NC}"
    echo ""

    read -p "Install R packages? (y/n): " install_packages
    if [ "$install_packages" = "y" ]; then
        submit_job "00_install_packages.sh" "Install R packages" ""
        echo -e "${YELLOW}Wait for package installation to complete before proceeding.${NC}"
        return
    fi

    read -p "Run CellRanger count? (y/n): " run_cellranger
    if [ "$run_cellranger" = "y" ]; then
        job1=$(submit_job "01_cellranger_count.sh" "CellRanger count" "")
    fi

    read -p "Run QC and filtering? (y/n): " run_qc
    if [ "$run_qc" = "y" ]; then
        if [ -n "${job1:-}" ]; then
            job2=$(submit_job "02_qc_filtering.sh" "QC and filtering" "$job1")
        else
            job2=$(submit_job "02_qc_filtering.sh" "QC and filtering" "")
        fi
    fi

    read -p "Run integration and clustering? (y/n): " run_integration
    if [ "$run_integration" = "y" ]; then
        if [ -n "${job2:-}" ]; then
            job3=$(submit_job "03_integration_clustering_no_harmony.sh" "Integration and clustering (No-Harmony)" "$job2")
        else
            job3=$(submit_job "03_integration_clustering_no_harmony.sh" "Integration and clustering (No-Harmony)" "")
        fi
    fi

    read -p "Run differential expression? (y/n): " run_de
    if [ "$run_de" = "y" ]; then
        if [ -n "${job3:-}" ]; then
            job4=$(submit_job "04_differential_expression.sh" "Differential expression" "$job3")
        else
            job4=$(submit_job "04_differential_expression.sh" "Differential expression" "")
        fi
    fi

    read -p "Run pathway enrichment? (y/n): " run_pathway
    if [ "$run_pathway" = "y" ]; then
        if [ -n "${job4:-}" ]; then
            job5=$(submit_job "05_pathway_enrichment.sh" "Pathway enrichment" "$job4")
        else
            job5=$(submit_job "05_pathway_enrichment.sh" "Pathway enrichment" "")
        fi
    fi
}

# Parse command line arguments
STEP=""
INTERACTIVE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            ;;
        -s|--step)
            STEP="$2"
            shift 2
            ;;
        -i|--interactive)
            INTERACTIVE=true
            shift
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            usage
            ;;
    esac
done

# Run interactive mode if requested
if [ "$INTERACTIVE" = true ]; then
    interactive_mode
    exit 0
fi

# Run specific step or all steps
case ${STEP} in
    0)
        echo -e "${BLUE}Installing R packages...${NC}"
        submit_job "00_install_packages.sh" "Install R packages" ""
        echo -e "${YELLOW}Note: Wait for this to complete before running analysis steps.${NC}"
        ;;
    1)
        echo -e "${BLUE}Running CellRanger count...${NC}"
        submit_job "01_cellranger_count.sh" "CellRanger count" ""
        ;;
    2)
        echo -e "${BLUE}Running QC and filtering...${NC}"
        submit_job "02_qc_filtering.sh" "QC and filtering" ""
        ;;
    3)
        echo -e "${BLUE}Running integration and clustering (No-Harmony)...${NC}"
        submit_job "03_integration_clustering_no_harmony.sh" "Integration and clustering (No-Harmony)" ""
        ;;
    4)
        echo -e "${BLUE}Running differential expression...${NC}"
        submit_job "04_differential_expression.sh" "Differential expression" ""
        ;;
    5)
        echo -e "${BLUE}Running pathway enrichment...${NC}"
        submit_job "05_pathway_enrichment.sh" "Pathway enrichment" ""
        ;;
    all)
        echo -e "${BLUE}Running complete pipeline...${NC}"
        echo -e "${YELLOW}Note: Make sure R packages are already installed (step 0)${NC}"
        echo ""

        job1=$(submit_job "01_cellranger_count.sh" "CellRanger count" "")
        job2=$(submit_job "02_qc_filtering.sh" "QC and filtering" "$job1")
        job3=$(submit_job "03_integration_clustering_no_harmony.sh" "Integration and clustering (No-Harmony)" "$job2")
        job4=$(submit_job "04_differential_expression.sh" "Differential expression" "$job3")
        job5=$(submit_job "05_pathway_enrichment.sh" "Pathway enrichment" "$job4")

        echo -e "${GREEN}All jobs submitted successfully!${NC}"
        echo ""
        echo "Job dependencies:"
        echo "  CellRanger: $job1"
        echo "  QC: $job2 (depends on $job1)"
        echo "  Clustering: $job3 (depends on $job2)"
        echo "  Diff Expression: $job4 (depends on $job3)"
        echo "  Pathway: $job5 (depends on $job4)"
        echo ""
        echo "Monitor jobs with: squeue -u $(whoami)"
        ;;
    "")
        echo -e "${RED}Error: No step specified${NC}"
        usage
        ;;
    *)
        echo -e "${RED}Error: Invalid step: ${STEP}${NC}"
        usage
        ;;
esac

echo ""
echo -e "${GREEN}Pipeline submission complete!${NC}"
echo -e "Monitor jobs with: ${BLUE}squeue -u $(whoami)${NC}"
echo -e "Check logs in: ${BLUE}logs/${NC}"
echo ""
