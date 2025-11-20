#!/bin/bash

################################################################################
# Interactive Method Selection Helper
################################################################################

WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
cd "${WORK_DIR}"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  Integration Method Selection Helper${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

# Check if both methods have been run
NO_HARMONY_FILE="results/clustering/seurat_integrated_no_harmony.rds"
HARMONY_FILE="results/clustering/seurat_integrated_harmony.rds"

if [ ! -f "${NO_HARMONY_FILE}" ] && [ ! -f "${HARMONY_FILE}" ]; then
    echo -e "${RED}Error: No integration results found!${NC}"
    echo ""
    echo "You need to run at least one integration method first:"
    echo "  sbatch scripts/03_integration_clustering_no_harmony.sh"
    echo "  OR"
    echo "  sbatch scripts/03_integration_clustering.sh"
    exit 1
fi

# Check what's available
echo -e "${YELLOW}Available Results:${NC}"
if [ -f "${NO_HARMONY_FILE}" ]; then
    echo -e "  ${GREEN}✓${NC} No Harmony: ${NO_HARMONY_FILE}"
fi
if [ -f "${HARMONY_FILE}" ]; then
    echo -e "  ${GREEN}✓${NC} With Harmony: ${HARMONY_FILE}"
fi
echo ""

# Check if comparison has been run
COMPARISON_FILE="results/clustering/comparison/recommendation.csv"

if [ -f "${COMPARISON_FILE}" ]; then
    echo -e "${GREEN}Comparison results found!${NC}"
    echo ""
    echo -e "${YELLOW}Automated Recommendation:${NC}"
    cat "${COMPARISON_FILE}"
    echo ""

    # Extract recommendation
    recommendation=$(tail -n 1 "${COMPARISON_FILE}" | cut -d',' -f1 | tr -d '"')
    mixing_no=$(tail -n 1 "${COMPARISON_FILE}" | cut -d',' -f2)
    mixing_harmony=$(tail -n 1 "${COMPARISON_FILE}" | cut -d',' -f3)

    echo -e "${YELLOW}Key Metrics:${NC}"
    echo "  Mixing score (No Harmony): ${mixing_no}"
    echo "  Mixing score (With Harmony): ${mixing_harmony}"
    echo ""
    echo -e "${GREEN}Recommended: ${recommendation}${NC}"
    echo ""
else
    echo -e "${YELLOW}No comparison results found yet.${NC}"
    echo ""

    if [ -f "${NO_HARMONY_FILE}" ] && [ -f "${HARMONY_FILE}" ]; then
        echo "Both methods have been run. Would you like to compare them now?"
        read -p "Run comparison script? (y/n): " run_comparison

        if [ "$run_comparison" = "y" ]; then
            echo ""
            echo "Running comparison..."
            Rscript scripts/compare_integration_methods.R
            echo ""
            echo "Comparison complete! Re-run this script to see results."
            exit 0
        fi
    else
        echo "Run both integration methods first, then compare:"
        echo "  1. sbatch scripts/03_integration_clustering_no_harmony.sh"
        echo "  2. sbatch scripts/03_integration_clustering.sh"
        echo "  3. sbatch scripts/compare_integration_methods.sh"
        echo "  4. Run this script again"
        exit 0
    fi
fi

# Ask user which method to use
echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  Which method do you want to use for downstream analysis?${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""
echo "1) No Harmony (standard Seurat workflow)"
echo "2) With Harmony (batch correction)"
echo "3) View comparison plots first"
echo "4) Exit"
echo ""
read -p "Enter your choice (1-4): " choice

case $choice in
    1)
        selected_method="No Harmony"
        selected_file="${NO_HARMONY_FILE}"
        ;;
    2)
        selected_method="With Harmony"
        selected_file="${HARMONY_FILE}"
        ;;
    3)
        echo ""
        echo "Opening comparison plots..."
        echo ""
        ls -lh results/clustering/comparison/*.png
        echo ""
        echo "View these files:"
        echo "  - 01_umap_comparison.png (most important)"
        echo "  - 02_cluster_composition_comparison.png"
        echo ""
        echo "Re-run this script when ready to choose."
        exit 0
        ;;
    4)
        echo "Exiting."
        exit 0
        ;;
    *)
        echo -e "${RED}Invalid choice${NC}"
        exit 1
        ;;
esac

# Verify file exists
if [ ! -f "${selected_file}" ]; then
    echo -e "${RED}Error: Selected file not found: ${selected_file}${NC}"
    echo "Run the corresponding integration script first."
    exit 1
fi

echo ""
echo -e "${GREEN}You selected: ${selected_method}${NC}"
echo -e "${GREEN}File: ${selected_file}${NC}"
echo ""

# Create symlink for downstream analysis
DEFAULT_FILE="results/clustering/seurat_integrated.rds"

echo "Creating symbolic link for downstream analysis..."
ln -sf "$(basename ${selected_file})" "${DEFAULT_FILE}"

if [ -L "${DEFAULT_FILE}" ]; then
    echo -e "${GREEN}✓ Symlink created: ${DEFAULT_FILE} → ${selected_file}${NC}"
else
    echo -e "${RED}✗ Failed to create symlink${NC}"
    exit 1
fi

echo ""
echo -e "${YELLOW}Summary:${NC}"
echo "  Method: ${selected_method}"
echo "  File: ${selected_file}"
echo "  Linked as: ${DEFAULT_FILE}"
echo ""
echo -e "${GREEN}Ready to proceed with downstream analysis!${NC}"
echo ""
echo "Next steps:"
echo "  1. Review clustering results:"
echo "     ls results/clustering/plots/"
echo ""
echo "  2. Continue with differential expression:"
echo "     ./run_pipeline.sh --step 4"
echo ""
echo -e "${BLUE}================================================================${NC}"
