#!/bin/bash

################################################################################
# Pipeline Status Checker
# Quick utility to check the status of pipeline execution
################################################################################

WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
cd "${WORK_DIR}"

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  Pipeline Status Check${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

# Check SLURM jobs
echo -e "${YELLOW}Current SLURM Jobs:${NC}"
echo "----------------------------------------------------------------"
squeue -u $(whoami) -o "%.18i %.20j %.8T %.10M %.6D %.20S" 2>/dev/null || echo "No jobs running"
echo ""

# Function to check if output exists
check_output() {
    local step=$1
    local file=$2
    local description=$3

    if [ -f "$file" ] || [ -d "$file" ]; then
        echo -e "  ${GREEN}✓${NC} $description"
        if [ -f "$file" ]; then
            size=$(du -h "$file" | cut -f1)
            echo "      Size: $size"
        fi
        return 0
    else
        echo -e "  ${RED}✗${NC} $description"
        return 1
    fi
}

# Check Step 0: R Packages
echo -e "${YELLOW}Step 0: R Package Installation${NC}"
echo "----------------------------------------------------------------"
if [ -f "logs/00_install_packages.out" ]; then
    if grep -q "SUCCESS" logs/00_install_packages.out 2>/dev/null; then
        echo -e "  ${GREEN}✓${NC} R packages installed successfully"
    elif grep -q "FAIL" logs/00_install_packages.out 2>/dev/null; then
        echo -e "  ${RED}✗${NC} Package installation had errors"
    else
        echo -e "  ${YELLOW}?${NC} Installation may be in progress"
    fi
else
    echo -e "  ${RED}✗${NC} Not started"
fi
echo ""

# Check Step 1: CellRanger
echo -e "${YELLOW}Step 1: CellRanger Count${NC}"
echo "----------------------------------------------------------------"
check_output "1" "results/cellranger/KO/outs/web_summary.html" "KO web summary"
check_output "1" "results/cellranger/WT/outs/web_summary.html" "WT web summary"
check_output "1" "results/cellranger/KO/outs/filtered_feature_bc_matrix" "KO count matrix"
check_output "1" "results/cellranger/WT/outs/filtered_feature_bc_matrix" "WT count matrix"
echo ""

# Check Step 2: QC
echo -e "${YELLOW}Step 2: Quality Control${NC}"
echo "----------------------------------------------------------------"
check_output "2" "results/qc/ko_filtered.rds" "KO filtered object"
check_output "2" "results/qc/wt_filtered.rds" "WT filtered object"
check_output "2" "results/qc/combined_filtered.rds" "Combined filtered object"
check_output "2" "results/qc/qc_summary.csv" "QC summary"
echo ""

# Check Step 3: Clustering
echo -e "${YELLOW}Step 3: Integration & Clustering${NC}"
echo "----------------------------------------------------------------"
check_output "3" "results/clustering/seurat_integrated.rds" "Integrated Seurat object"
check_output "3" "results/clustering/cluster_markers_all.csv" "Cluster markers"
check_output "3" "results/clustering/cluster_composition.csv" "Cluster composition"
echo ""

# Check Step 4: Differential Expression
echo -e "${YELLOW}Step 4: Differential Expression${NC}"
echo "----------------------------------------------------------------"
check_output "4" "results/diff_expression/tables/global_degs_KO_vs_WT.csv" "Global DEGs"
check_output "4" "results/diff_expression/tables/cluster_specific_degs.csv" "Cluster DEGs"
check_output "4" "results/diff_expression/summary_report.csv" "DE summary"
echo ""

# Check Step 5: Pathway Enrichment
echo -e "${YELLOW}Step 5: Pathway Enrichment${NC}"
echo "----------------------------------------------------------------"
check_output "5" "results/pathway_analysis/tables/GO_BP_all_degs.csv" "GO:BP results"
check_output "5" "results/pathway_analysis/tables/KEGG_all_degs.csv" "KEGG results"
check_output "5" "results/pathway_analysis/enrichment_summary.csv" "Enrichment summary"
echo ""

# Check logs for errors
echo -e "${YELLOW}Recent Log Activity:${NC}"
echo "----------------------------------------------------------------"
if [ -d "logs" ]; then
    latest_log=$(ls -t logs/*.out 2>/dev/null | head -1)
    if [ -n "$latest_log" ]; then
        echo "Latest log: $latest_log"
        echo "Last 5 lines:"
        tail -5 "$latest_log" | sed 's/^/  /'
    else
        echo "No log files found"
    fi
else
    echo "No logs directory"
fi
echo ""

# Summary
echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  Summary${NC}"
echo -e "${BLUE}================================================================${NC}"

completed=0
total=5

[ -f "results/cellranger/KO/outs/web_summary.html" ] && ((completed++))
[ -f "results/qc/combined_filtered.rds" ] && ((completed++))
[ -f "results/clustering/seurat_integrated.rds" ] && ((completed++))
[ -f "results/diff_expression/tables/global_degs_KO_vs_WT.csv" ] && ((completed++))
[ -f "results/pathway_analysis/enrichment_summary.csv" ] && ((completed++))

echo "Pipeline progress: ${completed}/${total} steps completed"

if [ $completed -eq $total ]; then
    echo -e "${GREEN}✓ Pipeline completed successfully!${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. Review results in results/ directory"
    echo "  2. Open CellRanger reports: results/cellranger/*/outs/web_summary.html"
    echo "  3. Examine DEGs: results/diff_expression/tables/global_degs_KO_vs_WT.csv"
elif [ $completed -eq 0 ]; then
    echo -e "${YELLOW}Pipeline not started yet${NC}"
    echo ""
    echo "To start:"
    echo "  1. Update reference path in scripts/01_cellranger_count.sh"
    echo "  2. Run: ./run_pipeline.sh --step 0 (install packages)"
    echo "  3. Run: ./run_pipeline.sh --step all"
else
    echo -e "${YELLOW}Pipeline in progress (${completed}/${total} steps done)${NC}"
    echo ""
    echo "Monitor with: watch -n 30 ./check_status.sh"
fi

echo ""
echo -e "${BLUE}================================================================${NC}"
