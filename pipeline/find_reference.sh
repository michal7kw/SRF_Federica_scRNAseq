#!/bin/bash

################################################################################
# Find Existing CellRanger Reference Genomes on HPC
################################################################################

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  Searching for CellRanger References on HPC${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

echo -e "${YELLOW}Searching common locations...${NC}"
echo "This may take a few minutes..."
echo ""

# Search for CellRanger references
echo "Searching for refdata-gex-GRCh38-2020-A..."
find /beegfs -type d -name "refdata-gex-GRCh38-2020-A" 2>/dev/null | while read -r ref_path; do
    echo -e "${GREEN}✓ Found:${NC} ${ref_path}"

    # Verify it's a valid reference
    if [ -f "${ref_path}/reference.json" ]; then
        echo -e "  ${GREEN}✓ Valid CellRanger reference${NC}"
        echo ""
        echo "  Use this path in your pipeline:"
        echo -e "  ${BLUE}${ref_path}${NC}"
        echo ""
    else
        echo -e "  ${RED}✗ Invalid or incomplete reference${NC}"
        echo ""
    fi
done

echo ""
echo "Searching for any CellRanger references (refdata-gex-*)..."
find /beegfs -type d -name "refdata-gex-*" 2>/dev/null | while read -r ref_path; do
    echo -e "${YELLOW}Found:${NC} ${ref_path}"

    if [ -f "${ref_path}/reference.json" ]; then
        echo -e "  ${GREEN}✓ Valid CellRanger reference${NC}"

        # Check if it's human
        if [[ "${ref_path}" == *"GRCh38"* ]] || [[ "${ref_path}" == *"hg38"* ]]; then
            echo -e "  ${GREEN}✓ Human genome (GRCh38/hg38)${NC}"
            echo ""
            echo "  Recommended for your pipeline:"
            echo -e "  ${BLUE}${ref_path}${NC}"
        elif [[ "${ref_path}" == *"GRCm"* ]] || [[ "${ref_path}" == *"mm"* ]]; then
            echo -e "  ${YELLOW}⚠ Mouse genome (not suitable for human organoids)${NC}"
        fi
    fi
    echo ""
done

echo -e "${BLUE}================================================================${NC}"
echo ""
echo -e "${YELLOW}If no references found:${NC}"
echo "  Run: ./download_reference.sh"
echo ""
echo -e "${YELLOW}If a valid reference was found:${NC}"
echo "  1. Copy the path shown above"
echo "  2. Update pipeline/scripts/01_cellranger_count.sh line 35"
echo "  3. Replace with the correct path"
echo ""
echo -e "${BLUE}================================================================${NC}"
