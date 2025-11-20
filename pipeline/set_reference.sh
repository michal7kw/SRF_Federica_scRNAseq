#!/bin/bash

################################################################################
# Update CellRanger Reference Path in Pipeline
################################################################################

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

SCRIPT_PATH="scripts/01_cellranger_count.sh"

echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  Update CellRanger Reference Path${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

# Check if script exists
if [ ! -f "${SCRIPT_PATH}" ]; then
    echo -e "${RED}Error: ${SCRIPT_PATH} not found!${NC}"
    echo "Make sure you're running this from the pipeline directory."
    exit 1
fi

# Get current path
current_path=$(grep "^TRANSCRIPTOME=" "${SCRIPT_PATH}" | cut -d'"' -f2)
echo -e "${YELLOW}Current reference path:${NC}"
echo "  ${current_path}"
echo ""

# Ask for new path
echo -e "${YELLOW}Enter the full path to your CellRanger reference:${NC}"
echo "(e.g., /beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/refdata-gex-GRCh38-2020-A)"
echo ""
read -p "Path: " new_path

# Validate path exists
if [ ! -d "${new_path}" ]; then
    echo -e "${RED}Error: Directory does not exist: ${new_path}${NC}"
    exit 1
fi

# Validate it's a CellRanger reference
if [ ! -f "${new_path}/reference.json" ]; then
    echo -e "${YELLOW}Warning: reference.json not found in this directory.${NC}"
    echo "This may not be a valid CellRanger reference."
    read -p "Continue anyway? (y/n): " continue_anyway
    if [ "$continue_anyway" != "y" ]; then
        echo "Aborted."
        exit 1
    fi
fi

# Update the script
echo ""
echo -e "${YELLOW}Updating ${SCRIPT_PATH}...${NC}"

# Create backup
cp "${SCRIPT_PATH}" "${SCRIPT_PATH}.backup"
echo -e "${GREEN}✓ Backup created: ${SCRIPT_PATH}.backup${NC}"

# Update the path
sed -i "s|^TRANSCRIPTOME=.*|TRANSCRIPTOME=\"${new_path}\"|" "${SCRIPT_PATH}"

# Verify update
updated_path=$(grep "^TRANSCRIPTOME=" "${SCRIPT_PATH}" | cut -d'"' -f2)

if [ "${updated_path}" = "${new_path}" ]; then
    echo -e "${GREEN}✓ Reference path updated successfully!${NC}"
    echo ""
    echo -e "${YELLOW}New path:${NC}"
    echo "  ${updated_path}"
    echo ""
    echo -e "${GREEN}Your pipeline is now ready to run!${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. ./run_pipeline.sh --step 0  (install R packages, if not done)"
    echo "  2. ./run_pipeline.sh --step all  (run full pipeline)"
else
    echo -e "${RED}Error: Update failed!${NC}"
    echo "Restoring backup..."
    mv "${SCRIPT_PATH}.backup" "${SCRIPT_PATH}"
    exit 1
fi

echo ""
echo -e "${BLUE}================================================================${NC}"
