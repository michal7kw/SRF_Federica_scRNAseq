#!/bin/bash

################################################################################
# Download CellRanger Reference Genome
# Human GRCh38 (2020-A release from 10x Genomics)
################################################################################

set -euo pipefail

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  CellRanger Reference Genome Download${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

# Set download location
GENOME_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome"
REF_NAME="refdata-gex-GRCh38-2020-A"
DOWNLOAD_URL="https://cf.10xgenomics.com/supp/cell-exp/${REF_NAME}.tar.gz"

echo -e "${YELLOW}Download location: ${GENOME_DIR}${NC}"
echo -e "${YELLOW}Reference: ${REF_NAME}${NC}"
echo ""

# Check if already exists
if [ -d "${GENOME_DIR}/${REF_NAME}" ]; then
    echo -e "${GREEN}Reference genome already exists at:${NC}"
    echo -e "${GREEN}${GENOME_DIR}/${REF_NAME}${NC}"
    echo ""
    echo -e "${YELLOW}Verifying reference structure...${NC}"

    # Check for required directories
    if [ -d "${GENOME_DIR}/${REF_NAME}/fasta" ] && \
       [ -d "${GENOME_DIR}/${REF_NAME}/genes" ] && \
       [ -d "${GENOME_DIR}/${REF_NAME}/star" ]; then
        echo -e "${GREEN}✓ Reference structure looks correct!${NC}"
        echo ""
        echo "You can use this path in your pipeline:"
        echo -e "${BLUE}${GENOME_DIR}/${REF_NAME}${NC}"
        echo ""
        echo "Update line 35 in scripts/01_cellranger_count.sh with this path."
        exit 0
    else
        echo -e "${RED}✗ Reference structure is incomplete or corrupted${NC}"
        echo "Recommendation: Delete and re-download"
        echo ""
        read -p "Delete and re-download? (y/n): " confirm
        if [ "$confirm" != "y" ]; then
            echo "Aborted."
            exit 1
        fi
        rm -rf "${GENOME_DIR}/${REF_NAME}"
    fi
fi

# Create genome directory if it doesn't exist
mkdir -p "${GENOME_DIR}"
cd "${GENOME_DIR}"

echo -e "${YELLOW}Step 1: Checking for existing download...${NC}"
echo "----------------------------------------------------------------"
if [ -f "${REF_NAME}.tar.gz" ]; then
    echo -e "${GREEN}Found existing tarball: ${REF_NAME}.tar.gz${NC}"
    echo "Skipping download..."
else
    echo "Tarball not found. Starting download..."
    echo ""
    echo -e "${YELLOW}Step 2: Downloading reference genome...${NC}"
    echo "----------------------------------------------------------------"
    echo "URL: ${DOWNLOAD_URL}"
    echo "Size: ~11 GB (this will take 15-30 minutes)"
    echo ""

    # Download with progress bar
    wget --continue --progress=bar:force "${DOWNLOAD_URL}"

    echo ""
    echo -e "${GREEN}✓ Download complete!${NC}"
fi

echo ""
echo -e "${YELLOW}Step 3: Extracting reference genome...${NC}"
echo "----------------------------------------------------------------"
echo "This will take 5-10 minutes..."

tar -xzf "${REF_NAME}.tar.gz"

echo -e "${GREEN}✓ Extraction complete!${NC}"

echo ""
echo -e "${YELLOW}Step 4: Verifying reference structure...${NC}"
echo "----------------------------------------------------------------"

# Verify directory structure
if [ -d "${REF_NAME}/fasta" ] && \
   [ -d "${REF_NAME}/genes" ] && \
   [ -d "${REF_NAME}/star" ] && \
   [ -f "${REF_NAME}/reference.json" ]; then
    echo -e "${GREEN}✓ Reference genome structure verified!${NC}"
else
    echo -e "${RED}✗ Reference structure verification failed!${NC}"
    echo "The extracted files may be incomplete."
    exit 1
fi

echo ""
echo -e "${YELLOW}Step 5: Checking file sizes...${NC}"
echo "----------------------------------------------------------------"
du -sh "${REF_NAME}"
du -sh "${REF_NAME}/fasta"
du -sh "${REF_NAME}/genes"
du -sh "${REF_NAME}/star"

echo ""
echo -e "${YELLOW}Cleaning up...${NC}"
echo "----------------------------------------------------------------"
read -p "Delete the tarball to save space? (y/n): " delete_tar
if [ "$delete_tar" = "y" ]; then
    rm "${REF_NAME}.tar.gz"
    echo -e "${GREEN}✓ Tarball deleted${NC}"
else
    echo "Tarball kept: ${REF_NAME}.tar.gz"
fi

echo ""
echo -e "${GREEN}================================================================${NC}"
echo -e "${GREEN}  Reference Genome Setup Complete!${NC}"
echo -e "${GREEN}================================================================${NC}"
echo ""
echo "Reference path:"
echo -e "${BLUE}${GENOME_DIR}/${REF_NAME}${NC}"
echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo "1. Copy the path above"
echo "2. Edit: pipeline/scripts/01_cellranger_count.sh"
echo "3. Update line 35 with the path above:"
echo ""
echo "   TRANSCRIPTOME=\"${GENOME_DIR}/${REF_NAME}\""
echo ""
echo "Or run this command to update automatically:"
echo ""
echo "sed -i 's|TRANSCRIPTOME=\"/path/to/refdata-gex-GRCh38-2020-A\"|TRANSCRIPTOME=\"${GENOME_DIR}/${REF_NAME}\"|' pipeline/scripts/01_cellranger_count.sh"
echo ""
echo -e "${GREEN}================================================================${NC}"
