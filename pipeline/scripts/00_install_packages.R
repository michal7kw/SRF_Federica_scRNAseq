################################################################################
# Script: 00_install_packages.R
# Description: Install required R packages for scRNA-seq analysis pipeline
# Author: Bioinformatics Pipeline
# Date: 2025-11-14
################################################################################

cat("\n================================================================\n")
cat("Installing Required R Packages\n")
cat("================================================================\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Function to install packages if not already installed
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  } else {
    cat("Already installed:", pkg, "\n")
  }
}

cat("Step 1: Installing CRAN packages...\n")
cat("----------------------------------------------------------------\n")

cran_packages <- c(
  "dplyr",
  "ggplot2",
  "patchwork",
  "Matrix",
  "RColorBrewer",
  "gridExtra",
  "reshape2",
  "devtools",
  "BiocManager"
)

for (pkg in cran_packages) {
  install_if_missing(pkg, bioc = FALSE)
}

cat("\n")
cat("Step 2: Installing Bioconductor packages...\n")
cat("----------------------------------------------------------------\n")

bioc_packages <- c(
  "BiocGenerics",
  "S4Vectors",
  "IRanges",
  "GenomeInfoDb",
  "GenomicRanges",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  "DOSE",
  "pathview"
)

for (pkg in bioc_packages) {
  install_if_missing(pkg, bioc = TRUE)
}

cat("\n")
cat("Step 3: Installing Seurat and related packages...\n")
cat("----------------------------------------------------------------\n")

# Install Seurat
if (!requireNamespace("Seurat", quietly = TRUE)) {
  cat("Installing Seurat (this may take a while)...\n")
  install.packages("Seurat", dependencies = TRUE)
} else {
  cat("Already installed: Seurat\n")
}

# Install SeuratObject
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  cat("Installing SeuratObject...\n")
  install.packages("SeuratObject", dependencies = TRUE)
} else {
  cat("Already installed: SeuratObject\n")
}

cat("\n")
cat("Step 4: Installing additional analysis packages...\n")
cat("----------------------------------------------------------------\n")

# Harmony for integration
if (!requireNamespace("harmony", quietly = TRUE)) {
  cat("Installing harmony...\n")
  install.packages("harmony", dependencies = TRUE)
} else {
  cat("Already installed: harmony\n")
}

# Clustree for cluster visualization
if (!requireNamespace("clustree", quietly = TRUE)) {
  cat("Installing clustree...\n")
  install.packages("clustree", dependencies = TRUE)
} else {
  cat("Already installed: clustree\n")
}

# EnhancedVolcano for volcano plots
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  cat("Installing EnhancedVolcano...\n")
  BiocManager::install("EnhancedVolcano", update = FALSE, ask = FALSE)
} else {
  cat("Already installed: EnhancedVolcano\n")
}

# pheatmap for heatmaps
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  cat("Installing pheatmap...\n")
  install.packages("pheatmap", dependencies = TRUE)
} else {
  cat("Already installed: pheatmap\n")
}

cat("\n================================================================\n")
cat("Package Installation Complete\n")
cat("================================================================\n\n")

cat("Verifying installation...\n")
cat("----------------------------------------------------------------\n")

# Verify key packages
key_packages <- c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "harmony", "clustree",
  "clusterProfiler", "enrichplot", "org.Hs.eg.db", "EnhancedVolcano"
)

all_installed <- TRUE
for (pkg in key_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "\n")
  } else {
    cat("✗", pkg, "- INSTALLATION FAILED\n")
    all_installed <- FALSE
  }
}

cat("\n")
if (all_installed) {
  cat("SUCCESS: All packages installed successfully!\n")
  cat("You are ready to run the scRNA-seq analysis pipeline.\n")
} else {
  cat("WARNING: Some packages failed to install.\n")
  cat("Please review error messages above and install manually.\n")
}

cat("\n================================================================\n\n")
