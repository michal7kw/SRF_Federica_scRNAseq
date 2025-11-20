################################################################################
# Script: 05_pathway_enrichment.R
# Description: Gene Set Enrichment Analysis and pathway analysis
# Input: Differential expression results from step 04
# Output: Enrichment results, pathway plots
# Author: Bioinformatics Pipeline
# Date: 2025-11-14
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(DOSE)
  library(pathview)
  library(ggrepel)
  library(RColorBrewer)
})

# Set working directory
work_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Federica_top/SRF_Federica/SRF_Federica_scRNAseq/pipeline"
setwd(work_dir)

# Create output directories
dir.create("results/pathway_analysis/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/pathway_analysis/tables", recursive = TRUE, showWarnings = FALSE)

cat("\n================================================================\n")
cat("Pathway Enrichment Analysis\n")
cat("================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

################################################################################
# 1. Load differential expression results
################################################################################
cat("Step 1: Loading differential expression results...\n")
cat("----------------------------------------------------------------\n")

global_degs <- read.csv("results/diff_expression/tables/global_degs_KO_vs_WT.csv")

cat("  Total genes tested:", nrow(global_degs), "\n")
cat("  Significant DEGs (padj < 0.05):", sum(global_degs$p_val_adj < 0.05), "\n\n")

################################################################################
# 2. Prepare gene lists
################################################################################
cat("Step 2: Preparing gene lists for enrichment...\n")
cat("----------------------------------------------------------------\n")

# Filter significant DEGs
sig_degs <- global_degs %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5)

cat("  Significant DEGs (padj < 0.05, |log2FC| > 0.5):", nrow(sig_degs), "\n")

# Separate up and down regulated genes
up_genes <- sig_degs %>%
  filter(avg_log2FC > 0.5) %>%
  pull(gene)

down_genes <- sig_degs %>%
  filter(avg_log2FC < -0.5) %>%
  pull(gene)

cat("  Up-regulated in KO:", length(up_genes), "\n")
cat("  Down-regulated in KO:", length(down_genes), "\n\n")

# Convert gene symbols to Entrez IDs
convert_genes <- function(genes) {
  result <- bitr(genes,
                fromType = "SYMBOL",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Hs.eg.db)
  return(result)
}

cat("  Converting gene symbols to Entrez IDs...\n")
all_sig_entrez <- convert_genes(sig_degs$gene)
up_entrez <- convert_genes(up_genes)
down_entrez <- convert_genes(down_genes)

cat("    All DEGs converted:", nrow(all_sig_entrez), "/", nrow(sig_degs), "\n")
cat("    Up-regulated converted:", nrow(up_entrez), "/", length(up_genes), "\n")
cat("    Down-regulated converted:", nrow(down_entrez), "/", length(down_genes), "\n\n")

# Prepare ranked gene list for GSEA
ranked_genes <- global_degs %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC))

ranked_genes_entrez <- convert_genes(ranked_genes$gene)
ranked_genes <- ranked_genes %>%
  filter(gene %in% ranked_genes_entrez$SYMBOL) %>%
  left_join(ranked_genes_entrez, by = c("gene" = "SYMBOL"))

gene_list <- ranked_genes$avg_log2FC
names(gene_list) <- ranked_genes$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat("  Ranked gene list prepared:", length(gene_list), "genes\n\n")

################################################################################
# 3. Gene Ontology (GO) Enrichment - All DEGs
################################################################################
cat("Step 3: GO enrichment analysis (all DEGs)...\n")
cat("----------------------------------------------------------------\n")

# GO Biological Process
cat("  Running GO:BP enrichment...\n")
go_bp <- enrichGO(gene = all_sig_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  write.csv(go_bp@result, "results/pathway_analysis/tables/GO_BP_all_degs.csv", row.names = FALSE)
  cat("    GO:BP terms enriched:", nrow(go_bp@result[go_bp@result$p.adjust < 0.05, ]), "\n")
} else {
  cat("    No significant GO:BP terms found\n")
}

# GO Molecular Function
cat("  Running GO:MF enrichment...\n")
go_mf <- enrichGO(gene = all_sig_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
  write.csv(go_mf@result, "results/pathway_analysis/tables/GO_MF_all_degs.csv", row.names = FALSE)
  cat("    GO:MF terms enriched:", nrow(go_mf@result[go_mf@result$p.adjust < 0.05, ]), "\n")
} else {
  cat("    No significant GO:MF terms found\n")
}

# GO Cellular Component
cat("  Running GO:CC enrichment...\n")
go_cc <- enrichGO(gene = all_sig_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
  write.csv(go_cc@result, "results/pathway_analysis/tables/GO_CC_all_degs.csv", row.names = FALSE)
  cat("    GO:CC terms enriched:", nrow(go_cc@result[go_cc@result$p.adjust < 0.05, ]), "\n")
} else {
  cat("    No significant GO:CC terms found\n")
}

cat("\n")

################################################################################
# 4. KEGG Pathway Enrichment
################################################################################
cat("Step 4: KEGG pathway enrichment...\n")
cat("----------------------------------------------------------------\n")

kegg <- enrichKEGG(gene = all_sig_entrez$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)

if (!is.null(kegg) && nrow(kegg@result) > 0) {
  kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(kegg_readable@result, "results/pathway_analysis/tables/KEGG_all_degs.csv", row.names = FALSE)
  cat("  KEGG pathways enriched:", nrow(kegg@result[kegg@result$p.adjust < 0.05, ]), "\n\n")
} else {
  cat("  No significant KEGG pathways found\n\n")
}

################################################################################
# 5. Separate enrichment for up/down regulated genes
################################################################################
cat("Step 5: Separate enrichment for up/down regulated genes...\n")
cat("----------------------------------------------------------------\n")

# UP-REGULATED in KO
if (length(up_entrez$ENTREZID) > 10) {
  cat("  UP-regulated genes:\n")

  go_bp_up <- enrichGO(gene = up_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

  if (!is.null(go_bp_up) && nrow(go_bp_up@result) > 0) {
    write.csv(go_bp_up@result, "results/pathway_analysis/tables/GO_BP_up_in_KO.csv", row.names = FALSE)
    cat("    GO:BP terms:", nrow(go_bp_up@result[go_bp_up@result$p.adjust < 0.05, ]), "\n")
  }

  kegg_up <- enrichKEGG(gene = up_entrez$ENTREZID,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)

  if (!is.null(kegg_up) && nrow(kegg_up@result) > 0) {
    kegg_up_readable <- setReadable(kegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(kegg_up_readable@result, "results/pathway_analysis/tables/KEGG_up_in_KO.csv", row.names = FALSE)
    cat("    KEGG pathways:", nrow(kegg_up@result[kegg_up@result$p.adjust < 0.05, ]), "\n")
  }
}

# DOWN-REGULATED in KO
if (length(down_entrez$ENTREZID) > 10) {
  cat("  DOWN-regulated genes:\n")

  go_bp_down <- enrichGO(gene = down_entrez$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         readable = TRUE)

  if (!is.null(go_bp_down) && nrow(go_bp_down@result) > 0) {
    write.csv(go_bp_down@result, "results/pathway_analysis/tables/GO_BP_down_in_KO.csv", row.names = FALSE)
    cat("    GO:BP terms:", nrow(go_bp_down@result[go_bp_down@result$p.adjust < 0.05, ]), "\n")
  }

  kegg_down <- enrichKEGG(gene = down_entrez$ENTREZID,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

  if (!is.null(kegg_down) && nrow(kegg_down@result) > 0) {
    kegg_down_readable <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(kegg_down_readable@result, "results/pathway_analysis/tables/KEGG_down_in_KO.csv", row.names = FALSE)
    cat("    KEGG pathways:", nrow(kegg_down@result[kegg_down@result$p.adjust < 0.05, ]), "\n")
  }
}

cat("\n")

################################################################################
# 6. Gene Set Enrichment Analysis (GSEA)
################################################################################
cat("Step 6: Running GSEA...\n")
cat("----------------------------------------------------------------\n")

# GSEA GO
cat("  GSEA GO:BP...\n")
gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  write.csv(gsea_go@result, "results/pathway_analysis/tables/GSEA_GO_BP.csv", row.names = FALSE)
  cat("    GSEA GO:BP terms:", nrow(gsea_go@result[gsea_go@result$p.adjust < 0.05, ]), "\n")
}

# GSEA KEGG
cat("  GSEA KEGG...\n")
gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = 'hsa',
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
  gsea_kegg_readable <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(gsea_kegg_readable@result, "results/pathway_analysis/tables/GSEA_KEGG.csv", row.names = FALSE)
  cat("    GSEA KEGG pathways:", nrow(gsea_kegg@result[gsea_kegg@result$p.adjust < 0.05, ]), "\n")
}

cat("\n")

################################################################################
# 7. Visualization plots
################################################################################
cat("Step 7: Generating enrichment visualization plots...\n")
cat("----------------------------------------------------------------\n")

# GO BP bar plot
if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  p1 <- barplot(go_bp, showCategory = 20, title = "GO Biological Process (All DEGs)")
  ggsave("results/pathway_analysis/plots/01_GO_BP_barplot.png",
         plot = p1, width = 12, height = 8, dpi = 300)

  p2 <- dotplot(go_bp, showCategory = 20, title = "GO Biological Process (All DEGs)")
  ggsave("results/pathway_analysis/plots/02_GO_BP_dotplot.png",
         plot = p2, width = 12, height = 8, dpi = 300)
}

# KEGG bar plot
if (!is.null(kegg) && nrow(kegg@result) > 0) {
  p3 <- barplot(kegg_readable, showCategory = 20, title = "KEGG Pathways (All DEGs)")
  ggsave("results/pathway_analysis/plots/03_KEGG_barplot.png",
         plot = p3, width = 12, height = 8, dpi = 300)

  p4 <- dotplot(kegg_readable, showCategory = 20, title = "KEGG Pathways (All DEGs)")
  ggsave("results/pathway_analysis/plots/04_KEGG_dotplot.png",
         plot = p4, width = 12, height = 8, dpi = 300)
}

# GSEA plots
if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  p5 <- dotplot(gsea_go, showCategory = 20, title = "GSEA GO:BP", split = ".sign") +
        facet_grid(.~.sign)
  ggsave("results/pathway_analysis/plots/05_GSEA_GO_dotplot.png",
         plot = p5, width = 14, height = 8, dpi = 300)

  # Ridge plot
  p6 <- ridgeplot(gsea_go, showCategory = 20)
  ggsave("results/pathway_analysis/plots/06_GSEA_GO_ridgeplot.png",
         plot = p6, width = 12, height = 10, dpi = 300)
}

if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
  p7 <- dotplot(gsea_kegg_readable, showCategory = 20, title = "GSEA KEGG", split = ".sign") +
        facet_grid(.~.sign)
  ggsave("results/pathway_analysis/plots/07_GSEA_KEGG_dotplot.png",
         plot = p7, width = 14, height = 8, dpi = 300)
}

# Enrichment map (network)
if (!is.null(go_bp) && nrow(go_bp@result) > 5) {
  go_bp_simple <- pairwise_termsim(go_bp)
  p8 <- emapplot(go_bp_simple, showCategory = 30, cex_label_category = 0.6)
  ggsave("results/pathway_analysis/plots/08_GO_enrichment_map.png",
         plot = p8, width = 14, height = 12, dpi = 300)
}

# GSEA enrichment plots for top pathways
if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  top_gsea <- head(gsea_go@result[order(gsea_go@result$p.adjust), ], 5)

  for (i in 1:nrow(top_gsea)) {
    pathway_id <- top_gsea$ID[i]
    pathway_name <- gsub("[/:]", "_", substr(top_gsea$Description[i], 1, 50))

    tryCatch({
      p <- gseaplot2(gsea_go, geneSetID = pathway_id, title = top_gsea$Description[i])
      filename <- paste0("results/pathway_analysis/plots/GSEA_plot_", i, "_", pathway_name, ".png")
      ggsave(filename, plot = p, width = 10, height = 6, dpi = 300)
    }, error = function(e) {
      cat("    Error plotting GSEA for:", pathway_id, "\n")
    })
  }
}

cat("  Visualization plots saved\n\n")

################################################################################
# 8. Generate summary report
################################################################################
cat("Step 8: Generating summary report...\n")
cat("----------------------------------------------------------------\n")

summary_data <- data.frame(
  Analysis = character(),
  Count = numeric(),
  stringsAsFactors = FALSE
)

# Add results counts
summary_data <- rbind(summary_data, data.frame(
  Analysis = "Total DEGs (padj<0.05, |log2FC|>0.5)",
  Count = nrow(sig_degs)
))

if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  summary_data <- rbind(summary_data, data.frame(
    Analysis = "GO:BP terms enriched",
    Count = nrow(go_bp@result[go_bp@result$p.adjust < 0.05, ])
  ))
}

if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
  summary_data <- rbind(summary_data, data.frame(
    Analysis = "GO:MF terms enriched",
    Count = nrow(go_mf@result[go_mf@result$p.adjust < 0.05, ])
  ))
}

if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
  summary_data <- rbind(summary_data, data.frame(
    Analysis = "GO:CC terms enriched",
    Count = nrow(go_cc@result[go_cc@result$p.adjust < 0.05, ])
  ))
}

if (!is.null(kegg) && nrow(kegg@result) > 0) {
  summary_data <- rbind(summary_data, data.frame(
    Analysis = "KEGG pathways enriched",
    Count = nrow(kegg@result[kegg@result$p.adjust < 0.05, ])
  ))
}

if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  summary_data <- rbind(summary_data, data.frame(
    Analysis = "GSEA GO terms",
    Count = nrow(gsea_go@result[gsea_go@result$p.adjust < 0.05, ])
  ))
}

if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
  summary_data <- rbind(summary_data, data.frame(
    Analysis = "GSEA KEGG pathways",
    Count = nrow(gsea_kegg@result[gsea_kegg@result$p.adjust < 0.05, ])
  ))
}

write.csv(summary_data, "results/pathway_analysis/enrichment_summary.csv", row.names = FALSE)

cat("\n  Enrichment Summary:\n")
print(summary_data)

cat("\n================================================================\n")
cat("Pathway Enrichment Analysis Complete\n")
cat("================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nOutputs:\n")
cat("  - Enrichment tables: results/pathway_analysis/tables/\n")
cat("  - Visualization plots: results/pathway_analysis/plots/\n")
cat("  - Summary: results/pathway_analysis/enrichment_summary.csv\n")
cat("================================================================\n\n")
