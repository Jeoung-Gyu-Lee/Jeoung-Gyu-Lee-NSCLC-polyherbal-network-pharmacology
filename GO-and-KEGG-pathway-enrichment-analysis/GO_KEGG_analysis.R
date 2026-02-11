# 1.Install and load Required Packages
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler", dependencies = TRUE)
BiocManager::install("enrichplot", dependencies = TRUE)
BiocManager::install("org.Hs.eg.db", dependencies = TRUE)
BiocManager::install("DOSE", dependencies = TRUE)

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)

# 2. GO Enrichment Function
GO_enrichment <- function(gene_list, ont_type = c("BP", "MF", "CC"), padj_cutoff = 0.05) {
  enrichGO(gene          = gene_list,
           OrgDb         = org.Hs.eg.db,
           keyType       = "ENTREZID",
           ont           = "ALL",
           pAdjustMethod = "BH",
           pvalueCutoff  = padj_cutoff,
           readable      = TRUE)
}

# 3. KEGG Enrichment Function
KEGG_enrichment <- function(gene_list, padj_cutoff = 0.05) {
  kegg <- enrichKEGG(gene          = gene_list,
                     organism      = "hsa",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = padj_cutoff)
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  return(kegg)
}

# 4. Run Enrichment Analysis
# Load gene list and translate to ENTREZID for analysis
file_path <- file.path(getwd(), "data", "Common_target.csv")
gene_list <- read.csv(file_path, header = FALSE)
colnames(gene_list)[1] <- "SYMBOL"
gene_list <- bitr(gene_list$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO (BP, MF, CC) and KEGG Enrichment
GO <- GO_enrichment(gene_list$ENTREZID, ont_type = "All", padj_cutoff = 0.05)
KEGG <- KEGG_enrichment(gene_list$ENTREZID, padj_cutoff = 0.05)

# 5. Save results
result_dir <- file.path(getwd(), "results")
setwd(result_dir)
write.csv(as.data.frame(GO), "GO_enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(KEGG), "KEGG_enrichment.csv", row.names = FALSE)

# 6. Visualization
# GO plot
GO_dot <- dotplot(
      GO,
      x            = "Count",
      color        = "p.adjust",
      showCategory = 10,
      split        = "ONTOLOGY" 
    )

pdf("Dotplot_GO.pdf", width  = 8, height = 10)
print(GO_dot +
    facet_grid(
      ONTOLOGY ~ .,
      scales = "free_y",
      space  = "free_y"
    ) +
    labs(
      title = "",
      color = "-log10(p.adjust)"
    ) +
    theme_minimal() +
    theme(
      axis.text.y      = element_text(colour = "black"),
      panel.border     = element_rect(colour    = "black",
                                      fill      = NA,
                                      linewidth = 0.8),
      panel.background = element_rect(fill   = "white",
                                      colour = NA),
      strip.background = element_rect(fill      = "grey90",
                                      colour    = "black",
                                      linewidth = 0.7),
      strip.text.y     = element_text(angle = 0,
                                      face  = "bold",
                                      size  = 12),
      strip.placement  = "outside",
      panel.spacing.y  = unit(0.5, "lines"),
      legend.position  = "right",
      plot.title       = element_text(hjust = 0.5,
                                      face  = "bold",
                                      size  = 1)
    )
)
dev.off()

# KEGG plot
pdf("Barplots_KEGG.pdf", width = 8, height = 7)
print(barplot(KEGG, showCategory = 20))
dev.off()

# Emapplot and Cnetplot
KEGG_sim <- pairwise_termsim(KEGG)

pdf("cnetplot_KEGG.pdf", width = 12, height = 10)
cnetplot(KEGG_sim, showCategory = 20, circular = FALSE, colorEdge = TRUE)
dev.off()

pdf("emapplot_KEGG.pdf", width = 10, height = 8)
emapplot(KEGG_sim, showCategory = 20, layout = "kk")
dev.off()