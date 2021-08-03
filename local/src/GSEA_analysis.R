### GSEA enrichment analysis

library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

gene_res_f <- snakemake@input[["gene_res_freq"]]
GSEA_r <- snakemake@output[["GSEA_r"]]
GSEA_ridgeplot <- snakemake@output[["GSEA_ridgeplot"]]
gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE)
###order
geneList <- gene_res_df$Freq
names(geneList) <- as.character(gene_res_df$gene)
geneList <- sort(geneList, decreasing = TRUE)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, human_gene_symbol) ### altrimenti chiede gli id numerici



em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)

#GSEA_r <- write.table(em@result, quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(em@result, file = GSEA_r, quote = TRUE, sep = " ", row.names = TRUE,
            col.names = TRUE)

ridgeplot(em, showCategory = 20)
ggsave(GSEA_ridgeplot)


