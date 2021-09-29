### GO enrichment analysis

library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

gene_res_f <- snakemake@input[["gene_list"]]
gene_univ_f <- snakemake@input[["gene_univ"]]
GO_r <- snakemake@output[["GO_r"]]
GO_barplot <- snakemake@output[["GO_barplot"]]
#threshold <- snakemake@wildcards[['thr']]

gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = FALSE)

#gene_univ_f <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/genes_residuals_universe.tsv"
gene_univ_df <- read.table(gene_univ_f, quote = "", sep = "\t", header = FALSE)

###order ##scegliere un treshold per discriminare quali geni tenere in base alle frequenze
#gene_freq_10 <- subset(gene_res_df, gene_res_df$Freq > threshold)
#geneList <- gene_freq_10$gene
geneList <- as.character(gene_res_df$V1)

### as.character universe
geneUni <- gene_univ_df$V1
geneUni <- as.character(geneUni)

#m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  #dplyr::select(gs_name, human_gene_symbol)

ego <- enrichGO(gene          = geneList,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)

write.table(ego@result, file = GO_r, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)
barplot(ego, showCategory = 20)
ggsave(GO_barplot)


