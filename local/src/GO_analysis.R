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
out_dir <- snakemake@output[["out_dir"]]

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

egocc <- enrichGO(gene          = geneList,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont           = "CC",
                pAdjustMethod = "BH",  
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)

egomf <- enrichGO(gene          = geneList,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont           = "MF",
                pAdjustMethod = "BH",  
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)

egobp <- enrichGO(gene          = geneList,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",  
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = FALSE)

dir.create(out_dir)
setwd(out_dir)

barplot(egocc, showCategory = 20)
ggsave("CC.pdf")
barplot(egomf, showCategory = 20)
ggsave("MF.pdf")
barplot(egobp, showCategory = 20)
ggsave("BP.pdf")

setwd("..")
# call also for MF, BP, produce three separated plots (either use a directory for output or three separate filenames, as you prefer :))
# rbind the ego@result in a single dataframe after having addead a column 'ontology' with CC, MF or BP , use p.adjust to correct the nominal pvalue in a single run (https://www.biostars.org/p/12182/)

egoall_df <- rbind(egocc@result, egobp@result, egomf@result)
egoall_df$p.adjust <- p.adjust(egoall_df$pvalue, method='BH')
write.table(egoall_df, file = GO_r, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

# then print a single dataframe with all the information together
# we should understand how to do the same for GSEA...