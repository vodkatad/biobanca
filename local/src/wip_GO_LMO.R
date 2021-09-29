### GO for missing genes in CMS LMO

library(tidyverse)
library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

## carico meda per selezionare la classe di interesse e vsd per l'espressione

meda <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_f <- filter(meda_f, grepl("LMO_BASALE", type))

vsd <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE)
genes <- rownames(vsd)
rownames(vsd) <- NULL
vsd <- cbind(genes,vsd)
vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
colnames(vsd)[colnames(vsd) == 'genes'] <- 'symbol'
vsd_subset <- vsd[, colnames(vsd) %in% c(meda_f$sample_id_R, "symbol")]

## vsd_genes contiene i geni in ENTREZ ID che servono a CMS

vsd_genes <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/VSD_genes.tsv"
vsd_genes <- read.table(vsd_genes, quote = "", sep = "\t", header = TRUE)
vsd_genes <- vsd_genes[!grepl(",", vsd_genes$description),,drop = FALSE]
vsd_genes <- vsd_genes %>% na.exclude

vsd_f <- merge(vsd_genes, vsd_subset, by = "symbol")
vsd_f$symbol <- NULL
vsd_f <- vsd_f %>% remove_rownames %>% column_to_rownames(var="description")


# calcolo la deviazione standard delle righe e segnarsi quelle uguali a uno
# apply, matrice, direzione (per chiamare righe e colonne) e funzione
sds <- apply(vsd_f, 1, sd)
zerosd <- as.data.frame(names(sds[sds==0]))
colnames(zerosd) <- "GENES"
vsd_f <- tibble::rownames_to_column(vsd_f, "GENES")

vsd_GO <- merge(vsd_f, zerosd, by = "GENES")

### Prepare lists of gene for GO analysis
vsd_GO_f <- dplyr::select(vsd_GO, c("GENES"))
vsd_f2 <- dplyr::select(vsd_f, c("GENES"))

gene <- as.character(vsd_GO_f$GENES)
geneUni <- as.character(vsd_f2$GENES)


ego <- enrichGO(gene          = gene,
                universe      = geneUni,
                OrgDb         = "org.Hs.eg.db",
                keyType = "ENTREZID",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = FALSE)

#write.table(ego@result, file = "GO_LMO_results.tsv", quote = FALSE, sep = "\t", row.names = TRUE,
#            col.names = TRUE)
barplot(ego, showCategory = 20)
ggsave("barplot_go_lmo.pdf")

### scrivo un tsv che contenga i geni non espressi, serviranno poi per raggruppare gli stessi 
### geni nelle altri classi

write.table(vsd_GO, file = "notexpressed_genes_lmo.tsv", quote = FALSE, sep ="\t", col.names = TRUE,
            row.names = FALSE)

