library(tidyverse)

cor_f <- snakemake@input[["simo_cor"]]
pdo_f <- snakemake@input[["pdo"]]
results <- snakemake@output[["res"]]

#cor_f <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/LMX-LMO_correlation_simo.tsv.gz"
cor <- read.table(cor_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#pdo_f <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
pdo <- read.table(pdo_f, quote = "",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pdo <- pdo %>% filter(buoni == TRUE)

cor <- cor[colnames(cor) %in% pdo$smodel]
cor <- cor[rownames(cor) %in% pdo$smodel,]

#if (all(colnames(cor)!=rownames(cor))) {
#  stop('Brutto llama!')
#}

write.table(cor, file = results, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")