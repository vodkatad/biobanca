## pdo buoni filter for correlation

library(tidyverse)

lmo_f <- snakemake@input[["lmo"]]
lmx_f <- snakemake@input[["lmx"]]
pdo_f <- snakemake@input[["pdo"]]
lmo_r <- snakemake@output[["lmo_res"]]
lmx_r <- snakemake@output[["lmx_res"]]

#lmo <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/LMO_BASALE_mean_gene_genealogyall.tsv.gz"
lmo <- read.table(lmo_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#lmx <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/LMX_BASALE_mean_gene_genealogyall.tsv.gz"
lmx <- read.table(lmx_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#pdo <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
pdo <- read.table(pdo_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pdo <- pdo %>% filter(buoni == "TRUE")

lmo <- lmo[rownames(lmo) %in% pdo$smodel,]
lmx <- lmx[rownames(lmx) %in% pdo$smodel,]

write.table(lmo, file = lmo_r, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(lmx, file = lmx_r, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)