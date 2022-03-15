library(tidyverse)

pdo <- snakemake@input[["buoni"]]
cet <- snakemake@input[["pdo_cet"]]
result <- snakemake@output[["res"]]

#pdo <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
pdo <- read.table(pdo, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pdo <- pdo %>% filter(buoni == TRUE)

#cet <- "/scratch/trcanmed/biobanca/dataset/V1/cetuximab_tgi/pdo_cetuxi.tsv"
cet <- read.table(cet, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cet <- cet[cet$case %in% pdo$smodel,]

write.table(cet, result, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)