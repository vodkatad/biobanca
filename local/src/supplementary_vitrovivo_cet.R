## supplementary vitro vivo cetuxi 

library(tidyverse)

ctgi_f <- snakemake@input[["pdo_cetuxi"]]
vivo_f <- snakemake@input[["cetuxi_w3"]]
res <- snakemake@output[["sup"]]

#ctgi_f <- "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/pdo_cetuxi_buoni.tsv"
ctgi <- read.table(ctgi_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#vivo_f <- "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/cetuxi_perc_w3_buoni.tsv"
vivo <- read.table(vivo_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

merged <- merge(ctgi, vivo, all.x = TRUE)
names(merged)[names(merged) == "perc"] <- "perc_cetux_3_weeks"

write.table(merged, file = "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/supplementary_vitro_vivo_cetuxi.tsv", quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = FALSE)


