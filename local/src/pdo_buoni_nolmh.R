library(tidyverse)

pdo <- snakemake@input[["inp"]]
res <- snakemake@output[["out"]]

#pdo <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
pdo <- read.table(pdo, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#pdo <- pdo %>% filter(buoni == "TRUE")

pdofromlmh <- c("CRC1337", "CRC1342", "CRC1563", "CRC1566")
pdo <- pdo[!pdo$smodel %in% pdofromlmh,]
buoni <- pdo

write.table(buoni, file = res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)