library(tidyverse)

pdx <- snakemake@input[["pdx_f"]]
pdo <- snakemake@input[["pdo_f"]]
output <- snakemake@output[[1]]

#pdx <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv"
pdx <- read.table(pdx, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(pdx) <- c("genealogy", "CMS_PDX")

#pdo <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMO_BASALE_nc_smodel.tsv"
pdo <- read.table(pdo, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(pdo) <- c("genealogy", "CMS_PDO")

merged <- merge(pdx, pdo, by = "genealogy")
merged <- merged %>% remove_rownames %>% column_to_rownames(var="genealogy")

write.table(merged, output, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

