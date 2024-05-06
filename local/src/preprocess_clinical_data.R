## control msi mss e preparazione clinical data

library(tidyverse)
library(readxl)

c_f <- snakemake@input[["data"]]
msi_f <- snakemake@input[["m"]]
done <- snakemake@output[["res"]]

#clinical <- read_excel("/scratch/trcanmed/biobanca/local/share/data/Extended_Data_Table_1_NEW_DEF.xlsx")
#msi <- "/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/eugy_msi_long.tsv"

clinical <- read_excel(c_f)
msi <- read.table(msi_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(msi) <- c("GENEALOGY", "MSI_MSS")
msi$CASE <- substr(msi$GENEALOGY, 0, 7)
msi$type <- substr(msi$GENEALOGY, 8, 10)


merged <- merge(clinical, msi, by = "CASE", all.x = TRUE)
casi_dupli <- merged[duplicated(merged$CASE),]
casi_dupligen <- casi_dupli$CASE
mergeddupli <- merged %>% filter(CASE %in% casi_dupligen)
duplicati <- as.data.frame(cbind(mergeddupli$CASE, mergeddupli$GENEALOGY, mergeddupli$MSI_MSS))

merged <- merged[!duplicated(merged$CASE),]
rownames(merged) <- merged$CASE
#inserisco crc1331 mss a mano poiché eliminando i duplicati sarebbe rimasto nt
#invece lmx è stato testato e ci teniamo quello
merged["CRC1331", "MSI_MSS"] <- "MSS"
merged$type <- NULL
merged$GENEALOGY <- NULL

merged["CRC1241", -c(1)] <- NA
merged["CRC1575", -c(1)] <- NA
merged["CRC0578", -c(1)] <- NA

write.table(merged, file = done, quote = FALSE, sep = "\t", 
            col.names = TRUE, row.names = FALSE)
