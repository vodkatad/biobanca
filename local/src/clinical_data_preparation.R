## control msi mss e preparazione clinical data

clinical <- read_excel("/scratch/trcanmed/biobanca/local/share/data/XENTURION_CLINICAL DATA_Eugenia_06092022_DEFdef.xlsx")
clinical$NAME <- NULL
msi <- "/scratch/trcanmed/pdxopedia/local/share/data/eugy_msi_long.tsv"
buoni <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"

msi <- read.table(msi, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(msi) <- c("GENEALOGY", "MSI_MSS")
msi$CASE <- substr(msi$GENEALOGY, 0, 7)
msi$type <- substr(msi$GENEALOGY, 8, 10)

buoni <- read.table(buoni, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

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

## elimina crc1374 e metti NA ai due non crc
merged <- merged %>% filter(!CASE == "CRC1874")
merged["CRC1241", -c(1)] <- NA
merged["CRC1575", -c(1)] <- NA
merged["CRC0578", -c(1)] <- NA

#y minuscola
merged["CRC1875", "THERAPY BEFORE (Y/N)"] <- "Y"

write.table(merged, file = "/scratch/trcanmed/biobanca/local/share/data/clinical_data_done.tsv", quote = FALSE, sep = "\t", 
            col.names = TRUE, row.names = FALSE)
