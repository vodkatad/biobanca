pdo <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
pdo <- read.table(pdo, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#pdo <- pdo %>% filter(buoni == TRUE)
colnames(pdo) <- c("CASE", "buoni")

pdx <- "/scratch/trcanmed/biobanca/local/share/data/List_matched_PDO-PDX_DNA_Revised_SimoMay2021.tsv"
pdx <- read.table(pdx, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

#setdiff(pdo$CASE, merged$CASE)
#[1] "CRC1634" "CRC1779" "CRC1782" "CRC1921"

##togliere gli lmh anche dagli arricchimenti?

merged <- merge(pdo, pdx, by = "CASE")
merged <- data.frame(merged$CASE, merged$buoni, merged$PDO.lineage, merged$MATCHED.LMX.LMO)
merged$passage <- substr(merged$merged.MATCHED.LMX.LMO, 13, 14)
merged$class <- substr(merged$merged.MATCHED.LMX.LMO, 8, 10)
merged <- merged %>% filter(class == "LMX")
merged$class <- NULL

table(merged$merged.buoni, merged$passage)
