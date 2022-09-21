library(tidyverse)

fra <- snakemake@input[["mut"]]
res <- snakemake@output[["results"]]

#fra <- "/scratch/trcanmed/pdxopedia/local/share/data/dtb_mutations_vlookupAndrea_expandedPRLM.tsv"
fra <- read.table(fra, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

fra2 <- fra[!(fra$X=="PR"),]
fra2$X <- NULL          
fra2$GenealogyID <- NULL
fra2$flag_H.X <- NULL

## removing these because hcc and not CRC
hcc <- c("CRC1379", "CRC1479", "CRC1737")
fra2 <- fra2 %>% filter(!CASE %in% hcc)

fra3 <- fra2
cases <- fra3$CASE              

for (i in seq(length(cases))) {
  fra3$KRAS2 <- ifelse(fra3$KRAS == "wt", "False", "True")
  fra3$NRAS2 <- ifelse(fra3$NRAS == "wt", "False", "True")
  fra3$BRAF2 <- ifelse(fra3$BRAF == "wt", "False", "True")
  fra3$PIK3CA2 <- ifelse(fra3$PIK3CA == "wt", "False", "True")
}

fra3$KRAS <- NULL
fra3$NRAS <- NULL
fra3$BRAF <- NULL
fra3$PIK3CA <- NULL
colnames(fra3) <- c('smodel', 'KRAS', 'NRAS', 'BRAF', 'PIK3CA')
rownames(fra3) <- fra3$smodel
fra3$smodel <- NULL
tfra <- as.data.frame(t(fra3))
tfra <- tibble::rownames_to_column(tfra, "genes")

write.table(tfra, res, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)