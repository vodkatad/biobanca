#### dati puntuali di cambiamento di espressione dopo cetuximab 
#### dei 13 target selezionati nei 5 PDXT usati per lo screening

library(tidyverse)
library(openxlsx)

target_f <- snakemake@input[["tar"]]
tmm_f <- snakemake@input[["expr"]]
status_f <- snakemake@input[["treat"]]
result <- snakemake@output[["tsv"]]

#target <- read.xlsx("/scratch/trcanmed/biobanca/local/share/data/targets_drug_screening.xlsx",colNames = FALSE)
target <- read.xlsx(target_f, colNames = FALSE)
colnames(target) <- c("drug", "target")
target <- target$target
target <- trimws(target)
target <- c(target, "FGFR1", "FGFR2", "HDAC5", "HDAC11", "HDAC6", "HDAC10", "STAT3")

#tmm <- read.table("/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/tmm.tsv.gz",
#                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tmm <- read.table(tmm_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- rownames(tmm)
rownames(tmm) <- NULL
tmm <- cbind(genes,tmm)
tmm <- tmm %>% mutate(genes = gsub("H_", "", genes))
colnames(tmm)[colnames(tmm) == 'genes'] <- 'symbol'
tmm <- tmm %>% filter(symbol %in% target)
rownames(tmm) <- tmm$symbol
tmm$symbol <- NULL

tmm <- as.data.frame(t(tmm))

# whoiswho <- "/scratch/trcanmed/biobanca/local/share/data/whoiswho_validation_xen_nolmh.tsv"
# ww <- read.table(whoiswho, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# ww <- ww %>% filter(type == "Validation successful")
# 
# casi <- substr(rownames(tmm), 1, 7)
# casi <- unique(casi)
# 
# length(unique(casi))
# length(intersect(casi, ww$CASE))
print("Check effettuato sono tutti validati")
tmm$casi <- substr(rownames(tmm), 1, 7)

casi_drug <- c("CRC0059", "CRC0322", "CRC0327", "CRC1331", "CRC0069", "CRC0076", "CRC0399",
               "CRC0534", "CRC0542", "CRC0743", "CRC1272", "CRC0116")
tmm <- tmm %>% filter(casi %in% casi_drug)
tmm$casi <- NULL

#status_f <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/samples_data"
status <- read.table(status_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
status <- status %>% filter(sample %in% casi_drug)

tmm$id <- rownames(tmm)

merged <- merge(tmm, status, by = "id")
rownames(merged) <- merged$id
merged$id <- NULL
merged$sample <- NULL
merged$batch <- NULL
merged$treat <- gsub("i", "", merged$treat)
merged$Genealogy_ID <- rownames(merged) 
names(merged)[names(merged) == "treat"] <- "Treatment"
merged <- merged[,c(19,18,1:17)]


write.table(merged, result,
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
