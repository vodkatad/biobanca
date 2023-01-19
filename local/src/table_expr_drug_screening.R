#### dati puntuali di cambiamento di espressione dopo cetuximab 
#### dei 13 target selezionati nei 5 PDXT usati per lo screening

target <- read.xlsx("/scratch/trcanmed/biobanca/local/share/data/targets_drug_screening.xlsx",colNames = FALSE)
colnames(target) <- c("drug", "target")
target <- target$target
target <- trimws(target)
target <- c(target, "FGFR1", "FGFR2", "HDAC5", "HDAC11", "HDAC6", "HDAC10", "STAT3")

tmm <- read.table("/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/tmm.tsv.gz",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- rownames(tmm)
rownames(tmm) <- NULL
tmm <- cbind(genes,tmm)
tmm <- tmm %>% mutate(genes = gsub("H_", "", genes))
colnames(tmm)[colnames(tmm) == 'genes'] <- 'symbol'
tmm <- tmm %>% filter(symbol %in% target)
rownames(tmm) <- tmm$symbol
tmm$symbol <- NULL

tmm <- as.data.frame(t(tmm))
tmm$casi <- substr(rownames(tmm), 1, 7)

casi_drug <- c("CRC0059", "CRC0322", "CRC0327", "CRC1272", "CRC1331")
tmm <- tmm %>% filter(casi %in% casi_drug)
tmm$casi <- NULL

status <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/samples_data"
status <- read.table(status, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
status <- status %>% filter(sample %in% casi_drug)

tmm$id <- rownames(tmm)

merged <- merge(tmm, status, by = "id")
rownames(merged) <- merged$id
merged$id <- NULL
merged$sample <- NULL
merged$batch <- NULL

write.table(merged, "/scratch/trcanmed/biobanca/dataset/V1/drug_screening/table_expr_drug_screening.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
