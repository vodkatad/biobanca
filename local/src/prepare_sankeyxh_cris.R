library(tidyverse)

cris_lmx <- snakemake@input[["lmx"]]
cris_lmh <- snakemake@input[["lmh"]]
# pdobuoni <- snakemake@input[["true"]]
results <- snakemake@output[["merged"]]

#crislmx <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd_cris_LMX_BASALE_prediction_result_nc.tsv"
crislmx <- read.table(cris_lmx, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
crislmx$model <- rownames(crislmx)
rownames(crislmx) <- NULL
crislmx <- select(crislmx, sample.names, predict.label2)
colnames(crislmx) <- c("model", "CRIS_PDX")

#crislmh <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd_cris_LMH_prediction_result_nc.tsv"
crislmh <- read.table(cris_lmh, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
crislmh$model <- substr(crislmh$sample.names, 1, 7)
crislmh <- select(crislmh, model, predict.label2)
colnames(crislmh) <- c("model", "CRIS_LMH")

merged <- merge(crislmx, crislmh, by = "model")
merged <- merged %>% remove_rownames %>% column_to_rownames(var="model")

#pdo_buoni <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
# pdo_buoni <- read.table(pdobuoni, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# pdo_buoni <- pdo_buoni %>% filter(buoni == "TRUE")
# colnames(pdo_buoni) <- c("model", "buoni")

# merged_buoni <- merge(merged, pdo_buoni, by = "model")
# merged_buoni$buoni <- NULL
# merged_buoni <- merged_buoni %>% remove_rownames %>% column_to_rownames(var="model")

write.table(merged, file = results, quote = FALSE, sep = "\t", col.names = TRUE)