library(tidyverse)

cms_lmx <- snakemake@input[["lmx"]]
cms_lmh <- snakemake@input[["lmh"]]
#pdobuoni <- snakemake@input[["true"]]
results <- snakemake@output[["merged"]]

#cmslmx <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/LMX_BASALE_CMScaller_mediated.tsv"
cmslmx <- read.table(cms_lmx, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmslmx$model <- rownames(cmslmx)
rownames(cmslmx) <- NULL
cmslmx <- select(cmslmx, model, prediction)
colnames(cmslmx) <- c("model", "CMS_PDX")

#cmslmh <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/LMH_CMScaller.tsv"
cmslmh <- read.table(cms_lmh, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmslmh$model <- substr(rownames(cmslmh), 1, 7)
rownames(cmslmh) <- NULL
cmslmh <- select(cmslmh, model, prediction)
colnames(cmslmh) <- c("model", "CMS_LMH")

merged <- merge(cmslmx, cmslmh, by = "model")
merged <- merged %>% remove_rownames %>% column_to_rownames(var="model")

#pdo_buoni <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
# pdo_buoni <- read.table(pdobuoni, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# pdo_buoni <- pdo_buoni %>% filter(buoni == "TRUE")
# colnames(pdo_buoni) <- c("model", "buoni")

# merged_buoni <- merge(merged, pdo_buoni, by = "model")
# merged_buoni$buoni <- NULL
# merged_buoni <- merged_buoni %>% remove_rownames %>% column_to_rownames(var="model")

write.table(merged, file = results, quote = FALSE, sep = "\t", col.names = TRUE)