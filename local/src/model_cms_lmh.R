library(tidyverse)

modelcms <- snakemake@input[["tsv"]]
cms_lmh <- snakemake@input[["lmh"]]
results <- snakemake@output[["res"]]

#model_cms <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/model_cms.tsv"
model_cms <- read.table(modelcms, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#cmslmh <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/LMH_CMScaller.tsv"
cmslmh <- read.table(cms_lmh, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmslmh$model <- substr(rownames(cmslmh), 1, 7)
cmslmh <- select(cmslmh, prediction, model)
colnames(cmslmh) <- c("CMS_LMH", "model")
mergedlmh <- merge(model_cms, cmslmh, by = "model")
mergedlmh$CMS_PDO <- NULL
mergedlmh <- mergedlmh %>% remove_rownames %>% column_to_rownames(var="model")

write.table(mergedlmh, file = results, quote = FALSE, sep = "\t", col.names = TRUE)