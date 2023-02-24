library(tidyverse)

cris_f <- snakemake@input[["cris"]]
res <- snakemake@output[["results"]]

#cris_f <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/vsd_model_cris-right_validated_parallelo.tsv"
cris <- read.table(cris_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cris <- cris[,c("model", "PDO", "PDX")]
cris <- cris[,c("model", "PDX", "PDO")]
colnames(cris) <- c("model", "CRIS_PDX", "CRIS_PDO")
write.table(cris, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)