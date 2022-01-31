library(tidyverse)

cms_f <- snakemake@input[["cms"]]
results <- snakemake@output[["reorder"]]

cms_model <- read.table(cms_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
res <- cms_model[, c("model", "CMS_PDX", "CMS_PDO")]

write.table(res, file = results, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)