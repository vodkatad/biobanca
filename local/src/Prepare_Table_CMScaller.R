### prepare VSD table for CMS

library(tidyverse)
vsd <- snakemake@input[["expr"]]
vsd_f <- snakemake@output[["VSD_f"]]

vsd_file <- read.table(vsd, quote = "", sep = "\t", header = TRUE)
genes <- rownames(vsd_file)
symbols <- gsub("H_", "", genes)
#rownames(vsd_file) <- str_remove(rownames(vsd_file), 'H_')

res <- data.frame(symbol = symbols)
write.table(res, file = vsd_f, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

