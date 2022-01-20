### prepare VSD table for CMS

library(tidyverse)
fpkm <- snakemake@input[["expr"]]
fpkm_f <- snakemake@output[["FPKM_f"]]

fpkm_file <- read.table(fpkm, quote = "", sep = "\t", header = TRUE)
symbols <- rownames(fpkm_file)
#symbols <- gsub("H_", "", genes)
#rownames(vsd_file) <- str_remove(rownames(vsd_file), 'H_')

res <- data.frame(symbol = symbols)
write.table(res, file = fpkm_f, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

