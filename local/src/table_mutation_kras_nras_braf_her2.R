## table mutation KRAS_NRAS_BRAF and ERBB2

library(tidyverse)
library(WriteXLS)


d <- read.table("/scratch/trcanmed/biobanca/dataset/V1/cetuximab/atp/targeted_annot.tsv",
                quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
d <- read.table(df, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(d) <- d$Row.names
d$Row.names <- NULL

write.xlsx(d, file = "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/atp/targeted_annot.xlsx")
d$CTG_5000 <- NULL
d$Cetuximab_dVw3 <- NULL
d$triple_wt <- NULL
d$alterations <- NULL
write.xlsx(d, file = "/scratch/trcanmed/biobanca/dataset/V1/cetuximab/atp/targeted_annot_clean.xlsx")
