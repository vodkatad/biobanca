library (tidyverse)
library(openxlsx)

max_f <- snakemake@input[['max_inh']]
anova_f <- snakemake@input[['anova']]
res <- snakemake@output[["results"]]
res2 <- snakemake@output[["results2"]]

#anova_f <- "/scratch/trcanmed/biobanca/dataset/V1/drug_screening/anova.tsv"
a <- read.table(anova_f, quote = "", sep = "\t", header = T, stringsAsFactors = F)
a$model_cond <- paste0(a$MODEL, ",", a$DRUG, ",", a$CONDITION)

#max_f <- "/scratch/trcanmed/biobanca/dataset/V1/drug_screening/max_drug_value.tsv"
max <- read.table(max_f, quote = "", sep = "\t", header = T, stringsAsFactors = F)
max$model_cond <- paste0(max$MODEL, ",", max$DRUG, ",", max$CONDITION)

merged <- merge(a, max, by = "model_cond")
merged$model_cond <- NULL
merged$DRUG.y <- NULL
merged$CONDITION.y <- NULL
merged$MODEL.y <- NULL
names(merged)[names(merged) == "MODEL.x"] <- "MODEL"
names(merged)[names(merged) == "DRUG.x"] <- "DRUG"
names(merged)[names(merged) == "CONDITION.x"] <- "CONDITION"

write.table(merged, file = res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.xlsx(merged, file = res2)

