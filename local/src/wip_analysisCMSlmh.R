### CMSclasses comparison

library(tidyverse)
library(data.table)

cmsall2 <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/CMScaller.tsv"
cmsall2 <- read.table(cmsall2, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmsall2 <- rownames_to_column(cmsall2)
names(cmsall2)[1] <- "Genealogy"
cmsalls2 <- cmsall2 %>% select(Genealogy, prediction)
names(cmsalls2)[2] <- "prediction_all"

cmslmh <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMH_CMScaller.tsv"
cmslmh <- read.table(cmslmh, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmslmh <- rownames_to_column(cmslmh)
names(cmslmh)[1] <- "Genealogy"
cmslmhs <- cmslmh %>% select(Genealogy, prediction)
names(cmslmhs)[2] <- "prediction_lmh"

merged2 <- merge(cmsalls2, cmslmhs, by = "Genealogy")
merged2 <- merged2 %>% add_column(new_col = NA)
names(merged2)[4] <- "results"
merged2[is.na(merged2$prediction_all), 'prediction_all'] <- "NC"
merged2[is.na(merged2$prediction_lmh), 'prediction_lmh'] <- "NC"
merged2$results <- ifelse(merged2$prediction_all == merged2$prediction_lmh, 'same', 'different')

tablemerged2 <- setDT(merged2)[, 100 * .N / nrow(merged2), by = results]

different2 <- subset(merged2, merged2$results == "different")

tdifferentall2 <- setDT(different2)[, 100 * .N / nrow(different2), by = prediction_all]
tdifferentlmh <- setDT(different2)[, 100 * .N / nrow(different2), by = prediction_lmh]
