### CMSclasses comparison

library(tidyverse)
library(data.table)

cmsall <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/CMScaller.tsv"
cmsall <- read.table(cmsall, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmsall <- rownames_to_column(cmsall)
names(cmsall)[1] <- "Genealogy"
cmsalls <- cmsall %>% select(Genealogy, prediction)
names(cmsalls)[2] <- "prediction_all"

cmslmo <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMO_BASALE_CMScaller.tsv"
cmslmo <- read.table(cmslmo, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmslmo <- rownames_to_column(cmslmo)
names(cmslmo)[1] <- "Genealogy"
cmslmos <- cmslmo %>% select(Genealogy, prediction)
names(cmslmos)[2] <- "prediction_lmo"

merged <- merge(cmsalls, cmslmos, by = "Genealogy")
merged <- merged %>% add_column(new_col = NA)
names(merged)[4] <- "results"
merged[is.na(merged$prediction_all), 'prediction_all'] <- "NC"
merged[is.na(merged$prediction_lmo), 'prediction_lmo'] <- "NC"
merged$results <- ifelse(merged$prediction_all == merged$prediction_lmo, 'same', 'different')

tablemerged <- setDT(merged)[, 100 * .N / nrow(merged), by = results]

different <- subset(merged, merged$results == "different")
### different classes CRC0177LMO0D04021002R01000, CRC1344LMO0A04007001R01000, CRC1390LMO0B03018001R01000
tdifferentall <- setDT(different)[, 100 * .N / nrow(different), by = prediction_all]
tdifferentlmo <- setDT(different)[, 100 * .N / nrow(different), by = prediction_lmo]
