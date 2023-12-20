library(tidyverse)

mixed_f <- snakemake@input[["tsv"]]
res <- snakemake@output[["loh"]]
res_check<-snakemake@output[['nosense_loh']]
res_common<-snakemake@output[['common_loh']]

#mixed_f <- "/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_subs/loh/mixed_bed/CRC1599PRX0A02002TUMD03000V2_CRC1599LMX0A02001TUMD03000V2.bed"
mixed <- read.table(mixed_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(mixed) <- c("chr_e", "start_e", "end_e", "A_e","B_e", "chr_l", "start_l", "end_l", "A_l","B_l", "intersection")
#mi salvo in un df diverso gli eventi strani==> B_e=0 and B_l!=0
mixed_check<-read.table(mixed_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(mixed_check) <- c("chr_e", "start_e", "end_e", "A_e","B_e", "chr_l", "start_l", "end_l", "A_l","B_l", "intersection")
#mi salvo in un df diverso gli eventi comuni==> B_e=0 and B_l=0
common<-read.table(mixed_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(common) <- c("chr_e", "start_e", "end_e", "A_e","B_e", "chr_l", "start_l", "end_l", "A_l","B_l", "intersection")


mixed$loh <- ifelse((mixed$B_e != 0 & mixed$B_l == 0), "LOH", "NONE")
mixed_check$nosense_loh <- ifelse((mixed$B_e == 0 & mixed$B_l != 0), "EARLY", "NONE")
common$common_loh <- ifelse((common$B_e == 0 & common$B_l == 0), "COMMON", "NONE")

loh <- mixed %>% filter(loh == "LOH")
nosense_loh<- mixed_check %>% filter(nosense_loh == "EARLY")
common<- common %>% filter(common_loh == "COMMON")

loh$start <- pmin(loh$start_e, loh$start_l)
loh$end <- pmax(loh$end_e, loh$end_l)

nosense_loh$start <- pmin(nosense_loh$start_e, nosense_loh$start_l)
nosense_loh$end <- pmax(nosense_loh$end_e, nosense_loh$end_l)

common$start <- pmin(common$start_e, common$start_l)
common$end <- pmax(common$end_e, common$end_l)

loh <- loh[,c(1, 11, 13, 14)]
loh <- loh[,c(1, 3, 4, 2)]
names(loh)[names(loh)== "chr_e"] <- "chr"
loh <- loh %>% filter(intersection > 100000)

nosense_loh <- nosense_loh[,c(1, 11, 13, 14)]
nosense_loh <- nosense_loh[,c(1, 3, 4, 2)]
names(nosense_loh)[names(nosense_loh)== "chr_e"] <- "chr"
nosense_loh <- nosense_loh %>% filter(intersection > 100000)

common <- common[,c(1, 11, 13, 14)]
common <- common[,c(1, 3, 4, 2)]
names(common)[names(common)== "chr_e"] <- "chr"
common <- common %>% filter(intersection > 100000)

write.table(loh, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(nosense_loh, file=res_check, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(common, file=res_common, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

