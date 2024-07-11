library(tidyverse)

mixed_f <- snakemake@input[["tsv"]]
res <- snakemake@output[["cnv"]]
#res_check<-snakemake@output[['nosense_loh']]
#res_common<-snakemake@output[['common_loh']]

#mixed_f <- "/scratch/trcanmed/biobanca/dataset/V1/LOH_earlylate/mixed_bed_cnv/CRC0095_cnv.bed"
mixed <- read.table(mixed_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(mixed) <- c("chr_e", "start_e", "end_e", "CNt_e", "chr_l", "start_l", "end_l", "CNt_l", "intersection")
#mi salvo in un df diverso gli eventi strani==> B_e=0 and B_l!=0
#mixed_check<-read.table(mixed_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#colnames(mixed_check) <- c("chr_e", "start_e", "end_e", "A_e","B_e", "chr_l", "start_l", "end_l", "A_l","B_l", "intersection")
#mi salvo in un df diverso gli eventi comuni==> B_e=0 and B_l=0
#common<-read.table(mixed_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#colnames(common) <- c("chr_e", "start_e", "end_e", "A_e","B_e", "chr_l", "start_l", "end_l", "A_l","B_l", "intersection")

mixed$cnv <- ifelse(mixed$CNt_e < 5 & mixed$CNt_l >= 5, "GAINED", ifelse(mixed$CNt_e >=2 & mixed$CNt_l <=1, "LOSS", "NONE"))

#mixed$loh <- ifelse((mixed$B_e != 0 & mixed$B_l == 0), "LOH", "NONE")
#mixed_check$nosense_loh <- ifelse((mixed$B_e == 0 & mixed$B_l != 0), "EARLY", "NONE")
#common$common_loh <- ifelse((common$B_e == 0 & common$B_l == 0), "COMMON", "NONE")
filtro <- c("GAINED", "LOSS")

gained <- mixed %>% filter(cnv %in% filtro)

gained$start <- pmax(gained$start_e, gained$start_l)
gained$end <- pmin(gained$end_e, gained$end_l) ## corretto
#loh$start <- pmax(loh$start_e, loh$start_l)
#loh$end <- pmin(loh$end_e, loh$end_l)

# nosense_loh$start <- pmax(nosense_loh$start_e, nosense_loh$start_l)
# nosense_loh$end <- pmin(nosense_loh$end_e, nosense_loh$end_l)
# 
# common$start <- pmax(common$start_e, common$start_l)
# common$end <- pmin(common$end_e, common$end_l)

gained <- gained[,c(1, 10, 11, 12)]
gained <- gained[,c(1, 3, 4, 2)]
names(gained)[names(gained)== "chr_e"] <- "chr"
#loh <- loh %>% filter(intersection > 100000)

# nosense_loh <- nosense_loh[,c(1, 11, 13, 14)]
# nosense_loh <- nosense_loh[,c(1, 3, 4, 2)]
# names(nosense_loh)[names(nosense_loh)== "chr_e"] <- "chr"
# #nosense_loh <- nosense_loh %>% filter(intersection > 100000)
# 
# common <- common[,c(1, 11, 13, 14)]
# common <- common[,c(1, 3, 4, 2)]
# names(common)[names(common)== "chr_e"] <- "chr"
#common <- common %>% filter(intersection > 100000)

write.table(gained, file=res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
#write.table(nosense_loh, file=res_check, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
#write.table(common, file=res_common, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

