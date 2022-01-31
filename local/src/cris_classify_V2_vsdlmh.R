library(CRISclassifier, quietly=TRUE)
library(tidyverse)


exprfile <- snakemake@input[["expr"]]
outPrefix <- snakemake@params[["prefix"]]

exprfile <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMH_gene_genealogyall.tsv"
ex <- read.table(exprfile, sep="\t", header=TRUE)
# ex <- as.data.frame(t(ex))
# 
# 
# SYMBOL <- rownames(ex)
# rownames(ex) <- NULL
# ex <- cbind(SYMBOL,ex)
# ex <- ex %>% mutate(SYMBOL = gsub("H_", "", SYMBOL))
# SYMBOL <- ex$SYMBOL
#ex <- ex %>% remove_rownames %>% column_to_rownames(var="SYMBOL")
res <- ex


#dims <- dim(ex)
#res <- ex[,c(dims[2],seq(1,(dims[2]-1)))]

write.table(res, file=paste0(outPrefix, ".tmp"), sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
cris_classifier(paste0(outPrefix, ".tmp"), outPrefix)
