library(CRISclassifier, quietly=TRUE)
library(tidyverse)

exprfile <- snakemake@input[["expr"]]
outPrefix <- snakemake@params[["prefix"]]

#exprfile <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/LMO_BASALE_mean_gene_genealogyall.tsv.gz"
ex <- read.table(exprfile, sep="\t", header=TRUE)

save.image("pippo.Rdata")
SYMBOL <- rownames(ex)
rownames(ex) <- NULL
ex <- cbind(SYMBOL,ex)
ex <- ex %>% mutate(SYMBOL = gsub("H_", "", SYMBOL))
# SYMBOL <- ex$SYMBOL
# ex <- ex %>% remove_rownames %>% column_to_rownames(var="SYMBOL")
res <- ex


#dims <- dim(ex)
#res <- ex[,c(dims[2],seq(1,(dims[2]-1)))]

write.table(res, file=paste0(outPrefix, ".tmp"), sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)
cris_classifier(paste0(outPrefix, ".tmp"), outPrefix)
