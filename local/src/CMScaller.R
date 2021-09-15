### CMScaller

library(tidyverse)
library(CMScaller)

vsd <- snakemake@input[["expr"]]
vsd_genes <- snakemake@input[["VSD_genes"]]
plot <- snakemake@output[["CMS_heatmap"]]
results <- snakemake@output[["RES"]]

vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE)
genes <- rownames(vsd)
rownames(vsd) <- NULL
vsd <- cbind(genes,vsd)
vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
colnames(vsd)[colnames(vsd) == 'genes'] <- 'symbol'

vsd_genes <- read.table(vsd_genes, quote = "", sep = "\t", header = TRUE)
vsd_genes <- vsd_genes[!grepl(",", vsd_genes$description),,drop = FALSE]
vsd_genes <- vsd_genes %>% na.exclude

vsd_f <- merge(vsd_genes, vsd, by = "symbol")
vsd_f$symbol <- NULL
vsd_f <- vsd_f %>% remove_rownames %>% column_to_rownames(var="description")

pdf(plot)
res <- CMScaller(emat=vsd_f, RNAseq=TRUE, FDR=0.05)
dev.off()

write.table(res, file = results, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)
