### CMS for classes 

library(tidyverse)
library(CMScaller)

meta_f <- snakemake@input[['metadata']]
vsd <- snakemake@input[["expr"]]
vsd_genes <- snakemake@input[["VSD_genes"]]
classes <- snakemake@wildcards[['sclass']]
plot <- snakemake@output[["CMS_heatmap"]]
results <- snakemake@output[["RES"]]

meda <- meta_f
meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_f <- filter(meda_f, grepl(classes, type))
meda_f$sample_id_new <- gsub('-','.', meda_f$sample_id_R, fixed=TRUE)

vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE)
genes <- rownames(vsd)
rownames(vsd) <- NULL
vsd <- cbind(genes,vsd)
vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
colnames(vsd)[colnames(vsd) == 'genes'] <- 'symbol'
vsd_subset <- vsd[, colnames(vsd) %in% c(meda_f$sample_id_new, "symbol")]

vsd_genes <- read.table(vsd_genes, quote = "", sep = "\t", header = TRUE)
vsd_genes <- vsd_genes[!grepl(",", vsd_genes$description),,drop = FALSE]
vsd_genes <- vsd_genes %>% na.exclude

vsd_f <- merge(vsd_genes, vsd_subset, by = "symbol")
vsd_f$symbol <- NULL
vsd_f <- vsd_f %>% remove_rownames %>% column_to_rownames(var="description")

pdf(plot)
res <- CMScaller(emat=vsd_f, FDR=0.05, RNAseq=TRUE)
graphics.off()

write.table(res, file = results, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)



