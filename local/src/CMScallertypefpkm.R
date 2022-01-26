### CMS for classes 

library(tidyverse)
library(CMScaller)

meta_f <- snakemake@input[['metadata']]
fpkm <- snakemake@input[["expr"]]
fpkm_genes <- snakemake@input[["FPKM_genes"]]
classes <- snakemake@wildcards[['sclass']]
plot <- snakemake@output[["CMS_heatmap"]]
results <- snakemake@output[["RES"]]

#meda <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda <- meta_f
meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_f <- filter(meda_f, grepl(classes, type))
meda_f$sample_id_new <- gsub('-','.', meda_f$sample_id_R, fixed=TRUE)

#fpkm <- "/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/fpkm_H.tsv.gz"
fpkm <- read.table(fpkm, quote = "", sep = "\t", header = TRUE)
genes <- rownames(fpkm)
rownames(fpkm) <- NULL
fpkm <- cbind(genes,fpkm)
colnames(fpkm)[colnames(fpkm) == 'genes'] <- 'symbol'
fpkm_subset <- fpkm[, colnames(fpkm) %in% c(meda_f$sample_id_new, "symbol")]

#fpkm_genes <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/fpkm_genes.tsv"
fpkm_genes <- read.table(fpkm_genes, quote = "", sep = "\t", header = TRUE)
fpkm_genes <- fpkm_genes[!grepl(",", fpkm_genes$description),,drop = FALSE]
fpkm_genes <- fpkm_genes %>% na.exclude

fpkm_f <- merge(fpkm_genes, fpkm_subset, by = "symbol")
fpkm_f$symbol <- NULL
fpkm_f <- fpkm_f %>% remove_rownames %>% column_to_rownames(var="description")

png(plot)
res <- CMScaller(emat=fpkm_f, FDR=0.05, RNAseq = TRUE)
graphics.off()

write.table(res, file = results, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)