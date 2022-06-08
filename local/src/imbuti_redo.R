try <- read.table(gzfile('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/big_merge.tsv.gz'), sep="\t", comment.char='', quote='')

mediane <- median(try$baseMean)
try$filter_PDX_S <- ifelse(!is.na(try$padj_PDX_S) & try$padj_PDX_S < 0.05 & try$log2FC_PDX_S > 0, TRUE, FALSE)
try$filter_PDX_R <- ifelse(!is.na(try$padj_PDX_R) & try$padj_PDX_R >= 0.05, TRUE, FALSE) 
try$filter_expr <- ifelse(try$baseMean > mediane, TRUE, FALSE)
try$allfilter <- apply(try[,grepl('filter_', colnames(try))], 1, all)
table(try$allfilter)

table(res$allfilter)

