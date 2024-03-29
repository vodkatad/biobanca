library(ggplot2)
library(pheatmap)
#setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S_no116')
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S')

deg <- read.table('treat_cutoff0.05-cetuxi.vs.NT.deseq2_sign_genedesc.tsv', comment.char="", quote='', header=T, sep="\t")
rownames(deg) <- deg$gene
deg$gene <- NULL

## heatmap of DEG ########
expr <- read.table(gzfile('fpkm.tsv.gz'), sep="\t", header=T)
rownames(expr) <- substr(rownames(expr), 3, nchar(rownames(expr)))
samples_data <- read.table('samples_data', sep="\t", header=T)
rownames(samples_data) <- samples_data$id
samples_data$id <- NULL

THR <- 0.5849625
THR <- 0
PC <- 1
deg_lf <- deg[deg$log2FoldChange > THR,]

deg_expr <- expr[rownames(expr) %in% rownames(deg_lf),]


# samples_data <- samples_data[order(samples_data$sample, samples_data$treat),]
# deg_expr <- deg_expr[,match(rownames(samples_data), colnames(deg_expr))]
# pheatmap(log(deg_expr+PC), show_rownames=TRUE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F)
# annot_row <- data.frame(row.names=rownames(deg), expr=log(deg$baseMean), log2FoldChange=deg$log2FoldChange)
# 
# annot_row <- annot_row[rownames(annot_row) %in% rownames(deg_expr),]
# pheatmap(log(deg_expr+PC), show_rownames=FALSE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F, annotation_row=annot_row)
# 
# annot_row <- annot_row[annot_row$expr>4.791111,]
# deg_expr <- deg_expr[rownames(deg_expr) %in% rownames(annot_row),]
# pheatmap(log(deg_expr+PC), show_rownames=FALSE, show_colnames=FALSE, annotation_col = samples_data, cluster_cols = F, annotation_row=annot_row)

#### merge
deg_up <- deg[deg$log2FoldChange > 0,]
deg_up$description <- as.character(deg_up$description)
pdx_r <- read.table('../Biodiversa_up5starOK_cetuxi_treat_PDX_R_no1241/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv', comment.char="", quote='', header=T, sep="\t") # saved as new_res
####pdx_r <- read.table('../Biodiversa_up5starOK_cetuxi_treat_PDX_R/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv', comment.char="", quote='', header=T, sep="\t")

pdx_s <- read.table('../Biodiversa_up5starOK_cetuxi_treat_PDX_S/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv', comment.char="", quote='', header=T, sep="\t")



prepare_df <- function(df, keep=c('log2FoldChange','padj')) {
  rownames(df) <- substr(rownames(df), 3, nchar(rownames(df)))
  df[, keep]
}



pdx_s <- prepare_df(pdx_s)
colnames(pdx_s) <- c('log2FC_PDX_S', 'padj_PDX_S')
res <- merge(deg_up, pdx_s, by="row.names", all.x=TRUE)
rownames(res) <- res$Row.names
res$Row.names <- NULL

pdx_r <- prepare_df(pdx_r)
colnames(pdx_r) <- c('log2FC_PDX_R', 'padj_PDX_R')
res <- merge(res, pdx_r, by="row.names", all.x=TRUE)
rownames(res) <- res$Row.names
res$Row.names <- NULL




##
res <- res[order(-res$log2FoldChange, -res$baseMean,res$padj),]

mediane <- median(res$baseMean)
res$filter_PDX_S <- ifelse(!is.na(res$padj_PDX_S) & res$padj_PDX_S < 0.05 & res$log2FC_PDX_S > 0, TRUE, FALSE)
res$filter_PDX_R <- ifelse(!is.na(res$padj_PDX_R) & res$padj_PDX_R >= 0.05, TRUE, FALSE) 
res$filter_expr <- ifelse(res$baseMean > mediane, TRUE, FALSE)
res$allfilter <- apply(res[,grepl('filter_', colnames(res))], 1, all)

nrow(deg)
nrow(res)
table(res$filter_PDX_S)
table(res$filter_PDX_R)
table(res$filter_PDX_S & res$filter_PDX_R)
table(res$filter_PDX_S & res$filter_PDX_R & res$filter_expr)
table(res$allfilter)

#write.table(res, gzfile('18v1_deg_filtersteps_tractability_group.tsv.gz'), sep="\t", quote=F, row.names=T)

res_new <- res



all(rownames(res)== rownames(res_new))

m <- merge(res, res_new, by='row.names')
m[m$allfilter.x != m$allfilter.y,]
table(m[m$allfilter.x != m$allfilter.y,'allfilter.y'])

buckets <- read.table("../../local/share/data/Supplementary_Table_9_list_bucket.txt", sep="\t", header=T)
buck <- aggregate(buckets$Cancer.Tissue.Type, buckets[, c(2,3,4)], FUN=function(x) {as.character(paste0(unique(x), collapse=","))} , simplify=T)
colnames(buck)[4] <- 'cancer_type'

mres <- merge(res_new, buck, by.x="row.names", by.y="Target", all.x=TRUE)
rownames(mres) <- mres$Row.names
mres$Row.names <- NULL
#write.table(res, gzfile('deg_filtersteps_tractability_group.tsv.gz'), sep="\t", quote=F, row.names=T)
selection <- mres[!is.na(mres$Tractability.Group) & mres$Tractability.Group <=1 & mres$allfilter,] 

rownames(res_new[res$allfilter != res_new$allfilter & !res_new$allfilter,])

