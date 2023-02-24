library(ggplot2)
library(pheatmap)
setwd('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S')

deg <- read.table('treat_cutoff0.05-cetuxi.vs.NT.deseq2_sign_genedesc.tsv', comment.char="", quote='', header=T, sep="\t")
rownames(deg) <- deg$gene
deg$gene <- NULL

## heatmap of DEG ########
#expr <- read.table(gzfile('fpkm.tsv.gz'), sep="\t", header=T)
# rownames(expr) <- substr(rownames(expr), 3, nchar(rownames(expr)))
# samples_data <- read.table('samples_data', sep="\t", header=T)
# rownames(samples_data) <- samples_data$id
# samples_data$id <- NULL
# 
# THR <- 0.5849625
THR <- 0
PC <- 1
deg_lf <- deg[deg$log2FoldChange > THR,]

#deg_expr <- expr[rownames(expr) %in% rownames(deg_lf),]


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
#pdx_r <- read.table('../Biodiversa_up5starOK_cetuxi_treat_PDX_R/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv', comment.char="", quote='', header=T, sep="\t")

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


res_new <- res

dgi <- read.table('/scratch/trcanmed/biobanca/local/share/data/dgidb_101.csv', sep="\t", header=TRUE, quote="", comment.char="")
dgig <- unique(dgi$search_term)
dgig <- as.character(dgig)
dgig <- c(dgig, "ABCA1")
dgig <- c(dgig, "MAOA") #added the 23/02/23 to reach his friend.
dgi_genes <- data.frame(gene.name=dgig, dgi=rep('druggable', length(dgig)))
# dgi_genes <- as.data.frame(dgig)
# names(dgi_genes)[names(dgi_genes)=="dgig"] <- "gene.name"
# dgi_genes$gene.name <- as.character(dgi_genes$gene.name)
# dgi_genes <- rbind(dgi_genes, "ABCA1")
# dgi_genes$dgi <- "DGI102"
# add ABCA1 as dgi


length(dgig)
length(intersect(row.names(res_new), dgig))
res_new_m <- merge(res_new, dgi_genes, by.x="row.names", by.y='gene.name', all.x=TRUE, sort=FALSE)
rownames(res_new_m) <- res_new_m$Row.names
res_new_m$Row.names <- NULL
# write.table(res_new_m, gzfile('18v1_deg_filtersteps_tractability_group.tsv.gz'), sep="\t", quote=F, row.names=T)
# 
# dim(res_new_m)
# dim(res_new)
# rownames(res_new_m) <- res_new_m$Row.names
# res_new_m$Row.names <- NULL
# res_new_m <- res_new_m[match(rownames(res_new), rownames(res_new_m)),]
# all(rownames(res)== rownames(res_new_m))
# nuovi <- res_new_m[res$allfilter != res_new_m$allfilter & res_new_m$allfilter,]
# 
# 
# 
# 
# write.table(nuovi, '~/stronzi_nodgi.tsv', sep="\t", quote=F)

buckets <- read.table("../../local/share/data/Supplementary_Table_9_list_bucket.txt", sep="\t", header=T)
buck <- aggregate(buckets$Cancer.Tissue.Type, buckets[, c(2,3,4)], FUN=function(x) {as.character(paste0(unique(x), collapse=","))} , simplify=T)
colnames(buck)[4] <- 'cancer_type'

mres <- merge(res_new_m, buck, by.x="row.names", by.y="Target", all.x=TRUE)
rownames(mres) <- mres$Row.names
mres$Row.names <- NULL
# keep only the new col Tractability.Group
mres$Tractability.Group <- NULL
mres$cancer_type <- NULL
## add the two manually selected cols
mres$genes <- rownames(mres)
bucket_chosen <- c("STAT3", "OGT", "FGFR1", "FGFR2")
mres["STAT3","Tractability.Bucket_manuallyselected"] <- "selected"
mres["OGT","Tractability.Bucket_manuallyselected"] <- "selected"
mres["FGFR1","Tractability.Bucket_manuallyselected"] <- "selected"
mres["FGFR2","Tractability.Bucket_manuallyselected"] <- "selected"

dgi_chosen <- c("NUAK2", "HDAC5", "HDAC6", "HDAC10", "HDAC11", "ULK1",
                "SMO", "DGAT1", "MERTK", "NOTCH1", "FGFR1", "FGFR2",
                "MKNK1", "CLK1", "DPEP1")
mres["NUAK2","DGIdb_manuallyselected"] <- "selected"
mres["HDAC5","DGIdb_manuallyselected"] <- "selected"
mres["HDAC6","DGIdb_manuallyselected"] <- "selected"
mres["HDAC10","DGIdb_manuallyselected"] <- "selected"
mres["HDAC11","DGIdb_manuallyselected"] <- "selected"
mres["ULK1","DGIdb_manuallyselected"] <- "selected"
mres["SMO","DGIdb_manuallyselected"] <- "selected"
mres["DGAT1","DGIdb_manuallyselected"] <- "selected"
mres["MERTK","DGIdb_manuallyselected"] <- "selected"
mres["NOTCH1","DGIdb_manuallyselected"] <- "selected"
mres["FGFR1","DGIdb_manuallyselected"] <- "selected"
mres["FGFR2","DGIdb_manuallyselected"] <- "selected"
mres["MKNK1","DGIdb_manuallyselected"] <- "selected"
mres["CLK1","DGIdb_manuallyselected"] <- "selected"
mres["DPEP1","DGIdb_manuallyselected"] <- "selected"

names(mres)[names(mres)=="dgi"] <- "DGIdb"
mres$genes <- NULL
#write.table(res, gzfile('deg_filtersteps_tractability_group.tsv.gz'), sep="\t", quote=F, row.names=T)
selection <- mres[!is.na(mres$Tractability.Group) & mres$Tractability.Group <=1 & mres$allfilter,] 

rownames(res_new[res$allfilter != res_new$allfilter & !res_new$allfilter,])

write.table(mres, file ="deg_filtersteps_tractability_group.tsv", quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = TRUE)

write.table(mres, file ="/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/deg_filtersteps_tractability_group.tsv", quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = TRUE)
