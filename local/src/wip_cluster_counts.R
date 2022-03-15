
load('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/ioodioilclustering.Rdata')
# run begin of cluster_fin.R
has_smt <- function(smodel, data, smt) {
  subset <- data[data$shortgen==smodel,]
  any(grepl(smt,  rownames(subset)))
}

all_pdo <- unique(as.character(clustered_data$shortgen))
table(sapply(all_pdo, has_smt, clustered_data, 'LMO'))
table(sapply(all_pdo, has_smt, clustered_data, 'LMH'))
table(sapply(all_pdo, has_smt, clustered_data, 'LMX'))

clustered_data$shortgen <- as.character(clustered_data$shortgen)
clustered_data <- clustered_data[clustered_data$shortgen!="CRC1451",]

only_one <- function(smodel, data) {
  subset <- data[data$shortgen == smodel,]
  return(length(unique(subset$cl_id)))
}

ncl <- sapply(all_pdo, only_one, clustered_data)
sum(ncl==1)
length(ncl)
sum(ncl==1)/length(ncl)
sum(ncl>1)/length(ncl)
sum(ncl>1)



indagine <- names(ncl[ncl>1])
ind <- clustered_data[clustered_data$shortgen %in% indagine, ]
ind[order(ind$shortgen),]


data <- data.frame(smodel=all_pdo)
data$switch <- ifelse(data$smodel %in% indagine, 'yes', 'no')

crisx <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_nc_smodel.tsv', sep="\t", header=T)
criso <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMO_BASALE_nc_smodel.tsv', sep="\t", header=T)
dim(data)
m1 <- merge(data, crisx, by.x="smodel",by.y="genealogy")
dim(m1)
m2 <- merge(m1, criso, by="smodel",by.y="genealogy")
dim(m2)

pdf <- as.matrix(table(m2$switch, m2$cris.x))
pdf <- pdf / rowSums(pdf)
pdf <- as.data.frame(pdf)
ggplot(data=pdf, aes(x=Var2, fill=Var1, y=Freq))+geom_col(position="dodge")+xlab('CRIS PDX')+ggtitle('cluster unstable')

pdf <- as.matrix(table(m2$switch, m2$cris.y))
pdf <- pdf / rowSums(pdf)
pdf <- as.data.frame(pdf)
ggplot(data=pdf, aes(x=Var2, fill=Var1, y=Freq))+geom_col(position="dodge")+xlab('CRIS PDO')+ggtitle('cluster unstable')

# TODO guardare vsd_smodel_cris-merged_validated.tsv
#########
library(pheatmap)
h_scores <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/h-scores_LMX_BASALE.tsv', sep="\t", header=T)

s <- read.table('~/dnmt3b_h_l.tsv', sep="\t", header=T)
#s <- s[]

scores <- h_scores
scores <- as.data.frame(t(scores))
scores$smodel <- as.character(substr(rownames(scores), 0, 7))
#indexes <- match(rasdep$smodel, scores$smodel) # do we keep replicates or not? #nope
scores <- merge(s, scores, by="smodel")
scores <- scores[order(scores$annot), ]
rownames(scores) <- make.unique(as.character(scores$smodel))

annot <- scores[,c('chosen', 'annot')]

data <- scores[, grepl('HALLMARK', colnames(scores))]

pheatmap(t(data), annotation_col=annot, cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames=F)

ggplot(data=scores, aes(y=HALLMARK_WNT_BETA_CATENIN_SIGNALING, x=chosen))+geom_boxplot(outliers.shape=NULL)+theme_bw()+geom_jitter()

library(ggsignif)
ggplot(data=scores, aes(y=HALLMARK_WNT_BETA_CATENIN_SIGNALING, x=chosen))+
  geom_boxplot(outlier.shape=NA)+theme_bw()+geom_jitter()+
  geom_signif(comparisons = list(c("yes", "no")))

ggplot(data=scores, aes(y=HALLMARK_WNT_BETA_CATENIN_SIGNALING, x=annot))+
  geom_boxplot(outlier.shape=NA)+theme_bw()+geom_jitter()

############ randomiz clu leaves

tot_boot <- 1000
res_tot_boot <- data.frame(matrix(ncol = tot_boot, nrow = length(rownames(res_true))), row.names = rownames(res_true))
# creazione dataframe vuoto tot_boot colonne tot_campioni righe
shuffle_labels_compute_spread <- function(my_clustered_data) {
  cl_id_var <- my_clustered_data$cl_id
  n_el <- length(cl_id_var)
  cl_id_var_sh <- cl_id_var[sample.int(n=n_el, size=n_el)]
  my_clustered_data$cl_id <- cl_id_var_sh
  res <- as.data.frame(sapply(smodels, count_clu_pt2, my_clustered_data))
  return(res)
}

for (i in seq(1, tot_boot))  {
  res_tot_boot[,i] <- shuffle_labels_compute_spread(clustered_data)
}

table(apply(res_tot_boot, 2, function(x){length(x[x==1])}))
res_tb_mean <- as.data.frame(apply(res_tot_boot, 1, mean))

summary(res_tb_mean)

#############3 trying new heatmap

smodels <- unique(substr(colnames(expr_selected),0,7))
n_samples <- length(smodels)

mycolors <- colors()
mycolors <- mycolors[!grepl('gr[ae]y\\d+', mycolors)]
mycolors <- mycolors[!grepl('dark', mycolors)]
mycolors <- mycolors[!grepl('light', mycolors)]
mycolors <- mycolors[!grepl('white', mycolors)]
mycolors <- mycolors[!grepl('cornsilk', mycolors)]
mycolors <- mycolors[!grepl('lavenderblush', mycolors)]
mycolors <- mycolors[!grepl('ivory', mycolors)]
mycolors <- mycolors[!grepl('bisque', mycolors)]
mycolors <- mycolors[!grepl('azure', mycolors)]
mycolors <- mycolors[!grepl('mistyrose', mycolors)]
mycolors <- mycolors[!grepl('snow', mycolors)]
mycolors <- mycolors[!grepl('seashell', mycolors)]
mycolors <- mycolors[!grepl('[a-z]+[1-3]', mycolors)]

mycolors <- c(mycolors, 'grey72', 'grey30')

df_rainbow <- c(mycolors[sample.int(n_samples)])
names(df_rainbow) <- smodels

ann_colors <- list(class=c(X='darkblue', O='firebrick1'), smodel=df_rainbow)

annot <- data.frame(row.names=rownames(cors))
annot$smodel <- substr(rownames(annot), 0, 7)

annot$class <- substr(rownames(annot), 10, 10)


dist_pearson <- as.dist(1-cor(expr_selected, method="pearson"))
hc1 <- hclust(dist_pearson, method = "complete" )


pheatmap(cors, cluster_rows = hc1, cluster_cols=hc1, 
         show_rownames = FALSE, show_colnames = FALSE, 
         annotation_col=annot, annotation_colors=ann_colors, treeheight_row=0)


criso <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMO_BASALE_prediction_result_nc.tsv', sep="\t", header=T)
crisx <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_prediction_result_nc.tsv', sep="\t", header=T)
cris <- rbind(criso, crisx)
cris <- data.frame(row.names=cris$sample.names, cris=cris$predict.label2)

mannot <- merge(annot, cris, by="row.names")
rownames(mannot) <- mannot$Row.names
mannot$Row.names <- NULL

ann_colors <- list(class=c(X='darkblue', O='firebrick1'), smodel=df_rainbow, 
                   cris=c(`CRIS-A`='darkorange2', `CRIS-B`='brown3', `CRIS-C`='darkblue', 
                          `CRIS-D`='darkgreen', `CRIS-E`='aquamarine2'))



pheatmap(cors, cluster_rows = hc1, cluster_cols=hc1, 
         show_rownames = FALSE, show_colnames = FALSE, 
         annotation_col=mannot, annotation_colors=ann_colors, treeheight_row=0)
