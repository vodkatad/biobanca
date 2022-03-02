
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
