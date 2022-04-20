library(ggplot2)
library(ggsignif)

counts <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/merged_hs_mm.tsv.gz'), sep="\t", header=TRUE)
meta <- read.table('/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final', sep="\t", header=TRUE, stringsAsFactors = FALSE)

meta$sample_id_R <- gsub('-2', '.2', meta$sample_id_R, fixed=TRUE)
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL

totreads <- colSums(counts)
dividetot <- t(counts)/totreads
dividetot <- t(dividetot)
cpm <- dividetot * 10**6

log2pc1cpm <- log2(cpm+1)

l2cpm_m <- log2pc1cpm[grepl('^M_', rownames(log2pc1cpm)),]
l2cpm_h <- log2pc1cpm[grepl('^H_', rownames(log2pc1cpm)),]

totH <- colSums(l2cpm_h)
totM <- colSums(l2cpm_m)
stopifnot(all(names(totH)==names(totM)))
total_reads <- data.frame(mouse=totM, human=totH, row.names = names(totH))

meta_total_reads <- merge(total_reads, meta, by.x="row.names", by.y='sample_id_R')
stopifnot(nrow(meta_total_reads) == nrow(total_reads))
stopifnot(nrow(meta_total_reads) == nrow(meta))
meta_total_reads$sample <- ifelse(grepl('H', meta_total_reads$type), 'human', ifelse(grepl('X', meta_total_reads$type), 'PDX', 'PDO'))

ggplot(data=meta_total_reads, aes(x=sample, y=mouse, color=sample))+geom_boxplot()+theme_bw()+
  ylab('total mouse log2CPM')+scale_color_manual(values=c("PDO"="darkblue", "PDX"= "firebrick1", "human"="bisque3"))

#geom_signif( comparisons = list(c("human", "xenograft"), c("organoid", 'xenograft')))+

ggplot(data=meta_total_reads, aes(x=sample, y=human, color=sample))+geom_boxplot()+theme_bw()+
  ylab('total human log2CPM')+scale_color_manual(values=c("PDO"="darkblue", "PDX"= "firebrick1", "human"="bisque3"))

ktm <- kruskal.test(formula=as.formula(mouse ~ sample), data=meta_total_reads)
kth <- kruskal.test(formula=as.formula(human ~ sample), data=meta_total_reads)
ktm$p.value
kth$p.value

head(meta_total_reads[order(meta_total_reads$human),])

#egrassi@godot:/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK$ grep CRC0164LMX0B02202TUMR02000 removed_samples_readsQC 
#CRC0164LMX0B02202TUMR02000
#egrassi@godot:/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK$ grep CRC0427LMX0A02004TUMR02000 removed_samples_readsQC 
#CRC0427LMX0A02004TUMR02000
#egrassi@godot:/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK$ grep CRC0771LMX0B02003TUMR01000 removed_samples_readsQC 
#CRC0771LMX0B02003TUMR01000
#egrassi@godot:/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK$ grep CRC1624PRX0B02001TUMR09000 removed_samples_readsQC 
#egrassi@godot:/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK$ grep CRC0427LMX0A02004TUMR02000 selected_metadata_annot_final_nolinfo_nooutlier_replisafe 
#egrassi@godot:/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK$ grep CRC0164LMX0B02202TUMR02000 selected_metadata_annot_final_nolinfo_nooutlier_replisafe
#egrassi@godot:/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK$ grep CRC0771LMX0B02003TUMR01000 selected_metadata_annot_final_nolinfo_nooutlier_replisafe
#egrassi@godot:/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK$ 
pdo <- meta_total_reads[meta_total_reads$sample == 'PDO',]
tail(pdo[order(pdo$mouse),])

ggplot(data=pdo, aes(x=sample, y=mouse, color=sample))+geom_boxplot()+theme_bw()+
  ylab('total mouse log2CPM')+ylim(c(0, 100))


ggplot(data=meta_total_reads, aes(x=sample, y=mouse, color=sample))+geom_boxplot()+theme_bw()+
  ylab('total mouse log2CPM')+scale_color_manual(values=c("PDO"="darkblue", "PDX"= "firebrick1", "human"="bisque3"))+
scale_y_log10()

ggplot(data=meta_total_reads, aes(x=type, y=mouse, color=sample))+geom_boxplot()+theme_bw()+
  ylab('total mouse log2CPM')+scale_color_manual(values=c("PDO"="darkblue", "PDX"= "firebrick1", "human"="bisque3"))+
  scale_y_log10()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data=meta_total_reads, aes(x=type, y=human, color=sample))+geom_boxplot()+theme_bw()+
  ylab('total mouse log2CPM')+scale_color_manual(values=c("PDO"="darkblue", "PDX"= "firebrick1", "human"="bisque3"))+
  scale_y_log10()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

