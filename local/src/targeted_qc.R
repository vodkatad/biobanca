library(ggplot2)
setwd('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_shallowseq_pdo')
qc_o <- read.table('all_aligned_dedup.tsv', sep="\t", header=F)
colnames(qc_o) <- c('id','nreads','cov')
all_qc_o <- read.table('all_aligned_dedup.tsv_old', sep="\t", header=F)
colnames(all_qc_o) <- c('id','nreads_all','cov_all')

qc <- merge(qc_o, all_qc_o, by="id")
qc$cov_all <- NULL

plot(qc$nreads, qc$nreads_all, main='PDO')
ggplot(data=qc, aes(x=nreads))+geom_histogram()+theme_bw()+theme(text=element_text(size=15))
ggplot(data=qc, aes(x=cov))+geom_histogram()+theme_bw()+theme(text=element_text(size=15))

setwd('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_shallowseq_pdx')
qc_x <- read.table('all_aligned_dedup.tsv', sep="\t", header=F)
colnames(qc_x) <- c('id','nreads','cov')
all_qc_x <- read.table('all_aligned_dedup.tsv_old', sep="\t", header=F)
colnames(all_qc_x) <- c('id','nreads','cov')
