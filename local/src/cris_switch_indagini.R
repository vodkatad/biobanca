d <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_model_cris-right.tsv', sep="\t", header=T)
head(d)
x <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_prediction_result_nc.tsv', sep="\t", header=T)
o <- read.table('/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMO_BASALE_prediction_result_nc.tsv', sep="\t", header=T)
wrong <- d[d$CRIS_PDO != d$CRIS_PDX,'model']
#x$wrong <- ifelse(x$sample.names)

x$model <- substr(x$sample.names,0,7)
x$wrong <- ifelse(x$model %in% wrong, 'switch','conserved')
library(ggsignif)
library(ggplot2)
ggplot(data=x, aes(x=wrong, y=dist.to.template))+geom_boxplot()+theme_bw()+geom_signif(comparisons = list(c("conserved", "switch")))


o$model <- substr(o$sample.names,0,7)
o$wrong <- ifelse(o$model %in% wrong, 'switch','conserved')
ggplot(data=o, aes(x=wrong, y=dist.to.template))+geom_boxplot()+theme_bw()+geom_signif(comparisons = list(c("conserved", "switch")))
