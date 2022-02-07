# TODO put in Snakemake with correct dependencies and input files and split xeno/pdo with wildcard
# get wanted list of models from input files or wildcard?
library(QDNAseq)
wanted <- c('CRC1449','CRC1917','CRC1337', 'CRC1342' ,'CRC1563', 'CRC1566')

setwd('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/cn_plots')

names <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/map_lmo_withS.tsv', sep="\t", header=F)
load('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/pdo_qdnaseq/cn_seg.Rdata')
colnames(names) <- c('id', 'genealogy')
names$smodel <- substr(names[,2], 0, 7)
wanted_long_gen <- names[names$smodel %in% wanted,]

#plot(copyNumbersSegmented[,colnames(copyNumbersSegmented) == 'CRC2113LMO0B01004001D02000'])
for (i in seq(1, nrow(wanted_long_gen))) {
  long_gen <- wanted_long_gen[i, 'genealogy']
  seen <- colnames(copyNumbersSegmented) == long_gen
  if(any(seen)) {
    pdf(paste0(long_gen, '.pdf'))
    plot(copyNumbersSegmented[,seen])
    graphics.off()
  } 
}


names <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/map_lmx_withS.tsv', sep="\t", header=F)
load('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/xeno_qdnaseq/cn_seg.Rdata')
colnames(names) <- c('id', 'genealogy')
names$smodel <- substr(names[,2], 0, 7)
wanted_long_gen <- names[names$smodel %in% wanted,]

#plot(copyNumbersSegmented[,colnames(copyNumbersSegmented) == 'CRC2113LMO0B01004001D02000'])
for (i in seq(1, nrow(wanted_long_gen))) {
  long_gen <- wanted_long_gen[i, 'genealogy']
  seen <- colnames(copyNumbersSegmented) == long_gen
  if(any(seen)) {
    pdf(paste0(long_gen, '.pdf'))
    plot(copyNumbersSegmented[,seen])
    graphics.off()
  } 
}

### Primo's MSI

library(QDNAseq)
wanted <- c('CRC1448')

setwd('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/cn_plots')

names <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/map_lmo_withS.tsv', sep="\t", header=F)
load('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/pdo_qdnaseq/cn_seg.Rdata')
colnames(names) <- c('id', 'genealogy')
names$smodel <- substr(names[,2], 0, 7)
wanted_long_gen <- names[names$smodel %in% wanted,]

  #plot(copyNumbersSegmented[,colnames(copyNumbersSegmented) == 'CRC2113LMO0B01004001D02000'])
  for (i in seq(1, nrow(wanted_long_gen))) {
    long_gen <- wanted_long_gen[i, 'genealogy']
    seen <- colnames(copyNumbersSegmented) == long_gen
    if(any(seen)) {
      jpeg(paste0(long_gen, '.jpg'))
      plot(copyNumbersSegmented[,seen])
      graphics.off()
    } 
  }
