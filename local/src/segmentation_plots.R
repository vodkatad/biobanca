# TODO put in Snakemake with correct dependencies and input files and split xeno/pdo with wildcard
# get wanted list of models from input files or wildcard?
library(QDNAseq)
#wanted <- c('CRC1449','CRC1917','CRC1337', 'CRC1342' ,'CRC1563', 'CRC1566')
#wanted <- c('CRC1888')
#wanted <- c('CRC0578')
#wanted <- c('CRC1575')
#wanted <- c('CRC1251', 'CRC0729', 'CRC0370')
wanted <- c('CRC1239')

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

### GS

library(QDNAseq)
#wanted <- c('CRC0123') #CRC0355LMO0B01006001D02000 CRC1620LMO0B01004001D02000
wanted <- c('CRC1432')

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

### switch in adjacent segments lfc
d <- read.table('/scratch/trcanmed/biobanca/dataset/V1/shallowseq/pdo_qdnaseq/cn_segm.tsv', sep='\t', header=T)
head(d)
rownames(d) <- d$feature
f$feature <- NULL
head(colnames(d))
d$feature <- NULL
d$chromosome <- NULL
d$start <- NULL
d$end <- NULL
rle(diff(c(1,1,0)))
diff(c(1,1,0))
diff(c(1,1,1,0))
rle(diff(c(1,1,1,0)))
tt <- rle(diff(d$CRC0123LMO0A01003001D02000))
tt
m$lengths[m$values==0]
tt$lengths[tt$values==0]
sum(tt$lengths[tt$values==0])
sum(tt$lengths[tt$values!=0])
rles <- apply(d, 2, function(x){rle(diff(x))})
head(rles)
changeLfc <- sapply(rles, function(tt) { sum(tt$lengths[tt$values!=0])}
)
changeLfc <- sapply(rles, function(tt) { sum(tt$lengths[tt$values!=0])})
sameLfc <- sapply(rles, function(tt) { sum(tt$lengths[tt$values==0])})
head(sameLfc)
summary(sameLfc)
summary(changeLfc)
changeLfc[order(changedFiles())]
changeLfc[order(changeLfc]
changeLfc[order(changeLfc)]
head(changeLfc[order(changeLfc)])





