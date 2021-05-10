
d_shallow <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/qdnaseq_multi/CRC0327-04-0.tsv.gz'), sep="\t", header=F)
d_30x <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset/CRC0327/sequenza/CRC0327-04-0.15000.tsv.gz'), sep="\t", header=F)

colnames(d_shallow) <- c('chr','b','e','cn')
colnames(d_30x) <- c('chr','b','e','signal')
d_30x$logfc <- log2(d_30x$signal)
