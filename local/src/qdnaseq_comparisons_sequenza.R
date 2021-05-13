
library(corrplot)
library(RColorBrewer)
samples <- c('CRC0327-04-0','CRC1078-02-0','CRC1078-02-1-C')

sfx <- '.tsv.gz'
sfx2 <- '.15000.tsv.gz'
p1 <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/qdnaseq_multi/'
p2 <- '/scratch/trcanmed/AF_spectra/dataset/CRC0327/sequenza/'
p3 <- '/scratch/trcanmed/AF_spectra/dataset/CRC1078/sequenza/'


load_sh <- function(x, pre, su) {
  d <- read.table(gzfile(paste0(pre, x, su)), sep="\t", header=F)
  colnames(d) <- c('chr','b','e','cn')
  d$n <- paste0(d$chr,"_", d$b)
  #d$cn <- log2(d$cn+1) # almost the same
  return(d[d$cn != 0,])
}


load_wgs <- function(x, pre, su) {
  d <- read.table(gzfile(paste0(pre, x, su)), sep="\t", header=F)
  colnames(d) <- c('chr','b','e','cn')
  d$n <- paste0(d$chr,"_", d$b)
  return(d)
}

l_sh <- lapply(samples, load_sh, p1, sfx)
l_wgs <- lapply(samples[2:3], load_wgs, p3, sfx2)
all <- c(l_sh,  list(load_wgs(samples[[1]], p2, sfx2)), l_wgs)

names(all) <- c(paste0('sh_', samples), paste0('wgs_',samples))

get_corr <- function(pair, all) {
  pair <- unlist(pair)
  n1 <- as.character(pair[[1]])
  n2 <- as.character(pair[[2]])
  m <- merge(all[[n1]], all[[n2]], by="n")
  cc <- cor.test(m$cn.x, m$cn.y)
  return(cc$estimate)
}


pairs <- expand.grid(names(all), names(all))

pears <- apply(pairs, 1, get_corr, all)
res <- matrix(pears, nrow=6)
colnames(res) <- names(all)
rownames(res) <- names(all)

corrplot.mixed(res)
corrplot.mixed(res,  lower.col = "black", upper.col= rev(brewer.pal(n = 8, name = "RdBu")), tl.pos="d", tl.cex=1.2, number.cex=1.8, cl.cex=1.8)

get_corr <- function(pair, all) {
  pair <- unlist(pair)
  n1 <- as.character(pair[[1]])
  n2 <- as.character(pair[[2]])
  m <- merge(all[[n1]], all[[n2]], by="n")
  cc <- cor.test(m$cn.x, m$cn.y)
  return(cc$p.value)
}


pairs <- expand.grid(names(all), names(all))

pears <- apply(pairs, 1, get_corr, all)
pres <- matrix(pears, nrow=6)
colnames(pres) <- names(all)
rownames(pres) <- names(all)


##########

d_shallow <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/qdnaseq_multi/CRC0327-04-0.tsv.gz'), sep="\t", header=F)
d_30x <- read.table(gzfile('/scratch/trcanmed/AF_spectra/dataset/CRC0327/sequenza/CRC0327-04-0.15000.tsv.gz'), sep="\t", header=F)

colnames(d_30x) <- c('chr','b','e','cn')
colnames(d_shallow) <- c('chr','b','e','signal')
d_30x$logfc <- log2(d_30x$signal)

d_shallow$n <- paste0(d_shallow$chr,"_", d_shallow$b)
d_30x$n <- paste0(d_30x$chr,"_", d_30x$b)

dd_shallow <- d_shallow[d_shallow$signal != 0,]

m <- merge(d_30x, dd_shallow, by="n")