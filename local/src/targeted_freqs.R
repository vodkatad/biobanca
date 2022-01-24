library(ggplot2)

xeno_af <- snakemake@input[['xenoTiers']]
pdo_af <- snakemake@input[['pdoTiers']]
xeno_genes <- snakemake@input[['xenoGenes']]
pdo_genes <- snakemake@input[['pdoGenes']]
treat_f <- snakemake@input[['treat']]
msk_f <- snakemake@input[['msk']]
tcga_f <- snakemake@input[['tcga']]

tcgamsk_xeno_f <- snakemake@output[['TCGAMSK_xeno']]
tcgamsk_pdo_f <- snakemake@output[['TCGAMSK_pdo']]
tcgamsk_xeno_zoom_f <- snakemake@output[['TCGAMSK_xeno_zoom']]
tcgamsk_pdo_zoom_f <- snakemake@output[['TCGAMSK_pdo_zoom']]
savedata <- snakemake@output[['preprocGeneAF']]

thr <- as.numeric(snakemake@wildcards[['AF']])

osnakemake <- snakemake
load(snakemake@input[['Rimage']])
#eval(parse(text=myriad))
snakemake <- osnakemake

pdx <- read.table(xeno_af, header=TRUE, sep="\t", row.names=1)
pdo <- read.table(pdo_af, header=TRUE, sep="\t", row.names=1)


gpdx <- read.table(xeno_genes, header=TRUE, sep="\t", row.names=1)
gpdo <- read.table(pdo_genes, header=TRUE, sep="\t", row.names=1)

#gpdx <- read.table('./mutect/merged.table_nomultiallele_wtiers', header=TRUE, sep="\t", row.names=1)
#gpdo <- read.table('../biobanca_targeted_pdo/mutect/merged.table_nomultiallele_wtiers', header=TRUE, sep="\t", row.names=1)

###

#all_genes <- intersect(unique(gpdo$genes), unique(gpdx$genes)) # y?
all_genes <- union(unique(gpdo$genes), unique(gpdx$genes))

get_genes_binary <- function(sample, AF, genes, thr, all_genes) {
    keep <- AF[,sample] > thr
    genes <- genes[keep,]
    return(all_genes %in% genes$genes)
}

gpdo <- gpdo[match(rownames(pdo), rownames(gpdo)),, drop=F] # we keep info for muts available in our AF tables
gpdx <- gpdx[match(rownames(pdx), rownames(gpdx)),, drop=F]
pdobing <- sapply(colnames(pdo), get_genes_binary, pdo, gpdo, thr, all_genes)
pdxbing <- sapply(colnames(pdx), get_genes_binary, pdx, gpdx, thr, all_genes)
rownames(pdobing) <- all_genes
rownames(pdxbing) <- all_genes # TRUE FALSE at the level of genes

pdo_genes <- apply(pdobing, 1, sum)
pdx_genes <- apply(pdxbing, 1, sum)
freqs <- data.frame(pdo_freq = pdo_genes/ncol(pdobing), pdx_freq=pdx_genes/ncol(pdxbing))
freqs$gene <- rownames(freqs)
tcga <- read.table(tcga_f, sep="\t")
msk <- read.table(msk_f, sep="\t")
colnames(tcga) <- c('gene', 'freq')
colnames(msk) <- c('gene', 'freq')

m1 <- merge(tcga, msk, all.x = TRUE, all.y = TRUE, by='gene')
colnames(m1) <- c('gene','tcga_freq','msk_freq')
m2 <- merge(m1, freqs, all.x = TRUE, all.y = TRUE, by='gene')

save.image('pippo.Rdata')
m2[is.na(m2$msk_freq),]$msk_freq <- 0
#m2[is.na(m2$tcga_freq),]$tcga_freq <- 0 # no missing here
m2[is.na(m2$pdx_freq),]$pdx_freq <- 0
m2[is.na(m2$pdo_freq),]$pdo_freq <- 0


cl <- c(rep('tcga',nrow(m2)), rep('msk',nrow(m2)))

m3 <- data.frame(x=c(m2$pdo_freq, m2$pdo_freq), 
                 y=c(m2$tcga_freq, m2$msk_freq), 
                 class=cl)

lmplot <- function(data, title, xlim=NULL, out) {
  fit1 <- data[data$class=='tcga',]
  fit2 <- data[data$class=='msk',]
  pi1 <- cor.test(fit1$x, fit1$y)
  pi2 <- cor.test(fit2$x, fit2$y)
  
  cap1 <- round(c(pi1$estimate, pi1$p.value, pi2$estimate, pi2$p.value),2)
  print(cap1)
  names(cap1) <- NULL
  cap <- paste(cap1, collapse= " ")
  if (is.null(xlim)) {
    ggplot(data=data, aes(x=x, y=y, color=class))+ geom_point()+ geom_smooth(method=lm, se=FALSE)+unmute_theme+labs(caption=cap)+ggtitle(title)
  } else {
    ggplot(data=data, aes(x=x, y=y, color=class))+ geom_point()+ geom_smooth(method=lm, se=FALSE)+unmute_theme+labs(caption=cap)+ggtitle(title)+xlim(xlim)
  }
  ggsave(out)
}

lmplot(m3, title='PDO', out=tcgamsk_pdo_f)
lmplot(m3, title='PDO', xlim=c(0, 0.2), out=tcgamsk_pdo_zoom_f)

m3 <- data.frame(x=c(m2$pdx_freq, m2$pdx_freq), 
                 y=c(m2$tcga_freq, m2$msk_freq), 
                 class=cl)
lmplot(m3, title='PDX', out=tcgamsk_xeno_f)
lmplot(m3, title='PDX', xlim=c(0, 0.2), out=tcgamsk_xeno_zoom_f)


save.image(savedata)