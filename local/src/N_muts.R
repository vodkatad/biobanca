library(ggplot2)

xeno_af <- snakemake@input[['xenoTiers']]
pdo_af <- snakemake@input[['pdoTiers']]
thr <- as.numeric(snakemake@wildcards[['AF']])

osnakemake <- snakemake
load(snakemake@input[['Rimage']])
#eval(parse(text=myriad))
snakemake <- osnakemake

pdx <- read.table(xeno_af, header=TRUE, sep="\t", row.names=1)
pdo <- read.table(pdo_af, header=TRUE, sep="\t", row.names=1)


pdxbin <- ifelse(pdx > thr, 1, 0)
n_muts_x <- colSums(pdxbin)
pdobin <- ifelse(pdo > thr, 1, 0)
n_muts_o <- colSums(pdobin)
#hist(n_muts_x, breaks=20)
#hist(n_muts_o, breaks=20)
pd <- data.frame(af=c(n_muts_o, n_muts_x), class=c(rep('pdo', length(n_muts_o)),rep('pdx',length(n_muts_x))))

ggplot(data=pd, aes(x=af, fill=class)) + geom_histogram(alpha=0.5, position='dodge')+xlab('nmuts')+unmute_theme+
    scale_fill_manual(values=c('darkblue', 'firebrick1'))
ggsave(snakemake@output[['histTiers']])
