library(ggplot2)

xeno_af <- snakemake@input[['xenoAF']]
pdo_af <- snakemake@input[['pdoAF']]
xeno_waf <- snakemake@input[['xenoWAF']]
pdo_waf <- snakemake@input[['pdoWAF']]
vs <- snakemake@params[['vs']]
log <- snakemake@log[['log']]

osnakemake <- snakemake
load(snakemake@input[['Rimage']])
#eval(parse(text=myriad))
snakemake <- osnakemake

xeno_df <- read.table(xeno_af, header=TRUE, sep="\t", row.names=1)
pdo_df <- read.table(pdo_af, header=TRUE, sep="\t", row.names=1)

sink(snakemake@log[[1]])
print(paste0('Starting with x= ', ncol(xeno_df), ' o= ', ncol(pdo_df)))
sink()

samples <- colnames(pdo_df)
if (sum(grepl('LMO', samples)) != ncol(pdo_df)) {
    stop('I was expecting only LMOs!')
}
passage <- as.numeric(substr(samples, 15, 17))
if (vs == "earlyxeno") {
    pdo_df <- pdo_df[, passage <= 3]
} else if (vs == "latexeno") {
    pdo_df <- pdo_df[, passage > 3]
} else {
    stop('Still to be implemented!')
}

sink(snakemake@log[[1]])
print(paste0(vs, ' Starting with x= ', ncol(xeno_df), ' o= ', ncol(pdo_df)))
sink()


filter_plot <- function(pdo, xeno, out) {
    # filter qc here if needed TODO

    sink(snakemake@log[[1]], append=TRUE)
    print(paste0('Filter pairs x=', ncol(xeno), ' o= ', ncol(pdo)))
    sink()

    models_xeno <- substr(colnames(xeno), 0, 7)
    #models_pdo <- substr(colnames(pdo_df), 0, 7)
    colnames(xeno) <- models_xeno
    colnames(pdo) <- substr(colnames(pdo), 0, 7)

    com <- length(intersect(models_xeno, colnames(pdo)))
    # selection of the good one only (validated PDOs)
    sink(snakemake@log[[1]], append=TRUE)
    print(paste0('After validation filtering we kept ', length(models_xeno), ' starting from x= ', ncol(xeno), 'and o= ', ncol(pdo), ' -> ', com))
    sink()

    xeno <- xeno[, models_xeno[order(models_xeno)]]
    pdo <- pdo[, models_xeno[order(models_xeno)]]

    afx <- unlist(xeno[xeno!=0])
    afo <- unlist(pdo[pdo!=0])

    ks <- ks.test(afx, afo)
    sink(snakemake@log[[1]], append=TRUE)
    print(paste0('ks test ', out,  ': ', ks$p.value))
    pd <- data.frame(af=c(afo, afx), class=c(rep('pdo', length(afo)),rep('pdx',length(afx))))

    #ggplot(data=pd, aes(x=af, color=class, y=..count..)) + geom_density()+scale_x_continuous(breaks=(seq(0, 1, by=0.05)))+unmute_theme
    ggplot(data=pd, aes(x=af, fill=class)) + geom_histogram(alpha=0.5, position='dodge')+
    scale_x_continuous(breaks=(seq(0, 1, by=0.05)))+unmute_theme+
    scale_fill_manual(values=c('darkblue', 'firebrick1'))
    ggsave(out)
    #ggplot(data=pd, aes(x=af, color=class)) + geom_density()+scale_x_continuous(breaks=(seq(0, 1, by=0.05)))+unmute_theme
    return(list(x=xeno, o=pdo))
}

af_all <- filter_plot(pdo_df, xeno_df, snakemake@output[['freqAll']])
##
xeno_df_w <- read.table(xeno_waf, header=TRUE, sep="\t", row.names=1)
pdo_df_w <- read.table(pdo_waf, header=TRUE, sep="\t", row.names=1)

samples <- colnames(pdo_df_w)
if (sum(grepl('LMO', samples)) != ncol(pdo_df_w)) {
    stop('I was expecting only LMOs!')
}
passage <- as.numeric(substr(samples, 15, 17))
if (vs == "earlyxeno") {    
    pdo_df_w <- pdo_df_w[, passage <= 3]
} else if (vs == "latexeno") {
    pdo_df_w <- pdo_df_w[, passage > 3]
} else {
    stop('Still to be implemented!')
}


af_tiers <- filter_plot(pdo_df_w, xeno_df_w, snakemake@output[['freqTiers']])

plot_n <- function(pdo, pdx, out) {
    pdxbin <- ifelse(pdx==0,0, 1)
    n_muts_x <- colSums(pdxbin)
    pdobin <- ifelse(pdo==0,0, 1)
    n_muts_o <- colSums(pdobin)
    #hist(n_muts_x, breaks=20)
    #hist(n_muts_o, breaks=20)

    ks <- ks.test(n_muts_o, n_muts_x)
    sink(snakemake@log[[1]], append=TRUE)
    print(paste0('ks test ', out,  ': ', ks$p.value))
    
    pd <- data.frame(af=c(n_muts_o, n_muts_x), class=c(rep('pdo', length(n_muts_o)),rep('pdx',length(n_muts_x))))
    #ggplot(data=pd, aes(x=af, color=class,y=..count..)) + geom_density()+xlab('nmuts')+unmute_theme
    ggplot(data=pd, aes(x=af, fill=class)) + geom_histogram(alpha=0.5, position='dodge')+xlab('nmuts')+unmute_theme+
    scale_fill_manual(values=c('darkblue', 'firebrick1'))
    ggsave(out)
}
plot_n(af_all[['o']], af_all[['x']], snakemake@output[['histAll']])
plot_n(af_tiers[['o']], af_tiers[['x']], snakemake@output[['histTiers']])

write.table(af_tiers[['o']], file=snakemake@output[['pdoTiers']], sep="\t", quote=FALSE)
write.table(af_tiers[['x']], file=snakemake@output[['xenoTiers']], sep="\t", quote=FALSE)

save.image('pippo.Rdata')