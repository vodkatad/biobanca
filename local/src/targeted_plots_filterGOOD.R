library(ggplot2)

xeno_af <- snakemake@input[['xenoAF']]
pdo_af <- snakemake@input[['pdoAF']]
xeno_waf <- snakemake@input[['xenoWAF']]
pdo_waf <- snakemake@input[['pdoWAF']]
good_f <- snakemake@input[['good']]
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

# list of validated models for PDO, the one we want to keep
good_df <- read.table(good_f, sep="\t", stringsAsFactors=FALSE, header=TRUE)
keep <- good_df[good_df$buoni, , drop=FALSE]

filter_plot <- function(pdo, xeno, keep, out) {
    pdo <- pdo[, colnames(pdo) != "CRC0177LMO0D04017001D02000"] # remove the mutated pdo -> keep the wt one (exception)
    xeno <- xeno[,colnames(xeno) != "CRC0177LMX0B05001TUMD07000"] # remove the wt xeno -> keep the mutated (normal line)
#     egrassi@godot:/scratch/trcanmed/biobanca/dataset/V1/targeted$ cat  /mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdx/mutect/merged_longformat_wtiers.tsv | grep CRC0177 | grep EGFR | grep chr7:55160233:G:A
# 21220   chr7:55160233:G:A       EGFR    missense_variant:ENST00000275493.7:c.1393G>A:exon12:p.G465R     CRC0177LMX0B05001TUMD07000      0
# 21828   chr7:55160233:G:A       EGFR    missense_variant:ENST00000275493.7:c.1393G>A:exon12:p.G465R     CRC0177LMX0B07002TUMD05000      0.195
# egrassi@godot:/scratch/trcanmed/biobanca/dataset/V1/targeted$ cat  /mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdo/mutect/merged_longformat_wtiers.tsv | grep CRC0177 | grep EGFR | grep chr7:55160233:G:A
# 21428   chr7:55160233:G:A       EGFR    missense_variant:ENST00000275493.7:c.1393G>A:exon12:p.G465R     CRC0177LMO0A04008002D02000      0
# 22042   chr7:55160233:G:A       EGFR    missense_variant:ENST00000275493.7:c.1393G>A:exon12:p.G465R     CRC0177LMO0D04017001D02000      0.276
    list_remove_xeno <- c("CRC1875") # for targeted we miss only CRC1875
    list_lmh <- colnames(pdo_df)[grepl('LMH', colnames(pdo_df))]
    model_lmh <- c(substr(list_lmh, 0,7), 'CRC0282')
    #if (which == "xo") {
    #expected_n <- expected_n - 1  # removal of CRC0177 brings us here
    pdo <- pdo[, !colnames(pdo) %in% list_lmh]
    pdo <- pdo[, !substr(colnames(pdo),0,7) %in% model_lmh]
    xeno <- xeno[, !substr(colnames(xeno),0,7) %in% model_lmh]

    #> setdiff(models_xeno, models_pdo)
    #[1] "CRC1870" "CRC1875" "CRC2041"
    # We know we were missing some PDOs TODO load from file if qc changes in the future
    xeno <- xeno[, !substr(colnames(xeno),0,7) %in% list_remove_xeno]

    sink(snakemake@log[[1]], append=TRUE)
    print(paste0('Filter pairs x=', ncol(xeno), ' o= ', ncol(pdo)))
    sink()

    models_xeno <- substr(colnames(xeno), 0, 7)
    #models_pdo <- substr(colnames(pdo_df), 0, 7)
    colnames(xeno) <- models_xeno
    colnames(pdo) <- substr(colnames(pdo), 0, 7)

    # selection of the good one only (validated PDOs)
    models_xeno <- models_xeno[models_xeno %in% keep$smodel]

    sink(snakemake@log[[1]], append=TRUE)
    print(paste0('After validation filtering we kept ', length(models_xeno), ' starting from x= ', ncol(xeno), 'and o= ', ncol(pdo)))
    sink()

    xeno <- xeno[, models_xeno[order(models_xeno)]]
    pdo <- pdo[, models_xeno[order(models_xeno)]]

    afx <- unlist(xeno[xeno!=0])
    afo <- unlist(pdo[pdo!=0])
    pd <- data.frame(af=c(afo, afx), class=c(rep('pdo', length(afo)),rep('pdx',length(afx))))
    #ggplot(data=pd, aes(x=af, color=class, y=..count..)) + geom_density()+scale_x_continuous(breaks=(seq(0, 1, by=0.05)))+unmute_theme
    ggplot(data=pd, aes(x=af, fill=class)) + geom_histogram(alpha=0.5, position='dodge')+
    scale_x_continuous(breaks=(seq(0, 1, by=0.05)))+unmute_theme+
    scale_fill_manual(values=c('darkblue', 'firebrick1'))
    ggsave(out)
    #ggplot(data=pd, aes(x=af, color=class)) + geom_density()+scale_x_continuous(breaks=(seq(0, 1, by=0.05)))+unmute_theme
    return(list(x=xeno, o=pdo))
}

af_all <- filter_plot(pdo_df, xeno_df, keep, snakemake@output[['freqAll']])
##
xeno_df_w <- read.table(xeno_waf, header=TRUE, sep="\t", row.names=1)
pdo_df_w <- read.table(pdo_waf, header=TRUE, sep="\t", row.names=1)

af_tiers <- filter_plot(pdo_df_w, xeno_df_w, keep, snakemake@output[['freqTiers']])

plot_n <- function(pdo, pdx, out) {
    pdxbin <- ifelse(pdx==0,0, 1)
    n_muts_x <- colSums(pdxbin)
    pdobin <- ifelse(pdo==0,0, 1)
    n_muts_o <- colSums(pdobin)
    #hist(n_muts_x, breaks=20)
    #hist(n_muts_o, breaks=20)
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