library(ggplot2)
library(pheatmap)
#xeno_f <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/xeno_qdnaseq/cn_log2.tsv'
#pdo_f <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/pdo_qdnaseq/cn_log2.tsv'

xeno_f <- snakemake@input[['xeno']]
pdo_f <- snakemake@input[['pdo']]
expected_n <- snakemake@params[['expected_pairs']] # 142
heatmap_f <- snakemake@output[['heatmap']]
density_f <- snakemake@output[['density']]
pearson_f <- snakemake@output[['pearson']]
rdata_f <- snakemake@output[['rdata']]
which <- snakemake@wildcards[['which']]
good_f <- snakemake@input[['good']]
save.image(rdata_f)


## for ggplot
osnakemake <- snakemake
load(snakemake@input[['Rimage']])
#eval(parse(text=myriad))
snakemake <- osnakemake

xeno_df <- read.table(xeno_f, quote="", sep="\t", header=TRUE, row.names = 1)
pdo_df <- read.table(pdo_f, quote="", sep="\t", header=TRUE, row.names = 1)
xeno_df[,c("chromosome","start","end")] <- NULL
pdo_df[,c("chromosome","start","end")] <- NULL

# We start keeping the EGFR wt PDO and the mutated organoid.  EGFR wt (the other one is mutant) with a different 'model'/short genealogy
#egrassi@godot:/scratch/trcanmed/biobanca/dataset/V1/shallowseq$ head -n1  xeno_qdnaseq/cn_segm.tsv  | tr "\t" "\n" | grep CRC0177
#CRC0177LMX0B05001TUMD06000
#CRC0177LMX0B07002TUMD04000
#egrassi@godot:/scratch/trcanmed/biobanca/dataset/V1/shallowseq$ head -n1  pdo_qdnaseq/cn_segm.tsv  | tr "\t" "\n" | grep CRC0177
#CRC0177LMO0A04008002D02000
#CRC0177LMO0D04017001D02000
# old code to keep 'right' matches
#colnames(pdo_df)[colnames(pdo_df)=="CRC0177LMO0A04008002D02000"] <- 'CRCE177LMO0A04008002D02000'
#colnames(xeno_df)[colnames(xeno_df)=="CRC0177LMX0B05001TUMD06000"] <- 'CRCE177LMX0B05001TUMD06000'
pdo_df <- pdo_df[, colnames(pdo_df) != "CRC0177LMO0D04017001D02000"] # remove the mutated pdo
xeno_df <- xeno_df[,colnames(xeno_df) != "CRC0177LMX0B05001TUMD06000"] # remove the wt xeno

# list of validated models for PDO, the one we want to keep
good_df <- read.table(good_f, sep="\t", stringsAsFactors=FALSE, header=TRUE)
keep <- good_df[good_df$buoni, , drop=FALSE]

list_remove_xeno <- c("CRC1870", "CRC1875","CRC2041")
list_lmh <- colnames(pdo_df)[grepl('LMH', colnames(pdo_df))]
model_lmh <- substr(list_lmh, 0,7)
if (which == "xo") {
  expected_n <- expected_n - 1  # removal of CRC0177 brings us here
  pdo_df <- pdo_df[, !colnames(pdo_df) %in% list_lmh]
  pdo_df <- pdo_df[, !substr(colnames(pdo_df),0,7) %in% model_lmh]
  xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% model_lmh]

  #> setdiff(models_xeno, models_pdo)
  #[1] "CRC1870" "CRC1875" "CRC2041"
  # We know we were missing some PDOs TODO load from file if qc changes in the future
  xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% list_remove_xeno]
  # TODO? make modular: first step to prepare matrix then general script for correlations
} else if (which == "xh") {
  pdo_df <- pdo_df[, colnames(pdo_df) %in% list_lmh]
  xeno_df <- xeno_df[, substr(colnames(xeno_df),0,7) %in% model_lmh]
} else if (which == "oh") {
  xeno_df <- pdo_df[, colnames(pdo_df) %in% list_lmh] # xeno now impersonate human samples
  pdo_df <- pdo_df[, !colnames(pdo_df) %in% list_lmh]
  pdo_df <- pdo_df[, substr(colnames(pdo_df),0,7) %in% model_lmh]
} else {
  stop(paste0("Dunno which type you requested! ", which))
}

models_xeno <- substr(colnames(xeno_df), 0, 7)
models_pdo <- substr(colnames(pdo_df), 0, 7)
if (length(intersect(models_xeno, models_pdo)) != expected_n) {
    stop('Wrong number of corresponding models!')
}

#xeno_df <- xeno_df[, order(models_xeno)] # eeeeh
#pdo_df <- pdo_df[, order(models_xeno)]

colnames(xeno_df) <- substr(colnames(xeno_df), 0, 7)
colnames(pdo_df) <- substr(colnames(pdo_df), 0, 7)

# selection of the good one only (validated PDOs)
models_xeno <- models_xeno[models_xeno %in% keep$smodel]

sink(snakemake@log[[1]])
print(paste0('After validation filtering we kept ', length(models_xeno), ' starting from ', expected_n))
sink()

xeno_df <- xeno_df[, models_xeno[order(models_xeno)]]
pdo_df <- pdo_df[, models_xeno[order(models_xeno)]]


if (!all(substr(colnames(xeno_df),0,7) == substr(colnames(pdo_df), 0,7))) {
  stop('Piciapirilla!')
}


xeno_df_bck <- xeno_df
pdo_df_bck <- pdo_df
xeno_df <- log2(xeno_df+1)
pdo_df <- log2(pdo_df+1)
pearson <- cor(xeno_df, pdo_df)

svg(heatmap_f)
pheatmap(pearson, cluster_rows=F , cluster_cols=F, labels_col="PDO", labels_row="PDX", fontsize.number=1.5)
graphics.off()
diag <- diag(pearson)

pearson2 <- pearson
diag(pearson2) <- rep(NA, length(diag))
all <- as.numeric(unlist(pearson2))
all <- all[!is.na(all)]
#all <- upper.tri(pearson, diag = FALSE) # this is not a simmetric matrix!
pdata <- data.frame(pearson=c(all, diag), type=c(rep('unmatched', length(all)),rep('matched', length(diag))))
#ggplot(data=pdata, aes(x=pearson, color=type))+geom_density()+unmute_theme_myriad
ggplot(data=pdata, aes(x=pearson, color=type))+geom_density(size=1.5)+unmute_theme
ggsave(density_f, height=31.7, width=31.7, units='cm')

ext <- substr(density_f, nchar(density_f)-3, nchar(density_f)) 
# this works only with 3 char extensions, TODO FIXME  https://stackoverflow.com/questions/29113973/get-filename-without-extension-in-r
density_f_mute <- paste0(substr(density_f, 0, nchar(density_f)-3), 'mute', ext)
print(density_f_mute)
ggplot(data=pdata, aes(x=pearson, color=type))+geom_density(size=1.5)+mute_theme
ggsave(density_f_mute, height=31.7, width=31.7, units='cm')

write.table(as.data.frame(pearson), file=pearson_f, quote=FALSE, sep="\t")
save.image(rdata_f)