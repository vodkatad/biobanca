# Compute delta of gistic scores for genes
library(pheatmap)

#pdo_gistic_f <- '/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/gistic/gistic_pdo/all_data_by_genes.txt'
#pdx_gistic_f <- '/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/gistic/gistic_xeno/all_data_by_genes.txt'
pdo_gistic_f <- snakemake@input[['pdo']]
pdx_gistic_f <- snakemake@input[['pdx']]
delta_f <- snakemake@output[['delta']]
heatmap_f <- snakemake@output[['delta_plot']]
expected_n <- snakemake@params[['expected_pairs']] # 142

load_prepare_df <- function(filename) {
  df <- read.table(filename, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  df$Gene.ID <- NULL
  df$Cytoband <- NULL
  df$Gene.Symbol <- sapply(strsplit(df$Gene.Symbol, '|', fixed=TRUE), function(x) {x[[1]][1]})
  #rownames(df) <- df$Gene.ID
  df
}

# right now we average rows with repeated gs
# TODO look up in GISTIC what's this
average_repe <- function(df) {
  freqs <- as.data.frame(table(df$Gene.Symbol))
  reps <- freqs[freqs$Freq > 1, 'Var1']
  uni <- df[!df$Gene.Symbol %in% reps,]
  repi <- df[df$Gene.Symbol %in% reps,]
  averaged <- t(sapply(reps, function(x) { d <- df[df$Gene.Symbol==x,]; colMeans(d[,-1]) } ))
  res <- rbind(averaged, uni[,-1])
  rownames(res) <- c(as.character(reps), uni$Gene.Symbol)
  res
}

pdo <- load_prepare_df(pdo_gistic_f)
pdx <- load_prepare_df(pdx_gistic_f)
pdo <- average_repe(pdo)
pdx <- average_repe(pdx)

# same order for rows/genes for both df
# all(rownames(pdo)==rownames(pdx))
# FALSE
# but intersect is ok...
nx <- nrow(pdx)
pdx <- pdx[rownames(pdo),]
nx2 <- nrow(pdx)
if (nx != nx2 | !all(rownames(pdo)==rownames(pdx))) {
  stop('Something off in gistic genes')
}


# we do only xo comparison now
colnames(pdo)[colnames(pdo)=="CRC0177LMO0A04008002D02000"] <- 'CRCE177LMO0A04008002D02000'
colnames(pdx)[colnames(pdx)=="CRC0177LMX0B05001TUMD06000"] <- 'CRCE177LMX0B05001TUMD06000'

list_remove_xeno <- c("CRC1870", "CRC1875","CRC2041")
list_lmh <- colnames(pdo)[grepl('LMH', colnames(pdo))]
model_lmh <- substr(list_lmh, 0,7)
pdo <- pdo[, !colnames(pdo) %in% list_lmh]
pdo <- pdo[, !substr(colnames(pdo),0,7) %in% model_lmh]
pdx <- pdx[, !substr(colnames(pdx),0,7) %in% model_lmh]
  
#[1] "CRC1870" "CRC1875" "CRC2041"
# We know we were missing some PDOs TODO load from file if qc changes in the future
pdx <- pdx[, !substr(colnames(pdx),0,7) %in% list_remove_xeno]

colnames(pdo) <- substr(colnames(pdo), 0, 7)
colnames(pdx) <- substr(colnames(pdx), 0, 7)
models_pdo <- colnames(pdo)  
models_pdx <- colnames(pdx)  

if (length(intersect(models_pdx, models_pdo)) != expected_n) {
  stop('Wrong number of corresponding models!')
}

delta <- pdo - pdx
write.table(delta, file=delta_f, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
pdf(heatmap_f)
pheatmap(delta, show_colnames = FALSE, show_rownames = FALSE)
graphics.off()