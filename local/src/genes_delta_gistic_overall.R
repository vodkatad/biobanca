# Compute delta of gistic scores for genes
library(pheatmap)

#pdo_gistic_f <- '/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/gistic/gistic_pdo/all_data_by_genes.txt'
#pdx_gistic_f <- '/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/gistic/gistic_xeno/all_data_by_genes.txt'
pdo_gistic_f <- snakemake@input[['pdo']]
pdx_gistic_f <- snakemake@input[['pdx']]
delta_f <- snakemake@output[['delta']]
plot_f <- snakemake@output[['delta_plot']]

load_prepare_df <- function(filename) {
  df <- read.table(filename, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c('gene','gistic')
  df
}

# right now we average rows with repeated gs
# TODO look up in GISTIC what's this
average_repe <- function(df) {
  freqs <- as.data.frame(table(df$gene))
  reps <- freqs[freqs$Freq > 1, 'Var1']
  uni <- df[!df$gene %in% reps,]
  repi <- df[df$gene %in% reps,]
  averaged <- as.data.frame(sapply(reps, function(x) { d <- df[df$gene==x,, drop=FALSE]; colMeans(d[,-1, drop=FALSE]) } ))
  #res <- rbind(averaged, uni[,-1, drop=F])
  #rownames(res) <- c(as.character(reps), uni$gene)
  res <- data.frame(row.names=c(as.character(reps), uni$gene), gistic=c(averaged[,1], uni$gistic))
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
#w_genes <- intersect(rownames(pdx), rownames(pdo))
#pdo <- pdo[w_genes,, drop=FALSE]
#pdx <- pdx[w_genes,, drop=FALSE]

m <- merge(pdo, pdx, by="row.names")
colnames(m) <- c('gene_symbol','pdo_gistic_score','pdx_gistic_score')
m$delta <- m$pdo_gistic_score - m$pdx_gistic_score

write.table(m, file=delta_f, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
pdf(plot_f)
hist(m$delta, breaks=100, main='PDO-PDX delta ggistic scores')
#lines(density(m$delta))
graphics.off()