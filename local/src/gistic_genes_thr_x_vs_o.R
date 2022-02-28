# Compute delta of gistic scores for genes
library(ggplot2)
library(ggrastr)

pdx_f <- snakemake@input[['pdx']]
pdo_f <- snakemake@input[['pdo']]
corrplot_f <- snakemake@output[['corrplot']]
corrs_f <- snakemake@output[['corrs']]
log <- snakemake@log[['log']]

osnakemake <- snakemake
load(snakemake@input[['Rimage']])
#eval(parse(text=myriad))
snakemake <- osnakemake

df_x <- read.table(pdx_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
df_o <- read.table(pdo_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
rownames(df_x) <- df_x$Gene.Symbol
rownames(df_o) <- df_o$Gene.Symbol
df_x$Locus.ID <- NULL
df_x$Gene.Symbol <- NULL
df_x$Cytoband <- NULL
df_o$Locus.ID <- NULL
df_o$Gene.Symbol <- NULL
df_o$Cytoband <- NULL

#save.image('pippo.Rdata')

# TODO FIXME
get_freq <- function(pdx, pdo, direction, kind) {
  if (direction == 1) {
    freq_x <- as.data.frame(apply(pdx, 1, function(x) { sum(x >= direction)/length(x) }))
    freq_o <- as.data.frame(apply(pdo, 1, function(x) { sum(x >= direction)/length(x) }))
  } else {
    freq_x <- as.data.frame(apply(pdx, 1, function(x) { sum(x <= direction)/length(x) }))
    freq_o <- as.data.frame(apply(pdo, 1, function(x) { sum(x <= direction)/length(x) }))
  }
  colnames(freq_x) <- kind
  colnames(freq_o) <- 'PDO'
  res <- merge(freq_x, freq_o, by='row.names')
  sink(log)
  print(nrow(res))
  print(nrow(freq_x))
  print(nrow(freq_o))
  sink()
  rownames(res) <- res$Row.names
  res$Row.names <- NULL
  return(res)
}

freqs_amp <- get_freq(df_x, df_o, 1, 'PDX')
freqs_del <- get_freq(df_x, df_o, -1, 'PDX')
freqs_amp$event <- 'amp'
freqs_del$event <- 'del'
freqs <- rbind(freqs_amp, freqs_del)

plot <- ggplot(data=freqs, aes_string(x='PDX', y='PDO', color='event'))+rasterise(geom_point(alpha=0.5, size=0.1), dpi=300)+geom_smooth(method='lm', size=0.2)+
        scale_color_manual(values=c('#c84440','#185492'))

save.image('pippo.Rdata')
plotbis <- function(plot, theme_unmute, theme_mute, name, h=2.5, w=2.5, units='in') {
  unmute <- plot + theme_unmute
  ggsave(filename=name, plot=unmute, height=h, width=w, units=units)
  ext <- substr(name, nchar(name)-3, nchar(name)) 
  # this works only with 3 char extensions, TODO FIXME  https://stackoverflow.com/questions/29113973/get-filename-without-extension-in-r
  name_mute <- paste0(substr(name, 0, nchar(name)-3), 'mute', ext)
  mute <- plot + theme_mute
  ggsave(filename=name_mute, plot=mute, height=h, width=w, units=units)
}

#plotbis(plot=plot, theme_unmute=unmute_theme, theme_mute=mute_theme, name=corrplot_f)
plotbis(plot=plot, theme_unmute=unmute_theme, theme_mute=mute_theme, name=corrplot_f, h=8, w=8)

cor_amp <- cor.test(freqs_amp$PDX, freqs_amp$PDO)
cor_del <- cor.test(freqs_del$PDX, freqs_del$PDO)

res <- data.frame(event=c('amp','del'), pearson=c(cor_amp$estimate, cor_del$estimate), pval=c(cor_amp$p.value, cor_del$p.value))
write.table(res, file=corrs_f, sep='\t', row.names=FALSE, quote=FALSE)