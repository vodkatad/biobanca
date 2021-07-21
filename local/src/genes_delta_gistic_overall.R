# Compute delta of gistic scores for genes
library(pheatmap)

library(ggplot2)
th <- function() {
  textSize <- 1.5
  theme <- theme_bw() +
    theme(
      strip.text = element_text(size = rel(textSize)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title = element_text(size = rel(1.8)),
      axis.text.x = element_text(size=rel(1.7)),
      axis.text.y = element_text(angle = 0,
                                 size = rel(1.7)),
      axis.line = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.length.y.left = unit(3,'mm'),
      
      legend.position = "top",
      legend.justification = "right",
      #legend.margin = margin(unit(0, "cm")),
      legend.title = element_text(size = rel(textSize), face = "bold"),
      legend.text = element_text(size = rel(1.2)),
      legend.background = element_rect(size=0.5, linetype="solid", color="black"),
      plot.title = element_text(
        face = "bold",
        size = rel(2),
        hjust = 0.5
      ),
      panel.border = element_blank(),
      plot.caption = element_text(size=rel(1))
    )
  theme
}
current_theme <- th()
#pdo_gistic_f <- 'genes_scores_pdo.tsv'
#pdx_gistic_f <- 'genes_scores_xeno.tsv'
pdo_gistic_f <- snakemake@input[['pdo']]
pdx_gistic_f <- snakemake@input[['pdx']]
delta_f <- snakemake@output[['delta']]
plot_f <- snakemake@output[['delta_plot']]
plot2_f <- snakemake@output[['scores_plot']]


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



gglm <- function(x, y, nx, ny, title, file) {
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y)
  ggplot(d, aes(x=x, y=y)) +geom_point(size=2)+current_theme+xlab(nx)+ylab(ny)+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle(title)+geom_smooth(method=lm)+current_theme
  ggsave(file, width=10, height=10, units="in")
}


gglm(m$pdo_gistic_score, m$pdx_gistic_score, 'pdo_ggistic','pdx_ggistic','Gene Gistic Scores', plot2_f)