#!/usr/bin/env Rscript
library(ggplot2)

cetuxi_f <- snakemake@input[['cetuxi']]
crispr_f <- snakemake@input[['crispr']]
annot_f <- snakemake@input[['annot']]
Rimage_f <- snakemake@input[['Rimage']]
log_f <- snakemake@log[['log']]
out_plot_f <- snakemake@output[['plot']]
out_plotlegend_f <- snakemake@output[['plot_legend']]
col_cetuxi <- snakemake@wildcards[['CTG']]
out_barplot_f <- snakemake@output[['barplot_legend']]
out_barplotnol_f <- snakemake@output[['barplot']]

load(Rimage_f)

crispr <- read.table(crispr_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cetuxi <- read.table(cetuxi_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(cetuxi)[1] <- 'smodel'
annot <- read.table(annot_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)

save.image('meh.Rdata')
merge_data <- merge(cetuxi, crispr, by="smodel")
merged_annot <- merge(merge_data, annot, by="smodel")

sink(log_f)
print('WARNING if rerun with more data from MarcoA check which non validated/double failed PDOs more than 1257 will be dragged in. We want 13 top right now:')
print('cetuxi:')
nrow(cetuxi)
print('crispr:')
nrow(crispr)
print('both:')
nrow(merge_data)
print('both annot:')
nrow(merged_annot)
ci <- cor.test(merged_annot[, col_cetuxi], merged_annot$ko_score)
ci
ci$p.value
sink()

gplot <- ggplot(data=merged_annot, aes_string(x=col_cetuxi, y='ko_score')) + 
      geom_point(size=0.1, aes(color=alterations)) +
      stat_smooth(method = "lm", col = "darkgrey", size=0.2) +
          unmute_theme +
          xlab('cetuximab viability') +
          ylab('EGFR avg ko score') +
          scale_colour_manual(values=c('darkorange', 'darkred', 'grey40')) +
          labs(color="Relevant somatic alterations")

ggsave(out_plotlegend_f, plot=gplot, height=2.5, width=2.5, units='in')

gplot_nol <- gplot + theme(legend.position = "none")

ggsave(out_plot_f, plot=gplot_nol, height=2.5, width=2.5, units='in')

merged_annot$sort <- merged_annot[,col_cetuxi]

gplot <- ggplot(data=merged_annot, aes(x=reorder(smodel, sort), y=ko_score, fill=alterations))+geom_col()+
  unmute_theme +
  xlab('cetuximab viability') +
  ylab('EGFR avg ko score') +
  scale_fill_manual(values=c('darkorange', 'darkred', 'grey40')) +
  labs(color="Relevant somatic alterations") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(out_barplot_f, plot=gplot, height=2.5, width=2.5, units='in')

gplot_nol <- gplot + theme(legend.position = "none")
ggsave(out_barplotnol_f, plot=gplot_nol, height=2.5, width=2.5, units='in')