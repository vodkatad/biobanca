#!/usr/bin/env Rscript
library(ggplot2)

g1_f <- snakemake@input[['g1']]
g2_f <- snakemake@input[['g2']]
Rimage_f <- snakemake@input[['Rimage']]
log_f <- snakemake@log[['log']]
out_plot_f <- snakemake@output[['plot']]
avg_f <- snakemake@output[['avg']]

load(Rimage_f)

data_1 <- read.table(g1_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
data_2 <- read.table(g2_f, sep="\t", header=TRUE, stringsAsFactors=FALSE)

merged_guides <- merge(data_1, data_2, by="smodel")

sink(log_f)
ci <- cor.test(merged_guides$EGFR_ds_1, merged_guides$EGFR_ds_2)
ci
ci$p.value
sink()

ggplot(data=merged_guides, aes(x=EGFR_ds_1, y=EGFR_ds_2)) + 
      geom_point(size=0.1) +
      stat_smooth(method = "lm", col = "darkgrey", size=0.2) +
          unmute_theme +
          xlab('EGFR guide n.1') +
          ylab('EGFR guide n.2')
          #theme(legend.position = "none"))

ggsave(out_plot_f, height=2.5, width=2.5, units='in')

res <- data.frame(smodel=merged_guides$smodel, ko_score=rowMeans(merged_guides[,c(2,3)]))

write.table(res, file=avg_f, quote=FALSE, sep="\t", row.names=FALSE)