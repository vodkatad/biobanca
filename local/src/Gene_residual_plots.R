### gene residual plots

library(tidyverse)

gene_residual_f <- snakemake@input[["residuals"]]
gene_res<- snakemake@output[["gene_plot"]]
gene_res_zoom <- snakemake@output [["gene_plot_zoom"]]
gene_residual_frequency <- snakemake@output[["gene_res_freq"]]

gene_residual_df <- read.table(gene_residual_f, quote = "", sep = "\t", header = FALSE)
colnames(gene_residual_df) <- c("gene")

gene_rep <- as.data.frame(table(gene_residual_df))
colnames(gene_rep) <- c("gene", "Freq")
write.table(gene_rep, file=gene_residual_frequency, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
gene_rep <- gene_rep[order(gene_rep$Freq), ]

ggplot(data = gene_rep, aes(x=reorder(gene,Freq), y=Freq)) + geom_col() + theme(axis.text.x = element_blank())
ggsave(gene_res)

##zoom top20

gene_rep_zoom <- as.data.frame(gene_rep %>% filter(Freq > 20))
ggplot(data = gene_rep_zoom, aes(x=reorder(gene,Freq), y=Freq)) + geom_col() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 7))
ggsave(gene_res_zoom)