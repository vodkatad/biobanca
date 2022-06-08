data <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/ls')
odata <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv')
pdata <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDX_S/treat_cutoff0.05-cetuxi.vs.NT.deseq2.tsv')
head(odata)
m <- merge(odata, pdata)
m <- merge(odata, pdata, by="row.names")
dim(m)
dim(odata)
dim(pdata)
library(ggplot2)
load('/scratch/trcanmed/biobanca/dataset/V1/theme_5.Rdata')
ggplot(data=m, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj.x)))+
  geom_point(size=0.1)+theme_bw()+xlab('logFC Organoids')+ylab('logFC Xenografts')+
  scale_color_distiller(palette = "YlOrRd", direction=1)+unmute_theme+theme(legend.position = "none")
#ggsave("~/test_biobanca_cor_ctx.svg", width=55, height=55, units="mm")
ggsave("~/test_biobanca_cor_ctx.pdf", width=55, height=55, units="mm")

ggplot(data=m, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj.x)))+
  geom_point(size=0.1)+theme_bw()+xlab('logFC Organoids')+ylab('logFC Xenografts')+
  scale_color_distiller(palette = "YlOrRd", direction=1)+unmute_theme
#ggsave("~/test_biobanca_cor_ctx.svg", width=55, height=55, units="mm")
ggsave("~/test_biobanca_cor_ctx_legend.pdf", width=55, height=55, units="mm")



ggplot(data=m, aes(x=log2FoldChange.x,y=log2FoldChange.y, color=-log10(padj.x)))+
  geom_point()+theme_bw()+xlab('logFC Organoids')+ylab('logFC Xenografts')+
  scale_color_distiller(palette = "YlOrRd", direction=1)