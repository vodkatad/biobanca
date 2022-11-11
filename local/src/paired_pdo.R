library(pheatmap)

setwd('/scratch/trcanmed/biobanca/dataset/V1/targeted_all/')
pairs <- read.table('/scratch/trcanmed/biobanca/local/share/data/paired.tsv', sep='\t', header=TRUE, row.names=1)
load('preprocGeneAF_0.05.Rdata')

subset <- pdobing[,colnames(pdobing) %in% rownames(pairs)]

dim(subset)
dim(pairs)
subset <- ifelse(subset, 1, 0)

pheatmap(subset, annotation_col = pairs)