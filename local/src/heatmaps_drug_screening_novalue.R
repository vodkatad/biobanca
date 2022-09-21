library(reshape)
library(pheatmap)

max_inh_f <- snakemake@input[['max_inh']]
anova_f <- snakemake@input[['anova']]
heatmap_f <- snakemake@output[['heatmap']]

#d <- read.table('/scratch/trcanmed/biobanca/dataset/V1/drug_screening/max_drug_value.tsv', sep="\t", header=TRUE)
d <- read.table(max_inh_f, sep="\t", header=TRUE)
dm <- cast(d, formula=DRUG ~ MODEL, value = "MAX_VALUE", add.missing=TRUE, fill=NA)
row.names(dm) <- dm$DRUG
dm$DRUG <- NULL
pheatmap(dm, cluster_cols=F, cluster_rows=F)

#da <- read.table('/scratch/trcanmed/biobanca/dataset/V1/drug_screening/anova.tsv', sep="\t", header=TRUE)
da <- read.table(anova_f, sep="\t", header=TRUE)

pheatmap(dm, cluster_cols=F, cluster_rows=F, fontsize=5, height=4, width=4, file=heatmap_f)