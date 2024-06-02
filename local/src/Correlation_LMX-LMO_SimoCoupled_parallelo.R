#!/usr/bin/env Rscript

library(getopt)
library(tidyverse)
library(tibble)

set.seed(42)

 opts <- matrix(c(
   'help', 'h', 0, 'logical',
   'simo', 's', 1, 'character',
   'vsd', 'v', 1, 'character',
   'meda', 'm', 1, 'character',
   'genes', 'g', 1, 'character',
   'right', 'r', 1, 'character',
   'output', 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

options(warn=1) # to print warnings as they occurr


meda <- read.table(gzfile(opt$meda), sep='\t', quote="", header=TRUE)
meda$sample_id_R <- gsub('-', '.', meda$sample_id_R, fixed = TRUE)
right <- read.table(gzfile(opt$right), sep='\t', quote="", header=FALSE, stringsAsFactors = FALSE)
### carica metadati e TIENI i basali
meda <- meda[meda$type %in% right$V1,]
simo <- read.table(gzfile(opt$simo), sep='\t', quote="", header=TRUE)

vsd <- read.table(gzfile(opt$vsd), sep='\t', quote="", row.names=1, header=TRUE)
names(vsd) <- gsub('-', '.', names(vsd), fixed = TRUE)
rownames(vsd) <- str_remove(rownames(vsd), 'H_')

d <- vsd[,names(vsd) %in% meda$sample_id_R]

#save.image("pippo83.Rdata")
### filtering expression data: we want high sd genes but not clear outliers / not expressed genes
### filter not expressed genes
means <- apply(d, 1, mean)
med <- median(means)

de <- d[means > med,]

sds <- apply(de, 1, sd)
### there are some very high sd?
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
### 0.1086  0.3540  0.4702  0.5313  0.6411  3.7110 
# noisy_genes_thr <- quantile(sds, 0.9) # tenere o no? non � una differenza cos� enorme come prima...toglilo

# de <- de[sds < noisy_genes_thr,]
# sds <- sds[sds < noisy_genes_thr]

### now we keep thet top 10% variable genes 
sds <- sds[order(-sds)]

n <- length(sds)
keep <- head(sds, round(0.10*n)) #prova con e senza
keep_genes <- names(keep)
desd <- de[rownames(de) %in% keep_genes,]

print(n)
print(length(keep))
write.table(as.data.frame(keep_genes), file=opt$genes, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# desd1 <- desd
# names(desd1) <- substr(names(desd1), 1, 12)
# lmo <- desd1[grepl('LMO', names(desd1),, fixed=TRUE)]
# lmx <- desd1[grepl('LMX', names(desd1),, fixed=TRUE)]
# lmo <- lmo[,names(lmo) %in% simo$LMO_lineage]
# lmx <- lmx[,names(lmx) %in% simo$LMX_lineage]

desd2 <- desd
lmo2 <- desd2[grepl("LMO", names(desd2), fixed = TRUE)]
#names(lmo2) <- gsub('-', '.', names(lmo2), fixed = TRUE)
lmx2 <- desd2[grepl("LMX", names(desd2), fixed = TRUE)]
lmo2 <- lmo2[,names(lmo2) %in% simo$PDO_long_lineage]
lmx2 <- lmx2[,names(lmx2) %in% simo$PDX_long_lineage]

# TODO check so we do not have replicates or they are sorted 'correctly' here?
#> length(simo$LMX_lineage)
#[1] 55
#> length(simo$LMO_lineage)
#[1] 55
#> length(unique(substr(simo$LMO_lineage,0,7)))
#[1] 55
#> length(unique(substr(simo$LMX_lineage,0,7)))
#[1] 55
# Since we do not have replicates at the smodel level in Simo's right pairs we
# can safely substr right now and know that we do not have repetitions that will
# mess things up.
# colnames(lmo) <- substr(colnames(lmo),0,7)
# colnames(lmx) <- substr(colnames(lmx),0,7)
# clmo <- lmo[, names(lmx)]
# clmx <- lmx[, names(clmo)]
# ### check also for genes
# all_genes <- intersect(rownames(lmo), rownames(lmx))
# clmo <- clmo[all_genes,]
# clmx <- clmx[all_genes,]
# 
# 
# if (all(colnames(clmo)!=colnames(clmx)) & all(rownames(clmo)!=rownames(clmx))) {
#   stop('Brutto llama!')
# }
# 
# res <- cor(clmo, clmx)

###desd2

colnames(lmo2) <- substr(colnames(lmo2),0,7)
colnames(lmx2) <- substr(colnames(lmx2),0,7)
clmo2 <- lmo2[, names(lmx2)]
clmx2 <- lmx2[, names(clmo2)]
### check also for genes
all_genes2 <- intersect(rownames(lmo2), rownames(lmx2))
clmo2 <- clmo2[all_genes2,]
clmx2 <- clmx2[all_genes2,]

if (all(colnames(clmo2)!=colnames(clmx2)) & all(rownames(clmo2)!=rownames(clmx2))) {
  stop('Brutto llama!')
}

res2 <- cor(clmo2, clmx2)


write.table(res2, gzfile(opt$output), sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
