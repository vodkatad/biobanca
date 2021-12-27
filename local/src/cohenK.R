#!/usr/bin/env Rscript

set.seed(42)

library(tidyverse)
library(psych)

df <- snakemake@input[[1]]
output <- snakemake@output[[1]]


df <- read.table(gzfile(df), sep='\t', quote="", header=TRUE, row.names=1)

# df[,c(3,4,6:ncol(df))] <- NULL

res <- as.data.frame(cohen.kappa(df)$kappa)

write.table(res, output, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
