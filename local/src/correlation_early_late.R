library(tidyverse)

early_f <- snakemake@input[["e"]]
late_f <- snakemake@input[["l"]]
res <- snakemake@output[["result"]]

#early_f <- "/scratch/trcanmed/biobanca/dataset/V1/LOH_earlylate/depth_ratio/heatmap_data/seg_dim_100000_early.csv"
early <- read.csv(early_f)
rownames(early) <- early$segment_id
early$segment_id <- NULL

#late_f <- "/scratch/trcanmed/biobanca/dataset/V1/LOH_earlylate/depth_ratio/heatmap_data/seg_dim_100000_late.csv"
late <- read.csv(late_f)
rownames(late) <- late$segment_id
late$segment_id <- NULL

#x <- cor(early, late)
cor <- as.data.frame(cor(early, late))

write.csv(cor, file = res, quote = FALSE, col.names = TRUE, row.names = TRUE)