library(tidyverse)

csv_f <- snakemake@input[["cor"]]
res <- snakemake@output[["result"]]
log_f <- snakemake@log[['log']]


#csv_f <- "/scratch/trcanmed/biobanca/dataset/V1/LOH_earlylate/depth_ratio/heatmap_data/seg_dim_100000.csv"
csv <- read.csv(csv_f)
rownames(csv) <- csv$X
csv$X <- NULL
cor <- as.matrix(csv)

diagonal_values <- diag(cor)
mean_diagonal <- mean(diagonal_values)
median_diagonal <- median(diagonal_values)

non_diagonal_values <- cor[lower.tri(cor) | upper.tri(cor)]
mean_non_diagonal <- mean(non_diagonal_values)
median_non_diagonal <- median(non_diagonal_values)

diag <- as.data.frame(matrix(ncol = 2, nrow = 2))
colnames(diag) <- c("mean", "median")
rownames(diag) <- c("diag", "non_diag")

diag["diag","mean"] <- mean_diagonal
diag["diag","median"] <- median_diagonal
diag["non_diag","mean"] <- mean_non_diagonal
diag["non_diag","median"] <- median_non_diagonal

sink(log_f, append=TRUE)
w <- wilcox.test(diagonal_values, non_diagonal_values)
w
sink()

write.table(diag, res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)