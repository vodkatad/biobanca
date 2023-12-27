#% alterazioni geni/casi

load("/scratch/trcanmed/biobanca/dataset/V1/targeted/oncoprint_0.05.RData")

rowSums(both)
diff_mat <- (td2o-both) + (td2x-both)
rowSums(diff_mat)

both_df <- as.data.frame(rowSums(both))
colnames(both_df) <- "n_mut_comuni"
both_df$genes <- row.names(both_df)

diff_mat_sum <- as.data.frame(rowSums(diff_mat))
colnames(diff_mat_sum) <- "n_mut_singole"
diff_mat_sum$genes <- rownames(diff_mat_sum)

merged <- merge(both_df, diff_mat_sum, by = "genes")
rownames(merged) <- merged$genes
merged$genes <- NULL
merged$percentage <- (merged$n_mut_singole/(merged$n_mut_comuni+merged$n_mut_singole))*100

merged <- merged[order(-merged$percentage),]
write.table(merged, file="/home/mferri/diff_mutation_genes.tsv", quote = FALSE, sep = "\t", col.names = TRUE)

colSums(both)
diff_mat_col <- (td2o-both) + (td2x-both)
colSums(diff_mat_col)

both_df_col <- as.data.frame(colSums(both))
colnames(both_df_col) <- "n_mut_comuni"
both_df_col$cases <- row.names(both_df_col)

diff_mat_sum_col <- as.data.frame(colSums(diff_mat_col))
colnames(diff_mat_sum_col) <- "n_mut_singole"
diff_mat_sum_col$cases <- rownames(diff_mat_sum_col)

merged_col <- merge(both_df_col, diff_mat_sum_col, by = "cases")
rownames(merged_col) <- merged_col$cases
merged_col$cases <- NULL
merged_col$percentage <- (merged_col$n_mut_singole/(merged_col$n_mut_comuni+merged_col$n_mut_singole))*100

merged_col <- merged_col[order(-merged_col$percentage),]
write.table(merged_col, file="/home/mferri/diff_mutation_cases.tsv", quote = FALSE, sep = "\t", col.names = TRUE)
