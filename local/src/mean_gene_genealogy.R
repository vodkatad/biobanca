### Script to calculate median value between model 

library(tidyverse)

### load files
### vsd_expr_f <- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/replicati_genealogy/small_vsd.tsv.gz"
### meta_f <- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/replicati_genealogy/selected_metadata_annot_final_nolinfo_nooutlier"
### rep_gen_f <- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/replicati_genealogy/LMO_BASALE_replicates.tsv"
rep_gen_f <- snakemake@input[["replicates"]]
meta_f <- snakemake@input[["metadata"]]
vsd_expr_f <- snakemake@input[["expr"]]
mean_gene_genealogy <- snakemake@output[["mean_btw_genealogy"]]
vsd_expr_df <- read.table(gzfile(vsd_expr_f), quote = "", sep = "\t", header=TRUE)
meta_df<-read.table(meta_f, quote = "", sep = "\t", header = TRUE)
rep_gen_df <- read.table(rep_gen_f, quote = "", sep = "\t", header=TRUE)
rownames(rep_gen_df) <- rep_gen_df$model
rep_gen_df$model <- NULL
rep_gen_df$genealogy <- gsub("-2", ".2", rep_gen_df$genealogy, fixed = TRUE)

### Styling vsd_expr_df
vsd_expr_df <- as.data.frame(t(vsd_expr_df))
vsd_expr_df <- tibble::rownames_to_column(vsd_expr_df)
colnames(vsd_expr_df)[1] <- "sample_id_R"

### Filter for LMO_BASALE
class <- snakemake@wildcards[["sclass"]]
meta_df_lmo_b <- filter(meta_df, grepl(class, type))
meta_df_lmo_b$sample_id_R <- gsub("-2", ".2", meta_df_lmo_b$sample_id_R, fixed = TRUE)
vsd_expr_df_lmo <- vsd_expr_df[vsd_expr_df$sample_id_R %in% meta_df_lmo_b$sample_id_R, ]

### Calculate the median value
mean_df <- data.frame(matrix(NA, nrow=nrow(rep_gen_df), ncol = ncol(vsd_expr_df_lmo)-1))
rownames(mean_df) <- rownames(rep_gen_df)
colnames(mean_df) <- colnames(vsd_expr_df_lmo)[2:ncol(vsd_expr_df_lmo)]
i=1
for (i in seq(1, nrow(rep_gen_df))) {
  rep_gen_df_split <- strsplit(rep_gen_df[i,"genealogy"], ",")
  rep_gen_df_split_unlist <- unlist(rep_gen_df_split)
  x <- rep_gen_df_split_unlist[1]
  y <- rep_gen_df_split_unlist[2]
  a <- as.numeric((vsd_expr_df_lmo[vsd_expr_df_lmo$sample_id_R==x, ])[-1]) ##metto il -1 perchè la prima colonna sono i genalogy
  b <- as.numeric((vsd_expr_df_lmo[vsd_expr_df_lmo$sample_id_R==y, ])[-1])
  mean_i <- (a+b)/2
  mean_df[i,] <- mean_i
}

write.table(mean_df, file = mean_gene_genealogy, sep='\t', col.names = TRUE, row.names= TRUE, quote = FALSE)




