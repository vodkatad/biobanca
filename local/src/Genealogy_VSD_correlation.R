### Correlation of genealogy and expression from vsd files 

library(tidyverse)

### load input files 
### rep_gen_f <- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/replicati_genealogy/LMO_BASALE_replicates.tsv"
### vsd_expr_f <- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/replicati_genealogy/small_vsd.tsv.gz"
### meta_f <- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/replicati_genealogy/selected_metadata_annot_final_nolinfo_nooutlier"
rep_gen_f <- snakemake@input[["replicates"]]
meta_f <- snakemake@input[["metadata"]]
vsd_expr_f <- snakemake@input[["expr"]]
summary_f <- snakemake@output[["summary"]]
genealogy_VSD_correlation <- snakemake@output[["replicates_correlation"]]

save.image("repliche.Rdata")
rep_gen_df <- read.table(rep_gen_f, quote = "", sep = "\t", header=TRUE)
rownames(rep_gen_df) <- rep_gen_df$model
rep_gen_df$model <- NULL
rep_gen_df$genealogy <- gsub("-2", ".2", rep_gen_df$genealogy, fixed = TRUE)
vsd_expr_df <- read.table(gzfile(vsd_expr_f), quote = "", sep = "\t", header=TRUE)
meta_df<-read.table(meta_f, quote = "", sep = "\t", header = TRUE)

### Filter for lmo_basali in meta dataset
class <- snakemake@wildcards[["sclass"]]
meta_df_lmo_basale<-filter(meta_df, grepl(class, type))

### filtering expression data: whigh sd genes but not clear outliers/not expressed genes
### filter not expressed genes
means <- apply(vsd_expr_df, 1, mean)
med <- median(means)

### I obtain a data.table that contain only the genes  mean expression over the median
vsd_expr_df_filtered <- as.data.frame(vsd_expr_df[means > med,])

### I want to calculate the standard deviation
sds <- apply(vsd_expr_df_filtered, 1, sd) 

### now we keep the top 10% variable genes 
sds <- sds[order(-sds)]
n <- length(sds)
keep <- head(sds, round(0.10*n))
keep_genes <- names(keep)
desd <- vsd_expr_df_filtered[rownames(vsd_expr_df_filtered) %in% keep_genes,]

### Transmute e filter desd for LMO_BASALI
desd <- as.data.frame(t(desd))
desd <- tibble::rownames_to_column(desd)
colnames(desd)[1] <- "sample_id_R"
meta_df_lmo_basale$sample_id_R <- gsub("-2", ".2", meta_df_lmo_basale$sample_id_R, fixed = TRUE)
vsd_expr_df_lmo=desd[desd$sample_id_R %in% meta_df_lmo_basale$sample_id_R, ]
rownames(vsd_expr_df_lmo) <- vsd_expr_df_lmo$sample_id_R
vsd_expr_df_lmo$sample_id_R <- NULL


### Split genealogy
pearson_df <- data.frame(row.names = rownames(rep_gen_df), pvalue = rep(0, length=nrow(rep_gen_df)),
                                                                        correlation = rep(0, length = nrow(rep_gen_df)))
                        
for (i in seq(1, nrow(rep_gen_df))) {
  rep_gen_df_split <- strsplit(rep_gen_df[i,"genealogy"], ",")
  rep_gen_df_split_unlist <- unlist(rep_gen_df_split)
  x <- rep_gen_df_split_unlist[1]
  y <- rep_gen_df_split_unlist[2]
  pearson_i <- cor.test(as.numeric(vsd_expr_df_lmo[rownames(vsd_expr_df_lmo)==x, ]), 
                                      as.numeric(vsd_expr_df_lmo[rownames(vsd_expr_df_lmo)==y, ]), method = "pearson")
  pearson_df[i, "pvalue"] <- pearson_i$p.value
  pearson_df[i, "correlation"] <- pearson_i$estimate
  }

write.table(pearson_df, file = genealogy_VSD_correlation, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

sink(summary_f) 

summary(pearson_df)

sink()