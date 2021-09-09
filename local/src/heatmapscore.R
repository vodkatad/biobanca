### heatmap for DELTAmax

library(pheatmap)
library(reshape)

res_models_drugs_f <- snakemake@input[["res_models_drugs_f"]]
out_dir <- snakemake@output[["out_dir"]]


dmax_f <- res_models_drugs_f
dmax_df <- read.table(dmax_f, quote = "", header = TRUE, sep="\t")
min_v <- min(dmax_df$SCORE)
max_v <- max(dmax_df$SCORE)
# do we need to create dir or does snakemake take care of it?
if (dir.exists(out_dir)) {
   unlink(out_dir)
} else {
   dir.create(out_dir)
}

setwd(out_dir)

exp <- c(1,2,3)
for (i in seq(1, length(exp))) {
  exp_time_i <- exp[i]
  subset_dmax_df <- subset(dmax_df, dmax_df$EXP == exp_time_i )
  dm <- cast(subset_dmax_df, formula = "MODEL_DRUG~DOSE", value = "SCORE")
  rownames(dm) <- dm$MODEL_DRUG
  dm$MODEL_DRUG <- NULL
  pheatmap(dm, cluster_rows = FALSE, cluster_cols = FALSE, 
           filename = paste0("Synergy_Score_EXP", exp_time_i, ".pdf"), 
           main = paste0("EXP", exp_time_i), display_numbers = TRUE, number_format = "%.2f")
  pheatmap(dm, cluster_rows = FALSE, cluster_cols = FALSE, 
           breaks = seq(min_v, max_v, by = 0.5), filename = paste0("Synergy_Score_EXP", exp_time_i, "_0.5.pdf"), 
           main = paste0("EXP", exp_time_i, "0.05"), display_numbers = TRUE, number_format = "%.2f")
}

