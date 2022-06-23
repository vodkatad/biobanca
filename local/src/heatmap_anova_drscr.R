### General heatmap for drug screening with anova values

library(reshape)
library(tidyverse)
library(pheatmap)

anova_f <- snakemake@input[["anova"]]
anova_l <- snakemake@input[["anova_long"]]
anova_n <- snakemake@input[["anova_norm"]]
first_drscr <- snakemake@output[["h_anova"]]
second_drscr <- snakemake@output[["h_long"]]
norm_drscr <- snakemake@output[["h_anova_norm"]]

#anova_f <- "/scratch/trcanmed/biobanca/dataset/V1/drug_screening/anova.tsv"
anova <- read.table(anova_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

anova_df <- cast(anova, formula = "MODEL+CONDITION~DRUG", value = "padj")

anova_df$MODEL_CONDITION <- paste(anova_df$MODEL, anova_df$CONDITION)
anova_df$MODEL <- NULL
anova_df$CONDITION <- NULL

anova_fin <- anova_df %>% remove_rownames %>% column_to_rownames(var="MODEL_CONDITION")
anova_hm <- as.data.frame(t(anova_fin))

pdf(first_drscr)
pheatmap(anova_hm, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

#anova_l <- "/scratch/trcanmed/biobanca/dataset/V1/drug_screening/anova_long.tsv"
anova_long <- read.table(anova_l, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

anova_df_long <- cast(anova_long, formula = "MODEL+TYPE~DRUG", value = "padj")

anova_df_long$MODEL_TYPE <- paste(anova_df_long$MODEL, anova_df_long$TYPE)
anova_df_long$MODEL <- NULL
anova_df_long$TYPE <- NULL

anova_fin_long <- anova_df_long %>% remove_rownames %>% column_to_rownames(var="MODEL_TYPE")
anova_hm_long <- as.data.frame(t(anova_fin_long))

pdf(second_drscr)
pheatmap(anova_hm_long, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

#anova_n <- "/scratch/trcanmed/biobanca/dataset/V1/drug_screening/anova_norm.tsv"
anova_norm <- read.table(anova_n, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

anova_norm_df <- cast(anova_norm, formula = "MODEL+CONDITION~DRUG", value = "padj")

anova_norm_df$MODEL_CONDITION <- paste(anova_norm_df$MODEL, anova_norm_df$CONDITION)
anova_norm_df$MODEL <- NULL
anova_norm_df$CONDITION <- NULL

anova_norm_fin <- anova_norm_df %>% remove_rownames %>% column_to_rownames(var="MODEL_CONDITION")
anova_norm_hm <- as.data.frame(t(anova_norm_fin))

pdf(norm_drscr)
pheatmap(anova_norm_hm, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

