### General heatmap for drug screening with anova values

library(reshape)
library(tidyverse)
library(pheatmap)

anova_f <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/drug_screening/anova.tsv"
anova <- read.table(anova_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

anova_df <- cast(anova, formula = "MODEL+CONDITION~DRUG", value = "padj")

anova_df$MODEL_CONDITION <- paste(anova_df$MODEL, anova_df$CONDITION)
anova_df$MODEL <- NULL
anova_df$CONDITION <- NULL

anova_fin <- anova_df %>% remove_rownames %>% column_to_rownames(var="MODEL_CONDITION")
anova_hm <- as.data.frame(t(anova_fin))

heatmap <- pheatmap(anova_hm, cluster_rows = FALSE, cluster_cols = FALSE)