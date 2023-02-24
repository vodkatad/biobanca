### CMS with simo coupled

library(tidyverse)

simo_f <- snakemake@input[["simo"]]
cmslmo_f <- snakemake@input[["pdo"]]
cmslmx_f <- snakemake@input[["pdx"]]
results <- snakemake@output[["results"]]

#simo_f <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/cris/SIMO_v2.tsv"
simo <- read.table(simo_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(simo)[names(simo) == "LMO_lineage"] <- "PDO_lineage"
names(simo)[names(simo) == "LMX_lineage"] <- "PDX_lineage"

#cmslmo_f <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMO_BASALE_CMScaller.tsv"
cmslmo <- read.table(cmslmo_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmslmo_df <- cmslmo %>% select(prediction)
cmslmo_df <- tibble::rownames_to_column(cmslmo_df, "PDO_long_lineage")
cmslmo_df$PDO_lineage <- substr(cmslmo_df$PDO_long_lineage, 0, 12)
names(cmslmo_df)[names(cmslmo_df) == 'prediction'] <- 'CMS_PDO'

#cmslmx_f <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMX_BASALE_CMScaller.tsv"
cmslmx <- read.table(cmslmx_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmslmx_df <- cmslmx %>% select(prediction)
cmslmx_df <- tibble::rownames_to_column(cmslmx_df, "PDX_long_lineage")
cmslmx_df$PDX_lineage <- substr(cmslmx_df$PDX_long_lineage, 0, 12)
names(cmslmx_df)[names(cmslmx_df) == 'prediction'] <- 'CMS_PDX'

simo_lmo <- merge(simo, cmslmo_df, by="PDO_long_lineage")
simo_lmo$model <- substr(simo_lmo$PDO_long_lineage, 0, 7)
simo_lmx <- merge(simo, cmslmx_df, by="PDX_long_lineage")
simo_lmx$model <- substr(simo_lmx$PDX_long_lineage, 0, 7)

# We can safely merge on model only because we removed the wrong lineages with the previous merge
# with simo pairs using the whole lineage.
simo_df <- merge(simo_lmo, simo_lmx, by = "model") 
model_cms <- simo_df %>% select(model, CMS_PDO, CMS_PDX, PDO_lineage.x, PDX_lineage.y)

# wrong_pdo <- c()
# wrong_pdx <- c()
# res <- data.frame(stringsAsFactors = FALSE)
# smodel <- unique(model_cms$model)
# 
# #i=1
# ## with this loop we keep in wrong pdo/pdx all the "misclassificated" CMS genealogy,
# ## also the NC
# ## we print in res all the cases single or if there is a replicate with the same PDO classification
# for (i in seq(1, length(smodel))) {
#   smodel_i <- smodel[i]
#   cms_i <- model_cms[model_cms$model == smodel_i,]
#   n_cms_pdo <- length(unique(cms_i$CMS_PDO))
#   n_cms_pdx <- length(unique(cms_i$CMS_PDX))
#   if (n_cms_pdo > 1) {
#     wrong_pdo <- c(wrong_pdo, smodel_i)
#   } else if (n_cms_pdx > 1) {
#     wrong_pdx <- c(wrong_pdx, smodel_i)
#   } else { #(n_cms_pdo == 1 && n_cms_pdx == 1) 
#     res <- rbind(res, c(smodel_i, unique(cms_i$CMS_PDO), unique(cms_i$CMS_PDX)), stringsAsFactors = FALSE)
#   } 
# }
# 
# pdo_cms_wrong <- filter(model_cms, model %in% wrong_pdo)
# pdx_cms_wrong <- filter(model_cms, model %in% wrong_pdx)

res <- model_cms[,1:3]
colnames(res) <- c("model", "CMS_PDO", "CMS_PDX")
res <- res[,c(1,3,2)]
#reorder the column in order to have PDXs that goes into PDOs for sankey
#res <- res[, c("model", "CMS_PDX", "CMS_PDO")]

#save.image("pippo.Rdata")

# sink(snakemake@log[[1]])
# print(paste0("wrong pdo: ", length(wrong_pdo)))
# print(paste0("wrong pdx: ", length(wrong_pdx)))
# sink()

write.table(res, file = results, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)

# write.table(pdo_cms_wrong, file = results_pdo, quote = FALSE, sep = "\t", 
#             row.names = FALSE, col.names = TRUE)