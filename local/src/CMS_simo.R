### CMS with simo coupled

library(tidyverse)

simo_f <- snakemake@input[["simo"]]
cmslmo_f <- snakemake@input[["pdo"]]
cmslmx_f <- snakemake@input[["pdx"]]
results <- snakemake@output[["results"]]

#simo_f <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/cris/SIMO_v2.tsv"
simo <- read.table(simo_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#cmslmo_f <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMO_BASALE_CMScaller.tsv"
cmslmo <- read.table(cmslmo_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmslmo_df <- cmslmo %>% select(prediction)
cmslmo_df <- tibble::rownames_to_column(cmslmo_df, "Genealogy")
cmslmo_df$PDO_lineage <- substr(cmslmo_df$Genealogy, 0, 12)
names(cmslmo_df)[names(cmslmo_df) == 'prediction'] <- 'CMS_PDO'
                                 
#cmslmx_f <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/LMX_BASALE_CMScaller.tsv"
cmslmx <- read.table(cmslmx_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cmslmx_df <- cmslmx %>% select(prediction)
cmslmx_df <- tibble::rownames_to_column(cmslmx_df, "Genealogy")
cmslmx_df$PDX_lineage <- substr(cmslmx_df$Genealogy, 0, 12)
names(cmslmx_df)[names(cmslmx_df) == 'prediction'] <- 'CMS_PDX'

simo_lmo <- merge(simo, cmslmo_df, by="PDO_lineage")
simo_lmo$model <- substr(simo_lmo$PDO_lineage, 0, 7)
simo_lmx <- merge(simo, cmslmx_df, by="PDX_lineage")
simo_lmx$model <- substr(simo_lmx$PDX_lineage, 0, 7)

simo_df <- merge(simo_lmo, simo_lmx, by = "model")
model_cms <- simo_df %>% select(model, CMS_PDO, CMS_PDX, PDO_lineage.x, PDX_lineage.y)

wrong_pdo <- c()
wrong_pdx <- c()
res <- data.frame(stringsAsFactors = FALSE)
smodel <- unique(model_cms$model)

#i=1
for (i in seq(1, length(smodel))) {
  smodel_i <- smodel[i]
  cms_i <- model_cms[model_cms$model == smodel_i,]
  n_cms_pdo <- length(unique(cms_i$CMS_PDO))
  n_cms_pdx <- length(unique(cms_i$CMS_PDX))
  if (n_cms_pdo == 1 & n_cms_pdx == 1) {
    res <- rbind(res, c(smodel_i, unique(cms_i$CMS_PDO), unique(cms_i$CMS_PDX)), stringsAsFactors = FALSE)
      }   
    else {
      if (n_cms_pdo > 1) {
        wrong_pdo <- c(wrong_pdo, smodel_i)
      }
      if (n_cms_pdx > 1) {
        wrong_pdx <- c(wrong_pdx, smodel_i)
      }  
    }
  }

pdo_cms_wrong <- filter(model_cms, model %in% wrong_pdo)
pdx_cms_wrong <- filter(model_cms, model %in% wrong_pdx)

colnames(res) <- c("model", "CMS_PDO", "CMS_PDX")

save.image("pippo.Rdata")

sink(snakemake@log[[1]])
print(paste0("wrong pdo: ", length(wrong_pdo)))
print(paste0("wrong pdx: ", length(wrong_pdx)))
sink()

write.table(res, file = results, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)