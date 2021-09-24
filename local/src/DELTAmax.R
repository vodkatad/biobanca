### Analysis DELTAmax between cetuximab e combo

library(readxl)
library(tidyverse)
library(stringr)
library(reshape)
library(pheatmap)

drugs_tables <- snakemake@input[["drugs_tables"]]
anova <- snakemake@input[["anova"]]
res_models_drugs_f <- snakemake@output[["res_models_drugs"]]

anova_d <- anova
anova_df <- read.table(anova_d, quote = "", sep = "\t", header = TRUE)
anova_df <- subset(anova_df, anova_df$padj < 0.05)
models <- c("CRC0322", "CRC0059", "CRC1272", "CRC1331", "CRC0327")

# Drug_Screening_Tables <- read_excel(("Drug_Screening_Tables.xlsx"), sheet = models[i])
# data <- as.data.frame(Drug_Screening_Tables)
# data$DRUG <- str_trim(data$DRUG, "both")
# data$DRUG <- str_replace_all(data$DRUG, " +", "_")
# res$model <- models[i]
 

library(rio)
data_list <- import_list(drugs_tables, setclass = "tsv", rbind = TRUE)
order_m <- data.frame(order=c("1", "2", "3", "4", "5"), models)
data_list <- merge(data_list, order_m, by.x = "_file", by.y = "order")
save.image('pippo.Rdata')
data_list$DRUG <- str_trim(data_list$DRUG, "both")
data_list$DRUG <- str_replace_all(data_list$DRUG, " +", "_")
data_list$model_drug_cond <- c(paste(data_list$models, data_list$DRUG, data_list$CONDITION))

anova_df$DRUG <- str_replace_all(anova_df$DRUG, "__", "_")
anova_df$model_drug_cond <- c(paste(anova_df$MODEL, anova_df$DRUG, anova_df$CONDITION))

total <- merge(data_list, anova_df, by="model_drug_cond")  
total$`_file` <- NULL
total$DRUG.y <- NULL
total$CONDITION.y <- NULL
total$models <- NULL
total$model_drug_cond <- NULL
total$PVALUE<- NULL
total$padj<- NULL
colnames(total)[colnames(total) == 'CONDITION.x'] <- 'CONDITION'
colnames(total)[colnames(total) == 'DRUG.x'] <- 'DRUG'

total$model_drug <- c(paste(total$MODEL, total$DRUG))

all_models_drugs <- unique(total$model_drug)

res_models_drugs <- data.frame()

compute_score_drug <- function(sd, DOSE){
  w0 <- sd[sd$CONDITION=="Combo", colnames(sd) == "DOSE_1"]
  w1 <- sd[sd$CONDITION=="Combo", colnames(sd) == DOSE]
  v0 <- sd[sd$CONDITION=="Mono", colnames(sd) == "DOSE_1"]
  v1 <- sd[sd$CONDITION=="Mono", colnames(sd) == DOSE]
  #print(paste(w0, w1, v0, v1))
  res_d <- ((v1-w1)/(v0-w0))-1
  return(res_d)
}

doses <- c("DOSE_2", "DOSE_3", "DOSE_4", "DOSE_5")
exp <- c(1, 2, 3)
for (i in seq(1, length(all_models_drugs))) {
  model_drug_i <- all_models_drugs[i]
  subset_drug <- subset(total, model_drug == model_drug_i)
    for (k in seq(1, length(exp))) {
      subset_drug_model <- subset(subset_drug, EXP == exp[k])
      print(k)
        for (j in seq(1, length(doses))){
          print(paste0('expk',exp[k]))
          if ((nrow(subset_drug_model)==2)) {
            score <- compute_score_drug(subset_drug_model, doses[j])
            res_f <- c(model_drug_i, doses[j], score, exp[k])
            res_models_drugs <-rbind(res_models_drugs, res_f, stringsAsFactors = FALSE)
           
              } else {
                  res_f2 <- c(model_drug_i, doses[j], exp[k])
                  #print(paste0(c("Skipping", res_f2), collapse = " "))
                  }
            }  
      }
}

colnames(res_models_drugs) <- c("MODEL_DRUG", "DOSE", "SCORE", "EXP")

write.table(res_models_drugs, file = res_models_drugs_f, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)
