##kras passage

library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(tidyverse)

passage <- snakemake@input[["pas"]]
fra_f <- snakemake@input[["mut"]]
buoni_f <- snakemake@input[["pdo"]]
loc_f <- snakemake@input[["side"]]
result <- snakemake@output[["res"]]
plot_fit <- snakemake@output[["fit_plot"]]
log_f <- snakemake@log[['log']]

#passage <- "/scratch/trcanmed/biobanca/local/share/data/passaggi_query_las_Simo_march2022.txt"
txt <- read.table(passage, quote= "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

txt$origin <- substr(txt$T.Cell.line.has.aliquot...Genealogy.ID_0, 8, 10)
lmx <- txt %>% filter(origin == "LMX")
lmx$model <- substr(lmx$Genealogy.ID, 1, 7)
lmx$passage <- as.numeric(substr(lmx$T.Cell.line.has.aliquot...Genealogy.ID_0, 13, 14))
#lmx <- lmx %>% distinct(model, .keep_all= TRUE)
lmx <- lmx[!(duplicated(lmx$model) | duplicated(lmx$model, fromLast = TRUE)), ]
#lmx <- lmx %>% filter(passage == 1)

casi <- as.data.frame(cbind(lmx$Genealogy.ID, lmx$model, lmx$T.Cell.line.has.aliquot...Genealogy.ID_0, lmx$origin, lmx$passage))
colnames(casi) <- c("Genealogy.ID", "model", "Genealogy.ID_0", "origin", "passage")

#fra_f <- "/scratch/trcanmed/biobanca/dataset/V1/enrichment/fra_mutational_annotation.tsv"
fra <- read.table(fra_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(fra) <- fra$genes
fra$genes <- NULL
tfra <- as.data.frame(t(fra))
tfra$NRAS <- NULL
tfra$BRAF <- NULL
tfra$PIK3CA <- NULL
tfra$model <- rownames(tfra)

merged <- merge(tfra, casi, by= "model")
#ggplot(data=merged, aes(x=KRAS))+geom_bar()+theme_bw()
merged$Genealogy.ID <- NULL
merged$Genealogy.ID_0 <- NULL
merged$origin <- NULL

#buoni_f <- "/scratch/trcanmed/biobanca/local/share/data/whoiswho_validation_xen_nolmh.tsv"
pdo <- read.table(buoni_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pdo <- pdo[pdo$type !=  'Validation not performed',]
for (i in rownames(pdo)) {
  if (pdo[i, "type"] == "Validation successful"){
    pdo [i, "type"] <- TRUE
  } else {
    pdo [i, "type"] <- FALSE
  }
}
colnames(pdo) <- c("model", "validated")

res <- merge(merged, pdo, by = "model")
res$KRAS <- as.character(res$KRAS)
res$passage <- as.numeric(res$passage)
for (i in seq(length(res$model))){
 res[i,2] <- ifelse(res[i,2] == "True", "Yes", "No")
 res[i,4] <- ifelse(res[i,4] == TRUE, "validated", "non_validated")
}
res$KRAS <- as.factor(res$KRAS)
res$validated <- as.factor(res$validated)

#loc_f <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/clinical_data_for_circos_revision.tsv"
loc <- read.table(loc_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(loc)[names(loc) == "CASE"] <- "model"

res2 <- merge(res, loc, by = "model")
res2 <- res2[, c("model", "KRAS.x", "passage", "validated", "SITE.OF.PRIMARY")]
colnames(res2) <- c("model", "KRAS", "Passage", "Validation", "Site.of.Primary")
rownames(res2) <- res2$model

# for (i in rownames(res2)) {
#   if (res2[i, "Passage"] == 1 || res2[i, "Passage"] == 2) {
#     res2[i, "name_Passage"] <- "Early Passage"
#   } else {
#     res2[i, "name_Passage"] <- "Late Passage"
#   }
# }
sink(log_f)
print(table(res2$Validation))
sink()
#fit.full <- glm(KRAS ~ validated + passage, data=res,family=binomial())
#fit.alt <- glm(KRAS ~ passage + validated, data = res, family = binomial())
fit.kras <- glm(Validation ~ Passage + Site.of.Primary + KRAS, data=res2,family=binomial())
#pdf("/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/multivariata_KRAS_passage_site.pdf",useDingbats=FALSE)
pdf(plot_fit, useDingbats=FALSE)
plot_model(fit.kras, axis.lim = c(0.1, 2), title = "Validated") 
dev.off()
fit <- as.data.frame(summary.glm(fit.kras)$coefficients)

conf_intervals <- confint(fit.kras)
fit2 <- cbind(fit, exp(conf_intervals))
## odds ratio
# Obtain coefficients
coefficients <- coef(fit.kras)
odds_ratios <- as.data.frame(exp(coefficients))
colnames(odds_ratios) <- "odds ratio"
fit3 <- cbind(fit2, odds_ratios)

rownames(fit3) <- c("Intercept", "Passage", "Site of primary (right colon)", "KRAS (mutant)")
write.table(fit3, file = result, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

## part for revision
for (i in rownames(res2)) {
   if (res2[i, "Passage"] <=3) {
     res2[i, "name_Passage"] <- "Early Passage"
   } else {
     res2[i, "name_Passage"] <- "Late Passage"
   }
 }
#colors = c("royalblue2", "firebrick2")
table_test_3 <- res2
table_test_3 <- table_test_3 %>% filter(name_Passage == "Early Passage")
table_test_3 <- as.data.frame(table(table_test_3$KRAS, table_test_3$name_Passage))
rownames(table_test_3) <- c("KRAS_wt", "KRAS_mut")
table_test_3$Var1 <- NULL
table_test_3$Var2 <- NULL
colnames(table_test_3) <- "Early Passage 1-3"
table_test_3 <- as.data.frame(t(table_test_3))

table_test_late <- res2
table_test_late <- table_test_late %>% filter(name_Passage == "Late Passage")
table_test_late <- as.data.frame(table(table_test_late$KRAS, table_test_late$name_Passage))
rownames(table_test_late) <- c("KRAS_wt", "KRAS_mut")
table_test_late$Var1 <- NULL
table_test_late$Var2 <- NULL
colnames(table_test_late) <- "Late Passage 4-9"
table_test_late <- as.data.frame(t(table_test_late))

table_test <- rbind(table_test_3, table_test_late)

# questo non viene significativo e riflette il fatto che se mettiamo un threshold secco early/late non vediamo l'effetto (succedeva anche in multivariata)
# la continua invece è borderline
fisher.test(table_test)

passage_kras_val <- res2
passage_kras_val <- passage_kras_val %>% filter(Passage == 3)
table3_kras_val <- as.data.frame(matrix(ncol = 2, nrow = 2))
rownames(table3_kras_val) <- c("KRASwt", "KRASmut")
colnames(table3_kras_val) <- c("not_validated", "validated")
contingency <- table(passage_kras_val$Validation, passage_kras_val$KRAS)
fisher.test(as.matrix(contingency))

status <- read.table("/scratch/trcanmed/biobanca/local/share/data/XENTURION_DEF_SML_12-10.tsv", quote = "",
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)
res2 <- res
names(res2)[names(res2) == 'model'] <- 'CASE'

merged <- merge(res2, status, by = "CASE")
engraf <- as.data.frame(cbind(merged$CASE, merged$KRAS, merged$passage, merged$DERIVATION..successful.failed.))
colnames(engraf) <- c("model", "KRAS", "passage", "derivation")
engraf$passage <- as.numeric(engraf$passage)
engraf$model <- as.character(engraf$model)

fit.engraf <- glm(derivation ~ passage + KRAS, data=engraf, family=binomial())
plot_model(fit.engraf)
