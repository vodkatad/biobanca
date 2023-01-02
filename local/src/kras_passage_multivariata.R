##kras passage

library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(tidyverse)


passage <- "/scratch/trcanmed/biobanca/local/share/data/passaggi_query_las_Simo_march2022.txt"
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

fra_f <- "/scratch/trcanmed/biobanca/dataset/V1/enrichment/fra_mutational_annotation.tsv"
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

pdo <- "/scratch/trcanmed/biobanca/local/share/data/whoiswho_validation_xen_nolmh.tsv"
pdo <- read.table(pdo, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
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

loc <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/clinical_data_for_circos.tsv"
loc <- read.table(loc, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(loc)[names(loc) == "CASE"] <- "model"

res2 <- merge(res, loc, by = "model")
res2 <- res2[, c("model", "KRAS.x", "passage", "validated", "SITE.OF.PRIMARY")]
colnames(res2) <- c("model", "KRAS", "Passage", "Validation", "Site.of.Primary")
rownames(res2) <- res2$model

for (i in rownames(res2)) {
  if (res2[i, "Passage"] == 1 || res2[i, "Passage"] == 2) {
    res2[i, "Passage"] <- "Early Passage"
  } else {
    res2[i, "Passage"] <- "Late Passage"
  }
}

fit.full <- glm(KRAS ~ validated + passage, data=res,family=binomial())
fit.alt <- glm(KRAS ~ passage + validated, data = res, family = binomial())
fit.kras <- glm(Validation ~ Passage + Site.of.Primary + KRAS, data=res2,family=binomial())
pdf("/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/multivariata_KRAS_passage_site.pdf",useDingbats=FALSE)
plot_model(fit.kras, axis.lim = c(0.1, 2), title = "Validated") 
dev.off()
fit <- as.data.frame(summary.glm(fit.kras)$coefficients)

write.table(res2, file = "passage_validation_site_kras.tsv", quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = FALSE)

#colors = c("royalblue2", "firebrick2")
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
