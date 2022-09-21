## clinical data multivariata

library(tidyverse)
library(sjPlot)


print("se clinical data non viene caricata controllare tramite visual studio code che 
alcune righe non contengano un invio non visibile su rstudio e cancellare l'invio
controllare che il tab sia corretto per mantere le colonne corrette")
## se clinical data non viene caricata controllare tramite visual studio code che 
## alcune righe non contengano un invio non visibile su rstudio e cancellare l'invio
## controllare che il tab sia corretto per mantere le colonne corrette

## carico i dati clinici, filtro per i pdo buoni e per le mutazioni annotate di fra
cl_f <- "/scratch/trcanmed/biobanca/local/share/data/clinical_data_done.tsv"
cl <- read.table(cl_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

buoni_f <- "/scratch/trcanmed/biobanca/local/share/data/biobanca_pdo_buoni.tsv"
buoni <- read.table(buoni_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(buoni)[names(buoni) == "smodel"] <- "CASE"

xen_f <- "/scratch/trcanmed/biobanca/local/share/data/XENTURION_DEF_SML_12-10.tsv"
xen <- read.table(xen_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
np <- xen %>% filter(VALIDATION..successful.failed.not.performed. == "Not performed")
np <- np$CASE
qc <- xen %>% filter(QC..passed.not.passed. == "Not passed")
qc <- qc$CASE
## che fare con gli LMH?

add <- c(np, qc)
added <- data.frame(matrix(ncol = 2, nrow = 52))
colnames(added) <- c("CASE", "buoni")
added$CASE <- add
added$buoni <- "Removed"

finalexen <- as.data.frame(rbind(buoni, added))

merged <- merge(cl, finalexen, by = "CASE")
## uno in meno perché già eliminata la metastasi sincrona

fra_f <- "/scratch/trcanmed/biobanca/dataset/V1/enrichment/fra_mutational_annotation.tsv"
fra <- read.table(fra_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(fra) <- fra$genes
fra$genes <- NULL
tfra <- as.data.frame(t(fra))
tfra$CASE <- rownames(tfra)
rownames(tfra) <- NULL

res <- merge(merged, tfra, by = "CASE", all.x = TRUE)
rownames(res) <- res$CASE
row_names_df_to_remove<-c("CRC1379", "CRC1479", "CRC1737")
res <- res[!(row.names(res) %in% row_names_df_to_remove),]

## creo un df con solo i dati necessari per la multivariata
res <- res[c("CASE", "AGE.AT.COLLECTION", "SEX", "SITE.OF.PRIMARY", "T", "N", "THERAPY.BEFORE..Y.N.",
             "KRAS", "BRAF", "NRAS", "MSI_MSS", "buoni")]
names(res)[names(res) == "T"] <- "Classification_T"
names(res)[names(res) == "N"] <- "Classification_N"

## sistemo il df per poter accorpare i siti del primario

res["CRC1241", "SITE.OF.PRIMARY"] <- "NONE"
res["CRC1575", "SITE.OF.PRIMARY"] <- "NONE"
res["CRC0578", "SITE.OF.PRIMARY"] <- "NONE"
for (i in seq(res$CASE)) {
  if (res[i, "SITE.OF.PRIMARY"] == "SIGMOID COLON" | res[i, "SITE.OF.PRIMARY"] == "ANO" | 
      res[i, "SITE.OF.PRIMARY"] == "LEFT COLON" | res[i, "SITE.OF.PRIMARY"] == "RECTUM" |
      res[i, "SITE.OF.PRIMARY"] == "SPLENIC FLEXURE"){ 
    res[i, "SITE.OF.PRIMARY"] <- "LEFT COLON"
  } else if (res[i, "SITE.OF.PRIMARY"] == "CAECUM" | res[i, "SITE.OF.PRIMARY"] == "TRANSVERSE COLON" | 
             res[i, "SITE.OF.PRIMARY"] == "HEPATIC FLEXURE" | res[i, "SITE.OF.PRIMARY"] == "RIGHT COLON" ) {
    res[i, "SITE.OF.PRIMARY"] <- "RIGHT COLON"
  } else {
    res[i, "SITE.OF.PRIMARY"] <- NA
  }
}

## sistemo T ed N per poterli rendere numerici senza perdere quelli con la lettera
res["CRC1241", "Classification_T"] <- "N"
res["CRC1575", "Classification_T"] <- "N"
res["CRC0578", "Classification_T"] <- "N"
res["CRC1390", "Classification_T"] <- "N"

res["CRC1241", "Classification_N"] <- "N"
res["CRC1575", "Classification_N"] <- "N"
res["CRC0578", "Classification_N"] <- "N"
res["CRC1390", "Classification_N"] <- "N"


for (i in seq(res$CASE)) {
  if (res[i, "Classification_T"] == "4A") {
    res[i, "Classification_T"] <- "4"
  } else if (res[i, "Classification_T"] == "4B") {
    res[i, "Classification_T"] <- "4"
  } else {
    res[i, "Classification_T"] <- res[i, "Classification_T"]
  }
}

for (i in seq(res$CASE)) {
  if (res[i, "Classification_N"] == "1A") {
    res[i, "Classification_N"] <- "1"
  } else if (res[i, "Classification_N"] == "1B") {
    res[i, "Classification_N"] <- "1"
  } else if (res[i, "Classification_N"] == "1C") {
    res[i, "Classification_N"] <- "1"
  } else if (res[i, "Classification_N"] == "2A") {
    res[i, "Classification_N"] <- "2"
  } else if (res[i, "Classification_N"] == "2B") {
    res[i, "Classification_N"] <- "2"
  } else {
    res[i, "Classification_N"] <- res[i, "Classification_N"]
  }
}

## rendo tutto numerico o factor a seconda delle necessità

res$AGE.AT.COLLECTION <-  gsub('N', NA, res$AGE.AT.COLLECTION)
res$AGE.AT.COLLECTION <- as.numeric(res$AGE.AT.COLLECTION)

res$SEX <-  gsub('N', NA, res$SEX)
res$SEX <- as.factor(res$SEX)

res$Classification_T <- gsub("N", NA, res$Classification_T)
res$Classification_T <- as.numeric(res$Classification_T)

res$Classification_N <- gsub("N", NA, res$Classification_N)
res$Classification_N <- as.numeric(res$Classification_N)

res$SITE.OF.PRIMARY <- gsub("NONE", NA, res$SITE.OF.PRIMARY)
res$SITE.OF.PRIMARY <- as.factor(res$SITE.OF.PRIMARY)

## stage non lo metto più perché sono definiti tutti 4
#res$STAGE <-  gsub('N', NA, res$STAGE)
#res$STAGE <- as.numeric(res$STAGE)

res$THERAPY.BEFORE..Y.N. <- as.factor(res$THERAPY.BEFORE..Y.N.)

is.na(res$MSI_MSS) <- NA
res$MSI_MSS <-  gsub('NT', NA, res$MSI_MSS)
res$MSI_MSS <- as.factor(res$MSI_MSS)

res$buoni <-  gsub('TRUE', "True", res$buoni)
res$buoni <- gsub("FALSE","False", res$buoni)
res$buoni <- as.factor(res$buoni)
## scrivo la tabella per i circos che comprenda anche i casi
write.table(res, file = "complete_data_for_circos.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

res$CASE <- NULL
#res_prova <- res 
#res_prova <- res_prova %>% filter(!KRAS == "True")

## faccio il fit

#fit.full <- glm(buoni ~ SEX + AGE.AT.COLLECTION, data=res, family=binomial())
#fit.full <- glm(buoni ~ SEX + AGE.AT.COLLECTION + STAGE + THERAPY.BEFORE..Y.N. , data=res, family=binomial())
fit.full <- glm(buoni ~ SEX + AGE.AT.COLLECTION + THERAPY.BEFORE..Y.N. + SITE.OF.PRIMARY + Classification_T + Classification_N+ KRAS + BRAF + NRAS, data=res, family=binomial())
pdf(plot_fit)
plot_model(fit.full)
dev.off()
#fit.full <- glm(buoni ~ SEX + AGE.AT.COLLECTION + THERAPY.BEFORE..Y.N. + SITE.OF.PRIMARY + BRAF + NRAS, data=res_prova, family=binomial())
fit <- as.data.frame(summary.glm(fit.full)$coefficients)
write.table(fit, file = res_fit, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)