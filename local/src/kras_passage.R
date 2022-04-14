library(ggplot2)
library(dplyr)
txt <- "/scratch/trcanmed/biobanca/local/share/data/passaggi_query_las_Simo_march2022.txt"
txt <- read.table(txt, quote= "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

txt$origin <- substr(txt$T.Cell.line.has.aliquot...Genealogy.ID_0, 8, 10)
lmx <- txt %>% filter(origin == "LMX")
lmx$model <- substr(lmx$Genealogy.ID, 1, 7)
lmx$passage <- as.numeric(substr(lmx$T.Cell.line.has.aliquot...Genealogy.ID_0, 13, 14))
#lmx <- lmx %>% distinct(model, .keep_all= TRUE)
lmx <- lmx[!(duplicated(lmx$model) | duplicated(lmx$model, fromLast = TRUE)), ]
lmx <- lmx %>% filter(passage == 1)

casi <- as.data.frame(cbind(lmx$Genealogy.ID, lmx$model, lmx$T.Cell.line.has.aliquot...Genealogy.ID_0, lmx$origin, lmx$passage))
colnames(casi) <- c("Genealogy.ID", "model", "Genealogy.ID_0", "origin", "passage")

fra <- "/scratch/trcanmed/biobanca/dataset/V1/enrichment/fra_mutational_annotation.tsv"
fra <- read.table(fra, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(fra) <- fra$genes
fra$genes <- NULL
tfra <- as.data.frame(t(fra))
tfra$NRAS <- NULL
tfra$BRAF <- NULL
tfra$PIK3CA <- NULL
tfra$model <- rownames(tfra)

merged <- merge(tfra, casi, by= "model")
ggplot(data=merged, aes(x=KRAS))+geom_bar()+theme_bw()