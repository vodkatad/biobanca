## jaccard's score with CMS and CRIS
library(tidyverse)

file_cms <- snakemake@input[["cms"]]
file_cris <- snakemake@input[["cris"]]
result_cms <- snakemake@output[["cms_t"]]
result_cris <- snakemake@output[["cris_t"]]

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

## CMS 
#cms <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/reorder_model_cms_validated.tsv"
cms <- read.table(file_cms, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cms$CMS_PDX <- cms$CMS_PDX %>% replace_na('NC')
cms$CMS_PDO <- cms$CMS_PDO %>% replace_na('NC')

cms1 <- cms %>% filter(CMS_PDX == "CMS1")
##mettere anche i cms1 che ci sono solo nei pdo
cms1pdo <- cms %>% filter(CMS_PDO == "CMS1")
allcms1 <- as.character(cms1$model)
cases1c <- as.character(cms1pdo$model)
cms1j <- jaccard(allcms1, cases1c)

cms2 <- cms %>% filter(CMS_PDX == "CMS2")
cms2pdo <- cms %>% filter(CMS_PDO == "CMS2")
allcms2 <- cms2$model
cases2c <- cms2pdo$model
cms2j <- jaccard(allcms2, cases2c)

cms3 <- cms %>% filter(CMS_PDX == "CMS3")
cms3pdo <- cms %>% filter(CMS_PDO == "CMS3")
allcms3 <- cms3$model
cases3c <- cms3pdo$model
cms3j <- jaccard(allcms3, cases3c)

cms4 <- cms %>% filter(CMS_PDX == "CMS4")
cms4pdo <- cms %>% filter(CMS_PDO == "CMS4")
allcms4 <- cms4$model
cases4c <- cms4pdo$model
cms4j <- jaccard(allcms4, cases4c)

cmsNC <- cms %>% filter(CMS_PDX == "NC")
cmsNCpdo <- cms %>% filter(CMS_PDO == "NC")
allcmsNC <- cmsNC$model
casesNCc <- cmsNCpdo$model
cmsNCj <- jaccard(allcmsNC, casesNCc)

#cris <- "/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_model_cris-right_validated.tsv"
cris <- read.table(file_cris, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

crisa <- cris %>% filter(CRIS_PDX == "CRIS-A")
crisapdo<- cris %>% filter(CRIS_PDO == "CRIS-A")
allcrisa <- as.character(crisa$model)
crisacon <- as.character(crisapdo$model)
crisaj <- jaccard(allcrisa, crisacon)

crisb <- cris %>% filter(CRIS_PDX == "CRIS-B")
crisbpdo <- cris %>% filter(CRIS_PDO == "CRIS-B")
allcrisb <- as.character(crisb$model)
crisbcon <- as.character(crisbpdo$model)
crisbj <- jaccard(allcrisb, crisbcon)

crisc <- cris %>% filter(CRIS_PDX == "CRIS-C")
criscpdo <- cris %>% filter(CRIS_PDO == "CRIS-C")
allcrisc <- as.character(crisc$model)
crisccon <- as.character(criscpdo$model)
criscj <- jaccard(allcrisc, crisccon)

crisd <- cris %>% filter(CRIS_PDX == "CRIS-D")
crisdpdo <- cris %>% filter(CRIS_PDO == "CRIS-D")
allcrisd <- as.character(crisd$model)
crisdcon <- as.character(crisdpdo$model)
crisdj <- jaccard(allcrisd, crisdcon)

crise <- cris %>% filter(CRIS_PDX == "CRIS-E")
crisepdo <- cris %>% filter(CRIS_PDO == "CRIS-E")
allcrise <- as.character(crise$model)
crisecon <- as.character(crisepdo$model)
crisej <- jaccard(allcrise, crisecon)

cms_class <- c("CMS1", "CMS2", "CMS3", "CMS4", "NC")
jaccard_cms <- c(cms1j, cms2j, cms3j, cms4j, cmsNCj)

jcms <- cbind(cms_class)
jcms <- as.data.frame(cbind(jcms, jaccard_cms))
rownames(jcms) <- jcms$cms_class
jcms$cms_class <- NULL

cris_class <- c("CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E")
jaccard_cris <- c(crisaj, crisbj, criscj, crisdj, crisej)

jcris <- cbind(cris_class)
jcris <- as.data.frame(cbind(jcris, jaccard_cris))
rownames(jcris) <- jcris$cris_class
jcris$cris_class <- NULL

write.table(jcms, file = result_cms, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(jcris, file = result_cris, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
