### CRIS with simo coupled

library(tidyverse)

simo_in <- snakemake@input[["simo"]]
lmx_in <- snakemake@input[["pdx"]]
lmo_in <- snakemake@input[["pdo"]]

save.image("linage.R")

simo <- read.table(simo_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmx <- read.table(lmx_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmo <- read.table(lmo_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

simo[,c(2,3)] <- NULL
lmx <- lmx[,c(1,2)]
names(lmx) <- c("PDX_lineage", "PDX")
lmo <- lmo[,c(1,2)]
names(lmo) <- c("PDO_lineage", "PDO")
lmx$PDX_lineage <- substr(lmx$PDX_lineage, 0, 12)
lmo$PDO_lineage <- substr(lmo$PDO_lineage, 0, 12)

lmx_m <- merge(simo, lmx, by="PDX_lineage")
lmx_m[,c(2)] <- NULL
lmo_m <- merge(simo, lmo, by="PDO_lineage")
lmo_m[,c(2)] <- NULL


lmx_m$model <- substr(lmx_m$PDX_lineage, 0, 7)
lmo_m$model <- substr(lmo_m$PDO_lineage, 0, 7)

simo_df <- merge(lmo_m, lmx_m, by = "model")
model_cris <- simo_df[,c(1,5,3)]

save.image("pippo.RData")


wrong_pdo <- c()
wrong_pdx <- c()
smodel <- unique(model_cris$model)
right <- data.frame(stringsAsFactors = FALSE)
wrong <- data.frame(stringsAsFactors = FALSE)


for (i in seq(1, length(smodel))) {
    tmp <- model_cris[model_cris$model == smodel[i],]
    if ( nrow(tmp) > 1 ) {
        if ( length(unique(tmp$PDO)) != 1 ) {
            wrong_pdo <- c(wrong_pdo, smodel[i])
            wrong <- rbind(wrong, tmp)
        }
        if ( length(unique(tmp$PDX)) != 1 ) {
            wrong_pdx <- c(wrong_pdx, smodel[i])
            wrong <- rbind(wrong, tmp)
        }
        if ( length(unique(tmp$PDO)) == 1 && length(unique(tmp$PDX)) == 1 ) {
            right <- rbind(right, unique(tmp))
        }
    } else  {
        right <- rbind(right, tmp)
    } 
}

colnames(right) <- c("model", "CRIS_PDX", "CRIS_PDO")
colnames(wrong) <- c("model", "CRIS_PDX", "CRIS_PDO")

sink(snakemake@log[[1]])
print(paste0("wrong pdo: ", length(wrong_pdo)))
print(paste0("wrong pdx: ", length(wrong_pdx)))
sink()

write.table(right, snakemake@output[["right"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(wrong, snakemake@output[["wrong"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
