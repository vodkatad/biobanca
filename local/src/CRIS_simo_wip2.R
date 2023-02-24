### CRIS with simo coupled

library(tidyverse)

simo_in <- snakemake@input[["simo"]] # /scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/SIMO_v2.tsv
lmx_in <- snakemake@input[["pdx"]] # /scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_prediction_result_nc
lmo_in <- snakemake@input[["pdo"]] # /scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMO_BASALE_prediction_result_nc

# for testing and evaluating the code, 22/02/23
#simo_in <- '/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/SIMO_v2.tsv'
#lmx_in <- '/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_prediction_result_nc'
#lmo_in <- '/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMO_BASALE_prediction_result_nc'

simo <- read.table(simo_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmx <- read.table(lmx_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmo <- read.table(lmo_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#simo[,c(2,3)] <- NULL
lmx <- lmx[,c(1,2)]
names(lmx) <- c("PDX_long_lineage", "PDX")
lmo <- lmo[,c(1,2)]
names(lmo) <- c("PDO_long_lineage", "PDO")
lmx$PDX_lineage <- substr(lmx$PDX_long_lineage, 0, 12)
lmo$PDO_lineage <- substr(lmo$PDO_long_lineage, 0, 12)

lmx_m <- merge(simo, lmx, by="PDX_lineage")
lmx_m[,c(2)] <- NULL
lmo_m <- merge(simo, lmo, by="PDO_lineage")
lmo_m[,c(2)] <- NULL


lmx_m$model <- substr(lmx_m$PDX_lineage, 0, 7)
lmo_m$model <- substr(lmo_m$PDO_lineage, 0, 7)

simo_df <- merge(lmo_m, lmx_m, by = "model")
#model_cris <- simo_df[,c(1,5,3)]

#save.image("pippo.RData")

wrong_pdo <- c()
wrong_pdx <- c()
smodel <- unique(simo_df$model)
right <- data.frame(stringsAsFactors = FALSE)
wrong <- data.frame(stringsAsFactors = FALSE)


for (i in seq(1, length(smodel))) {
  tmp <- simo_df[simo_df$model == smodel[i],]
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


choose_lmo <- function(lineages) {
  passages <- as.numeric(substr(lineages, 13, 14))
  ordered <- lineages[order(passages)]
  return(ordered[1])
}


get_delta <- function(vec1, vec2) {
  # two vec have the same length
  same <- TRUE
  i <- 1
  while (same && i <= length(vec1)) {
    if (vec1[i] != vec2[i]) {
      same <- FALSE
    }
    i <- i + 1
  }
  if (same) {
    return(length(vec1)+1)
  } else {
    return(i-1)
  }
}

choose_lmx <- function(lineages, simo_asked) {
  lwanted <- unlist(strsplit(simo_asked, ''))
  llineages <- strsplit(lineages, '')
  distances <- sapply(llineages, get_delta, lwanted) 
  names(distances) <- lineages
  distances <- distances[order(-distances)]
  return(names(distances[1]))
}


smodel_right <-  unique(right$model)
right_right <- data.frame(stringsAsFactors = FALSE)
for (i in seq(1, length(smodel_right))) {
  tmp <- right[right$model == smodel_right[i],]
  if ( nrow(tmp) > 1 ) {
    if ( length(unique(tmp$PDO)) == 1 && length(unique(tmp$PDX)) == 1 ) { # we have a set of replicates with the same predicted CRIS class
      # we need to choose between them
      # if we have more than 1 PDO...the lowest freeze/thaw cycle
      
      #if (length(unique(tmp$PDO_long_lineage)) > 1 && length(unique(tmp$PDX_long_lineage)) == 1 ) {
      #  chosen_lmo <- choose_lmo(tmp$PDO_long_lineage)
      #  chosen <- c(chosen_lmo, unique(tmp$PDX_long_lineage))
      # if we have more than 1 PDX
      #} else if (length(unique(tmp$PDX_long_lineage)) > 1 && length(unique(tmp$PDO_long_lineage)) == 1) {
      #  print('never never?')
      #  chosen_lmx <- choose_lmx(tmp$PDX_long_lineage, unique(tmp$MATCHED_LMX.LMO.y))
      #  chosen <- c(unique(tmp$PDO_long_lineage), chosen_lmx)
      #} else {
        # dramma vero!
      #print('never?')
      chosen_lmo <- choose_lmo(tmp$PDO_long_lineage)
      chosen_lmx <- choose_lmx(tmp$PDX_long_lineage, unique(tmp$MATCHED_LMX.LMO.y))
      chosen <- c(chosen_lmo, chosen_lmx)
      chosen_df <- data.frame(PDO_long_lineage=chosen[1], PDX_long_lineage=chosen[2])
      right_right <- rbind(right_right, chosen_df)
    }
  } else {
    right_right <- rbind(right_right, tmp[, c('PDO_long_lineage', 'PDX_long_lineage')])
  }
}

## prova di codice non funzionale perché non sono tutte coppie, quindi gli indici sbarellano
# keep <- data.frame()
# single <- data.frame()
# i = 2
# for (i in seq(2, nrow(right))) {
#   if (right[i,'model'] == right[(i-1), "model"]) { 
#     if (right[i, "PDO_passage"] <= right[(i-1), "PDO_passage"]) {
#       keep <- rbind(keep, right[i,])
#     } else {
#       keep <- rbind(keep, right[(i-1),])
#     }
#   } else {
#     single <- rbind(single, right[i,])
#   }
# }

# if (!(right[i, "model"]) %in% keep[i,"model"]) {
#   keep <- rbind(keep, right[i,])
# } else {
#   double <- rbind(double, right[i,])
# else if (right[i,'model'] != right[(i-1), "model"] && !(right[i, "model"] %in% keep[i,])) {
#   single <- rbind(single, right[i,])

colnames(right) <- c("model", "CRIS_PDX", "CRIS_PDO")
colnames(wrong) <- c("model", "CRIS_PDX", "CRIS_PDO")
right_right$LMO_lineage <- substr(right_right$PDO_long_lineage, 0, 12)
right_right$LMX_lineage <- substr(right_right$PDX_long_lineage, 0, 12)

sink(snakemake@log[[1]])
print(paste0("wrong pdo: ", length(wrong_pdo)))
print(paste0("wrong pdx: ", length(wrong_pdx)))
sink()

write.table(right, snakemake@output[["right"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(wrong, snakemake@output[["wrong"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(right_right, snakemake@output[["lineage"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)