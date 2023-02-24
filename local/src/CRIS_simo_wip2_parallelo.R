### CRIS with simo coupled

library(tidyverse)

simo_in <- snakemake@input[["simo"]] # /scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/SIMO_v2.tsv
lmx_in <- snakemake@input[["pdx"]] # /scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMX_BASALE_prediction_result_nc.tsvq
lmo_in <- snakemake@input[["pdo"]] # /scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris/vsd_cris_LMO_BASALE_prediction_result_nc.tsv
log <- snakemake@log[[1]]

sink(log)
print('#Fasta like file with the info on all the model with replicates and the choices that we made for them')
sink()

save.image('lamelmainsegue.Rdata')

simo <- read.table(simo_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmx <- read.table(lmx_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmo <- read.table(lmo_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Non passando mai ad avere i genealogy in colonne non faccio gsub e tengo i .2 come sono nei file di cris di partenza.

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

# We can safely merge on model only because we removed the wrong lineages with the previous merge
# with simo pairs using the whole lineage.
lmx_m$model <- substr(lmx_m$PDX_lineage, 0, 7)
lmo_m$model <- substr(lmo_m$PDO_lineage, 0, 7)

simo_df <- merge(lmo_m, lmx_m, by = "model")

smodel <- unique(simo_df$model)

choose_lmo <- function(lineages, log) {
  passages <- as.numeric(substr(lineages, 15, 17)) # CRC0051LMO0A02[006]001R01000
  ordered <- lineages[order(passages)]
  smodel <- unique(substr(lineages, 0, 7))
  sink(log, append=TRUE)
  print(paste0('>LMO ', smodel))
  print(ordered)
  print(passages)
  print(length(unique(passages)))
  sink()
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

choose_lmx <- function(lineages, simo_asked, log) {
  lwanted <- unlist(strsplit(simo_asked, ''))
  llineages <- strsplit(lineages, '')
  distances <- sapply(llineages, get_delta, lwanted) 
  names(distances) <- lineages
  distances <- distances[order(-distances)]
  smodel <- unique(substr(lineages, 0, 7))
  sink(log, append=TRUE)
  print(paste0('>LMX ', smodel))
  print(distances)
  sink()
  return(names(distances[1]))
}

count_right_gen <- function(lineages, simo_asked) {
  gen <- as.numeric(substr(unique(lineages), 13, 14))
  genw <- as.numeric(substr(simo_asked, 13, 14))
  if (gen == genw) {
    return(1)
  } else {
    return(0)
  }
}

right_right <- data.frame(stringsAsFactors = FALSE)
right <- data.frame(stringsAsFactors = FALSE)
right_gen <- 0 # how many xeno do we have with the same gen as the one chosen for DNA?
for (i in seq(1, length(smodel))) {
  tmp <- simo_df[simo_df$model == smodel[i],]
  if (nrow(tmp) > 1) {
    # we have a set of replicates that we need to evaluate
    # we need to choose between them
    # if we have more than 1 PDO...the _passage_ : CRC0051LMO0A02[006]001R01000 -> queste tre cifre sono il passaggio (NON LO SCONGELAEMNTO DIOMAJALE)
    if (length(unique(tmp$PDO_long_lineage)) != 1) {
      chosen_lmo <- choose_lmo(tmp$PDO_long_lineage, log)
    } else {
      chosen_lmo <- tmp$PDO_long_lineage[1]
    }
    checklmx <- unique(tmp$MATCHED_LMX.LMO.y) 
    if (length(checklmx) != 1) { # un controllo in più, lo xeno di riferimento del file di simo per ogni modello deve cmq essere uno solo
      print(tmp)
      stop('Problems lmx lineage!')
    }
    right_gen <- right_gen + count_right_gen(tmp$PDX_long_lineage, checklmx) # it's always only 1 the unique so I can do this right now
    if (length(unique(tmp$PDX_long_lineage)) != 1) {
      chosen_lmx <- choose_lmx(tmp$PDX_long_lineage, checklmx, log)
    } else  {
      chosen_lmx <- tmp$PDX_long_lineage[1]
    }
    chosen <- c(chosen_lmo, chosen_lmx)
    chosen_df <- data.frame(PDO_long_lineage=chosen[1], PDX_long_lineage=chosen[2])
    right_right <- rbind(right_right, chosen_df)
    forcrisright <- tmp[tmp$PDO_long_lineage==chosen[1] & tmp$PDX_long_lineage==chosen[2],]
    if (nrow(forcrisright) != 1) {
      print(tmp)
      print(chosen)
      stop('Problems chosen lineages!')
    }
    right <- rbind(right, forcrisright)
    # to keep right output == to the previous one also for the 12  "wrong" that now we keep I select from tmp the line with the chosen genealogy
  } else {
    right_right <- rbind(right_right, tmp[, c('PDO_long_lineage', 'PDX_long_lineage')])
    right <- rbind(right, tmp)
  }
}

right_right$LMO_lineage <- substr(right_right$PDO_long_lineage, 0, 12)
right_right$LMX_lineage <- substr(right_right$PDX_long_lineage, 0, 12)

write.table(right_right, snakemake@output[["lineage"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(right, snakemake@output[["right"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

sink(log, append=TRUE)
print('Xeno right gen')
print(right_gen)
sink()
#egrassi@godot:/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris$ wc  -l vsd_model_cris-right.tsv 
#114 vsd_model_cris-right.tsv
#egrassi@godot:/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris$ wc -l vsd_model_cris-right_parallelo.tsv 
#85 vsd_model_cris-right_parallelo.tsv

# this apparent problem is due to duplicates that were put the the not parallelo script in right, if I reduce them to unique smodel cris cris we seem to be fine:
#egrassi@godot:/scratch/trcanmed/biobanca/dataset/V1/trans_sign/cris$ sort vsd_model_cris-right.tsv | cut -f 1,6,11 |sort | uniq  | wc -l
#73

# 73+12 = 85, 12 was the number of 'wrong' previously

# controlli che non ci sia cmq mai scelta per gli xeno:
#check_x <- function(model, data) {
#  sub <- data[data$model == model,]
#  length(unique(sub$PDX_long_lineage))
#}

#sss <- sapply(smodel, check_x, simo_df)
#sss[sss!=1]