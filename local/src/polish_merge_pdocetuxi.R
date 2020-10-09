#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = T)
m1 <- args[1]
m2 <- args[2]
output <- args[3]
prefix <- args[4]

pdo_ctg <- read.table(m1, header=TRUE, sep="\t", na.strings="", stringsAsFactors = FALSE, quote="")
pdo_imaging <- read.table(m2, header=TRUE, sep="\t", na.strings="", stringsAsFactors = FALSE, quote="")
# rename CRC0177
# we keep the mutated one because PDX data are on the mutated one.
# Simo's mail 17 March 2020
#    - se sai per il CRC0177 quale sia il dato che c'e` in vivo, dal tuo excel direi per l'EGFR mut, e` corretto?
#Yep. PD il mutato, no info sul wt. Ricordo ai tempi di aver chiesto check al Galimi,, che sa chi è mutato dove.
remove_ctg <- which(pdo_ctg$CASE=="CRC0177 EGFR wt") 
if (length(remove_ctg) != 1) {
    stop('I wound more than one CRC0177 WT?')
}
pdo_ctg <- pdo_ctg[-remove_ctg,]
pdo_ctg[pdo_ctg$CASE=="CRC0177 EGFR mut",'CASE'] <- "CRC0177"

remove_imaging <- which(pdo_imaging$CASE=="CRC0177" & pdo_imaging$xeno.mut=='4PLE')
if (length(remove_imaging) != 1) {
    stop('I wound more than one CRC0177 WT?')
}
pdo_imaging <- pdo_imaging[-remove_imaging,]
pdo_imaging[pdo_imaging$CASE=="CRC0177" & pdo_imaging$xeno.mut=='EGFR','CASE'] <- "CRC0177"


pdo_ctg <- pdo_ctg[,c(1,4,5,6,7)]
pdo_imaging <- pdo_imaging[,c(1,4,5,6,7)]
merge_pdo <- merge(pdo_ctg, pdo_imaging, by='CASE')
colnames(merge_pdo) <- c('case','CTG_5000',"CTG_1250","CTG_10000","CTG_20000","imaging_5000","imaging_1250","imaging_10000","imaging_20000")

# the -1 is to remove the case columnm definitely not nice to work in this way with R, better to remove it, use as rownames and put back as a column only when calling write.table
tab <- as.data.frame(apply(merge_pdo[,-1], 2, function(x) {sum(is.na(x))}))
colnames(tab) <- ('n_NA')
write.table(data.frame('exp'=rownames(tab), tab), file=paste0('missingdata_', prefix, '.tsv'), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

png(paste0('highlevelcor_', prefix, '.png'))
plot(merge_pdo[,-1])
graphics.off()

#For the first models we remove 10000 cells, to many NAs.
merge_pdo$CTG_10000 <- NULL
merge_pdo$imaging_10000 <- NULL
write.table(merge_pdo, file=output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)