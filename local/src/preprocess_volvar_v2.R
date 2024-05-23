library(ggplot2)
library(openxlsx)
library(reshape)

tgi_xlsx <- snakemake@input[['xlsx']]
cet_f <- snakemake@input[['cet']]
volvar_f <-snakemake@output[['volvar_ave']]
log_f <- snakemake@log[['log']]
classes <- snakemake@wildcards[['sclass']]

# for the live session to begin with
#setwd('/scratch/trcanmed/biobanca/local/share/data/cetuxi')
#tgi_xlsx <-'aggregati_TGI_fixed0.xlsx'

d <- read.xlsx(tgi_xlsx, sheet="aggregati_all", colNames = TRUE, rowNames = FALSE)
names(d)[names(d) == 'TGI%(Median)'] <- 'TGI_Median'
names(d)[names(d) == 'TGI%(Average)'] <- 'TGI_Average'
names(d)[names(d) == 'Vol.Var.%.PHI'] <- 'volvarNT'
d$GenID <- ifelse(is.na(d$GenID),"",d$GenID)

WCOL = 'volvarNT'


# 776-799 are repeated entries?
# 416-427
d <- d[-seq(763, 774),]# repeated CRC0177 friends 403 414
d <- d[-seq(403, 414),] # repeated CRC0177 friends
d <- d[!d$GenID == "CRC0177LMX0C03",]
d <- d[!d$GenID == "CRC0177LMX0D03",]
d <- d[!d$GenID == "CRC0177LMX0E03",]
d <- d[!d$GenID == "CRC0177LMX0B",]

## eliminate M040SX_CRC0030LMX0A, CRC0166LMXB, CRC1723LMXA
d <- d[d$GenID != "CRC0030LMX0A",]
d <- d[!d$GenID == "CRC0166LMX0B",]
d <- d[!d$Nome.Gruppo.LAS == "CRC1723LMX0A.2018-09-14.0",]

## need to fill in missing GenID in the middle of their blocks
#d[is.na(d$GenID),'GenID'] <- "" # previous loading from tsv had "", with read.xlsx NAs 
previous <- ''
stretch <- c()
for (i in seq(1, nrow(d))) {
  if (d[i,'GenID']=="") {
    if (previous == "") {
      previous <- d[i-1,'GenID']
    }
    stretch <- c(stretch, i)
  } else if (d[i,'GenID']!="") {
    if (previous != "") {
      if (previous == d[i,'GenID']) {
        for (j in seq(1, length(stretch))) {
          d[stretch[j], 'GenID'] <- previous 
        }
      } else {
        print(paste0('in doubt for ', i))
      }
    }
    previous <- ""
    stretch <- c()
  }
}

stopifnot(all(d$GenID != ""))

# measures are in blocks for each model, we put 1 ... n to identify third week
nMeas <- lapply(unique(d$GenID), function(x) { dd <- d[d$GenID == x,]; return(seq(1, nrow(dd))) })
d$nMeas <- unlist(nMeas)

# check with other method that everything is ok for this assignment of days
d$nMeas2 <- rep(NA, nrow(d))
prev <- ""
seqj <- 1
for (i in seq(1, nrow(d))) {
  if (prev != "") {
    if (prev != d[i,'GenID']) {
      seqj <- 1
    } else {
      seqj <- seqj+1
    }
  }
  d$nMeas2[i] <- seqj
  prev <- d[i,'GenID']
}


stopifnot(all(d$nMeas == d$nMeas2))

# checking that we have 6 measures for everybody
df <- as.data.frame(table(d$GenID))
stopifnot(all(df$Freq == 6))

has_meas <- function(gen, wantedMeas, data) {
  d <- data[data$GenID == gen,]
  wanted <- d[d$nMeas==wantedMeas,]
  if (nrow(wanted) > 1) {
    print(paste0("Error! ", gen))
    return(NA)
  }
  if (!is.na(wanted[,WCOL])) {
    return(wanted[,WCOL]) 
  } else {
    return(NA)
  }
}

sink(file=log_f)
print("Not NA and all week 3:")
has3 <- sapply(unique(d$GenID), has_meas, 3, d)
length(has3[!is.na(has3)])
length(has3)
print("Not NA and all week 2:")
has2 <- sapply(unique(d$GenID), has_meas, 2, d)
length(has2[!is.na(has2)])
sink()

df3 <- data.frame(smodel=names(has3), vw3=has3)
df2 <- data.frame(smodel=names(has2), vw2=has2)

mdf <- merge(df2, df3, by="smodel", all.x=T)

sink(log_f, append=TRUE)
cor.test(mdf$vw2, mdf$vw3)
sink()

mdf$ssmodel <- substr(mdf$smodel, 0, 7)
getaverage <- function(ssmodel, data, col) {
  subset <- data[data$ssmodel==ssmodel,]
  subset <- subset[!is.na(subset[,col]),, drop=FALSE]
  return(mean(subset[,col]))
}

ave3w <- sapply(unique(mdf$ssmodel), getaverage, mdf, 'vw3')

res <- as.data.frame(ave3w)

sink(file=log_f, append=TRUE)
print("Not NA and all week 3 for avg models:")
nrow(res)
table(is.na(res[,1]))
sink()

#res <- res[rownames(res) != "CRC0066",, drop=FALSE]
res <- res[!is.na(res[,1]),, drop=FALSE]
write.table(res, volvar_f, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

save.image('cojoni.Rdata')
q(save='no')