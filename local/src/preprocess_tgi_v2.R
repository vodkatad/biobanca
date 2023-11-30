library(ggplot2)
library(openxlsx)
library(reshape)

tgi_xlsx <- snakemake@input[['xlsx']]
cet_f <- snakemake@input[['cet']]
tgi_nodup_f <-snakemake@output[['tgi_nodup']]
tgi_ave_f <-snakemake@output[['tgi_ave']]
log_f <- snakemake@log[['log']]
scattereci_f <- snakemake@output[['scattereci']]
delta_f <- snakemake@output[['delta']]
classes <- snakemake@wildcards[['sclass']]

save.image("pippo.Rdata")

# for the live session to begin with
#setwd('/scratch/trcanmed/biobanca/local/share/data/cetuxi')
#tgi_xlsx <-'aggregati_TGI_fixed0.xlsx'

d <- read.xlsx(tgi_xlsx, sheet="aggregati_all", colNames = TRUE, rowNames = FALSE)
names(d)[names(d) == 'TGI%(Median)'] <- 'TGI_Median'
names(d)[names(d) == 'TGI%(Average)'] <- 'TGI_Average'
d$GenID <- ifelse(is.na(d$GenID),"",d$GenID)

WCOL = classes


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

# salvare file intermedio filled per poter filtrare i dati di crescita dei NT
# passati a Francesca con query LAS secondo data di inizio terapia segnata qui
# del 'giusto' exp group
write.table(d, file=tgi_nodup_f, sep="\t", quote=FALSE, row.names=FALSE)

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

df3 <- data.frame(smodel=names(has3), TGI3w=has3)
df2 <- data.frame(smodel=names(has2), TGI2w=has2)

mdf <- merge(df2, df3, by="smodel", all.x=T)

#ggplot(data=mdf, aes(x=reorder(-TGI2w,smodel), y=TGI2w))+geom_col()+theme_bw()+theme(axis.text.x=element_blank())
#ggplot(data=mdf, aes(x=reorder(-TGI3w,smodel), y=TGI3w))+geom_col()+theme_bw()+theme(axis.text.x=element_blank())

#ggplot(data=mdf, aes(x=TGI3w,y=TGI2w))+geom_point()+theme_bw()
#mdf2 <- mdf[mdf$TGI3w > -1000 & mdf$TGI3w < 2000,]
#ggplot(data=mdf2, aes(x=TGI3w,y=TGI2w))+geom_point()+theme_bw()
sink(log_f, append=TRUE)
cor.test(mdf$TGI2w, mdf$TGI3w)
sink()

#
#cet <- read.table("/scratch/trcanmed/biobanca/dataset/V1/cetuximab/cetuxi_perc_w3.tsv", sep="\t", header=TRUE)
cet <- read.table(cet_f, sep="\t", header=TRUE)
colnames(cet) <- c('ssmodel', 'w3')

mdf$ssmodel <- substr(mdf$smodel, 0, 7)
mm <- merge(mdf, cet, by="ssmodel")
mm$recist3w <- ifelse(mm$w3 < -50, 'OR', ifelse(mm$w3 > 35, 'PD', 'SD'))
ggplot(data=mm, aes(x=TGI3w,y=TGI2w, color=recist3w))+geom_point()+theme_bw()#+xlim(c(-500, 500))+ylim(c(-500, 500))
ggsave(scattereci_f)
#ggplot(data=mm, aes(x=TGI3w,y=TGI2w, color=recist3w))+geom_point()+theme_bw()

## TODO possible to write.table here to investigate exceptions recist-TGI

ssmodels <- unique(mdf$ssmodel)

getdelta <- function(ssmodel, data, col) {
  subset <- data[data$ssmodel==ssmodel,]
  return(abs(subset[1,col]-subset[2,col]))
}

smodelsta <- as.data.frame(table(mdf$ssmodel))
dup <- smodelsta[smodelsta$Freq == 2, 'Var1']
# 177 and CRC1432 are the only > 2, we'll manage them separately TODO FIXME
deltas2w <- sapply(dup, getdelta, mdf, 'TGI2w')
deltas3w <- sapply(dup, getdelta, mdf, 'TGI3w')

delta <- data.frame(ssmodel=dup, delta3w=deltas3w, delta2w=deltas2w)

m2 <- merge(delta, mdf, by="ssmodel")
#m2 <- m2[m2$ssmodel != "CRC0066",]
m2$arm <- as.factor(substr(m2$smodel, 12,12))
ggplot(data=m2, aes(x=reorder(ssmodel, -delta3w), y=TGI3w, color=arm))+
  geom_point()+
  theme_bw()+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(delta_f)


# discuss vs 2w?
#ggplot(data=m2, aes(x=reorder(ssmodel, -delta2w), y=TGI2w, color=arm))+
#  geom_point()+
#  theme_bw()+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



getaverage <- function(ssmodel, data, col) {
  subset <- data[data$ssmodel==ssmodel,]
  subset <- subset[!is.na(subset[,col]),, drop=FALSE]
  return(mean(subset[,col]))
}

ave3w <- sapply(unique(mdf$ssmodel), getaverage, mdf, 'TGI3w')

res <- as.data.frame(ave3w)

sink(file=log_f, append=TRUE)
print("Not NA and all week 3 for avg models:")
nrow(res)
table(is.na(res[,1]))
sink()

#res <- res[rownames(res) != "CRC0066",, drop=FALSE]
res <- res[!is.na(res[,1]),, drop=FALSE]
write.table(res, tgi_ave_f, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

q(save='no')


getAB <- function(ssmodel, data, col) {
  subset <- data[data$ssmodel==ssmodel,]
  subset <- subset[!is.na(subset[,col]),, drop=FALSE]
  return(as.data.frame(table(substr(subset$smodel, 12, 13))))
}

s <- sapply(unique(mdf$ssmodel), getAB, mdf, 'TGI3w')
