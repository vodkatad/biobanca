setwd('~/work/biobanca/tgi/')
d <- read.table("GTI Biobanca_Eugy.txt", sep="\t", header=TRUE, comment.char="", stringsAsFactors = FALSE)
dbck <- d

WCOL = 'TGI..Median.'
#vw3 <- read.table("")

# 776-799 are repeated entries?
# 416-427
## need to fill in missing GenID

d <- d[-seq(775, 798),] # repeated CRC0177 friends
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

all(d$GenID != "")

# something is ko here (because we had repeated values!!!)
nMeas <- sapply(unique(d$GenID), function(x) { dd <- d[d$GenID == x,]; return(seq(1, nrow(dd))) })
d$nMeas <- unlist(nMeas)

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


all(d$nMeas == d$nMeas2)

df <- as.data.frame(table(d$GenID))
df[df$Freq != 6,]

# salvare file intermedio filled per poter filtrare i dati di crescita dei NT
# passati a Francesca con query LAS secondo data di inizio terapia segnata qui
# del 'giusto' exp group
# TODO
write.table(d, file="TGI_filled_GenID_removeddup.tsv", sep="\t", quote=FALSE, row.names=FALSE)

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

has3 <- sapply(unique(d$GenID), has_meas, 3, d)
length(has3[!is.na(has3)])
length(has3)

has2 <- sapply(unique(d$GenID), has_meas, 2, d)
length(has2[!is.na(has2)])

library(ggplot2)
df3 <- data.frame(smodel=names(has3), TGI3w=has3)
df2 <- data.frame(smodel=names(has2), TGI2w=has2)

mdf <- merge(df2, df3, by="smodel", all.x=T)

ggplot(data=mdf, aes(x=reorder(-TGI2w,smodel), y=TGI2w))+geom_col()+theme_bw()+theme(axis.text.x=element_blank())
ggplot(data=mdf, aes(x=reorder(-TGI3w,smodel), y=TGI3w))+geom_col()+theme_bw()+theme(axis.text.x=element_blank())

ggplot(data=mdf, aes(x=TGI3w,y=TGI2w))+geom_point()+theme_bw()

cet <- read.table("cetuxi_3w.tsv", sep="\t")
colnames(cet) <- c('ssmodel', 'w3')

mdf$ssmodel <- substr(mdf$smodel, 0, 7)
mm <- merge(mdf, cet, by="ssmodel")
mm$recist3w <- ifelse(mm$w3 < -50, 'OR', ifelse(mm$w3 > 35, 'PD', 'SD'))
ggplot(data=mm, aes(x=TGI3w,y=TGI2w, color=recist3w))+geom_point()+theme_bw()+xlim(c(-250, 250))+ylim(c(-250, 250))
ggplot(data=mm, aes(x=TGI3w,y=TGI2w, color=recist3w))+geom_point()+theme_bw()
