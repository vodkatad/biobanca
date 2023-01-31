# normal del
PD <- c(17,10)
PR <- c(5,1)
SD <- c(18,6)

data <- rbind(PD, PR, SD)
colnames(data) <- c('Normal ploidy', 'Deletion')

fisher.test(data)

data2 <- data[-3,]
fisher.test(data2)

tdata <- t(data)
dataf <- t(t(tdata)/colSums(tdata))
colSums(dataf)# check all 1
pd <- melt(dataf)
colnames(pd) <-c('chr19q13.31','recist', 'fraction')
pdm <- melt(data)
colnames(pdm) <- c('recist', 'chr19q13.31', 'n')
pd$id <- paste0(pd$recist, pd$chr19q13.31)
pdm$id <- paste0(pdm$recist, pdm$chr19q13.31)
pdm$recist <- NULL
pdm$chr19q13.31 <- NULL
pdd <- merge(pd, pdm, by="id")

pdd$recist <- factor(pdd$recist, levels=c('PD', 'SD', 'PR'))
ggplot(data=pdd, aes(y=fraction, x=factor(1),fill=as.factor(chr19q13.31)))+geom_col()+theme_bw()+facet_wrap(~recist)+
  coord_polar(theta = "y")+scale_fill_manual(values=c('blue', 'darkgrey'))

pdd2 <- pdd[pdd$recist != "SD",]
ggplot(data=pdd2, aes(y=fraction, x=factor(1),fill=chr19q13.31))+geom_col()+theme_void()+facet_wrap(~recist)+
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))+
  coord_polar(theta = "y", start=0)+scale_fill_manual(values=c('blue', 'darkgrey'))+xlab("")+
  theme(text = element_text(size = 20))     
ggsave('~/chr19_chemio.pdf')

PD <- c(5,17,5)
PR <- c(0,3,3)
SD <- c(6,8,10)
data <- rbind(PD, PR, SD)
colnames(data) <- c('Deletion', 'WT', 'Amplification')
fisher.test(data)

 
data2 <- data[-3,]
fisher.test(data2)

library(ggplot2)
library(reshape)
tdata <- t(data)
dataf <- t(t(tdata)/colSums(tdata))
colSums(dataf)# check all 1
pd <- melt(dataf)

colnames(pd) <-c('POLD1_status','recist', 'fraction')
pd$recist <- factor(pd$recist, levels=c('PD', 'SD', 'PR'))
ggplot(data=pd, aes(y=fraction, x=factor(1),fill=as.factor(POLD1_status)))+geom_col()+theme_bw()+facet_wrap(~recist)+
  coord_polar(theta = "y")+scale_fill_manual(values=c('red', 'blue', 'darkgrey')

                                                                                                                                                                