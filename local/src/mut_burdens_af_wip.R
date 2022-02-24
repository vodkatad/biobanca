library(ggplot2)
setwd('/scratch/trcanmed/biobanca/dataset/V1/targeted')
load('preprocGeneAF_0.02.Rdata')

pdo_mutb <- apply(pdo, 2, function(x){sum(x!=0)})
pdx_mutb <- apply(pdx, 2, function(x){sum(x!=0)})

all(names(pdo_mutb)==names(pdx_mutb))

mutb <- data.frame(row.names=names(pdo_mutb), pdo=pdo_mutb, pdx=pdx_mutb)

ggplot(data=mutb, aes(x=pdx, y=pdo))+geom_point()+geom_smooth(method="lm")+theme_bw()

m1 <- merge(pdx, pdo, all.x = TRUE, all.y = TRUE, by='row.names')
rownames(m1) <- m1$Row.names
m1$Row.names <- NULL
m1[is.na(m1)] <- 0
m1[m1 < thr] <- 0


#
models <- colnames(pdo)

get_common <- function(smodel, data) {
  data <- data.frame(pdx=m1[, paste0(smodel, '.x')], pdo=m1[, paste0(smodel, '.y')], row.names=rownames(m1))
  data <- data[apply(data,1, function(x) {all(x!=0)}),]
  return(data.frame(pdo=data$pdo, pdx=data$pdx, smodel=rep(smodel, nrow(data))))
}

common <- lapply(models, get_common, m1)

commonall <- do.call(rbind, common)

ex <- commonall[commonall$smodel == "CRC1272",]
ex$mut <- seq(1,nrow(ex))

exx <- melt(ex, id.vars="mut", measure.vars=c('pdo','pdx'))
ggplot(data=exx, aes(x=variable, y=value))+geom_point()+geom_line(aes(group=mut))+theme_bw()

### kras only
kras <- m1[rownames(m1)=="chr12:25245350:C:T", ]
n <- colnames(m1)
kk <- kras!=0
kras <- kras[kk]
names(kras) <- n[kk]

dd <- data.frame(af=kras)
rownames(dd) <- names(kras)
dd$smodel <- substr(rownames(dd), 0, 7)
dd$type <-  substr(rownames(dd), 9, 9)
dd$class <- ifelse(dd$type == 'x', 'pdx', 'pdo')

ggplot(data=dd, aes(x=class, y=af))+geom_point()+geom_line(aes(group=smodel))+theme_bw()

library(reshape)
long <- melt(commonall, id.vars="smodel")
ggplot(data=long, aes(x=smodel, fill=variable, y=value))+geom_boxplot()
## /mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdo/mutect/merged.table_nomultiallele
all <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdo/mutect/merged.table_nomultiallele', sep= "\t", header=T, row.names=1)