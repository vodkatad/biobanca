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
allo <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdo/mutect/merged.table_nomultiallele', sep= "\t", header=T, row.names=1)
allx <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdx/mutect/merged.table_nomultiallele', sep= "\t", header=T, row.names=1)

m <- merge(allo, allx, all.x=TRUE, all.y=TRUE, by="row.names")
m$Row.names <- NULL
m[is.na(m)] <- 0
#m[m < 0.3] <- 0


data <- data.frame(o=m$CRC1888LMO0A04009001D02000, x=m$CRC1888LMX0B01001TUMD05000)
data[apply(data, 1, function(x) {any(x!=0)}),]

## compare CRC1888 and CRC1307
data <- data.frame(o=m$CRC1888LMO0A04009001D02000, x=m$CRC1888LMX0B01001TUMD05000)
d <- data[apply(data, 1, function(x) {any(x!=0)}),]

jac <- proxy::simil(d, by_rows = FALSE, method = "Jaccard")


data <- data.frame(o=m$CRC1307LMO0A02008004D02000, x=m$CRC1307LMX0A01001TUMD05000)
d <- data[apply(data, 1, function(x) {any(x!=0)}),]

jac <- proxy::simil(d, by_rows = FALSE, method = "Jaccard")



####


mx <- m[,grepl('LMX', colnames(m))]
mo <- m[,grepl('LMO', colnames(m))]
colnames(mx) <- substr(colnames(mx), 0, 7)
colnames(mo) <- substr(colnames(mo), 0, 7)

com <- intersect(colnames(mx), colnames(mo))
mx <- mx[, com]
mo <- mo[, com]

jac <- proxy::simil(mx, mo, by_rows = FALSE, method = "Jaccard")

which.max(jac[,colnames(jac)=='CRC1888'])

pearson <- jac
diag <- diag(pearson)
pearson2 <- pearson
diag(pearson2) <- rep(NA, length(diag))
all <- as.numeric(unlist(pearson2))
all <- all[!is.na(all)]
#all <- upper.tri(pearson, diag = FALSE) # this is not a simmetric matrix!
pdata <- data.frame(pearson=c(all, diag), type=c(rep('unmatched', length(all)),rep('matched', length(diag))))

ggplot(data=pdata, aes(x=pearson, color=type))+geom_density()+geom_density(size=1.5)+scale_color_manual(values=c('#004D40','#FFC107'))+xlab('jaccard')+theme_bw()+theme(text=element_text(size=10))

ggplot(data=pdata[pdata$type=="matched",], aes(x=pearson, fill=type))+geom_histogram(alpha=0.7, position="dodge")+scale_fill_manual(values=c('#004D40','#FFC107'))+xlab('jaccard')+theme_bw()+theme(text=element_text(size=20))
ggplot(data=pdata[pdata$type=="unmatched",], aes(x=pearson, fill=type))+geom_histogram(alpha=0.7, position="dodge")+scale_fill_manual(values=c('#004D40','#FFC107'))+xlab('jaccard')+theme_bw()+theme(text=element_text(size=20))


########### gene level comparisons
load('/scratch/trcanmed/biobanca/dataset/V1/targeted/preprocGeneAF_0.05.Rdata')

# select mutations present in at least 5 samples for both o and x
muts_info <- unique(rbind(gpdo, gpdx))

m1 <- merge(pdx, pdo, all.x = TRUE, all.y = TRUE, by='row.names')
rownames(m1) <- m1$Row.names
m1$Row.names <- NULL
m1[is.na(m1)] <- 0
m1[m1 < thr] <- 0

count_muts <- function(muts, models, n_thr) {
  x <- grepl('.x', models, fixed=T)
  o <- grepl('.y', models, fixed=T)
  x_n <- sum(muts[x]!=0)
  o_n <- sum(muts[o]!=0)
  return (x_n>=n_thr && o_n >= n_thr)
}

wanted_muts <- apply(m1, 1, count_muts, colnames(m1), 5)
wanted_muts_id <- rownames(m1)[wanted_muts]
muts_info[rownames(muts_info) %in% wanted_muts_id,]


## ok we need to collapse on genes cause we have too few on APC/TP53/KRAS instead and working
# separately does not make too much sense. iterate on all genes and select the ones with muts common o-x in 5 models
# iterate on genes
 # for all its muts select the common ones
  # if more than 5 proceed
compare_o_x <- function(gene, mut_table, gene_muts, n_thr) {
  my_muts <- gene_muts[gene_muts$genes == gene, , drop=FALSE]
  common_x_o <- sapply(unique(rownames(my_muts)), compare_o_x_one_mut, mut_table)
  # we sum over columns to get the n of models mutated both in x o and if the overall sum if >= 5 we can go on
  total_muts <- colSums(common_x_o)
  if (sum(total_muts) >= n_thr) {
    wanted <- names(total_muts[total_muts!=0]) # we are still not selecting for both here TODO FIXME
    # only the muts same in x-o
    data <- mut_table[rownames(mut_table) %in% wanted, ]
    # only the mutated samples
    data[, apply(data, 2, function(x){any(x!=0)})]
    data$id <- rownames(data)
    long <- melt(data, id.vars="id")
    long <- long[long$value != 0,]
    long$type <-  substr(long$variable, 9, 9)
    long$class <- ifelse(long$type == 'x', 'pdx', 'pdo')
    long$smodel <- substr(long$variable, 0, 7)
    long$group <- paste0(long$id, long$smodel)
    #ggplot(data=long, aes(x=class, y=value))+geom_point()+geom_line(aes(group=group))+theme_bw()
     
    # we can remove what's not in both here with table on group - ggplot does that himself?
    
    ddd <- as.data.frame(table(long$group))
    w <- ddd[ddd$Freq==2, 'Var1']
    long <- long[long$group %in% w,]
    #ggplot(data=long, aes(x=class, y=value))+geom_point()+geom_line(aes(group=group))+theme_bw()
    
    #2287 chr12:25245350:C:A CRC1090.y 0.448    y   pdo CRC1090 chr12:25245350:C:ACRC1090 ??
    
    longx <- long[long$class == "pdx",]
    longo <- long[long$class == "pdo",]
    longx <- longx[order(longx$id, longx$smodel),]
    longo <- longo[order(longo$id, longo$smodel),]
    all(longx$group==longo$group)
    ttest <- t.test(longx$value, longo$value, paired=TRUE)
    #oo <- long$value > 0.9 & long$class=="pdo"
    #o <- long$value < 0.7 & long$class=="pdx"
    #ooo <- intersect(long[o,'group'], long[oo,'group'])
    #long[long$group %in% ooo,]
    return(ttest$p.value)
  } else {
    return(NA)
  }
}

compare_o_x_one_mut <- function(mut, mut_table) {
  x <- mut_table[rownames(mut_table)==mut, grepl('.x', colnames(mut_table), fixed=TRUE)]
  y <- mut_table[rownames(mut_table)==mut, grepl('.y', colnames(mut_table), fixed=TRUE)]
  common <- x != 0 & y != 0
  return(common)
}
muts_info$genes <- as.character(muts_info$genes)
try <- sapply(unique(muts_info$genes), compare_o_x, m1, muts_info, 5)



### work on one mut only only
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

