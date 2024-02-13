setwd('/scratch/trcanmed/biobanca/dataset/V1/shallowseq')
load('validated_pearson_xo.Rdata')

long <- read.table('validated_heatmap_cors_xo.tsv', header=T, stringsAsFactors = F)
long$smodel <- substr(long$gen, 0, 7)
long <- long[grepl('LMO', long$gen),]

stopifnot(ncol(pdo_df)==125)

long <- long[long$smodel %in% colnames(pdo_df),]
stopifnot(nrow(long)==125)

long$passage <- as.numeric(substr(long$gen, 15, 17))
table(long$passage)

long$group <- ifelse(long$passage <= 3, 'early', ifelse(long$passage <=8, 'mid', 'late'))

dfc <- data.frame(cor=diag, smodel=names(diag))
m <- merge(dfc, long, by="smodel")

m$group <- factor(m$group, levels=c('early', 'mid', 'late'))
ggplot(data=m, aes(y=cor, x=group, fill=group))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw()+ylab('Pearson')

summary(aov(formula=as.formula("cor~group"), data=m))

kruskal.test(formula=as.formula("cor~group"), data=m)
#### jaccards muts early late
setwd('/scratch/trcanmed/biobanca/dataset/V1/targeted')
jac <- read.table('jac_matrix_0.05.tsv', header=T, stringsAsFactors = F)

long <- read.table('pdo_longgen_after_xo_filtering.tsv', header=F, stringsAsFactors = F)
colnames(long) <- 'gen'
long$smodel <- substr(long$gen, 0, 7)

stopifnot(ncol(jac)==124) # - 1 cause 1 without any muts

long <- long[long$smodel %in% colnames(jac),]
stopifnot(nrow(long)==124)

long$passage <- as.numeric(substr(long$gen, 15, 17))
table(long$passage)

long$group <- ifelse(long$passage <= 3, 'early', ifelse(long$passage <=8, 'mid', 'late'))
long$group <- factor(long$group, levels=c('early', 'mid', 'late'))
table(long$group)

diag <- diag(as.matrix(jac))
dfc <- data.frame(cor=diag, smodel=names(diag))
m <- merge(dfc, long, by="smodel")
ggplot(data=m, aes(y=cor, x=group, fill=group))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw()+ylab('Pearson')

summary(aov(formula=as.formula("cor~group"), data=m))

kruskal.test(formula=as.formula("cor~group"), data=m)
### jaccards
all <- read.table('/scratch/trcanmed/biobanca/dataset/V1/WES_early_ox_all/jac_matrix_0.05.tsv' , sep="\t", stringsAsFactors = F, header=T)
w <- read.table('/scratch/trcanmed/biobanca/dataset/V1/WES_early_ox/jac_matrix_0.05.tsv' , sep="\t", stringsAsFactors = F, header=T)

da <- diag(as.matrix(all))
dw <- diag(as.matrix(w))

ma <-data.frame(a=da)
mw <- data.frame(w=dw)

mm <- merge(ma, mw, by="row.names")


ggplot(data=mm, aes(x=a, y=w))+geom_point()+theme_bw(base_size=20)+xlab('All tiers')+ylab('Tiers 1-2-3')+geom_abline(intercept=0, slope=1)+xlim(0,1)+ylim(0,1)

mlo <- melt(mm, id="Row.names")
pd <- position_dodge(width=0.2)
ggplot(data=mlo, aes(x=variable, y=value))+geom_boxplot(outlier.shape=NA)+geom_jitter(aes(group=`Row.names`), position=pd)+
  geom_line(aes(group=`Row.names`), position=pd)+
  theme_bw(base_size=20)
  
