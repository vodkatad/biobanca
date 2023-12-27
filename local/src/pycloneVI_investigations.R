dnobin <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_subs/miao4', sep="\t", header=FALSE, stringsAsFactors = F)
# thr 0.05 miao 4, miao 2 no thr
#d <- read.table('bestbet_bin_0.12_0.24.tsv', sep="\t", header=TRUE)
colnames(dnobin) <- c('sample', 'nclone')
rownames(dnobin) <- sapply(strsplit(dnobin$sample, '/'), '[[', 4)
#d$lmodel <- substr(rownames(d), 0, 10)
dnobin$lmodel <- substr(rownames(dnobin), 0, 10)

#mm <- merge(d, dnobin, by="lmodel")

d <- dnobin

d$smodel <- substr(rownames(d), 0, 7)
d$mp <- substr(rownames(d), 8, 10)

d$mp <- factor(d$mp, levels=c('PRX', 'LMX'))


dd <- as.data.frame(table(d$smodel))
d2 <- d[d$smodel %in% dd[dd$Freq == 2,'Var1'],]

ggplot(data=d2, aes(x=smodel, y=nclone, fill=mp))+
  geom_col(position="dodge")+theme_bw()+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Patient')+ylab('nclone')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title=""))

lar <- function(x, data) {
  sub <- data[data$smodel == x,]
  sub[sub$mp=="LMX",'nclone'] - sub[sub$mp=="PRX",'nclone']
}


meslo <- function(x, data) {
  sub <- data[data$smodel == x,]
  sub[sub$mp=="LMX",'nclone']
}


dd2 <- as.data.frame(sapply(unique(d2$smodel), lar, d2))
dd2 <- dd2[order(dd2[,1]),, drop=FALSE]

dd3 <- as.data.frame(sapply(unique(d2$smodel), meslo, d2))
dd3 <- dd3[order(dd3[,1]),, drop=FALSE]



d3 <- merge(dd3, d2, by.x="row.names", by.y="smodel")
colnames(d3)[2] <- 'o'
colnames(d3)[1] <- 'smodel'
ggplot(data=d3, aes(x=reorder(smodel,o), y=nclone, fill=mp))+
  geom_col(position="dodge")+theme_bw()+
  theme(text=element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Patient')+ylab('N.Clones')+
  scale_fill_manual(values=c('#adacac', '#595959'))+
  guides(fill=guide_legend(title=""))

pp <- d3[d3$mp=="LMX",]
pp$passage <- as.numeric(substr(pp$sample, 39, 40))
table(pp$passage, pp$nclone)

ggplot(data=pp, aes(x=nclone))+geom_histogram()+theme_bw()+facet_wrap(~passage)

pp <- d3[d3$mp=="PRX",]
pp$passage <- as.numeric(substr(pp$sample, 39, 40))
table(pp$passage, pp$nclone)
ggplot(data=pp, aes(x=nclone))+geom_histogram()+theme_bw()+facet_wrap(~passage)

# bck is dd3 from MR
m <- merge(pp, bck, by="smodel")
ggplot(data=m, aes(x=intercept, y=nclone))+geom_point()+theme_bw()+xlab('effectiveMR')


#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ zcat pyclone/C*.tsv.gz   |head -n1 | cut -f3,4 > o 
#egrassi@godot:/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot$ for f in pyclone/C*.tsv.gz ;  do zcat $f | sed 1d | cut -f 3,4 |sort |uniq  ; done >> o
o <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_subs/o', sep="\t", header=T)
ggplot(data=o, aes(cellular_prevalence))+geom_histogram()+theme_bw()


m <- merge(dd3, bck, by="row.names")
colnames(m) <-  c('s', 'THR005', 'noTHR')
ggplot(data=m, aes(x=noTHR, y=THR005))+geom_point()+theme_bw()
m$d <- m$noTHR-m$THR005

ggplot(data=m, aes(d))+geom_histogram()+theme_bw()

