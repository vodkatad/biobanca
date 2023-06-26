wanted <- c('CRC1251','CRC0729','CRC0370')
d <- read.table('/scratch/trcanmed/biobanca/local/share/data/shallowseq/gistic/gistic_pdo/all_thresholded.by_genes.txt', sep="\t", header=TRUE)

all <- d[d$Gene.Symbol=="DNMT3B",]


table(unlist(all[, c(-1, -2, -3)]))

wantedlong <- sapply(wanted, function(x) {colnames(d)[grepl(x, colnames(d))]})

all[, wantedlong]

names(all)[all==2]

# load the segmented log2fc
segm <- read.table(gzfile('/scratch/trcanmed/biobanca/dataset/V1/shallowseq/pdo_segm_l2fc.tsv.gz'), sep="\t", header=TRUE)

dnmt3b_expr <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMO_BASALE-DNMT3B_tmm_ave.tsv', sep="\t", header=T)

segm$chr <- as.numeric(sapply(strsplit(rownames(segm), ":"), function(x) {x[[1]][1]}))
segm$coord <- sapply(strsplit(rownames(segm), ":"), function(x) {x[[2]][1]})
segm$b <- as.numeric(sapply(strsplit(segm$coord, "-"), function(x) {x[[1]][1]}))
segm$e <- as.numeric(sapply(strsplit(segm$coord, "-"), function(x) {x[[2]][1]}))

w_coord_b <- 32762424
w_coord_e <- 32809356

mintersect <- function(wb, we, wchr, df) {
  sel <- df[df$chr==wchr, ]
  sel <- sel[pmax(wb, sel$b)<pmin(we, sel$e), ]
  return(sel)
}

#g1 <- GRanges(seqnames="I", ranges=IRanges(start=c(10, 30), end=c(15,35)))
#2 <- GRanges(seqnames="I", ranges=IRanges(start=c(10, 37), end=c(15,40)))
#intersect(g1, g2)

onDNMT3b <- mintersect(w_coord_b, w_coord_e, 20, segm)
onDNMT3b$size <- onDNMT3b$e - onDNMT3b$b
        
lfcDNMT3b <- onDNMT3b[,grepl('CRC', colnames(onDNMT3b))]

avelfc <- as.data.frame(colMeans(lfcDNMT3b))
colnames(avelfc) <- c('avelfc')
m <- merge(avelfc, dnmt3b_expr,by.x="row.names", by.y="model")

m$sel <- ifelse(m$Row.names %in% wanted, 'selected', 'other')

colnames(m)[1] <- c('model')

ggplot(data=m, aes(x=avelfc, y=expr))+geom_point(aes(color=sel))+geom_smooth(method='lm')+theme_bw()
cor.test(m$expr, m$avelfc)

# 20q
#http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz

w_coord_b <- 28100000
w_coord_e <- 64444167

onDNMT3b <- mintersect(w_coord_b, w_coord_e, 20, segm)
onDNMT3b$size <- onDNMT3b$e - onDNMT3b$b

lfcDNMT3b <- onDNMT3b[,grepl('CRC', colnames(onDNMT3b))]

avelfc <- as.data.frame(colMeans(lfcDNMT3b))
colnames(avelfc) <- c('avelfc')
m2 <- merge(avelfc, dnmt3b_expr,by.x="row.names", by.y="model")

m2$sel <- ifelse(m2$Row.names %in% wanted, 'selected', 'other')

colnames(m2)[1] <- c('model')

ggplot(data=m2, aes(x=avelfc, y=expr))+geom_point(aes(color=sel))+geom_smooth(method='lm')+theme_bw()
cor.test(m2$expr, m2$avelfc)


#######
# another nearby gene

#w_coord_b <- 34268678
#w_coord_e <- 34314019
w_coord_b <- 62228706
w_coord_e <- 62279677

onDNMT3b <- mintersect(w_coord_b, w_coord_e, 20, segm)
onDNMT3b$size <- onDNMT3b$e - onDNMT3b$b

lfcDNMT3b <- onDNMT3b[,grepl('CRC', colnames(onDNMT3b))]

avelfc <- as.data.frame(colMeans(lfcDNMT3b))
colnames(avelfc) <- c('avelfc')
m2 <- merge(avelfc, dnmt3b_expr,by.x="row.names", by.y="model")

m$sel <- ifelse(m2$Row.names %in% wanted, 'selected', 'other')

colnames(m2)[1] <- c('model')

ggplot(data=m2, aes(x=avelfc, y=expr))+geom_point(aes(color=sel))+geom_smooth(method='lm')+theme_bw()
cor.test(m2$expr, m2$avelfc)

m3 <- merge(m, m2, by="model")
ggplot(data=m3, aes(x=avelfc.x, y=avelfc.y))+geom_point(aes(color=sel))+geom_smooth(method='lm')+theme_bw()
cor.test(m3$avelfc.x, m3$avelfc.y)
ggplot(data=m3, aes(x=avelfc.x, y=avelfc.y))+geom_point(aes(color=expr.x))+geom_smooth(method='lm')+theme_bw()

       