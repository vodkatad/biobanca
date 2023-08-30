data <- read.table('/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_replisafe', sep="\t", stringsAsFactors = F, header=T)
head(data)
lmo <- data[grepl('LMO_BASALE', data$type),]
lmo$scong <- as.numeric(substr(lmo$sample_id_R, 13, 14))
lmo$passaggio <-as.numeric(substr(lmo$sample_id_R, 15, 17))

x_axis_labels <- min(lmo[,'scong']):max(lmo[,'scong'])
ggplot(data=lmo, aes(x=scong))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)

ggplot(data=lmo, aes(x=passaggio))+geom_bar()+theme_bw()


targ_buoni <- read.table('/scratch/trcanmed/biobanca/local/share/data/targeted/map_all_validation.tsv',  sep="\t", stringsAsFactors = F, header=T)
lmo_t <- targ_buoni[!is.na(targ_buoni$validation_status) & targ_buoni$validation_status & grepl('LMO', targ_buoni$genalogy),]

lmo_t$scong <- as.numeric(substr(lmo_t$genalogy, 13, 14))
lmo_t$passaggio <- as.numeric(substr(lmo_t$genalogy, 15, 17))
ggplot(data=lmo_t, aes(x=passaggio))+geom_bar()+theme_bw()
x_axis_labels <- min(lmo_t[,'scong']):max(lmo_t[,'scong'])
ggplot(data=lmo_t, aes(x=scong))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)
x_axis_labels <- min(lmo_t[,'passaggio']):max(lmo_t[,'passaggio'])
ggplot(data=lmo_t, aes(x=passaggio))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)



shallow_buoni <- read.table('/scratch/trcanmed/biobanca/local/share/data/shallowseq/map_all_validation.tsv',  sep="\t", stringsAsFactors = F, header=T)
lmo_t <- shallow_buoni[!is.na(shallow_buoni$validation_status) & shallow_buoni$validation_status & grepl('LMO', shallow_buoni$genalogy),]

lmo_t$scong <- as.numeric(substr(lmo_t$genalogy, 13, 14))
lmo_t$passaggio <- as.numeric(substr(lmo_t$genalogy, 15, 17))
ggplot(data=lmo_t, aes(x=passaggio))+geom_bar()+theme_bw()
x_axis_labels <- min(lmo_t[,'scong']):max(lmo_t[,'scong'])
ggplot(data=lmo_t, aes(x=scong))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)


reps <- as.data.frame(table(substr(lmo$sample_id_R, 0,7)))
table(reps$Freq)

has_rep <- reps[reps$Freq==2,'Var1']

lmo$rep <- ifelse(substr(lmo$sample_id_R, 0,7) %in% has_rep, 'yes', 'no')

lmo$model <- substr(lmo$sample_id_R, 0,7)
delta_col <- function(model, data, col) {
  subd <- data[data$model == model,]
  if (nrow(subd)==2) {
    p1 <- min(subd[1, col],subd[2, col])
    p2 <- max(subd[1, col],subd[2, col])
    return(c(abs(subd[1, col]-subd[2, col]), p1, p2))
  } else {
    return(NA)
  }
}

delta_scong <- sapply(has_rep, delta_col, lmo, 'scong')
delta_passaggio <- sapply(has_rep, delta_col, lmo, 'passaggio')

delta_passaggiopd <- as.data.frame(t(delta_passaggio))
colnames(delta_passaggiopd) <- c('delta_passaggio', 'p1', 'p2')
delta_passaggiopd$model <- has_rep

pd <- data.frame(model=has_rep, delta_scong=delta_scong, delta_passaggio=delta_passaggiopd$delta_passaggio)
ggplot(data=pd, aes(x=delta_passaggio))+geom_bar()+theme_bw()
ggplot(data=pd, aes(x=delta_scong))+geom_bar()+theme_bw()

#
orderdp <- delta_passaggiopd[order(delta_passaggiopd$delta_passaggio),]
orderdp1 <- orderdp
library(reshape)
orderdp$delta_passaggio <- NULL
mm <- melt(orderdp)
mmm <- merge(mm, orderdp1, by="model")
mmm$larger5 <- ifelse(mmm$delta_passaggio > 5, 'yes', 'no')
ggplot(data=mmm, aes(x=reorder(model,delta_passaggio), y=value, color=larger5))+geom_point()+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

lmo <- lmo[lmo$rep == "yes",]
x_axis_labels <- min(lmo[,'passaggio']):max(lmo[,'passaggio'])
ggplot(data=lmo, aes(x=passaggio))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)

early_late <- function(model, data, col, low=3, high=8) {
  subd <- data[data$model == model,]
  subd <- subd[order(subd[,col]),]
  if (nrow(subd)==2) {
    if (subd[1, col] <= low && subd[2, col] >= high) {
      return(TRUE)
    } else{
      return(FALSE)
    }
  } else {
    return(NA)
  }
}

el <- sapply(has_rep, early_late, lmo, 'passaggio')


el2 <- sapply(has_rep, early_late, lmo, 'passaggio', low=10, high=15)
el3 <- sapply(has_rep, early_late, lmo, 'passaggio', low=5, high=15)


#### expr geni LMO_S + R
genes <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/local/share/data/xenturion_threesome', header=FALSE)

d <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/fpkm.tsv.gz', sep="\t")
ann <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/samples_data', sep="\t", header=TRUE)
rownames(ann) <- ann$id
ann$id <- NULL

dr <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_R/fpkm.tsv.gz', sep="\t")
annr <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_R/samples_data', sep="\t", header=TRUE)
rownames(annr) <- annr$id
annr$id <- NULL

ann$class <- 'S'
annr$class <-  'R'
ann <- rbind(ann, annr)
d <- merge(d, dr, by="row.names")
rownames(d) <- d$Row.names
d$Row.names <- NULL

expr <- d[rownames(d) %in% paste0('H_', genes$V1),]


expr <- expr[, rownames(ann)]

ann <- ann[order(ann$sample, ann$treat),]
expr <- expr[, match(rownames(ann), colnames(expr))]
expr <- log2(expr+1)
pheatmap(expr, annotation_col=ann, cluster_rows = F, cluster_cols=F)

# expr levels excel
toexc <- merge(t(expr), ann, by="row.names")
colnames(toexc)[1] <- 'genealogy'
library(WriteXLS)
WriteXLS(toexc, ExcelFileName = "cetuxi_pdo.xlsx")

fc <- function(model, data, sign, ann) {
  d <- data[rownames(data)==sign,]
  dt <- t(d)
  m <- merge(dt, ann, by="row.names")
  d <- m[m$sample==model, ]
  d <- d[order(d$Row.names),] # we resort to using gen id to identify treat/nt, but we check
  if (nrow(d) == 4) {
    fc <- mean(c(d[2, sign] - d[1, sign], d[4, sign] - d[3, sign])) 
    stopifnot(d[1, 'treat'] == "NT" && d[2, 'treat'] == "cetuxi" && d[3, 'treat'] == "NT" && d[4, 'treat'] == "cetuxi")
  } else {
    fc <- d[2, sign] - d[1, sign]
    stopifnot(d[1, 'treat'] == "NT" && d[2, 'treat'] == "cetuxi")
  }
  return(fc)
}

callfc <- function(model, data, ann) {
  genes <- rownames(data)
  res <- c()
  for (g in genes) {
    res <- c(res, fc(model, data, g, ann))
  }
  return(res)
}

fcs <- sapply(unique(ann$sample), callfc, expr, ann)
fcs <- as.data.frame(fcs)
rownames(fcs) <- rownames(expr)
colnames(fcs) <- unique(ann$sample)
#fcs$gene <- rownames(fcs)

td <- t(fcs)
response <- unique(ann[c('sample', 'class')])
fcs_re <- merge(td, response, by.x="row.names", by.y="sample")
setwd('/home/egrassi')
colnames(fcs_re)[1] <- 'model'
WriteXLS(as.data.frame(fcs), ExcelFileName = "cetuxi_pdo_fc.xlsx")


