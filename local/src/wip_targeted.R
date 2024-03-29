setwd('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdx/')
pdx <- read.table('./mutect/merged.table_nomultiallele', header=TRUE, sep="\t", row.names=1)
pdo <- read.table('../biobanca_targeted_pdo/mutect/merged.table_nomultiallele', header=TRUE, sep="\t", row.names=1)
afx <- unlist(pdx[pdx!=0])
afo <- unlist(pdo[pdo!=0])
pd <- data.frame(af=c(afo, afx), class=c(rep('pdo', length(afo)),rep('pdx',length(afx))))
ggplot(data=pd, aes(x=af, color=class, y=..count..)) + geom_density()+theme_bw()
ggplot(data=pd, aes(x=af, color=class)) + geom_density()+theme_bw()+scale_x_continuous(breaks=(seq(0, 1, by=0.05)))



##
pdx <- read.table('./mutect/merged.table_nomultiallele_wtiers', header=TRUE, sep="\t", row.names=1)
pdo <- read.table('../biobanca_targeted_pdo/mutect/merged.table_nomultiallele_wtiers', header=TRUE, sep="\t", row.names=1)

pdxbin <- ifelse(pdx==0,0, 1)
n_muts_x <- colSums(pdxbin)
pdobin <- ifelse(pdo==0,0, 1)
n_muts_o <- colSums(pdobin)

hist(n_muts_x, breaks=20)
hist(n_muts_o, breaks=20)


pd <- data.frame(af=c(n_muts_o, n_muts_x), class=c(rep('pdo', length(n_muts_o)),rep('pdx',length(n_muts_x))))

ggplot(data=pd, aes(x=af, color=class,y=..count..)) + geom_density()+xlab('nmuts')+current_theme


sanger <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/dataset/sanger_targeted_v2_genealogy/nmuts.tsv',sep="\t", header=F)
hist(sanger$V1, breaks=20)

sanger <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/dataset/sanger_targeted_v2_genealogy/n_muts_driver.tsv',sep="\t", header=F)
hist(sanger$V1, breaks=20)

gpdx <- read.table('./mutect/merged.table_nomultiallele_genes', header=TRUE, sep="\t", row.names=1)
gpdo <- read.table('../biobanca_targeted_pdo/mutect/merged.table_nomultiallele_genes', header=TRUE, sep="\t", row.names=1)

#gpdx <- read.table('./mutect/merged.table_nomultiallele_wtiers', header=TRUE, sep="\t", row.names=1)
#gpdo <- read.table('../biobanca_targeted_pdo/mutect/merged.table_nomultiallele_wtiers', header=TRUE, sep="\t", row.names=1)

###
AFt <- 0.2

all_genes <- intersect(unique(gpdo$gene), unique(gpdx$gene))

get_genes_binary <- function(sample, AF, genes, thr, all_genes) {
    keep <- AF[,sample] > thr
    genes <- genes[keep,]
    return(all_genes %in% genes)
}

gpdo <- gpdo[match(rownames(pdo), rownames(gpdo)),, drop=F] # ASSUMPTIONS ASSUMPTIONS FOTTITI ELENA
gpdx <- gpdx[match(rownames(pdx), rownames(gpdx)),, drop=F]
pdobing <- sapply(colnames(pdo), get_genes_binary, pdo, gpdo, AFt, all_genes)
pdxbing <- sapply(colnames(pdx), get_genes_binary, pdx, gpdx, AFt, all_genes)

rownames(pdobing) <- all_genes
rownames(pdxbing) <- all_genes

pdo_genes <- apply(pdobing, 1, sum)
pdx_genes <- apply(pdxbing, 1, sum)
freqs <- data.frame(pdo_freq = pdo_genes/ncol(pdobing), pdx_freq=pdx_genes/ncol(pdxbing))
freqs$gene <- rownames(freqs)
tcga <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/tcga_panel_freqs.tsv', sep="\t")
msk <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/msk_panel_freqs.tsv', sep="\t")
colnames(tcga) <- c('gene', 'freq')
colnames(msk) <- c('gene', 'freq')

m1 <- merge(tcga, msk, all.x = TRUE, by='gene')
colnames(m1) <- c('gene','tcga_freq','msk_freq')
m2 <- merge(m1, freqs, all.x = TRUE, by='gene')

m2[is.na(m2$msk_freq),]$msk_freq <- 0
m2[is.na(m2$pdx_freq),]$pdx_freq <- 0
m2[is.na(m2$pdo_freq),]$pdo_freq <- 0

#pd <- melt(m2)

ggplot(data=m2, aes(x=pdx_freq, y=tcga_freq))+ geom_point()+ geom_smooth(method=lm, se=FALSE)+current_theme

cl <- c(rep('tcga',nrow(m2)), rep('msk',nrow(m2)))
#cl2 <- c(rep('pdo',nrow(m2)), rep('pdx',nrow(m2)), rep('pdo',nrow(m2)), rep('pdx',nrow(m2)))

m3 <- data.frame(x=c(m2$pdo_freq, m2$pdo_freq), 
                 y=c(m2$tcga_freq, m2$msk_freq), 
                 class=cl)

lmplot <- function(data, title, xlim=NULL) {
  fit1 <- data[data$class=='tcga',]
  fit2 <- data[data$class=='msk',]
  pi1 <- cor.test(fit1$x, fit1$y)
  pi2 <- cor.test(fit2$x, fit2$y)
  
  cap1 <- round(c(pi1$estimate, pi1$p.value, pi2$estimate, pi2$p.value),2)
  print(cap1)
  names(cap1) <- NULL
  cap <- paste(cap1, collapse= " ")
  if (is.null(xlim)) {
    ggplot(data=data, aes(x=x, y=y, color=class))+ geom_point()+ geom_smooth(method=lm, se=FALSE)+current_theme+labs(caption=cap)+ggtitle(title)
  } else {
    ggplot(data=data, aes(x=x, y=y, color=class))+ geom_point()+ geom_smooth(method=lm, se=FALSE)+current_theme+labs(caption=cap)+ggtitle(title)+xlim(xlim)
  }
}

lmplot(m3, 'PDO')
lmplot(m3, 'PDO', c(0, 0.2))

m3 <- data.frame(x=c(m2$pdx_freq, m2$pdx_freq), 
                 y=c(m2$tcga_freq, m2$msk_freq), 
                 class=cl)
lmplot(m3, 'PDX')
lmplot(m3, 'PDX', c(0, 0.2))
###
pdo_df <- pdobing
xeno_df <- pdxbing

##### UPTO


pdo_df <- pdo_df[, colnames(pdo_df) != "CRC0282.01.1.A"]
#> colnames(pdobing)[grepl('CRC0177', colnames(pdobing))]
#[1] "CRC0177LMO0A04008002D02000" "CRC0177LMO0D04017001D02000"
#> colnames(pdxbing)[grepl('CRC0177', colnames(pdxbing))]
#[1] "CRC0177LMX0B05001TUMD07000" "CRC0177LMX0B07002TUMD05000"

colnames(pdo_df)[colnames(pdo_df)=="CRC0177LMO0A04008002D02000"] <- 'CRCE177LMO0A04008002D02000'
colnames(xeno_df)[colnames(xeno_df)=="CRC0177LMX0B05001TUMD07000"] <- 'CRCE177LMX0B05001TUMD06000'

# TODO? make modular: first step to prepare matrix then general script for correlations
#CRC1875
list_remove_xeno <- c("CRC1875") 
list_lmh <- colnames(pdo_df)[grepl('LMH', colnames(pdo_df))]
model_lmh <- substr(list_lmh, 0,7)

pdo_df <- pdo_df[, !colnames(pdo_df) %in% list_lmh]
pdo_df <- pdo_df[, !substr(colnames(pdo_df),0,7) %in% model_lmh]
xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% model_lmh]
  
#> setdiff(models_xeno, models_pdo)
#[1] "CRC1870" "CRC1875" "CRC2041"
# We know we were missing some PDOs TODO load from file if qc changes in the future
xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% list_remove_xeno]


expected_n <- 144 # 2 more than shallow cause we miss less pdx

models_xeno <- substr(colnames(xeno_df), 0, 7)
models_pdo <- substr(colnames(pdo_df), 0, 7)
if (length(intersect(models_xeno, models_pdo)) != expected_n) {
  stop('Wrong number of corresponding models!')
}

colnames(xeno_df) <- substr(colnames(xeno_df), 0, 7)
colnames(pdo_df) <- substr(colnames(pdo_df), 0, 7)

xeno_df <- xeno_df[, models_xeno[order(models_xeno)]]
pdo_df <- pdo_df[, models_xeno[order(models_xeno)]]

if (!all(substr(colnames(xeno_df),0,7) == substr(colnames(pdo_df), 0,7))) {
  stop('Piciapirilla!')
}

library(proxy)

jac <- proxy::simil(xeno_df, pdo_df, by_rows = FALSE, method = "Jaccard")

pheatmap(jac, cluster_rows=F , cluster_cols=F, labels_col="PDO", labels_row="PDX", fontsize.number=1.5)

pearson <- jac
diag <- diag(pearson)
pearson2 <- pearson
diag(pearson2) <- rep(NA, length(diag))
all <- as.numeric(unlist(pearson2))
all <- all[!is.na(all)]
#all <- upper.tri(pearson, diag = FALSE) # this is not a simmetric matrix!
pdata <- data.frame(pearson=c(all, diag), type=c(rep('unmatched', length(all)),rep('matched', length(diag))))
ggplot(data=pdata, aes(x=pearson, color=type))+geom_density()+theme_bw()+theme(text=element_text(size=20))+xlab('jaccard')

### work at the muts level
AFt <- 0.2

all_muts <- intersect(unique(rownames(pdo)), rownames(pdx))

get_muts_binary <- function(sample, AF, thr, all_muts) {
  keep <- AF[,sample] > thr
  AF <- AF[keep,]
  return(all_muts %in% rownames(AF))
}

pdobing <- sapply(colnames(pdo), get_muts_binary, pdo, AFt, all_muts)
pdxbing <- sapply(colnames(pdx), get_muts_binary, pdx, AFt, all_muts)

###
pdo_df <- pdobing
xeno_df <- pdxbing

pdo_df <- pdo_df[, colnames(pdo_df) != "CRC0282.01.1.A"]
#> colnames(pdobing)[grepl('CRC0177', colnames(pdobing))]
#[1] "CRC0177LMO0A04008002D02000" "CRC0177LMO0D04017001D02000"
#> colnames(pdxbing)[grepl('CRC0177', colnames(pdxbing))]
#[1] "CRC0177LMX0B05001TUMD07000" "CRC0177LMX0B07002TUMD05000"

colnames(pdo_df)[colnames(pdo_df)=="CRC0177LMO0A04008002D02000"] <- 'CRCE177LMO0A04008002D02000'
colnames(xeno_df)[colnames(xeno_df)=="CRC0177LMX0B05001TUMD07000"] <- 'CRCE177LMX0B05001TUMD06000'

# TODO? make modular: first step to prepare matrix then general script for correlations
#CRC1875
list_remove_xeno <- c("CRC1875") 
list_lmh <- colnames(pdo_df)[grepl('LMH', colnames(pdo_df))]
model_lmh <- substr(list_lmh, 0,7)

pdo_df <- pdo_df[, !colnames(pdo_df) %in% list_lmh]
pdo_df <- pdo_df[, !substr(colnames(pdo_df),0,7) %in% model_lmh]
xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% model_lmh]

#> setdiff(models_xeno, models_pdo)
#[1] "CRC1870" "CRC1875" "CRC2041"
# We know we were missing some PDOs TODO load from file if qc changes in the future
xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% list_remove_xeno]


expected_n <- 144 # 2 more than shallow cause we miss less pdx

models_xeno <- substr(colnames(xeno_df), 0, 7)
models_pdo <- substr(colnames(pdo_df), 0, 7)
if (length(intersect(models_xeno, models_pdo)) != expected_n) {
  stop('Wrong number of corresponding models!')
}

colnames(xeno_df) <- substr(colnames(xeno_df), 0, 7)
colnames(pdo_df) <- substr(colnames(pdo_df), 0, 7)

xeno_df <- xeno_df[, models_xeno[order(models_xeno)]]
pdo_df <- pdo_df[, models_xeno[order(models_xeno)]]


xeno_df <- xeno_df[, order(models_xeno)]
pdo_df <- pdo_df[, order(models_xeno)]
library(proxy)

jac <- proxy::simil(xeno_df, pdo_df, by_rows = FALSE, method = "Jaccard")

#Jaccard Similarity = (number of observations in both sets) / (number in either set)

pheatmap(jac, cluster_rows=F , cluster_cols=F, labels_col="PDO", labels_row="PDX", fontsize.number=1.5)

pearson <- jac
diag <- diag(pearson)
pearson2 <- pearson
diag(pearson2) <- rep(NA, length(diag))
all <- as.numeric(unlist(pearson2))
all <- all[!is.na(all)]
#all <- upper.tri(pearson, diag = FALSE) # this is not a simmetric matrix!
pdata <- data.frame(pearson=c(all, diag), type=c(rep('unmatched', length(all)),rep('matched', length(diag))))
ggplot(data=pdata, aes(x=pearson, color=type))+geom_density()+theme_bw()+theme(text=element_text(size=20))+xlab('jaccard')

#### oncoprint
library(ComplexHeatmap)
all_genes <- intersect(unique(gdpo$gene), unique(gpdx$gene))

get_genes_binary <- function(sample, AF, genes, thr, all_genes) {
  keep <- AF[,sample] > thr
  genes <- genes[keep,]
  return(all_genes %in% genes)
}

pdobing <- sapply(colnames(pdo), get_genes_binary, pdo, gpdo, AFt, all_genes)
pdxbing <- sapply(colnames(pdx), get_genes_binary, pdx, gpdx, AFt, all_genes)

rownames(pdobing) <- all_genes
rownames(pdxbing) <- all_genes

pdo_df <- pdobing
xeno_df <- pdxbing

pdo_df <- pdo_df[, colnames(pdo_df) != "CRC0282.01.1.A"]
#> colnames(pdobing)[grepl('CRC0177', colnames(pdobing))]
#[1] "CRC0177LMO0A04008002D02000" "CRC0177LMO0D04017001D02000"
#> colnames(pdxbing)[grepl('CRC0177', colnames(pdxbing))]
#[1] "CRC0177LMX0B05001TUMD07000" "CRC0177LMX0B07002TUMD05000"

colnames(pdo_df)[colnames(pdo_df)=="CRC0177LMO0A04008002D02000"] <- 'CRCE177LMO0A04008002D02000'
colnames(xeno_df)[colnames(xeno_df)=="CRC0177LMX0B05001TUMD07000"] <- 'CRCE177LMX0B05001TUMD06000'

# TODO? make modular: first step to prepare matrix then general script for correlations
#CRC1875
list_remove_xeno <- c("CRC1875") 
list_lmh <- colnames(pdo_df)[grepl('LMH', colnames(pdo_df))]
model_lmh <- substr(list_lmh, 0,7)

pdo_df <- pdo_df[, !colnames(pdo_df) %in% list_lmh]
pdo_df <- pdo_df[, !substr(colnames(pdo_df),0,7) %in% model_lmh]
xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% model_lmh]

#> setdiff(models_xeno, models_pdo)
#[1] "CRC1870" "CRC1875" "CRC2041"
# We know we were missing some PDOs TODO load from file if qc changes in the future
xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% list_remove_xeno]


colnames(xeno_df) <- substr(colnames(xeno_df), 0, 7)
colnames(pdo_df) <- substr(colnames(pdo_df), 0, 7)

###

d2 <- t(pdo_df)
d2 <- ifelse(d2, 1, 0)
d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w_CRC0078_PR.tsv', sep="\t", header=FALSE, stringsAsFactors = TRUE, row.names=1)
rownames(d) <- paste0(rownames(d),'LMO')
#also pdo then
cases <- rownames(d2)

#rimuovere i non mutati
#> sum(rowSums(d2) == 0)

d3 <- data.frame(matrix(NA, ncol=ncol(d2), nrow=nrow(d2)))

for (i in seq(1, ncol(d2))) {
  f <- as.factor(d2[,i])
  levels(f) <- c("","SNV")    
  d3[,i] <- f
}
rownames(d3) <- rownames(d2)
colnames(d3) <- colnames(d2)
#d3 <- apply(d2, 2, function(x) { f <- as.factor(x); levels(f) <- c("","SNV") })
#col = c("HOMDEL" = "blue", "AMP" = "red", "GAIN"= "#e86f0a", "MISSENSE" = "#008000", "STOP" = "#91f903", "INDEL" = "#01c601", "OTHER"="#9a4af9")
col = c("HOMDEL" = "blue", "AMP" = "red", "GAIN"= "#e86f0a", "SNV" = "#008000", "OTHER"="#9a4af9")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # bug red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  GAIN = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["GAIN"], col = NA))
  },
  # small green
  # MISSENSE = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["MISSENSE"], col = NA))
  # },
  # STOP = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["STOP"], col = NA))
  # },
  # INDEL = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["INDEL"], col = NA))
  # },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["SNV"], col = NA))
  },
  OTHER = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["OTHER"], col = NA))
  }
)

column_title = "OncoPrint for PDO"
#heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "GAIN", "MISSENSE","STOP","INDEL","OTHER"), 
#                            labels = c("Deep deletion", "High/Focal amplification", "Non Focal amplification","Missense SVN","Stop SNV","Indel","Non-coding SNV"))

heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "GAIN", "SNV","OTHER"), 
                            labels = c("Deep deletion", "High/Focal amplification", "Non Focal amplification","Coding SVN","Non-coding SNV"))

rownames(d3) <- substr(rownames(d3),0,10)
dd <- d[rownames(d) %in% rownames(d3),, drop=FALSE] # volumes only of sequenced cases
sample_order <- rownames(dd) # the order of volumes (they are sorted when we load them)
d4 <- d3[rownames(d3)%in% sample_order,] # muts only for cases with volumes
td4 <- t(d4)
td4 <- td4[, sample_order] # we order them correctly
#library(circlize)
#col_fun = colorRamp2(c(35, -50, -100), c("red", "white", "blue"))

#oncoPrint(td4,
#          alter_fun = alter_fun, col = col, 
#          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
#          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
#                                             Irinotecan = dd$perc, col=list(Irinotecan=col_fun)),
#          row_names_gp = gpar(fontsize=8), column_order=sample_order)

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50, 4, ifelse(x > 35, 2, 8))
  #res[is.na(x)] <- 'black'
  return(as.numeric(res))
}
library(RColorBrewer)
#colOR <- c("PD"="red", "SD"="blue","OR"="black")

# 
colnames(dd) <- 'perc'
dd$perc <- dd$perc * 100
dd <- dd[order(-dd$perc), , drop=FALSE]

sample_order <- rownames(dd)
td4 <- td4[,sample_order]
oncoPrint(td4,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                             Cetuximab = anno_barplot(dd$perc, gp = gpar(fill = get_recist(dd$perc))), height = unit(4, "cm")),
          row_names_gp = gpar(fontsize=8), column_order=sample_order)
#[1] 4

######


d2 <- t(xeno_df)
d2 <- ifelse(d2, 1, 0)
d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/april2020/cetuxi_w3.txt', sep="\t", header=FALSE, stringsAsFactors = TRUE, row.names=1)
#also pdo then
cases <- rownames(d2)

#rimuovere i non mutati
#> sum(rowSums(d2) == 0)

d3 <- data.frame(matrix(NA, ncol=ncol(d2), nrow=nrow(d2)))

for (i in seq(1, ncol(d2))) {
  f <- as.factor(d2[,i])
  levels(f) <- c("","SNV")    
  d3[,i] <- f
}
rownames(d3) <- rownames(d2)
colnames(d3) <- colnames(d2)
#d3 <- apply(d2, 2, function(x) { f <- as.factor(x); levels(f) <- c("","SNV") })
#col = c("HOMDEL" = "blue", "AMP" = "red", "GAIN"= "#e86f0a", "MISSENSE" = "#008000", "STOP" = "#91f903", "INDEL" = "#01c601", "OTHER"="#9a4af9")
col = c("HOMDEL" = "blue", "AMP" = "red", "GAIN"= "#e86f0a", "SNV" = "#008000", "OTHER"="#9a4af9")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # bug red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  GAIN = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["GAIN"], col = NA))
  },
  # small green
  # MISSENSE = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["MISSENSE"], col = NA))
  # },
  # STOP = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["STOP"], col = NA))
  # },
  # INDEL = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["INDEL"], col = NA))
  # },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["SNV"], col = NA))
  },
  OTHER = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["OTHER"], col = NA))
  }
)

column_title = "OncoPrint for PDX"
#heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "GAIN", "MISSENSE","STOP","INDEL","OTHER"), 
#                            labels = c("Deep deletion", "High/Focal amplification", "Non Focal amplification","Missense SVN","Stop SNV","Indel","Non-coding SNV"))

heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "GAIN", "SNV","OTHER"), 
                            labels = c("Deep deletion", "High/Focal amplification", "Non Focal amplification","Coding SVN","Non-coding SNV"))


dd <- d[rownames(d) %in% rownames(d3),, drop=FALSE] # volumes only of sequenced cases
sample_order <- rownames(dd) # the order of volumes (they are sorted when we load them)
d4 <- d3[rownames(d3)%in% sample_order,] # muts only for cases with volumes
td4 <- t(d4)
td4 <- td4[, sample_order] # we order them correctly
#library(circlize)
#col_fun = colorRamp2(c(35, -50, -100), c("red", "white", "blue"))

#oncoPrint(td4,
#          alter_fun = alter_fun, col = col, 
#          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
#          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
#                                             Irinotecan = dd$perc, col=list(Irinotecan=col_fun)),
#          row_names_gp = gpar(fontsize=8), column_order=sample_order)

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50, 4, ifelse(x > 35, 2, 8))
  #res[is.na(x)] <- 'black'
  return(as.numeric(res))
}
library(RColorBrewer)
#colOR <- c("PD"="red", "SD"="blue","OR"="black")

# 
colnames(dd) <- 'perc'
dd$perc <- dd$perc * 100
dd <- dd[order(-dd$perc), , drop=FALSE]

sample_order <- rownames(dd)
td4 <- td4[,sample_order]
oncoPrint(td4,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                             Cetuximab = anno_barplot(dd$perc, gp = gpar(fill = get_recist(dd$perc))), height = unit(4, "cm")),
          row_names_gp = gpar(fontsize=8), column_order=sample_order)

#### together #################################################

mergemut <- merge(pdo_df, xeno_df, by="row.names")#, all.x=TRUE, all.y=TRUE, fill=FALSE) # no need all the same genes:
#> setdiff(rownames(pdo_df), rownames(xeno_df))
#character(0)
rownames(mergemut) <- mergemut$Row.names
mergemut$Row.names <- NULL
mergemut <- mergemut[, colnames(mergemut) != "CRC0282.01.1.A"]

d2 <- t(mergemut)
d2 <- ifelse(d2, 1, 0)
d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w_CRC0078_PR.tsv', sep="\t", header=FALSE, stringsAsFactors = TRUE, row.names=1)
#also pdo then
cases <- rownames(d2)

#rimuovere i non mutati
#> sum(rowSums(d2) == 0)


d2pdo <- d2[grepl('LMO', rownames(d2)),]
# d3pdo <- data.frame(matrix(NA, ncol=ncol(d2pdo), nrow=nrow(d2pdo)))
# 
# for (i in seq(1, ncol(d2pdo))) {
#   f <- as.factor(d2pdo[,i])
#   levels(f) <- c("","OSNV")    
#   d3pdo[,i] <- f
# }

d2pdx <- d2[grepl('LMX', rownames(d2)),]
# d3pdx <- data.frame(matrix(NA, ncol=ncol(d2pdx), nrow=nrow(d2pdx)))
# 
# for (i in seq(1, ncol(d2pdx))) {
#   f <- as.factor(d2pdx[,i])
#   levels(f) <- c("","XSNV")    
#   d3pdx[,i] <- f
# }

td2x <- t(d2pdx)
td2o <- t(d2pdo)
colnames(td2x) <- substr(colnames(td2x), 0, 7)
colnames(td2o) <- substr(colnames(td2o), 0, 7)

td2x <- td2x[, colnames(td2x) != "CRC1875"]
#mat_list <- list(OSNV=td2o, XSNV=td2x)
both <- td2o & td2x
both2 <- t(apply(both, 1, as.numeric))
rownames(both2) <- rownames(both)
mat_list <- list(Both=both, PDO= td2o-both, PDX=td2x-both  ) 

#rownames(d3pdo) <- rownames(d2pdo)
#colnames(d3pdo) <- colnames(d2pdo)

#rownames(d3pdx) <- rownames(d2pdx)
#colnames(d3pdx) <- colnames(d2pdx)

#if (!all(colnames(d3pdo)==colnames(d3pdo))) {
#  print("LLLAMA!!!")
#}

#d3 <- rbind(d3pdx, d3pdo)
#d2 <- rbind(d2pdx, d3pdo)

#d3 <- apply(d2, 2, function(x) { f <- as.factor(x); levels(f) <- c("","SNV") })
#col = c("HOMDEL" = "blue", "AMP" = "red", "GAIN"= "#e86f0a", "MISSENSE" = "#008000", "STOP" = "#91f903", "INDEL" = "#01c601", "OTHER"="#9a4af9")
col = c("Both" = "blue", "PDO" = "red", "PDX"= "#0b7015")#, "XSNV" = "#ffcc00", "OTHER"="#9a4af9", "OSNV" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  #background = function(x, y, w, h) {
  #  grid.polygon(
  #    unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
  #    unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
  #    gp = gpar(fill = "grey", col = "white"))
  #  grid.polygon(
  #    unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
  #    unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
  #    gp = gpar(fill = "grey", col = "white"))
  #},
  # big blue
  Both = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Both"], col = NA))
  },
  # bug red
  PDO = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["PDO"], col = NA))
  },
  PDX = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["PDX"], col = NA))
  },
  # small green
  # MISSENSE = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["MISSENSE"], col = NA))
  # },
  # STOP = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["STOP"], col = NA))
  # },
  # INDEL = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["INDEL"], col = NA))
  # },
  # OSNV = function(x, y, w, h) {
  #   grid.polygon(
  #     unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
  #     unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
  #     gp = gpar(fill = col["OSNV"], col = "white"))
  #   
  # },
  # XSNV = function(x, y, w, h) {
  #   grid.polygon(
  #     unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
  #     unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
  #     gp = gpar(fill = col["XSNV"], col = "white"))
  # },
  OTHER = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["OTHER"], col = NA))
  }
)

column_title = "OncoPrint for PDXO"
#heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "GAIN", "MISSENSE","STOP","INDEL","OTHER"), 
#                            labels = c("Deep deletion", "High/Focal amplification", "Non Focal amplification","Missense SVN","Stop SNV","Indel","Non-coding SNV"))

heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "GAIN", "SNV","OTHER"), 
                            labels = c("Deep deletion", "High/Focal amplification", "Non Focal amplification","Coding SVN","Non-coding SNV"))


dd <- d[rownames(d) %in% rownames(d3),, drop=FALSE] # volumes only of sequenced cases
sample_order <- rownames(dd) # the order of volumes (they are sorted when we load them)
#d4 <- d3[rownames(d3)%in% sample_order,] # muts only for cases with volumes
td4 <- t(d3)
#td4 <- td4[, sample_order] # we order them correctly
#library(circlize)
#col_fun = colorRamp2(c(35, -50, -100), c("red", "white", "blue"))

#oncoPrint(td4,
#          alter_fun = alter_fun, col = col, 
#          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
#          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
#                                             Irinotecan = dd$perc, col=list(Irinotecan=col_fun)),
#          row_names_gp = gpar(fontsize=8), column_order=sample_order)

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50, 4, ifelse(x > 35, 2, 8))
  #res[is.na(x)] <- 'black'
  return(as.numeric(res))
}
library(RColorBrewer)
#colOR <- c("PD"="red", "SD"="blue","OR"="black")

# 
colnames(dd) <- 'perc'
dd$perc <- dd$perc * 100
dd <- dd[order(-dd$perc), , drop=FALSE]

sample_order <- rownames(dd)
td4 <- td4[,sample_order]
oncoPrint(td4,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                             Cetuximab = anno_barplot(dd$perc, gp = gpar(fill = get_recist(dd$perc))), height = unit(4, "cm")),
          row_names_gp = gpar(fontsize=8), column_order=sample_order)
#[1] 4

mat_list2 <- list()
#mat_list2[[1]] <- mat_list[[1]][,seq(1,30)] # y
#mat_list2[[2]] <- mat_list[[2]][,seq(1,30)]
#mat_list2[[3]] <- mat_list[[3]][,seq(1,30)]
mat_list2[[1]] <- mat_list[[1]]
mat_list2[[2]] <- mat_list[[2]]
mat_list2[[3]] <- mat_list[[3]]

names(mat_list2) <- names(mat_list)

m <- mat_list[[1]]

d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w_CRC0078_PR.tsv', sep="\t", header=FALSE, stringsAsFactors = TRUE, row.names=1)

length(intersect(rownames(d), colnames(m)))

s <- mat_list2[[1]]+mat_list2[[2]]+mat_list2[[3]]
su <- colSums(s)
su <- su[order(-su)]

oncoPrint(mat_list2, alter_fun = alter_fun, col = col)
oncoPrint(mat_list2, alter_fun = alter_fun, col = col, column_order = names(su))

oncoPrint(mat_list2,
          alter_fun = alter_fun, col = col, top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                                               Cetuximab = anno_barplot(dd$perc, gp = gpar(fill = get_recist(dd$perc))), height = unit(4, "cm")),
          row_names_gp = gpar(fontsize=8))

pp <- oncoPrint(mat_list, alter_fun = alter_fun, col = col);

png('oncoprint_test.png',width = 880, height = 880, units = "px"); print(pp); dev.off(
