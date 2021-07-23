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

ggplot(data=pd, aes(x=af, color=class,y=..count..)) + geom_density()+theme_bw()


sanger <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/dataset/sanger_targeted_v2_genealogy/nmuts.tsv',sep="\t", header=F)
hist(sanger$V1, breaks=20)

sanger <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/dataset/sanger_targeted_v2_genealogy/n_muts_driver.tsv',sep="\t", header=F)
hist(sanger$V1, breaks=20)

gpdx <- read.table('./mutect/merged.table_nomultiallele_genes', header=TRUE, sep="\t", row.names=1)
gpdo <- read.table('../biobanca_targeted_pdo/mutect/merged.table_nomultiallele_genes', header=TRUE, sep="\t", row.names=1)

###
AFt <- 0.2

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

column_title = "OncoPrint for PDO"
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