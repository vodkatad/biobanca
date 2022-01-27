#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
x <- args[1]
y <- args[2]
output <- args[3]
prefix <- args[4]
merge_col <- args[5]
yname <- args[6]
mut_f <- args[7]
cn_f <- args[8]
Rimage <- args[9]

load(Rimage)

save.image('pippo.Rdata')

setwd('/scratch/trcanmed/biobanca/dataset/V1/cetuximab/atp')
load('pippo.Rdata')

data_x <- read.table(x, sep="\t", header=TRUE)
# y should be a two column file
data_y <- read.table(y, sep="\t", header=TRUE)
if (ncol(data_y) != 2 ) {
    stop(paste0(y, ' needs to be a two columns file!'))
}
colnames(data_y)[2] <- yname
data_x <- data_x[, c('case','CTG_5000')]
mdata <- merge(data_x, data_y, by=merge_col)
rownames(mdata) <- mdata[, merge_col]
mdata[,merge_col] <- NULL

mut <- read.table(mut_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
cn <- read.table(cn_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)

# cn has some duplicates, we first set all missing as NA then remove them
cn2 <- as.data.frame(t(apply(cn, 1, function(x) { xx <- as.character(x[-1]); suppressWarnings(as.numeric(ifelse(xx =='-' | xx == '', NA, xx))) })))
colnames(cn2) <- colnames(cn)[-1]
cn2$case <- cn$CASE

# we keep for replicates the row without missing value, we print in stdout exceptions
duplicate_info_cn <- function(case, data) {
  subset <- data[data$case == case,, drop=FALSE]
  subset$case <- NULL
  if (nrow(subset)==1) {
    return(as.matrix(subset))
  } else {
    sums <- apply(subset, 1, sum)
    chosenrow <- names(sums[!is.na(sums)])
    if (length(chosenrow) != 1) {
      print(subset)
      chosenrow <- chosenrow[1]
      #     MET EGFR HER2
      #1056 6.0  0.9  0.8
      #1057 6.6  0.1  0.0
    }
    return(as.matrix(subset[rownames(subset) == chosenrow, , drop=FALSE]))
  }
}

uniqued_cn <- as.data.frame(t(sapply(unique(cn2$case), duplicate_info_cn, cn2)), stringsAsFactors=FALSE)
colnames(uniqued_cn) <- colnames(cn2)[c(1,2,3)]

## binarize for Umberto
bin_cn <- as.data.frame(apply(uniqued_cn, 2, function(x) { ifelse(x < 3 | is.na(x), 'WT', 'AMPL') } ))
# add back in CRC1185 for HER2
bin_cn[rownames(bin_cn)=="CRC1185", 'HER2'] <- 'AMPL'
#write.table(bin_cn, "/home/egrassi/binarized_cn_210127_scatter_annotated_manual.tsv", sep="\t", quote=FALSE)

duplicate_info_mut <- function(case, data) {
  subset <- data[data$CASE == case,]
  if (nrow(subset)==1) {
    return(as.matrix(subset))
  } else {
    # for duplicates in muts we keep the LM
    return(as.matrix(subset[subset$X == "LM",]))
  }
}

uniqued_mut <- as.data.frame(t(sapply(unique(mut$CASE), duplicate_info_mut, mut)), stringsAsFactors=FALSE)
colnames(uniqued_mut) <- colnames(mut)
uniqued_mut$X <- NULL
uniqued_mut$GenealogyID <- NULL
uniqued_mut$CASE <- NULL
uniqued_mut$flag_H.X <- NULL

#uniqued_mut2 <- as.data.frame(apply(uniqued_mut, 2, function(x) { ifelse(is.na(x),'wt', x) } ))
## there is a NA for PI3KCA but it has a KRAS mutation so that's fine for fourwt info
merge_pdata <- merge(mdata, bin_cn, by="row.names", all.x=TRUE)
rownames(merge_pdata) <- merge_pdata$Row.names
merge_pdata$Row.names <- NULL
merge_pdata2 <- merge(merge_pdata, uniqued_mut, by="row.names", all.x=TRUE)
# apply(merge_pdata2, 2, class) ?? character why
infofour <- apply(merge_pdata2[,c('KRAS', 'NRAS', 'BRAF', 'PIK3CA')], 1, function(x) {any(x!="wt")})
merge_pdata2$four_wt <- ifelse(infofour, '4_KO', '4_WT') 
# there is a NA for PI3KCA but it has a KRAS mutation so that's fine

lmp <- function(model) {
  #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  #f <- summary(modelobject)$fstatistic
  f <- model$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


ggplotRegression <- function(fit, data, name) {
  print(ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
          geom_point(data=data, aes(fill=four_wt, size=as.factor(HER2), shape=as.factor(EGFR),
                     colour=as.factor(MET))) +
          stat_smooth(method = "lm", col = "darkgrey") +
          labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                             " P =",signif(summary(fit)$coef[2,4], 5)))+
          theme_bw()+theme(text = element_text(size=15))+
          scale_fill_manual(values=c('red', 'blue'))+
          scale_colour_manual(values=c('black', 'orange'))+
          scale_shape_manual(values=c(22, 24))+ ## 22 24 with borders
          scale_size_manual(values=c(2, 4)))
  #ggsave(paste0(name, ".svg"))
}



compute_lm <- function(xcol, data, ycol) {
    model <- lm(data=data, formula=as.formula(paste0(ycol, "~", xcol)))
    ggplotRegression(model, data, xcol)
    sm <- summary(model)
    r2 <- sm$r.squared
    ar2 <- sm$adj.r.squared
    pval <- lmp(sm) 
    ll <- logLik(model)
    return(c(pval, r2, ar2,ll))
}

compute_lm('CTG_5000', merge_pdata2, 'Cetuximab_dVw3')
#res <- t(sapply(xcolumns, compute_lm, mdata, yname))
colnames(res) <- c('pvalue','R2','adjR2','logLik')
write.table(data.frame('exp'=rownames(res), res), file=output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
