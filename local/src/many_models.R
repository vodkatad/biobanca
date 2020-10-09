#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = T)
x <- args[1]
y <- args[2]
output <- args[3]
prefix <- args[4]
wanted_x <- args[5]
merge_col <- args[6]
yname <- args[7]
# for each of the columns in wanted_x (comma separated) computes a linear model y ~ col after merging x and y on their merge_col
data_x <- read.table(x, sep="\t", header=TRUE)
# y should be a two column file
data_y <- read.table(y, sep="\t", header=TRUE)
if (ncol(data_y) != 2 ) {
    stop(paste0(y, ' needs to be a two columns file!'))
}
colnames(data_y)[2] <- yname
mdata <- merge(data_x, data_y, by=merge_col)
rownames(mdata) <- mdata[, merge_col]
mdata[,merge_col] <- NULL

xcolumns <- unlist(strsplit(wanted_x, ','))

lmp <- function(model) {
  #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  #f <- summary(modelobject)$fstatistic
  f <- model$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

compute_lm <- function(xcol, data, ycol) {
    model <- lm(data=data, formula=as.formula(paste0(ycol, "~", xcol)))
    sm <- summary(model)
    r2 <- sm$r.squared
    ar2 <- sm$adj.r.squared
    pval <- lmp(sm) 
    ll <- logLik(model)
    return(c(pval, r2, ar2,ll))
}

save.image('forplot.Rdata') # I stop here to troubleshoot the ggplot for the model with a live Rstudio
res <- t(sapply(xcolumns, compute_lm, mdata, yname))
colnames(res) <- c('pvalue','R2','adjR2','logLik')
write.table(data.frame('exp'=rownames(res), res), file=output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

