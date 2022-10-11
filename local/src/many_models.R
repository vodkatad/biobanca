#!/usr/bin/env Rscript
library(ggplot2)
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


ggplotRegression <- function(fit, name) {
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
          geom_point() +
          stat_smooth(method = "lm", col = "blue") +
          labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                             " P =",signif(summary(fit)$coef[2,4], 5)))+
          theme_bw()+theme(text = element_text(size=15))
  ggsave(paste0(name, ".svg"))
}



compute_lm <- function(xcol, data, ycol) {
    model <- lm(data=data, formula=as.formula(paste0(ycol, "~", xcol)))
    n <- sum(is.na(data[,xcol]))
    ggplotRegression(model, xcol)
    sm <- summary(model)
    r2 <- sm$r.squared
    ar2 <- sm$adj.r.squared
    pval <- lmp(sm) 
    ll <- logLik(model)
    return(c(pval, r2, ar2,ll,n))
}

res <- t(sapply(xcolumns, compute_lm, mdata, yname))
colnames(res) <- c('pvalue','R2','adjR2','logLik', 'n')
write.table(data.frame('exp'=rownames(res), res), file=output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
