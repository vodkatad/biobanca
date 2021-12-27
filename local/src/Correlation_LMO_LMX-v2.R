#!/usr/bin/env Rscript
### Correlation LMO-LMX
### Script to calculate correlation btw LMO and LMX, produce correlation and scatter plot
library(tidyverse)
library(ggplot2)
library(tibble)
library(getopt)
library(pheatmap)

options(warn=1) # to print warnings as they occurr
#  'thr', 't', 1, 'numeric'
opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'xeno', 'x', 1, 'character',
  'pdo', 'o', 1, 'character',
  'outrsq', 'r', 1, 'character',
  'outdir', 'd', 1, 'character',
  'matrix', 'm', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (!is.null(opt$help) | is.null(opt$xeno) | is.null(opt$pdo) | is.null(opt$outrsq) | is.null(opt$outdir)
| is.null(opt$matrix)) {
  cat(getopt(opts, usage=TRUE))
  stop('-x -o -r -d -m are all mandatory arguments!')
}

### load file
#lmo_f <- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/correlation_LMO_LMX/LMO_means_test2.tsv"
#lmx_f <- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/correlation_LMO_LMX/LMX_means_test2.tsv"
lmo_f <- opt$pdo
lmx_f <- opt$xeno
lmo_df_p <- read.table(gzfile(lmo_f), quote = "", sep = "\t", header = TRUE, row.names = 1)
lmx_df_p <- read.table(gzfile(lmx_f), quote = "", sep = "\t", header = TRUE, row.names = 1)
lmo_df <- as.data.frame(t(lmo_df_p))
lmx_df <- as.data.frame(t(lmx_df_p))

### filtering expression data: which sd genes but not clear outliers/not expressed genes
### filter not expressed genes
means_lmo <- apply(lmo_df, 1, mean)
means_lmx <- apply(lmx_df, 1, mean)
med_lmo <- median(means_lmo)
med_lmx <- median(means_lmx)

### I obtain a data.table that contain only the genes  mean expression over the median
lmo_filtered <- as.data.frame(lmo_df[means_lmo > med_lmo,])
lmx_filtered <- as.data.frame(lmx_df[means_lmx > med_lmx,])

### I want to calculate the standard deviation
sds_lmo <- apply(lmo_filtered, 1, sd) 
sds_lmx <- apply(lmx_filtered, 1, sd) 

### now we keep the top 10% variable genes 
sds_lmo <- sds_lmo[order(-sds_lmo)]
n_lmo <- length(sds_lmo)
keep_lmo <- head(sds_lmo, round(0.10*n_lmo))
keep_genes_lmo <- names(keep_lmo)
desd_lmo <- lmo_filtered[rownames(lmo_filtered) %in% keep_genes_lmo,]

sds_lmx <- sds_lmx[order(-sds_lmx)]
n_lmx <- length(sds_lmx)
keep_lmx <- head(sds_lmx, round(0.10*n_lmx))
keep_genes_lmx <- names(keep_lmx)
desd_lmx <- lmx_filtered[rownames(lmx_filtered) %in% keep_genes_lmx,]

###models <- c(models_o, models_x) modifica con desd lmo e desd lmx
###df <- as.data.frame(table(models))
###df$models <- as.character(df$models)

colnames(desd_lmo) <- paste0(colnames(desd_lmo), "_LMO")
colnames(desd_lmx) <- paste0(colnames(desd_lmx), "_LMX")

### merging the dataset

merged_lmo_lmx <- merge(desd_lmo, desd_lmx, by = "row.names")
merged_lmo_lmx <- column_to_rownames(merged_lmo_lmx, var = "Row.names")


## finding the models with both LMO and LMX
samples <- substr(colnames(merged_lmo_lmx),0,7)
freqdf <- as.data.frame(table(samples), stringsAsFactors=FALSE)
paired <- freqdf[freqdf$Freq==2,'samples']


results_rsquared <- data.frame(row.names=paired, rsquared=rep(NA, length(paired)), slope=rep(NA, length(paired)),
                               intercept=rep(NA, length(paired)), pvalue=rep(NA, length(paired)))

setwd(opt$outdir)

save.image("residuals.Rdata")

### for loop
for (i in seq(1, length(paired))) {
  wanted_column_o <- paste0(paired[i],"_LMO")
  wanted_column_x <- paste0(paired[i],"_LMX")
  #expr_o <- merged_lmo_lmx[,colnames(merged_lmo_lmx) == wanted_column_o]
  #expr_x <- merged_lmo_lmx[,colnames(merged_lmo_lmx) == wanted_column_x]
  #df_single_loop <- data.frame(row.names=rownames(merged_lmo_lmx), x=expr_x, o=expr_o)
  #fit <- lm(data=df_single_loop, formula="o~x")
  fit <- lm(data=merged_lmo_lmx, formula=paste0(wanted_column_o, "~", wanted_column_x))
  sfit <- summary(fit)
  #results_rsquared[rownames(results_rsquared) == paired[i],'rsquared'] <- sfit$r.squared
  results_rsquared[i,'rsquared'] <- sfit$r.squared
  results_rsquared[i, 'slope'] <- sfit$coefficients[2,1]
  results_rsquared[i, 'intercept'] <- sfit$coefficients[1,1]
  #results_rsquared[rownames(results_rsquared) == paired[i], 'slope'] <- sfit$coefficients[2,1]
  #results_rsquared[rownames(results_rsquared) == paired[i], 'intercept'] <- sfit$coefficients[1,1]
  f <- sfit$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  results_rsquared[i, 'pvalue'] <- p
  residuals <- rstudent(fit)
  merged_lmo_lmx$outliers <- ifelse(abs(residuals)>2, "yes", "no") 
  ggplot(data=merged_lmo_lmx, aes_string(y=wanted_column_o, x=wanted_column_x))+geom_point(aes(color=outliers))+theme_bw()+stat_smooth(method="lm")+scale_color_manual(values = c("black", "red"))
  ggsave(paste0(paired[i], ".pdf"))
  residual <- as.data.frame(residuals, rownames=merged_lmo_lmx) ### row.names?
  write.table(residual, file = paste0(paired[i],'.tsv'), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
}

results_correlations <- as.data.frame(results_rsquared, quote="", sep ="\t", header = TRUE)
setwd('..')
write.table(results_correlations, file=opt$outrsq, 
            sep='\t', col.names = TRUE, row.names= TRUE, quote = FALSE)



# We can use cor if we are only interested in the pearson estimate (not pvalue or other things)
# We need to remove the not 'matched' pairs before and be sure that the columns are ordered in the same way
# check that the ordering works

colnames(desd_lmo) <- substr(colnames(desd_lmo),0,7)
colnames(desd_lmx) <- substr(colnames(desd_lmx),0,7)
clmo <- desd_lmo[, paired]
clmx <- desd_lmx[, paired]

### check also for genes
all_genes <- intersect(rownames(clmo), rownames(clmx))
clmo <- clmo[all_genes,]
clmx <- clmx[all_genes,]


if (all(colnames(clmo)!=colnames(clmx)) & all(rownames(clmo)!=rownames(clmx))) {
  stop('Brutto llama!')
}

res1 <- cor(clmo, clmx)
#res <- data.frame(pearson=diag(res1), row.names=rownames(res1))
write.table(res1, file=opt$matrix, 
            sep='\t', col.names = TRUE, row.names= TRUE, quote = FALSE)