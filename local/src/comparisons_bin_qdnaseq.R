library(ggplot2)
#library(corrplot)
xeno_f <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/xeno_qdnaseq/cn_log2.tsv'
pdo_f <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/pdo_qdnaseq/cn_log2.tsv'
expected_n <- 142

xeno_df <- read.table(xeno_f, quote="", sep="\t", header=TRUE, row.names = 1)
pdo_df <- read.table(pdo_f, quote="", sep="\t", header=TRUE, row.names = 1)
xeno_df[,c("chromosome","start","end")] <- NULL
pdo_df[,c("chromosome","start","end")] <- NULL

# We no longer have Vale's clones
#qc <- qc[!grepl('-', qc$id, fixed=T),]
list_lmh <- colnames(pdo_df)[grepl('LMH', colnames(pdo_df))]
model_lmh <- substr(list_lmh, 0,7)

#TODO save them and compare also them
pdo_df <- pdo_df[, !colnames(pdo_df) %in% list_lmh]
pdo_df <- pdo_df[, !substr(colnames(pdo_df),0,7) %in% model_lmh]
xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% model_lmh]


# We mark the EGFR mutant with a different 'model'/short genealogy
colnames(pdo_df)[colnames(pdo_df)=="CRC0177LMO0A04008002D02000"] <- 'CRCE177LMO0A04008002D02000'
colnames(xeno_df)[colnames(xeno_df)=="CRC0177LMX0B05001TUMD06000"] <- 'CRCE177LMX0B05001TUMD06000'

#> setdiff(models_xeno, models_pdo)
#[1] "CRC1870" "CRC1875" "CRC2041"
# We know we were missing some PDOs TODO load from file if qc changes in the future
list_remove_xeno <- c("CRC1870", "CRC1875","CRC2041")
xeno_df <- xeno_df[, !substr(colnames(xeno_df),0,7) %in% list_remove_xeno]

models_xeno <- substr(colnames(xeno_df), 0, 7)
models_pdo <- substr(colnames(pdo_df), 0, 7)

if (length(intersect(models_xeno, models_pdo)) != expected_n) {
  stop('Wrong number of corresponding models!')
}

xeno_df <- xeno_df[, order(models_xeno)]
pdo_df <- pdo_df[, order(models_xeno)]

pdo_df <- 2**pdo_df
xeno_df <- 2**xeno_df
pearson <- cor(xeno_df, pdo_df)
#corrplot(pearson)

pheatmap(pearson, cluster_rows=F , cluster_cols=F, labels_col="PDO", labels_row="PDX", fontsize.number=1.5)

#>summary(diag(pearson))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01106 0.24362 0.41652 0.40662 0.57976 1.00000 
diag <- diag(pearson)

pearson2 <- pearson
diag(pearson2) <- rep(NA, length(diag))
all <- as.numeric(unlist(pearson2))
all <- all[!is.na(all)]
pdata <- data.frame(pearson=c(all, diag), type=c(rep('unmatched', length(all)),rep('matched', length(diag))))
ggplot(data=pdata, aes(x=pearson, color=type))+geom_density()+theme_bw()+theme(text=element_text(size=20))
ggplot(data=pdata, aes(y=pearson, color=type))+geom_boxplot()+theme_bw()+theme(text=element_text(size=20))

summary(diag)
summary(all)

sd_x <- apply(xeno_df, 2, sd)
sd_o <- apply(pdo_df, 2, sd)
head(sd_o[order(sd_o)])
head(sd_x[order(sd_x)])

interpdelta <- function(x) {
  qq <- quantile(x,  c(0.5, 0.95))
  return(qq[2]-qq[1])
}

sd_x <- apply(xeno_df, 2, interpdelta)
sd_o <- apply(pdo_df, 2, interpdelta)

#########

ff <- function (obj, fnames = NULL) 
{
  calls <- assayDataElement(obj, "segmented")
  segments <- assayDataElement(obj, "segmented")
  fd <- fData(obj)
  pd <- pData(obj)
  if (is.null(fnames)) 
    fnames <- pd$name
  if (length(fnames) != length(pd$name)) {
    stop("Length of 'fnames' is too short: ", length(fnames), 
         " != ", length(pd$name))
  }
  oopts2 <- options(scipen = 100)
  on.exit(options(scipen = oopts2), add = TRUE)
  for (i in 1:ncol(calls)) {
    d <- cbind(fd[, 1:3], calls[, i], segments[, i])
    sel <- d[, 4] != 0 & !is.na(d[, 4])
    dsel <- d[sel, ]
    rleD <- rle(paste(d[sel, 1], d[sel, 4], sep = ":"))
    endI <- cumsum(rleD$lengths)
    posI <- c(1, endI[-length(endI)] + 1)
    chr <- dsel[posI, 1]
    pos <- dsel[posI, 2]
    end <- dsel[endI, 3]
    score <- dsel[posI, 4]
    segVal <- round(dsel[posI, 5], digits = 2)
    bins <- rleD$lengths
    out <- cbind(fnames[i], chr, pos, end, bins, segVal)
    colnames(out) <- c("SAMPLE_NAME", "CHROMOSOME", "START", 
                       "STOP", "DATAPOINTS", "LOG2_RATIO_MEAN")
    fname <- paste(fnames[i], ".seg", sep = "")
    write.table(out, fname, quote = FALSE, sep = "\t", append = FALSE, 
                col.names = TRUE, row.names = FALSE)
  }
}

