#!/usr/bin/env Rscript
#options(error = function() traceback(3))

# Given an RData with qdnaseq analyses done up to normalizeSegmentedBins it produces seg file
# not with calls but with log2FC

library(QDNAseq)
library(getopt)
library(Biobase)

#testing variables
#rdataf <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/pdo_qdnaseq/cn_seg.Rdata'
#outf <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/qdnaseq/qdnaseq.seg'

#opts management
opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'rdata', 'r', 1, 'character',
  'outseg', 's', 1, 'character'
  ), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$rdata) | !is.null(opt$help) | is.null(opt$outseg)) {
  cat(getopt(opts, usage=TRUE))
  stop('-r is mandatory')
}

rdataf <- opt$rdata
outsegf <- opt$outseg # loading rdata is always tricky...
load(rdataf)

my_segm_nocall <- function (obj, outfile) 
{
  segments <- assayDataElement(obj, "segmented")
  segments <- log2(segments+1)
  fd <- fData(obj)
  pd <- pData(obj)
  fnames <- pd$name
  for (i in 1:ncol(segments)) {
    d <- cbind(fd[, 1:3], (segments[, i]))
    sel <- !is.na(d[, 4]) # removed != 0 since we have segments (and did log2 before)
    dsel <- d[sel, ]
    rleD <- rle(paste(d[sel, 1], d[sel, 4], sep = ":"))
    endI <- cumsum(rleD$lengths)
    posI <- c(1, endI[-length(endI)] + 1)
    chr <- dsel[posI, 1]
    pos <- dsel[posI, 2]
    end <- dsel[endI, 3]
    #score <- dsel[posI, 4]
    segVal <- dsel[posI, 4]
    #segVal <- round(dsel[posI, 4], digits = 2)
    bins <- rleD$lengths
    out <- cbind(fnames[i], chr, pos, end, bins, segVal)
    colnames(out) <- c("SAMPLE_NAME", "CHROMOSOME", "START", 
                       "STOP", "DATAPOINTS", "LOG2_RATIO_MEAN")
    if (i == 1) {
      write.table(out, outfile, quote = FALSE, sep = "\t", append = FALSE, 
                col.names = TRUE, row.names = FALSE)
    } else {
      write.table(out, outfile, quote = FALSE, sep = "\t", append = TRUE, 
                  col.names = FALSE, row.names = FALSE)
    }
  }
}
my_segm_nocall(copyNumbersSegmented, outsegf)

