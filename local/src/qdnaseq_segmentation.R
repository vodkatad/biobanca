#!/usr/bin/env Rscript
#options(error = function() traceback(3))

# Given an RData with qdnaseq analyses done up to smoothOutlierBins it performs default segmentation
# with DNACopy
# and calls with CGHcall ?


library(QDNAseq)
library(getopt)

#testing variables
#rdataf <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/qdnaseq/qdnaseq.RData'
#outf <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/qdnaseq/qdnaseq.seg'

#opts management
opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'rdata', 'r', 1, 'character',
  'outsegm', 's', 1, 'character',
  'outcalls', 'a', 2, 'character',
  'cores', 'c', 2, 'integer'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$rdata) | !is.null(opt$help) | is.null(opt$outsegm)) {
  cat(getopt(opts, usage=TRUE))
  stop('-r is mandatory')
}

rdataf <- opt$rdata
outsegmf <- opt$outsegm
outcallsf <- opt$outcalls

# Do we go multicore?
if (!is.null(opt$cores)) {
  future::plan("multiprocess", workers=opt$cores)
}

load(rdataf)

# using sqrt cause I do not like pseudocounts so small (.Machine$double.xmin)
# no wordS:
#https://github.com/ccagc/QDNAseq/issues/49

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
#prints NA NA NA if sample names have '-', transformed to . in the previous qdnaseq step.

copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
#plot(copyNumbersSegmented)

exportBins(copyNumbersSegmented, file=outsegmf, type="segments", format="tsv", logTransform=FALSE, digits="no") # filter bins on what?

if (!is.null(outcallsf)) {
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  exportBins(copyNumbersCalled, file=outcallsf, type="calls", format="tsv", logTransform=FALSE, digits="no")
}
# the warning does not seem bad:
#> warnings()
#Warning messages:
#  1: In allprior/tot :
#  Recycling array of length 1 in vector-array arithmetic is deprecated.
#Use c() or as.vector() instead.

#dat <- assayData(copyNumbersCalled)[["segmented"]]
#basically it segments then reproject to bins...
#> table(apply(dat,1, function(x){!is.na(x[1])}))
#FALSE   TRUE 
#37422 168475 
# reproject to bins...how/when does it do it?