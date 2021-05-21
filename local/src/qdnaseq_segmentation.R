#!/usr/bin/env Rscript
#options(error = function() traceback(3))

# Given an RData with qdnaseq analyses done up to smoothOutlierBins it performs default segmentation
# with DNACopy and calls with CGHcall

#testing variables
rdataf <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/qdnaseq/qdnaseq.RData'
outf <- '/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/qdnaseq/qdnaseq.seg'

library(QDNAseq)
library(getopt)

#opts management
opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'rdata', 'r', 1, 'character',
  'cores', 'c', 2, 'integer'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$rdata) | !is.null(opt$help)) {
  cat(getopt(opts, usage=TRUE))
  stop('-r is mandatory')
}

# Do we go multicore?
if (!is.null(opt$cores)) {
  future::plan("multiprocess", workers=opt$cores)
}

load(rdataf)

# using sqrt cause I do not like pseudocounts so small (.Machine$double.xmin)
# no wordS:
#https://github.com/ccagc/QDNAseq/issues/49

# if 
#sampleNames(copyNumbersSmooth) <- gsub('-',"_",sampleNames(copyNumbersSmooth)) 
#

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
#prints NA NA NA if sample names have '-'.

copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)

exportBins(copyNumbersSegmented, file=outf, type="segments")
# ? does not work, prints the same as exporting copyNumber signal

# TODO sink
copyNumbersCalled <- callBins(copyNumbersSegmented)

# the warning does not seem bad:
#> warnings()
#Warning messages:
#  1: In allprior/tot :
#  Recycling array of length 1 in vector-array arithmetic is deprecated.
#Use c() or as.vector() instead.

###
copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
copyNumbersCalled <- callBins(copyNumbersSegmented)

dat <- assayData(copyNumbersCalled)[["segmented"]]

basically it's not a real segmentation...
> table(apply(dat,1, function(x){!is.na(x[1])}))

FALSE   TRUE 
37422 168475 
