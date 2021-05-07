#!/usr/bin/env Rscript
options(error = function() traceback(3))
# Given a text file with a list of bam files produces diagnostic plots and bed files with log2 binned values.

#testing variables
#bamf <- '/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_shallowseq_pdo/align/markedDup_CRC0022LMO0A04012002D02000.sorted.bam'
#binsize <- 15

library(QDNAseq)
library(QDNAseq.hg38)
library(getopt)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'bam_files', 'f', 1, 'character',
  'cores', 'c', 2, 'integer',
  'bin_size', 'b', 1, 'integer'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

if (is.null(opt$bam_files) | !is.null(opt$help) | is.null(opt$bin_size)) {
    #stop('Usage: qdnaseq -f bamfiles -b binsize [-c cores] / -f -b are mandatory')
    cat(getopt(opts, usage=TRUE))
    stop('-f and -b are mandatory')
}

binsize <- opt$bin_size
list_bams <- opt$bam_files

# We load the list of bam files, single column no header
bams <- read.table(list_bams, header=FALSE, stringsAsFactors=FALSE)
bamf <- bams[,1]

# We define the binning that we want to use, limited by the one available in QDNAseq.hg38 and chosen considering the sequencing depth.
bins <- getBinAnnotations(binSize=binsize, genome="hg38")

# Do we go multicore?
if (!is.null(opt$cores)) {
    future::plan("multiprocess", workers=opt$cores)
}

# Read bam files, by default removes q<37, and marked as duplicates, we determined depth with the same parameters.
# FUTURE param for pairedEnd
sample_names <- basename(bamf)
sample_names <- gsub("markedDup_","", sample_names, fixed=TRUE)
sample_names <- gsub(".sorted.bam", "", sample_names, fixed=TRUE)

readCounts <- binReadCounts(bins, bamfiles=bamf, pairedEnds=TRUE, bamnames=sample_names)

#plot
#png("bin_counts_isobars.png")
#par(mfrow=c(1,2))
#plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
# default is not filtering bins on mappability and g/c, but then correction is done.
highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

#plot
#png("isobars.png")
#isobarPlot(readCountsFiltered)
#graphics.off()

# Here correction on mappability and g/c
readCountsFiltered <- estimateCorrection(readCountsFiltered)

#plot
png("noise_filtered.png")
noisePlot(readCountsFiltered)
graphics.off()

copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

#plot
#plot(copyNumbersSmooth)
#graphics.off()

#logTransform 	
#If TRUE (default), data will be log2-transformed.
#exportBins(copyNumbersSmooth, file="CRC0022LMO0A04012002D02000.txt")
#exportBins(copyNumbersSmooth, file="LGG150.igv", format="igv")

#e <- assayData(copyNumbersSmooth)

#e$copynumber

#> summary(e$copynumber[,1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.00    0.80    1.00    1.04    1.19    7.74   37422 

# -1000? come fai il log2 allora?

#exportBins(copyNumbersSmooth, file="%s_cn_log2.bed", format="bed")
exportBins(copyNumbersSmooth, file="cn_log2.tsv", format="tsv", type="copynumber", digits="no_rounding")


#Values below –1000 in each chromosome were floored to the first value greater than –1000 in the same chromosome. 

#Raw 306
#10log2ratio values were then segmented using the ASCAT22algorithm implemented in the ASCAT 307
#R package   v2.0.7.Whole-genome   sequence   reads  from EuroPDX BRCA   tumors  and 308corresponding tumo
save.image("qdnaseq.RData")