bamf <- '/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_shallowseq_pdo/align/markedDup_CRC0022LMO0A04012002D02000.sorted.bam'
binsize <- 15

library(QDNAseq)
library(QDNAseq.hg38)

bins <- getBinAnnotations(binSize=binsize, genome="hg38")
##future::plan("multiprocess", workers=4)
# by default removes q<37, and marked as duplicates
readCounts <- binReadCounts(bins, bamfiles=bamf, pairedEnds=TRUE)

#plot
plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
# default not filtering bins on mappability and bases!!
highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

#plot
isobarPlot(readCountsFiltered)

# but here it's corrected
readCountsFiltered <- estimateCorrection(readCountsFiltered)

#plot
noisePlot(readCountsFiltered)

copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
#plot
plot(copyNumbersSmooth)

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

exportBins(copyNumbersSmooth, file="CRC0022LMO0A04012002D02000.bed", format="bed")


#Values below –1000 in each chromosome were floored to the first value greater than –1000 in the same chromosome. 

#Raw 306
#10log2ratio values were then segmented using the ASCAT22algorithm implemented in the ASCAT 307
#R package   v2.0.7.Whole-genome   sequence   reads  from EuroPDX BRCA   tumors  and 308corresponding tumo