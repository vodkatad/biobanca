#!/opt/R/R-3.6.3/bin/Rscript

library(getopt)
library(ggsignif)
library(tidyverse)
library(ggrepel)
library(reshape)

set.seed(42)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'vsd_in', 'v', 1, 'character',
  'meda_in', 'm', 1, 'character',
  'simo_in', 's', 1, 'character',
  'output', 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)


meda <- read.table(gzfile(opt$meda_in), sep='\t', quote="", header=TRUE)
meda$sample_id <- gsub('.', '-', meda$sample_id, fixed = TRUE)

simo <- read.table(gzfile(opt$simo_in), sep='\t', quote="", header=TRUE)

vsd <- read.table(gzfile(opt$vsd_in), sep='\t', quote="", row.names=1, header=TRUE)
names(vsd) <- gsub('.', '-', names(vsd), fixed = TRUE)
rownames(vsd) <- str_remove(rownames(vsd), 'H_')
d <- vsd[,names(vsd) %in% meda$sample_id]

### filtering expression data: we want high sd genes but not clear outliers / not expressed genes
### filter not expressed genes
means <- apply(d, 1, mean)
med <- median(means)

de <- d[means > med,]

sds <- apply(de, 1, sd)

### now we keep thet top 10% variable genes 
sds <- sds[order(-sds)]

n <- length(sds)
keep <- head(sds, round(0.10*n)) #prova con e senza
keep_genes <- names(keep)
desd <- de[rownames(de) %in% keep_genes,]

samples <- names(desd)
### these are all the replicates (for the future)
replicates <- samples[grepl('2$', samples)] # 13, all LMO ## df con solo replicati e density plot delle due espressioni
replicates <- substr(replicates,1,12)

correlations_good <- as.data.frame(matrix(nrow=nrow(simo)*10, ncol=3)) #, stringsAsFactors = FALSE
correlations_bad <- as.data.frame(matrix(nrow=nrow(simo)*10, ncol=3))
correlations_repl <- as.data.frame(matrix(nrow=nrow(simo)*10, ncol=3))

colnames(correlations_good) <- c('pearson','pval','sample')
colnames(correlations_bad) <- c('pearson','pval','sample')
colnames(correlations_repl) <- c('pearson','pval','sample')

fpdom <- data.frame(row.names = rownames(desd))

added <- 1

desd1 <- desd
names(desd1) <- substr(names(desd1), 1, 12)

pdo_indexes <- grep('LMO', colnames(desd1), fixed=TRUE)
pdx_indexes <- grep('LMX', colnames(desd1), fixed=TRUE)

for (i in 1:nrow(simo)) {
    fpdo <- desd1[,colnames(desd1)==simo[i,'LMO_lineage'], drop=FALSE]
    fpdx <- desd1[,colnames(desd1)==simo[i,'LMX_lineage'], drop=FALSE]
    if ( ncol(fpdo) > 0 & ncol(fpdx) > 0 ) {
        if ( ncol(fpdo) > 1 ) {
            pe <- cor.test(fpdo[,1], fpdo[,2], method = "pearson")
            correlations_repl[added,] <- c(pe$estimate, pe$p.value,substr(names(fpdx)[1],1,7))
            # correlations_repl[added,'sample'] <- substr(names(fpdo)[1],1,7)
            fpdom[,1] <- rowMeans(fpdo)
            fpdo <- fpdom
        }
        for (i in 1:ncol(fpdo)) {
            for (j in 1:ncol(fpdx)) {
                pe <- cor.test(fpdo[,i], fpdx[,j], method = "pearson")
                correlations_good[added,] <- c(pe$estimate, pe$p.value,substr(names(fpdx)[1],1,7))
                # correlations_good[added,'sample'] <- substr(names(fpdo)[1],1,7)
                ### here we pick a wrong association and check that the two are not the good couple
                done <- FALSE
                while (!done) {
                    ir <- sample(pdo_indexes, 1, replace=FALSE)
                    jr <- sample(pdx_indexes, 1, replace=FALSE)
                    foundpdo <- colnames(desd1)[ir]
                    foundpdx <- colnames(desd1)[jr]
                    if (substr(foundpdo, 1, 7) != substr(foundpdx,1,7)) { # add lista di coppie già pescate e checkare non sia già stata presa
                        done <- TRUE
                    }
                }
                ### pdo random per ogni pdx? checkare che non si peschi due volte la stessa coppia sbagliata? fare più giri rispetto ai good?
                pe <- cor.test(desd1[,ir], desd1[,jr])
                correlations_bad[added,] <- c(pe$estimate, pe$p.value,substr(names(fpdx)[1],1,7))
                # correlations_good[added,'sample'] <- substr(names(fpdo)[1],1,7)
                added <- added + 1
            }
            
        }
    }
} # check che siano tutti sign i pval
# consistenza dei replicati dei pdo

correlations_good <- correlations_good[complete.cases(correlations_good),]
correlations_good$type <- ifelse(correlations_good$pearson < 0.82, "meh", "good") # thr deciso guardando il grafico
write.table(correlations_good, gzfile(opt$output), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
# correlations_good$type <- 'good'
correlations_bad <- correlations_bad[complete.cases(correlations_bad),]
correlations_bad$type <- 'bad'

labels <- correlations_good[correlations_good$pearson < 0.82,]$sample

correlations_repl <- correlations_repl[complete.cases(correlations_repl),]

correlations_good[,c('pval','sample')] <- NULL
correlations_bad[,c('pval','sample')] <- NULL

good <- melt(correlations_good)
bad <- melt(correlations_bad)
good$pearson <- as.numeric(good$pearson)
bad$pearson <- as.numeric(bad$pearson)

all <- rbind(good, bad)

long <- melt(all, id.vars = 'type')

p <- ggplot(data=long, aes(x=value, fill=type)) +
    geom_density(alpha=0.4) +
    xlab("Pearson") + ylab("Density") + 
    ggtitle("PDX-PDO genes expression correlations") +
    theme_bw() #+
    # geom_text(aes(label = bad), size = 5)

# png("Side/PDX-PDO_expression_correlations.png", width=14, height=7, units="in", type="cairo", res=300)
# p
# dev.off()
