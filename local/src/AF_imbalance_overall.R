library(ggplot2)
library(ggpubr)
library(reshape)

# Inital investigation code to see if we could work on single muts instead than collapsing for whole genes
# in mut_burdes_af_wip.R, from lines 129 onwards.
ttest_f <- snakemake@output[['ttest']]
plots_d <- snakemake@params[['plots']]
n_thr <- as.numeric(snakemake@params[['muts_thr']])
genes_f <- snakemake@input[['genes']]

load(snakemake@input[['Rimage']])

# We get info on all muts (id -> gene)
muts_info <- unique(rbind(gpdo, gpdx))

# We merge x and o info and apply our AF thr (inherithed from Rimage), we have 0 for under AF thr or absent mut and the AF otherwise
# in a matrix with smodel.x as colnames for xenos and smodel.y for pdos.
# We use all.x/y=TRUE to keep also muts present in only x or o and fill with 0s there (this is redundant for this code but
# is useful to start with a reasonable data structure that do not lose information if in the future we want to add checks not only
# on common muts but consider and plot 0 for missing/under thr ones).
m1 <- merge(pdx, pdo, all.x = TRUE, all.y = TRUE, by='row.names')
rownames(m1) <- m1$Row.names
m1$Row.names <- NULL
m1[is.na(m1)] <- 0
m1[m1 < thr] <- 0


# High level logic:
#iterate on all genes and select the ones with muts common o-x in at least n_thr models
# iterate on genes
# for all its muts select the x-o common ones
# if more than 5 models are commonly mutated overall for the gene proceed

# function that will be called with a single gene and compute t.test, print plot if n_thr is reached
compare_o_x <- function(gene, mut_table, gene_muts, n_thr, log) {
  my_muts <- gene_muts[gene_muts$genes == gene, , drop=FALSE]
  common_x_o <- sapply(unique(rownames(my_muts)), compare_o_x_one_mut, mut_table)
  # we sum over columns to get the n of models mutated both in x o and if the overall sum if >= n_thr we can go on
  total_muts <- colSums(common_x_o)
  #if (sum(total_muts) >= n_thr) { # commented out cause here we look only at not NA genes from previous rule / AF_imbalance.R
    wanted <- names(total_muts[total_muts!=0]) # we are still not selecting for "both" here but just mutated
    # only the selected muts:
    data <- mut_table[rownames(mut_table) %in% wanted, ]
    # only the mutated samples:
    data <- data[, apply(data, 2, function(x){any(x!=0)})]
    # we now have a table with only the models with muts in the selected muts (those mutated in both x and o in at least 5 cumulative
    # models). We need to switch to the long format for ggplot and we still need to remove the samples that are mutated only in
    # x or o.
    
    data$id <- rownames(data)
    long <- melt(data, id.vars="id")
    # we can remove not mutated samples now # FIXME if we want to consider missing muts this needs to be changed
    long <- long[long$value != 0,]
    # and get better columns for x-o etc
    long$type <-  substr(long$variable, 9, 9)
    long$class <- ifelse(long$type == 'x', 'PDX', 'PDXT')
    long$smodel <- substr(long$variable, 0, 7)
    # group is a variable with mut id / model useful for pairing things (e.g. draw lines pairing dots in ggplot)
    long$group <- paste0(long$id, long$smodel)

    # logging purposes of the number of muts that are mutated in a single x or o and that right now we are removing
    # we can remove what's not in both x-o here with table on group
    #2287 chr12:25245350:C:A CRC1090.y 0.448    y   pdo CRC1090 chr12:25245350:C:ACRC1090 
    # was an example of these rogue cases.
    paired <- as.data.frame(table(long$group))
    withpairs <- paired[paired$Freq==2, 'Var1']
    long <- long[long$group %in% withpairs,]

    # todo fix y axis interval if we want to use in supplementary
    plot <- ggplot(data=long, aes(x=class, y=value))+geom_point()+
    geom_line(aes(group=group))+theme_bw()+ggtitle(gene)+unmute_theme+xlab('')+ylab('AF')
    #ggsave(paste0(gene, ".png"))

    # setup of data for the t-test, we need two vectors of AF with ordered by mut/smodel (available in group here)
    longx <- long[long$class == "PDX",]
    longo <- long[long$class == "PDXT",]
    longx <- longx[order(longx$group),]
    longo <- longo[order(longo$group),]
    if (! all(longx$group==longo$group)) { 
      # this should be taken care of by withpairs but I could have forgot some corner cases
      stop(paste0('Something messy in putting together pairs of muts! ', gene))
    }
    res <- data.frame(x=longx$value, o=longo$value)
    #oo <- long$value > 0.9 & long$class=="pdo"
    #o <- long$value < 0.7 & long$class=="pdx"
    #ooo <- intersect(long[o,'group'], long[oo,'group'])
    #long[long$group %in% ooo,]
    return(list(plot, res))
  #} else {
  #  return(NA)
  #}
}

# function that returns a T/F vector of n_pairs_models that is TRUE only if the given mut is mutated both in x and o
compare_o_x_one_mut <- function(mut, mut_table) {
  x <- mut_table[rownames(mut_table)==mut, grepl('.x', colnames(mut_table), fixed=TRUE)]
  y <- mut_table[rownames(mut_table)==mut, grepl('.y', colnames(mut_table), fixed=TRUE)]
  common <- x != 0 & y != 0
  return(common)
}

#muts_info$genes <- as.character(muts_info$genes)
genes <- read.table(genes_f, sep="\t", header=TRUE)
genes <- genes[!is.na(genes$pval),]
genes <- genes[order(genes$pval),]
ttests <- lapply(unique(genes$gene), compare_o_x, m1, muts_info, n_thr)

plots <- lapply(ttests, '[[', 1)
afs <- lapply(ttests, '[[', 2)

afs <- do.call(rbind, afs)
ttest <- t.test(afs$x, afs$o, paired=TRUE)
sink(ttest_f)
print(ttest$p.value)
sink()

j <- 1
setwd(plots_d) # created by previous rule
for (i in seq(1, length(plots), by=12)) {
  end <- i + 11
  if (end > length(plots)) {
    end <- length(plots)
  }
  pdf(paste0('set_', j, '.pdf'), width=11, height=16) # paper="a4" was not working, need to specify sizes explicitly to use all the space
  print(ggarrange(plotlist=plots[i:end], ncol=2, nrow=6)+ theme(plot.margin = margin(0,0,0,0, "cm"))) # margin does not do anything at all
  graphics.off()
  j <- j + 1
  #ggsave(file=paste0('set_', i, '.pdf'), g) 
  # messy and I'm bored https://stackoverflow.com/questions/17059099/saving-grid-arrange-plot-to-file
}
