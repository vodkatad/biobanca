library(reshape)
library(tidyverse)

groups_f <- snakemake@input[['groups']]
mutmat_f <- snakemake@input[['mutmat']]
res_f <- snakemake@output[['fisher']]
log_f <- snakemake@log[['log']]
outfile <- snakemake@output[['fisher']]


save.image('pippo.Rdata')

muts <- read.table(mutmat_f, sep="\t", quote="", header=TRUE, stringsAsFactors = FALSE, row.names = 1)
muts <- ifelse(muts=="True", TRUE, FALSE)

groups <- read.table(groups_f, sep="\t", quote="", header=TRUE, stringsAsFactors = FALSE)
groups$group <- ifelse(groups$buoni == "TRUE", TRUE, FALSE)

common <- intersect(colnames(muts), groups$smodel)

sink(log_f)
print('With muts:')
ncol(muts)
print('With groups:')
nrow(groups)
print('With both:')
length(common)
sink()

muts <- muts[, common, drop=FALSE]
groups <- groups[groups$smodel %in% common, , drop=FALSE]
groups <- groups[match(common, groups$smodel), ,drop=FALSE]

for (i in seq(length(rownames(muts)))) {
  counts <- as.data.frame((rowSums(muts, na.rm = TRUE)))
}

colnames(counts) <- "n_true"
counts <- counts %>% filter(n_true > 4)
save_genes <- rownames(counts)
muts <- muts[rownames(muts) %in% save_genes,]

sink(log_f, append=TRUE)
print('Comparing: ')
table(groups$group)
sink()

if (! all(groups$smodel == colnames(muts))) {
  stop('Something fishy in matching mutation info with groups...')
}
                 
fisher_one_gene <- function(gene_tf, group_tf) {
  n_group <- sum(group_tf)
  n_not_group <- sum(!group_tf)
  n_mutated_in_group <- sum(gene_tf[group_tf])
  n_mutated_not_in_group <- sum(gene_tf[!group_tf])
  frac_mut_group <- n_mutated_in_group / n_group
  frac_mut_out_group <- n_mutated_not_in_group / n_not_group
  fi <- fisher.test(gene_tf, group_tf)
  return(c(fi$p.value, frac_mut_group, frac_mut_out_group, n_mutated_in_group, n_mutated_not_in_group, n_group, n_not_group))
}

at_least_one_muts <- muts[apply(muts, 1, function(x) any(x)),]
all_fisher_fracs <- apply(at_least_one_muts, 1, fisher_one_gene, groups$group)

res <- as.data.frame(t(all_fisher_fracs))
colnames(res) <- c('pvalue', 'frac_muts_group', 'frac_muts_outgroup', 'n_muts_group', 'n_muts_outgroup', 'n_group', 'n_outgroup')
res$p.adjust <- p.adjust(res$pvalue, method="BH")
res <- res[order(res$p.adjust),]

res$gene <- rownames(res)
res <- res[, c(ncol(res), seq(1, (ncol(res)-1)))]
write.table(res, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)