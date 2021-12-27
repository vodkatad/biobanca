library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(reshape)
library(tidyverse)

expression <- snakemake@input[["expr"]]
bootstrapimage <- snakemake@output[["boot"]]
dedrogramimage <- snakemake@output[["dendro"]]
resultscluster <- snakemake@output[["clust_res"]]
resultspvals <- snakemake@output[["pvals"]]

#save.image('ioodioilclustering.Rdata')

set.seed(500)

###dove c'è mille scrivere tot_boot

expr_selected <- read.table(expression, sep="\t", header=TRUE, stringsAsFactors = FALSE)
# read vsd_clustering.tsv.gz from /scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca
# expr_selected

dist_pearson <- as.dist(1-cor(expr_selected, method="pearson"))

#d <- dist(df, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(dist_pearson, method = "complete" )

smodels <- unique(substr(colnames(expr_selected),0,7))
n_samples <- length(smodels)

### blocco da rendere funzione
count_clu_pt1 <- function(hc, n) {
  clusters <- as.data.frame(cutree(hc, k = n))
  colnames(clusters) <- 'cl_id'
  clusters$shortgen <- substr(rownames(clusters), 0,7)
  return(clusters)
}
clustered_data <- count_clu_pt1(hc1, n_samples)
colnames(clustered_data) <- c('cl_id', "shortgen")

write.table(clustered_data, file = resultscluster, quote = FALSE, col.names = TRUE, row.names = FALSE)
# writetalble

count_clu_pt2 <- function(smodel, clustering) {
  clu <- clustering[clustering$shortgen == smodel,]
  return(length(unique(clu$cl_id)))
}

res_true <- as.data.frame(sapply(smodels, count_clu_pt2, clustered_data))
colnames(res_true) <- c("count_clu")
#write.table(res, file = "count_clustering_model.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE)
### res conta in quanti cluster appare un caso, salvo anche questo
#res$type <- "TRUE"


# bootstrap randomizzando colonna smodel di clustering
tot_boot <- 1000
res_tot_boot <- data.frame(matrix(ncol = tot_boot, nrow = length(rownames(res_true))), row.names = rownames(res_true))
# creazione dataframe vuoto tot_boot colonne tot_campioni righe
shuffle_labels_compute_spread <- function(my_clustered_data) {
  cl_id_var <- my_clustered_data$cl_id
  n_el <- length(cl_id_var)
  cl_id_var_sh <- cl_id_var[sample.int(n=n_el, size=n_el)]
  my_clustered_data$cl_id <- cl_id_var_sh
  res <- as.data.frame(sapply(smodels, count_clu_pt2, my_clustered_data))
  return(res)
}

for (i in seq(1, tot_boot))  {
  res_tot_boot[,i] <- shuffle_labels_compute_spread(clustered_data)
}
# chiamata funzione che calcola in quanti cluster sta ogni smodel
# aggiunta di res a dataframe vuoto creato prima
# }

res_tb_mean <- as.data.frame(apply(res_tot_boot, 1, mean))
colnames(res_tb_mean) <- c("count_clu")
#res_tb_mean$type <- "FALSE"

#res_fin <- rbind(res, res_tb_mean)

##rbind non potendomi mettere gli stessi nomi sulle righe mi aggiunge un uno, tento un altro
## approccio trovato online
res1 <- res_true
#res1$`sapply(smodels, count_clu_pt2, clustered_data)` = as.numeric(as.character(res1$`sapply(smodels, count_clu_pt2, clustered_data)`))
res1 <- tibble::rownames_to_column(res1, "model")
#res1$type <- NULL
colnames(res1) <- c("model", "count_clu_true")
res_tb2 <- res_tb_mean
#res_tb2$type <- NULL
res_tb2 <- tibble::rownames_to_column(res_tb2, "model")
colnames(res_tb2) <- c("model", "count_clu_bootstrap")
merged <- merge(res1, res_tb2, by = "model")


# merged2 <- merged
# merged2 <- merged2 %>% remove_rownames %>% column_to_rownames(var = "model")
#merged$count_clu_bootstrap <- NULL
#res_tot_boot <- tibble::rownames_to_column(res_tot_boot, "model")

# 
# res_t_f <- data.frame(row.names = rownames(merged2))
# res_t_f$T_F <- NA
# #i=1
# for (i in seq(length(rownames(merged)))){
#   res_t_f$T_F <- ifelse(merged2$count_clu_true < merged2$count_clu_bootstrap, FALSE, TRUE)
# }

tot_boot_count <- merge(res_tot_boot, res_true, by = "row.names")
tot_boot_count <- tot_boot_count %>% remove_rownames %>% column_to_rownames(var = "Row.names")
counts <- ncol(tot_boot_count)

pvalue_boot <- data.frame(matrix(ncol = tot_boot, nrow = length(rownames(tot_boot_count))), row.names = rownames(tot_boot_count))

for (i in seq(1, (counts-1))) {
  pvalue_boot[,i] <- (tot_boot_count[,i] < tot_boot_count[,counts])
}
##usiamo il minore pecrcè vogliamo vedere quando la randomizzazione ha lavorato meglio
##del nostro clustering per ottenere un pvalue del bootstrap 

sumtf <- apply(pvalue_boot, 1, sum)
pvalue_counts <- data.frame(row.names = names(sumtf), pvalue = sumtf/tot_boot)
write.table(pvalue_counts, file = resultspvals, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)


gdpMelt <- reshape::melt(merged,id="model",
                measure.vars=c("count_clu_true","count_clu_bootstrap"))
#gdpMelt <- gdpMelt[order(gdpMelt$value, decreasing = TRUE), ]


#save.image('pippo.Rdata')
gdpMelt$type <- as.factor(substr(gdpMelt$variable,11,14))

gdpMelt <- gdpMelt[order(gdpMelt$variable),]
gdpMelt$ordervalue <- rep(gdpMelt[gdpMelt$variable=="count_clu_true",'value'],2)

ggplot(gdpMelt, aes(x = reorder(model, -ordervalue), value, fill = type)) + 
  geom_col(position=position_dodge(0.8), width=0.8) + 
  scale_fill_brewer(palette = "Set1") + theme(axis.text.x = element_blank())

ggsave(bootstrapimage, width=135, height=68, dpi=300, units="mm")

# merged <- merged %>% remove_rownames %>% column_to_rownames(var="model")
# fare una colonna type per true e false con tutto true e tutto FALSE
# poi rbind  i robi a cui ho aggiunto le etichette density aes(x=value, fill=etichetta) fa due curve
# geomcol in cui x rownames y value fill etichetta e geomcol mi fa due colori

# merge di dataframenonopiuvuoto e res by=row.names

### colori - https://cran.r-project.org/web/packages/dendextend/vignettes/FAQ.html
dend <- as.dendrogram(hc1)
dendro_ordered_samples <- order.dendrogram(dend)

order_df <- data.frame(pos=colnames(expr_selected)[dendro_ordered_samples])
order_df$specie <- substr(order_df$pos, 8, 10)

dendro_colors <- as.factor(order_df$specie)
labels_colors_species <- as.character(unique(as.character(dendro_colors))) # sorted like the leaf order

nspecies <- length(levels(dendro_colors))
if (nspecies == 2) { 
  # O X
  levels(dendro_colors) <- list(black = "LMO", red = "LMX")
} else {
  # H O X
  levels(dendro_colors) <- list(black = "LMO", red = "LMX", blue = "LMH")
}
labels_colors <- as.character(unique(as.character(dendro_colors))) # sorted like the leaf order


plot(dend, leaflab="none")
#labels_colors(dend) <- dendro_colors
p <- color_branches(
  dend=dend,
  k = NULL,
  h = NULL,
  col=dendro_colors)
png(dedrogramimage)
plot(p, leaflab="none")
dev.off()
# https://rdrr.io/cran/dendextend/man/color_branches.html better to try to color branches

save.image("qua.Rdata")