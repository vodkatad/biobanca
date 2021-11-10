library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

# read vsd_clustering.tsv.gz from /scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca
# expr_selected

dist_pearson <- as.dist(1-cor(expr_selected, method="pearson"))

#d <- dist(df, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(dist_pearson, method = "complete" )

smodels <- unique(substr(colnames(expr_selected),0,7))
n_samples <- length(smodels)

clusters <- as.data.frame(cutree(hc1, k = n_samples))
colnames(clusters) <- 'cl_id'
clusters$smodel <- substr(rownames(clusters), 0,7)

count_clu <- function(smodel, clustering) {
  clu <- clustering[clustering$smodel == smodel,]
  return(length(unique(clu$cl_id)))
}

res <- as.data.frame(sapply(smodels, count_clu, clusters)) # baseline

# bootstrap randomizzando colonna smodel di clustering
# ricalcolo con count_clu

### colori - https://cran.r-project.org/web/packages/dendextend/vignettes/FAQ.html
dend <- as.dendrogram(hc1)
dendro_ordered_samples <- order.dendrogram(dend)

order_df <- data.frame(pos=colnames(expr_selected)[dendro_ordered_samples])
order_df$specie <- substr(order_df$pos, 8, 10)

dendro_colors <- as.numeric(as.factor(order_df$specie))

plot(dend, leaflab="none")
#labels_colors(dend) <- dendro_colors
p <- color_branches(
  dend=dend,
  k = NULL,
  h = NULL,
  col=dendro_colors)
plot(p, leaflab="none")

# https://rdrr.io/cran/dendextend/man/color_branches.html better to try to color branches