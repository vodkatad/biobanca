library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(reshape)
library(tidyverse)
library(pheatmap)

expression <- snakemake@input[["expr"]]
bootstrapimage <- snakemake@output[["boot"]]
dendrogramimage <- snakemake@output[["dendro"]]
circulardendro <- snakemake@output[["circleden"]]
heatmap <- snakemake@output[["pheat"]]
resultscluster <- snakemake@output[["clust_res"]]
resultspvals <- snakemake@output[["pvals"]]
log_f <- snakemake@log[['log']]


save.image('ioodioilclustering.Rdata')

set.seed(500)

###dove c'è mille scrivere tot_boot

expr_selected <- read.table(expression, sep="\t", header=TRUE, stringsAsFactors = FALSE)
# read vsd_clustering.tsv.gz from /scratch/trcanmed/DE_RNASeq/dataset/Class_comparison_biobanca
# expr_selected
model <- substr(colnames(expr_selected), 1, 7)
class <- substr(colnames(expr_selected), 8, 10)
lineage <- substr(colnames(expr_selected), 1, 10)

ex <- as.data.frame(cbind(colnames(expr_selected), model, class, lineage))
colnames(ex) <- c("lineage", "model", "class", "shortlin")
pdoex <- ex
pdoex <- pdoex %>% filter(class == "LMO")
o <- pdoex$model
pdxex <- ex
pdxex <- pdxex %>% filter(class == "LMX")
x <- pdxex$model

difference <- setdiff(o, x)
pdoex <- pdoex %>% filter(model %in% difference)
crc <- pdoex$lineage

expr_selected <- expr_selected[,!colnames(expr_selected) %in% crc]
save.image('table_expression_clustering.Rdata')
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

df_cluster_colors <- clustered_data[match(order_df$pos, rownames(clustered_data)),]

#df_rainbow <- data.frame(id=seq(1, n_samples), color=rainbow(n_samples))

mycolors <- colors()
mycolors <- mycolors[!grepl('gr[ae]y\\d+', mycolors)]
mycolors <- mycolors[!grepl('dark', mycolors)]
mycolors <- mycolors[!grepl('light', mycolors)]
mycolors <- mycolors[!grepl('white', mycolors)]
mycolors <- mycolors[!grepl('cornsilk', mycolors)]
mycolors <- mycolors[!grepl('lavenderblush', mycolors)]
mycolors <- mycolors[!grepl('ivory', mycolors)]
mycolors <- mycolors[!grepl('bisque', mycolors)]
mycolors <- mycolors[!grepl('azure', mycolors)]
mycolors <- mycolors[!grepl('mistyrose', mycolors)]
mycolors <- mycolors[!grepl('snow', mycolors)]
mycolors <- mycolors[!grepl('seashell', mycolors)]
mycolors <- mycolors[!grepl('[a-z]+[1-3]', mycolors)]

mycolors <- c(mycolors, 'grey72', 'grey30')

df_rainbow <- data.frame(id=seq(1, n_samples), color=mycolors[sample.int(n_samples)])
df_rainbow <- df_rainbow[match(df_cluster_colors$cl_id, df_rainbow$id),]
#df_rainbow$color_noff <- substr(df_rainbow$color, 1,7)

dendo <- dend %>% set("labels_colors", c("white")) %>%
  set("leaves_pch", c(17)) %>%  # node point type
  set("leaves_cex", 2) %>%  # node point size
  set("leaves_col", as.character(df_rainbow$color)) # %>% #node point colors
  #plot(type="rectangle") %>% legend(x = "topleft", legend = labels_colors_species, fill = labels_colors)

dendo <- color_branches(dendo, col=dendro_colors) 

png(dendrogramimage) 
  plot(dendo, type="rectangle") %>% legend(x = "topleft", legend = labels_colors_species, fill = labels_colors)
graphics.off()

#dend %>% set("branches_k_color", value = seq(1, n_samples), k = n_samples) %>% 
  #plot(main = "Customized colors")

ggd1 <- as.ggdend(dendo)
ggplot(ggd1, labels = FALSE) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x") 
ggsave(circulardendro)

### new heatmap 

smodels <- unique(substr(colnames(expr_selected),0,7))
n_samples <- length(smodels)

mycolors <- colors()
mycolors <- mycolors[!grepl('gr[ae]y\\d+', mycolors)]
mycolors <- mycolors[!grepl('dark', mycolors)]
mycolors <- mycolors[!grepl('light', mycolors)]
mycolors <- mycolors[!grepl('white', mycolors)]
mycolors <- mycolors[!grepl('cornsilk', mycolors)]
mycolors <- mycolors[!grepl('lavenderblush', mycolors)]
mycolors <- mycolors[!grepl('ivory', mycolors)]
mycolors <- mycolors[!grepl('bisque', mycolors)]
mycolors <- mycolors[!grepl('azure', mycolors)]
mycolors <- mycolors[!grepl('mistyrose', mycolors)]
mycolors <- mycolors[!grepl('snow', mycolors)]
mycolors <- mycolors[!grepl('seashell', mycolors)]
mycolors <- mycolors[!grepl('[a-z]+[1-3]', mycolors)]

mycolors <- c(mycolors, 'grey72', 'grey30')

df_rainbow <- c(mycolors[sample.int(n_samples)])
names(df_rainbow) <- smodels

ann_colors <- list(class=c(X='firebrick1', O='darkblue'), smodel=df_rainbow)

cors <- cor(expr_selected, method="pearson")

annot <- data.frame(row.names=rownames(cors))
annot$smodel <- substr(rownames(annot), 0, 7)

annot$class <- substr(rownames(annot), 10, 10)

# Commentate perchè è calcolato in cima
#dist_pearson <- as.dist(1-cor(expr_selected, method="pearson"))
#hc1 <- hclust(dist_pearson, method = "complete" )

graphics.off()
pdf(heatmap)
p1 <- pheatmap(cors, cluster_rows = hc1, cluster_cols=hc1, 
         show_rownames = FALSE, show_colnames = FALSE, 
         annotation_col=annot, annotation_colors=ann_colors, treeheight_row=0)
p1
graphics.off()

#saveRDS(list(cors=cors, annot_col=annot, annotation_colors=ann_colors, hc1=hc1), file="test.rds")

## info for log

has_smt <- function(smodel, data, smt) {
  subset <- data[data$shortgen==smodel,]
  any(grepl(smt,  rownames(subset)))
}

all_pdo <- unique(as.character(clustered_data$shortgen))
table(sapply(all_pdo, has_smt, clustered_data, 'LMO'))
table(sapply(all_pdo, has_smt, clustered_data, 'LMH'))
table(sapply(all_pdo, has_smt, clustered_data, 'LMX'))

clustered_data$shortgen <- as.character(clustered_data$shortgen)
#clustered_data <- clustered_data[clustered_data$shortgen!="CRC1451",]
all_pdo <- unique(as.character(clustered_data$shortgen))

only_one <- function(smodel, data) {
  subset <- data[data$shortgen == smodel,]
  return(length(unique(subset$cl_id)))
}

ncl <- sapply(all_pdo, only_one, clustered_data)
sum(ncl==1)
length(ncl)
sum(ncl==1)/length(ncl)
sum(ncl>1)/length(ncl)
sum(ncl>1)

## majority vote for cluster
majority <- function(smodel, data) {
  subset <- data[data$shortgen == smodel,]
  df <- as.data.frame(table(subset$cl_id))
  m <- max(df$Freq)
  return(sum(df$Freq == m))
  #return(m)
}

singlemax <- sapply(all_pdo, majority, clustered_data)
ss <- singlemax[ncl!=1]
print('models appearing in two different clusters but with a clear majority vote')
length(ss[ss==1])
head(ss[ss==1])
print('models appearing in two different clusters without a clear majority vote')
length(ss[ss!=1])
head(ss[ss!=1])


majority_vote <- function(smodel, data) {
  subset <- data[data$shortgen == smodel,]
  df <- as.data.frame(table(subset$cl_id), stringsAsFactors=FALSE)
  df[order(-df$Freq),]
  return(df[1,'Var1'])
}

right_one <- lapply(all_pdo, majority_vote, clustered_data)
names(right_one) <- all_pdo
right_s <- 0
for (i in seq(1, nrow(clustered_data))) {
  if (clustered_data[i, 'cl_id'] == right_one[[clustered_data[i, 'shortgen']]]) {
    right_s <- right_s + 1
  }
}

sink(log_f, append=TRUE)
print(paste0("N sample in correct cluster: ", right_s))
print(paste0("N sample: ", nrow(clustered_data)))
sink()

## quali modelli correlano poco con gli altri

cors2 <- cors
meancors <- as.data.frame(rowMeans(cors2))
meancors$genealogy <- rownames(meancors)
meancors <- as.data.frame(meancors[order(meancors$`rowMeans(cors2)`, decreasing = FALSE),])
meancors$genealogy <- NULL
write.table(meancors, file = "/home/mferri/medie_correlazioni.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
