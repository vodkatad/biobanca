#load("/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/table_expression_clustering.Rdata")
outf <- snakemake@output[['outf']]
load(snakemake@input[['rdata']])

dist_pearson <- as.data.frame(cor(expr_selected, method="pearson"))
dist_pearson$type <- substr(rownames(dist_pearson), 8, 10)
#dist_pearson <- dist_pearson %>% filter(!type == "LMX")
dist_pearson$type <- NULL
dist_pearson <- as.data.frame(t(dist_pearson))
dist_pearson$type <- substr(rownames(dist_pearson), 8, 10)
#dist_pearson <- dist_pearson %>% filter(type == "LMX")
dist_pearson$type <- NULL
#prova <- dist_pearson
## to find the coupled means
means <- c()
for (i in unique(substr(rownames(dist_pearson), 1, 7))) {
  d <- dist_pearson[,grepl(i, colnames(dist_pearson)), drop=FALSE]
  d <- d[grepl(i, rownames(d)),]
  d <- d[grepl('LMX', rownames(d)),grepl('LMO', colnames(d))]
  m <- mean(as.matrix(d)) 
  means <- c(means, m)
  #print(i)
}

# ma così teniamo x-x!
### cercare nelle colonne lo stesso caso 
# crea un vettore con la media tra replicati
# questa media si aggiunge un altro vettore che poi subità la media (media della media)
sink(outf, append=TRUE)

meansofmeanscouple <- mean(means)
summary(means)
sink()
# per ogni modello prendiamo delle sue righe tutte le colonne in cui non appare
# e poi segue uguale

means_lmxlmo_notcouple <- c()
for (i in unique(substr(rownames(dist_pearson), 1, 7))) {
  d <- dist_pearson[,!grepl(i, colnames(dist_pearson)), drop=FALSE]
  d <- d[grepl(i, rownames(d)),]
  d <- d[grepl('LMX', rownames(d)),grepl('LMO', colnames(d))]
  
  m <- mean(as.matrix(d)) 
  means_lmxlmo_notcouple <- c(means_lmxlmo_notcouple, m)
  #print(i)
}

sink(outf, append=TRUE)

summary(means_lmxlmo_notcouple)

wi <- wilcox.test(means,means_lmxlmo_notcouple)
wi$p.value
sink()