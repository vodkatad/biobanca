load("/scratch/trcanmed/biobanca/dataset/V1/trans_sign/expr/table_expression_clustering.Rdata")

dist_pearson <- as.data.frame(cor(expr_selected, method="pearson"))
dist_pearson$type <- substr(rownames(dist_pearson), 8, 10)
dist_pearson <- dist_pearson %>% filter(!type == "LMX")
dist_pearson$type <- NULL
dist_pearson <- as.data.frame(t(dist_pearson))
dist_pearson$type <- substr(rownames(dist_pearson), 8, 10)
dist_pearson <- dist_pearson %>% filter(type == "LMX")
dist_pearson$type <- NULL

#prova <- dist_pearson
## to find the coupled means
means <- c()
for (i in unique(substr(rownames(dist_pearson), 1, 7))) {
 d <- dist_pearson[,grepl(i, colnames(dist_pearson)), drop(FALSE)]
 d <- d[grepl(i, rownames(d)),]
 m <- mean(as.matrix(d)) 
 means <- c(means, m)
 #print(i)
}
### cercare nelle colonne lo stesso caso 
# crea un vettore con la media tra replicati
# questa media si aggiunge un altro vettore che poi subità la media (media della media)

meansofmeanscouple <- mean(means)

# per ogni modello prendiamo delle sue righe tutte le colonne in cui non appare
# e poi segue uguale

means_lmxlmo_notcouple <- c()
for (i in unique(substr(rownames(dist_pearson), 1, 7))) {
  d <- dist_pearson[,!grepl(i, colnames(dist_pearson)), drop(FALSE)]
  d <- d[grepl(i, rownames(d)),]
  m <- mean(as.matrix(d)) 
  means_lmxlmo_notcouple <- c(means_lmxlmo_notcouple, m)
  #print(i)
}

tdist_pearson <- as.data.frame(t(dist_pearson))
means_lmolmx_notcouple <- c()
for (i in unique(substr(rownames(tdist_pearson), 1, 7))) {
  d <- tdist_pearson[,!grepl(i, colnames(tdist_pearson)), drop(FALSE)]
  d <- d[grepl(i, rownames(d)),]
  m <- mean(as.matrix(d)) 
  means_lmolmx_notcouple <- c(means_lmolmx_notcouple, m)
  #print(i)
}

meansnotcouple <- mean(c(means_lmxlmo_notcouple, means_lmolmx_notcouple))
