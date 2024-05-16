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

meansofmeanscouple <- mean(means)
summary(means)

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


summary(means_lmxlmo_notcouple)

wi <- wilcox.test(means,means_lmxlmo_notcouple)
wi$p.value
# We had to add filtering on lmx on rows/lmo on cols to avoid getting someone cor with himself and the same pairs in a different
# order. Apparently I already did so but lost the diary with this details :()
############################# STOP HERE #########################

tdist_pearson <- as.data.frame(t(dist_pearson))
means_lmolmx_notcouple <- c()
for (i in unique(substr(rownames(tdist_pearson), 1, 7))) {
  d <- tdist_pearson[,!grepl(i, colnames(tdist_pearson)), drop=FALSE]
  d <- d[grepl(i, rownames(d)),]
  m <- mean(as.matrix(d)) 
  means_lmolmx_notcouple <- c(means_lmolmx_notcouple, m)
  #print(i)
}

meansnotcouple <- mean(c(means_lmxlmo_notcouple, means_lmolmx_notcouple))

# test media insieme lo stesso - credo di no avendo nrow!=ncol in dist pearson, quindi mettere insieme
# tutti i pearson spaiati (su righe e colonne) e farne la media insieme piuttosto che mediare gli spaiati su righe,
# poi gli spaiati su colonne e poi mediarli insieme da un risulato leggermente diverso.

means_lmxlmo_notcouple2 <- c()
m1 <- c()
m2 <- c()
for (i in unique(substr(rownames(dist_pearson), 1, 7))) {
  d <- dist_pearson[,!grepl(i, colnames(dist_pearson)), drop=FALSE]
  d <- d[grepl(i, rownames(d)),]
  others_on_row <- as.numeric(as.matrix(d))
  d2 <- dist_pearson[, grepl(i, colnames(dist_pearson)), drop=FALSE]
  d2 <- d2[!grepl(i, rownames(d2)),]
  others_on_col <- as.numeric(as.matrix(d2))
  means_lmxlmo_notcouple2 <- c(means_lmxlmo_notcouple2, mean(c(others_on_row, others_on_col)))
  m1 <- c(m1, mean(others_on_row))  
  m2 <- c(m2, mean(others_on_col))  
  #print(i)
}

meansnotcouple2 <- mean(c(means_lmxlmo_notcouple2))
meansnotcouple3 <- mean(c(m1, m2))

summary(meansnotcouple)
summary(meansnotcouple3)
summary(meansnotcouple2)