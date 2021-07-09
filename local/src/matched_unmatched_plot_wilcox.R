### Draw densit plot
library(ggplot2)
library(pheatmap)
res_f<- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/correlation_LMO_LMX/res.tsv"
res_df <- read.table(res_f, quote = "", sep = "\t", header = TRUE)

### prendo la diagonale, quindi i campioni matched con la funzione diag mentre gli unmatched li prendo con upper.tri per la
### parte sopra la diagonale e lower.tri per la parte inferiore alla diagonale

matched <- diag(as.matrix((res_df))) ### diag vuole la matrice perchè non converte da solo il data frame
unmatched <- c(res_df[upper.tri(res_df)], res_df[lower.tri(res_df)])

res_data <- data.frame(res_df=c(unmatched, matched), type=c(rep('unmatched', length(unmatched)), rep('matched', length(matched))))
ggplot(data=res_data, aes(x=res_df, color=type))+geom_density()+theme_bw()+theme(text=element_text(size=20))
ggsave(density_f)

#write.table(as.data.frame(pearson), file=pearson_f, quote=FALSE, sep="\t")
#save.image(rdata_f)

### test wilkcoxson
x <- matched
y <- unmatched
w <- wilcox.test(x, y)

wilcox <- data.frame(row.names = "wilcox.test", pvalue=w$p.value)
wilcox <- write.table(wilcox, quote = "", sep = "\t", row.names = TRUE, col.names = TRUE)


