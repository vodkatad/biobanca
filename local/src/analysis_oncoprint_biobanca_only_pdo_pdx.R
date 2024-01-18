## plot mutazioni singole

library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(ggplot2)

r <- snakemake@input[["rdata"]]
plot_pdf <- snakemake@output[["pdf"]]
res <- snakemake@output[["chisq"]]

#load("/scratch/trcanmed/biobanca/dataset/V1/targeted/oncoprint_0.05.RData")
load(r)

pdo <- as.data.frame(rowSums(td2o-both))
colnames(pdo) <- "O"
pdo$genes <- rownames(pdo)
#pdo <- pdo %>% filter(O > 1)

pdx <- as.data.frame(rowSums(td2x-both))
colnames(pdx) <- "X"
pdx$genes <- rownames(pdx)
#pdx <- pdx %>% filter(X > 3)

merged <- merge(pdo, pdx, by ="genes")
merged$sum <- merged$O + merged$X
merged <- merged[order(-merged$sum),]
merged$X_negative <- -merged$X
merged <- merged %>% filter(sum > 1)

colours <- setNames(c("darkblue", "firebrick1"), 
                    c("PDXTs", "PDX"))

p <- ggplot(merged, aes(x = factor(genes, levels = merged$genes))) +
  geom_bar(aes(y = O, fill = "PDXTs"), stat = "identity", position = "identity") +
  geom_bar(aes(y = X_negative, fill = "PDX"), stat = "identity", position = "identity")+
  xlab("Genes")+ylab("Number of mutations") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = colours) + labs(fill="")

pdf(plot_pdf)
print(p)
dev.off()

#"PDO" = "darkblue", "PDX"= "firebrick1"
statistic <- merged[,c(1,2,3)]
rownames(statistic) <- statistic$genes
statistic$genes <- NULL

x <- c()
y <- c()
dato <- c()

for (i in rownames(statistic)){
  x <- statistic[i,] 
  y <- chisq.test(x)
  dato <- rbind(dato, y$p.value)
}  

statistic$pvalue <- dato
statistic <- statistic[order(statistic$pvalue),]
colnames(statistic) <- c("PDXTs", "PDXs", "pvalue")

write.table(statistic, res, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
