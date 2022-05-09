library(onewaytests)
library(reshape)
library(ggplot2)
load('/scratch/trcanmed/biobanca/dataset/V1/drug_screening/anova_long.Rdata')
# bf.test(formula = value~ variable, data = long)
# bartlett.test(value ~ variable, data = long_ave)
data <- read.table('/scratch/trcanmed/biobanca/local/share/data/CRC0322_SBI_long.txt', sep = "\t", header=TRUE)
data2 <- data
long2 <- melt(data2)
ggplot(data=long2, aes(x=variable, y=value))+geom_point()

data$NT <- NULL
long <- melt(data)


bartlett.test(value ~ variable, data = long)


bf <- bf.test(formula = value~ variable, data = long)
paircomp(bf, adjust.method = "bonferroni")

library(DescTools)
DunnettTest(formula = value~ variable, data = long)

kruskal.test(formula = value~ variable, data = long)
pairwise.wilcox.test(x=long$value, g=long$variable)

pairwise.wilcox.test(x=long$value, g=long$variable, p.adjust.method = 'bonferroni')