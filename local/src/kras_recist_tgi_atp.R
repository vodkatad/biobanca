library(ggplot2)

tgi_f <- '/scratch/trcanmed/biobanca/dataset/V1/cetuximab_tgi/cetuxi_perc_w3_buoni.tsv'
volchange_f <- '/scratch/trcanmed/biobanca/dataset/V1/cetuximab/cetuxi_perc_w3_buoni.tsv'
vitro_f <- tgi <- '/scratch/trcanmed/biobanca/dataset/V1/cetuximab/pdo_cetuxi_buoni.tsv'
mut_f <- '/scratch/trcanmed/biobanca/dataset/V1/cetuximab/atp/threewt.tsv'

tgi <- read.table(tgi_f, sep='\t', header=TRUE)
volchange <- read.table(volchange_f, sep='\t', header=TRUE)
vitro <- read.table(vitro_f, sep='\t', header=TRUE)
mut <- read.table(mut_f, sep='\t', header=TRUE)

colnames(tgi)[2] <- 'tgi'
colnames(volchange)[2] <- 'dVperc'
nrow(tgi)
nrow(volchange)
nrow(vitro)
m1 <- merge(tgi, volchange, by='case')
nrow(m1)
m2 <- merge(m1, vitro, by='case')
nrow(m2)
m3 <- merge(m2, mut, by.x="case", by.y="smodel")
nrow(m3)

mid <- mean(m3$CTG_5000)
ggplot(data=m3, aes(x=dVperc, y=tgi, color=CTG_5000))+geom_point()+theme_bw()+scale_color_gradient(low="blue",
                                                                                                    high="red")

ggplot(data=m3, aes(x=dVperc, y=tgi, color=CTG_5000, shape=triple_wt))+geom_point()+theme_bw()+scale_color_gradient2(mid="grey", low="blue",
                                                                                                   high="red", midpoint=mid, space ="Lab")

mutonly <- m3[m3$triple_wt =="MUT",]


mid <- mean(mutonly$CTG_5000)

ggplot(data=mutonly, aes(x=dVperc, y=tgi, color=CTG_5000))+geom_point()+theme_bw()+scale_color_gradient2(mid="grey", low="blue",
                                                                                                    high="red", midpoint=mid, space ="Lab")