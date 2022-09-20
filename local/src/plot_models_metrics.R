#!/usr/bin/env Rscript
library(ggplot2)
library(reshape)
library(ggpubr)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
theme <- args[3]

load(theme)

data <- read.table(infile, sep="\t", header=TRUE)
data$R2 <- NULL
data$logLik <- NULL
data$pvalue <- -log10(data$pvalue)
m <- melt(data)
m$exp <- as.character(m$exp)
m$n <- as.numeric(sapply(strsplit(m$exp, "_"), function(x) {x[2]}))
m$expn <- sapply(strsplit(m$exp, "_"), function(x) {x[1]})
m <- m[order(m$expn, m$n),]

mp <- m[m$variable =="pvalue",]
mr <- m[m$variable =="adjR2",]

pplot <- ggplot(data=mp, aes(x=n, y=value, color=expn))+geom_point()+geom_line()+theme_bw()+theme(legend.position="none", text=element_text(size=15))+scale_x_log10()+scale_x_continuous(breaks=c(1250, 5000, 20000))+xlab('N.cells')+ylab("-Log10(pvalue)")+scale_color_manual(values=c('darkgoldenrod','darkgreen'))
Rplot <- ggplot(data=mr, aes(x=n, y=value, color=expn))+geom_point()+geom_line()+theme_bw()+theme(text=element_text(size=15), legend.position="bottom")+scale_x_log10()+scale_x_continuous(breaks=c(1250, 5000, 20000))+xlab('N.cells')+ylab("adjusted R2")+scale_color_manual(values=c('darkgoldenrod','darkgreen'))+labs(color="Measurement")

ggarrange(pplot, Rplot, labels = c("A", "B"), ncol = 1, nrow = 2)
ggsave(outfile)