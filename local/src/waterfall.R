#!/usr/bin/env Rscript
library(ggplot2)
args <- commandArgs(trailingOnly = T)
input <- args[1]
output <- args[2]
library(ggplot2)

data <- read.table(input, sep="\t", header=TRUE)

textSize <- 1
current_theme <-
  theme_bw() +
  theme(
    strip.text = element_text(size = rel(textSize)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = 12),
    #axis.text.x = element_text(size=9, angle = 90, hjust = 1, vjust=0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0,
                               size = 9),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y.left = unit(3,'mm'),
    axis.line.x = element_blank(),
    axis.line.y = element_line(),
    legend.position = "top",
    legend.justification = "right",
    #legend.margin = margin(unit(0, "cm")),
    legend.title = element_text(size = rel(textSize), face = "bold"),
    legend.text = element_text(size = rel(textSize)),
    legend.background = element_rect(size=0.5, linetype="solid", color="black"),
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    panel.border = element_blank()
  )
ggplot(data, aes(y=perc,x=reorder(case, -perc), fill=CTG_5000))+geom_col()+ylab('(3w-0w)/0w')+xlab("")+current_theme+scale_fill_distiller(palette="RdYlBu")
ggsave(output)