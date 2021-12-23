library(ggplot2)
set.seed(42)


#unmute_theme <- readRDS(file=snakemake@input[['unmute']])
#mute_theme <- readRDS(file=snakemake@input[['mute']])
unmute_f <- snakemake@output[['unmute']]
mute_f <- snakemake@output[['mute']]
load(snakemake@input[['Rimage']])
eval(parse(text=myriad))
##############################
#library(showtext)
#size <- as.numeric(snakemake@wildcards[['size']])
#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()
#############################

pdata <- data.frame(x=rnorm(100), y=rnorm(100))

ggplot(data=pdata, aes(x, y))+geom_point()+ggtitle('unmute')+unmute_theme
ggsave(unmute_f)
ggplot(data=pdata, aes(x, y))+geom_point()+ggtitle('mute')+mute_theme
ggsave(mute_f)