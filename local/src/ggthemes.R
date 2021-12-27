library(ggplot2)
#library(showtext)
size <- as.numeric(snakemake@wildcards[['size']])

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

textSize <- size
largerSize <- textSize + 5
unmute_theme_myriad <- theme_bw() +
theme(
	text = element_text(size = textSize, family = "myriad"),
	axis.title = element_text(size = largerSize),
	axis.text.x = element_text(size = textSize),#, angle = 90, vjust = 0.5, hjust=1)
	axis.text.y = element_text(size = textSize),
	plot.title = element_text(size = largerSize, hjust = 0.5)
)

#saveRDS(unmute_theme, file=snakemake@output[['unmute']])
unmute_theme <- theme_bw() +
theme(
	text = element_text(size = textSize, family='sans'),
	axis.title = element_text(size = largerSize),
	axis.text.x = element_text(size = textSize),#, angle = 90, vjust = 0.5, hjust=1)
	axis.text.y = element_text(size = textSize),
	plot.title = element_text(size = largerSize, hjust = 0.5)
)

#saveRDS(unmute_theme, file=snakemake@output[['unmute']])


mute_theme <- theme_bw() +
theme(
	text = element_blank()
)

myriad <- 'library(showtext); size <- as.numeric(snakemake@wildcards[["size"]]); font_add(family = "myriad", regular = snakemake@input[["myriad"]]);showtext_auto()'
save.image(snakemake@output[['Rimage']])
#saveRDS(mute_theme, file=snakemake@output[['mute']])

#pdata <- data.frame(x=rnorm(100), y=rnorm(100))
#ggplot(data=pdata, aes(x, y))+ggtitle('example')+geom_point()+unmute_theme
#ggsave('test.png')

