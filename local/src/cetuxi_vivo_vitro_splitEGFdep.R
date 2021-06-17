th <- function() {
  textSize <- 1.5
  current_theme <-
    theme_bw() +
    theme(
      strip.text = element_text(size = rel(textSize)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title = element_text(size = rel(1.8)),
      axis.text.x = element_text(size=rel(1.7)),
      axis.text.y = element_text(angle = 0,
                                 size = rel(1.7)),
      axis.line = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.length.y.left = unit(3,'mm'),
      
      legend.position = "top",
      legend.justification = "right",
      #legend.margin = margin(unit(0, "cm")),
      legend.title = element_text(size = rel(textSize), face = "bold"),
      legend.text = element_text(size = rel(1.2)),
      legend.background = element_rect(size=0.5, linetype="solid", color="black"),
      plot.title = element_text(
        face = "bold",
        size = rel(2),
        hjust = 0.5
      ),
      panel.border = element_blank(),
      plot.caption = element_text(size=rel(1))
    )
  current_theme
}
current_theme <- th()
library(ggplot2)

data <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/cetuximab/atp/merge.tsv', sep="\t", header=TRUE)
edep <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/local/share/data/EGF_dep', sep="\t", header=F)
data$slower <- 'No'
data[data$case %in% edep$V1,'slower'] <- "Yes"
ggplot(data=data, aes(x=CTG_5000, y=perc, color=slower))+geom_point()+current_theme+stat_smooth(method="lm")+scale_color_manual(values=c("blue","red"))

ds <- data[data$slower=="Yes",]
dn <- data[data$slower=="No",]

summary(lm(data=ds, formula="perc~CTG_5000"))
summary(lm(data=dn, formula="perc~CTG_5000"))