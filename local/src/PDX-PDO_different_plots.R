#!/usr/bin/env Rscript

library(getopt)
library(plotly)
library(networkD3)
library(psych)


set.seed(42)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'df_in', 'i', 1, 'character',
  'sankey_out', 's', 1, 'character',
  'classes_out', 'c', 1, 'character',
  'switch_out', 'w', 1, 'character',
  'switched_out', 'o', 1, 'character',
  'kappa_out', 'k', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)


df <- read.table(gzfile(opt$df_in), sep='\t', quote="", header=TRUE) # *at the end there is the cohen.kappa()

lmx_freq <- as.data.frame(table(df$prediction_LMX))
lmo_freq <- as.data.frame(table(df$prediction_LMO))

colnames(lmx_freq)[1] <- "CRIS"
colnames(lmx_freq)[2] <- "Freq"
lmx_freq$type <- rep("LMX",5)
colnames(lmo_freq)[1] <- "CRIS"
colnames(lmo_freq)[2] <- "Freq"
lmo_freq$type <- rep("LMO",5)

lmx_freq[,"%CRIS"] <- lmx_freq$Freq/nrow(df)*100 
lmo_freq[,"%CRIS"] <- lmo_freq$Freq/nrow(df)*100 

df2 <- rbind(lmx_freq,lmo_freq)
df2$`%CRIS` <- round(df2$`%CRIS`, 2)
df2$type <- factor(df2$type, levels=c('LMX','LMO'))

switched <- as.data.frame(table(df$switched))
colnames(switched) <- c("switch","freq")
switched$perc <- switched$freq/nrow(df)*100
switched$perc <- round(switched$perc, 2)


### plot for LMX-LMO frequency
png(opt$classes_out, width=8, height=5, units="in", type="cairo", res=300)
ggplot(data=df2, aes(x=CRIS, y=Freq, fill=type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=`%CRIS`), vjust=-0.3, position = position_dodge(0.9), size=3.5) +
  labs(title="LMX-LMO CRIS class", x="CRIS", y = "Freq") +
  scale_fill_manual(values=c("firebrick4","dodgerblue4")) +
  theme_minimal()
dev.off() 

### plot for switcher number
png(opt$switch_out, width=8, height=5, units="in", type="cairo", res=300)
ggplot(data=switched, aes(x=switch, y=freq, fill=switch)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=perc), vjust=-0.3, size=3.5) +
  labs(title="Switched CRIS LMX-LMO", x="Switched", y = "Freq") +
  scale_fill_manual(values=c("orange","navyblue")) +
  theme_minimal()
dev.off()

### plot for bars switching
switchtype <- as.data.frame(table(df$switch_type))
colnames(switchtype) <- c("switch_type","freq")
switchtype$perc <- switchtype$freq/nrow(df)*100
switchtype$perc <- round(switchtype$perc, 2)

png(opt$switched_out, width=8, height=5, units="in", type="cairo", res=300)
ggplot(data=switchtype, aes(x=switch_type, y=freq, fill=switch_type)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=freq), vjust=0.3, hjust=-1, size=3.5) +
  labs(title="Switch type between LMX and LMO", x="Switch type", y = "Freq") +
  coord_flip() +
  # scale_x_continuous(limits=order(factor(levels=switchtype$switch_type))) +
  theme_minimal()
dev.off()


### Sankey diagram preparation:
switchtype$source <- sub(" >.*", "", switchtype$switch_type)
switchtype$target <- sub(".*> ", "", switchtype$switch_type)
switchtype$target <- sub("-", ".", switchtype$target)

### From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(switchtype$source), as.character(switchtype$target)) %>% unique())
nodes$label <- nodes$name
nodes$label <- sub(".", "-", nodes$label, fixed=TRUE)
nodes$label <- factor(nodes$label, levels=c("CRIS-A","CRIS-B","CRIS-C","CRIS-D","CRIS-E"))
nodes <- nodes[order(nodes$label),]
nodes$color <- c("darkorange","darkorange","firebrick","firebrick","darkblue","darkblue","forestgreen","forestgreen","cyan","cyan")
nodes <- nodes[order(nodes$name),]
# we order them to have source then target to have consistent colors (before used labels to have easy way do define colors

### With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
switchtype$IDsource=match(switchtype$source, nodes$name)-1 
switchtype$IDtarget=match(switchtype$target, nodes$name)-1

# png(opt$sankey_out)
 fig <- plot_ly(
   type = "sankey",
   orientation = "h",
  
   node = list(
     label = nodes$label,
     color = nodes$color,
     pad = 15,
     thickness = 20,
     line = list(
       color = "black",
       width = 0.5
     ),
     x = c(rep(1,5), rep(0,5)),
     y = c(rep(0,5), rep(1,5))
   ),
  
   link = list(
     source = switchtype$IDsource,
     target = switchtype$IDtarget,
     value =  switchtype$perc
   )
 )
 fig <- fig %>% layout(
   title = "<b>CRIS class switch between LMX and LMO</b>",
   font = list(
     size = 13
   )
 )
# graphics.off()
# export(fig, file = opt$sankey_out)
###TODO: capire come generare l'immagine senza passare dall'html, o python libreria selenium
htmlwidgets::saveWidget(fig, file = opt$sankey_out)


### cohen.kappa() computation
df2 <- data.frame(row.names=df$LMX_lineage, PDX=df$prediction_LMX, PDO=df$prediction_LMO, stringsAsFactors=FALSE)
df3 <- data.frame(kappa=cohen.kappa(df2)$kappa, stringsAsFactors=FALSE)
write.table(df3, opt$kappa_out, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)