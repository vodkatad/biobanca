load("qua.Rdata")

library(dendextend)

df_cluster_colors <- clustered_data[match(order_df$pos, rownames(clustered_data)),]

#df_rainbow <- data.frame(id=seq(1, n_samples), color=rainbow(n_samples))

mycolors <- colors()
mycolors <- mycolors[!grepl('gr[ae]y\\d+', mycolors)]
mycolors <- mycolors[!grepl('dark', mycolors)]
mycolors <- mycolors[!grepl('light', mycolors)]
mycolors <- mycolors[!grepl('white', mycolors)]
mycolors <- mycolors[!grepl('cornsilk', mycolors)]
mycolors <- mycolors[!grepl('lavenderblush', mycolors)]
mycolors <- mycolors[!grepl('ivory', mycolors)]
mycolors <- mycolors[!grepl('bisque', mycolors)]
mycolors <- mycolors[!grepl('azure', mycolors)]
mycolors <- mycolors[!grepl('mistyrose', mycolors)]
mycolors <- mycolors[!grepl('snow', mycolors)]
mycolors <- mycolors[!grepl('seashell', mycolors)]
mycolors <- mycolors[!grepl('[a-z]+[1-3]', mycolors)]

mycolors <- c(mycolors, 'grey72', 'grey30')

df_rainbow <- data.frame(id=seq(1, n_samples), color=mycolors[sample.int(n_samples)])
df_rainbow <- df_rainbow[match(df_cluster_colors$cl_id, df_rainbow$id),]
#df_rainbow$color_noff <- substr(df_rainbow$color, 1,7)



#save.image("~/egrassi/plotti.Rdata")
dendo <- dend %>% set("labels_colors", c("white")) %>%
  set("leaves_pch", c(17)) %>%  # node point type
  set("leaves_cex", 2) %>%  # node point size
  set("leaves_col", as.character(df_rainbow$color))  #%>% #node point colors
  #plot(type="rectangle") %>% legend(x = "topleft", legend = labels_colors_species, fill = labels_colors)


dendo <- color_branches(dendo, col=dendro_colors) 
#%>%
#  plot(type="rectangle") %>% legend(x = "topleft", legend = labels_colors_species, fill = labels_colors)


#dend %>% set("branches_k_color", value = seq(1, n_samples), k = n_samples) %>% 
  #plot(main = "Customized colors")




ggd1 <- as.ggdend(dendo)
ggplot(ggd1, labels = FALSE) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x") 

#scale 
