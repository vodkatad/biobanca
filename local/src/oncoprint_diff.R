library(ggplot2)

xeno_af <- snakemake@input[['xenoTiers']]
pdo_af <- snakemake@input[['pdoTiers']]
xeno_genes <- snakemake@input[['xenoGenes']]
pdo_genes <- snakemake@input[['pdoGenes']]

thr <- as.numeric(snakemake@wildcards[['AF']])

osnakemake <- snakemake
load(snakemake@input[['Rimage']])
#eval(parse(text=myriad))
snakemake <- osnakemake

pdx <- read.table(xeno_af, header=TRUE, sep="\t", row.names=1)
pdo <- read.table(pdo_af, header=TRUE, sep="\t", row.names=1)


pdxbin <- ifelse(pdx > thr, 1, 0)
n_muts_x <- colSums(pdxbin)
pdobin <- ifelse(pdo > thr, 1, 0)




mergemut <- merge(pdo_df, xeno_df, by="row.names")#, all.x=TRUE, all.y=TRUE, fill=FALSE) # no need all the same genes:
#> setdiff(rownames(pdo_df), rownames(xeno_df))
#character(0)
rownames(mergemut) <- mergemut$Row.names
mergemut$Row.names <- NULL
mergemut <- mergemut[, colnames(mergemut) != "CRC0282.01.1.A"]

d2 <- t(mergemut)
d2 <- ifelse(d2, 1, 0)
d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w_CRC0078_PR.tsv', sep="\t", header=FALSE, stringsAsFactors = TRUE, row.names=1)
#also pdo then
cases <- rownames(d2)

#rimuovere i non mutati
#> sum(rowSums(d2) == 0)


d2pdo <- d2[grepl('LMO', rownames(d2)),]
# d3pdo <- data.frame(matrix(NA, ncol=ncol(d2pdo), nrow=nrow(d2pdo)))
# 
# for (i in seq(1, ncol(d2pdo))) {
#   f <- as.factor(d2pdo[,i])
#   levels(f) <- c("","OSNV")    
#   d3pdo[,i] <- f
# }

d2pdx <- d2[grepl('LMX', rownames(d2)),]
# d3pdx <- data.frame(matrix(NA, ncol=ncol(d2pdx), nrow=nrow(d2pdx)))
# 
# for (i in seq(1, ncol(d2pdx))) {
#   f <- as.factor(d2pdx[,i])
#   levels(f) <- c("","XSNV")    
#   d3pdx[,i] <- f
# }

td2x <- t(d2pdx)
td2o <- t(d2pdo)
colnames(td2x) <- substr(colnames(td2x), 0, 7)
colnames(td2o) <- substr(colnames(td2o), 0, 7)

td2x <- td2x[, colnames(td2x) != "CRC1875"]
#mat_list <- list(OSNV=td2o, XSNV=td2x)
both <- td2o & td2x
both2 <- t(apply(both, 1, as.numeric))
rownames(both2) <- rownames(both)
mat_list <- list(Both=both, PDO= td2o-both, PDX=td2x-both  ) 

#rownames(d3pdo) <- rownames(d2pdo)
#colnames(d3pdo) <- colnames(d2pdo)

#rownames(d3pdx) <- rownames(d2pdx)
#colnames(d3pdx) <- colnames(d2pdx)

#if (!all(colnames(d3pdo)==colnames(d3pdo))) {
#  print("LLLAMA!!!")
#}

#d3 <- rbind(d3pdx, d3pdo)
#d2 <- rbind(d2pdx, d3pdo)

#d3 <- apply(d2, 2, function(x) { f <- as.factor(x); levels(f) <- c("","SNV") })
#col = c("HOMDEL" = "blue", "AMP" = "red", "GAIN"= "#e86f0a", "MISSENSE" = "#008000", "STOP" = "#91f903", "INDEL" = "#01c601", "OTHER"="#9a4af9")
col = c("Both" = "blue", "PDO" = "red", "PDX"= "#0b7015")#, "XSNV" = "#ffcc00", "OTHER"="#9a4af9", "OSNV" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  #background = function(x, y, w, h) {
  #  grid.polygon(
  #    unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
  #    unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
  #    gp = gpar(fill = "grey", col = "white"))
  #  grid.polygon(
  #    unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
  #    unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
  #    gp = gpar(fill = "grey", col = "white"))
  #},
  # big blue
  Both = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Both"], col = NA))
  },
  # bug red
  PDO = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["PDO"], col = NA))
  },
  PDX = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["PDX"], col = NA))
  },
  # small green
  # MISSENSE = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["MISSENSE"], col = NA))
  # },
  # STOP = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["STOP"], col = NA))
  # },
  # INDEL = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
  #             gp = gpar(fill = col["INDEL"], col = NA))
  # },
  # OSNV = function(x, y, w, h) {
  #   grid.polygon(
  #     unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
  #     unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
  #     gp = gpar(fill = col["OSNV"], col = "white"))
  #   
  # },
  # XSNV = function(x, y, w, h) {
  #   grid.polygon(
  #     unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
  #     unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
  #     gp = gpar(fill = col["XSNV"], col = "white"))
  # },
  OTHER = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["OTHER"], col = NA))
  }
)

column_title = "OncoPrint for PDXO"
#heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "GAIN", "MISSENSE","STOP","INDEL","OTHER"), 
#                            labels = c("Deep deletion", "High/Focal amplification", "Non Focal amplification","Missense SVN","Stop SNV","Indel","Non-coding SNV"))

heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "GAIN", "SNV","OTHER"), 
                            labels = c("Deep deletion", "High/Focal amplification", "Non Focal amplification","Coding SVN","Non-coding SNV"))


dd <- d[rownames(d) %in% rownames(d3),, drop=FALSE] # volumes only of sequenced cases
sample_order <- rownames(dd) # the order of volumes (they are sorted when we load them)
#d4 <- d3[rownames(d3)%in% sample_order,] # muts only for cases with volumes
td4 <- t(d3)
#td4 <- td4[, sample_order] # we order them correctly
#library(circlize)
#col_fun = colorRamp2(c(35, -50, -100), c("red", "white", "blue"))

#oncoPrint(td4,
#          alter_fun = alter_fun, col = col, 
#          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
#          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
#                                             Irinotecan = dd$perc, col=list(Irinotecan=col_fun)),
#          row_names_gp = gpar(fontsize=8), column_order=sample_order)

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50, 4, ifelse(x > 35, 2, 8))
  #res[is.na(x)] <- 'black'
  return(as.numeric(res))
}
library(RColorBrewer)
#colOR <- c("PD"="red", "SD"="blue","OR"="black")

# 
colnames(dd) <- 'perc'
dd$perc <- dd$perc * 100
dd <- dd[order(-dd$perc), , drop=FALSE]

sample_order <- rownames(dd)
td4 <- td4[,sample_order]
oncoPrint(td4,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param, show_pct = FALSE,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                             Cetuximab = anno_barplot(dd$perc, gp = gpar(fill = get_recist(dd$perc))), height = unit(4, "cm")),
          row_names_gp = gpar(fontsize=8), column_order=sample_order)
#[1] 4

mat_list2 <- list()
#mat_list2[[1]] <- mat_list[[1]][,seq(1,30)] # y
#mat_list2[[2]] <- mat_list[[2]][,seq(1,30)]
#mat_list2[[3]] <- mat_list[[3]][,seq(1,30)]
mat_list2[[1]] <- mat_list[[1]]
mat_list2[[2]] <- mat_list[[2]]
mat_list2[[3]] <- mat_list[[3]]

names(mat_list2) <- names(mat_list)

m <- mat_list[[1]]

d <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/august2020/Treatments_Eugy_Ele_fix0cetuxi_201005_cetuxi3w_CRC0078_PR.tsv', sep="\t", header=FALSE, stringsAsFactors = TRUE, row.names=1)

length(intersect(rownames(d), colnames(m)))

s <- mat_list2[[1]]+mat_list2[[2]]+mat_list2[[3]]
su <- colSums(s)
su <- su[order(-su)]

oncoPrint(mat_list2, alter_fun = alter_fun, col = col)
oncoPrint(mat_list2, alter_fun = alter_fun, col = col, column_order = names(su))

oncoPrint(mat_list2,
          alter_fun = alter_fun, col = col, top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                                               Cetuximab = anno_barplot(dd$perc, gp = gpar(fill = get_recist(dd$perc))), height = unit(4, "cm")),
          row_names_gp = gpar(fontsize=8))

pp <- oncoPrint(mat_list, alter_fun = alter_fun, col = col);

png('oncoprint_test.png',width = 880, height = 880, units = "px"); print(pp); dev.off(