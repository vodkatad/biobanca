library(ComplexHeatmap)
library(ggplot2)

op_f <- snakemake@output[['op']]
op_data_f <- snakemake@output[['op_data']]
pie_f <- snakemake@output[['pie']]
thr <- as.numeric(snakemake@wildcards[['top']])
pdo_f <- snakemake@input[['pdo']]
xeno_f <- snakemake@input[['xeno']]

load(snakemake@input[['Rimage']])
#eval(parse(text=myriad))


pdobing <- read.table(pdo_f, sep="\t", row.names=1, stringsAsFactors=FALSE, header=TRUE)
pdxbing <- read.table(xeno_f, sep="\t", row.names=1, stringsAsFactors=FALSE, header=TRUE)
pdobing[,1] <- NULL
pdobing[,1] <- NULL
pdxbing[,1] <- NULL
pdxbing[,1] <- NULL


pdobing <- pdobing[, colnames(pdobing) != "CRC0177LMO0D04017001D02000"] # remove the mutated pdo
pdxbing <- pdxbing[,colnames(pdxbing) != "CRC0177LMX0B05001TUMD06000"] # remove the wt xeno

list_remove_xeno <- c("CRC1870", "CRC1875","CRC2041")
list_lmh <- colnames(pdobing)[grepl('LMH', colnames(pdobing))]
model_lmh <- substr(list_lmh, 0,7)

pdobing <- pdobing[, !colnames(pdobing) %in% list_lmh]
pdobing <- pdobing[, !substr(colnames(pdobing),0,7) %in% model_lmh]
pdxbing <- pdxbing[, !substr(colnames(pdxbing),0,7) %in% model_lmh]

#> setdiff(models_xeno, models_pdo)
#[1] "CRC1870" "CRC1875" "CRC2041"
# We know we were missing some PDOs TODO load from file if qc changes in the future
pdxbing <- pdxbing[, !substr(colnames(pdxbing),0,7) %in% list_remove_xeno]

colnames(pdobing) <- substr(colnames(pdobing), 0, 7)
colnames(pdxbing) <- substr(colnames(pdxbing), 0, 7)

mergemut <- merge(pdobing, pdxbing, by="row.names") ## TODO CHECK THEY HAVE ALL SAME GENES TODO
rownames(mergemut) <- mergemut$Row.names
mergemut$Row.names <- NULL
save.image('pippo.Rdata')
mergemut[mergemut==-1] <- 0
nampl <- apply(mergemut, 1, sum)
nampl <- nampl[order(-nampl)]
top <- head(nampl, n=thr)
#mergemut <- mergemut[names(top),]
mergemut <- mergemut[c('ERBB2','EGFR','MET','MYC', 'IGF2','FGFR1','FGFR2','IRS2','PTEN'),]
#> setdiff(rownames(pdo_df), rownames(xeno_df))
#character(0)

d2 <- t(mergemut)
#d2 <- ifelse(d2, 1, 0)
#d <- read.table(treat_f, sep="\t", header=FALSE, stringsAsFactors = TRUE, row.names=1)

#rimuovere i non mutati
#> sum(rowSums(d2) == 0)


d2pdo <- d2[grepl('.x', rownames(d2)),]
d2pdx <- d2[grepl('.y', rownames(d2)),]

td2x <- t(d2pdx)
td2o <- t(d2pdo)
colnames(td2x) <- substr(colnames(td2x), 0, 7)
colnames(td2o) <- substr(colnames(td2o), 0, 7)

both <- td2o & td2x
both2 <- t(apply(both, 1, as.numeric))
rownames(both2) <- rownames(both)
mat_list <- list(Both=both, PDO= td2o-both, PDX=td2x-both  ) 

col = c("Both" = "blue", "PDO" = "red", "PDX"= "#0b7015")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Both = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Both"], col = NA))
  },
  PDO = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["PDO"], col = NA))
  },
  PDX = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["PDX"], col = NA))
  }
)

column_title = "Xenturion OncoPrint"

##dd <- d[rownames(d) %in% rownames(d3),, drop=FALSE] # volumes only of sequenced cases
##sample_order <- rownames(dd) # the order of volumes
#d4 <- d3[rownames(d3)%in% sample_order,] # muts only for cases with volumes
#td4 <- t(d3)
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
#library(RColorBrewer)
#colOR <- c("PD"="red", "SD"="blue","OR"="black")

# 
#colnames(dd) <- 'perc'
#dd$perc <- dd$perc * 100
#dd <- dd[order(-dd$perc), , drop=FALSE]

s <- mat_list[[1]]+mat_list[[2]]+mat_list[[3]]
su <- colSums(s)
su <- su[order(-su)]

#oncoPrint(mat_list2, alter_fun = alter_fun, col = col)
op <- oncoPrint(mat_list, alter_fun = alter_fun, col = col, column_order = names(su),
                remove_empty_columns = TRUE, remove_empty_rows = TRUE)

#oncoPrint(mat_list2,
#          alter_fun = alter_fun, col = col, top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
#                                                                               Cetuximab = anno_barplot(dd$perc, gp = gpar(fill = get_recist(dd$perc))), height = unit(4, "cm")),
#          row_names_gp = gpar(fontsize=8))

#pp <- oncoPrint(mat_list, alter_fun = alter_fun, col = col)

#pdf(op_f)
#png(op_f, width = 880, height = 880, units = "px")
svg(op_f, pointsize=7)
print(op)
graphics.off()

save.image(op_data_f)


pd <- as.data.frame(sapply(mat_list, sum))
colnames(pd) <- "alterations"
pd$class <- rownames(pd)
ggplot(pd, aes(x="", y=alterations, fill=class))+ geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+unmute_theme+scale_fill_manual(values=c("blue", "red", "#0b7015"))
ggsave(pie_f)