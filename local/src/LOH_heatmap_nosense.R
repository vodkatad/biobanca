
library(tidyr)
library(pheatmap)
library(dplyr)
library(ggplot2)

res_path<-snakemake@input[['d']]
output_plot<-snakemake@output[['out']]
files <- list.files(path=res_path, pattern="*.tsv", full.names=TRUE, recursive=FALSE)

i=0
df<-read.table(files[1],quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df1_f<-read.table(files[1],quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

model<-strsplit(strsplit(files[1],split='/')[[1]][3],split = '_')[[1]][1]
df$model<-model
df1_f$model<-model
for (el in files){
  if(i!=0){
  model<-strsplit(strsplit(el,split='/')[[1]][3],split = '_')[[1]][1]
  tmp<-read.table(el,quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  tmp$model<-model
  df<-rbind(df,tmp)
  }
  i=1+1
}

df$start <- NULL
df$end <- NULL
df$event <- as.factor(df$event)

dfnochr <- df
dfnochr$chr <- NULL

result <- as.data.frame(dfnochr %>%
                          pivot_wider(names_from = segment_id, values_from = event))
rownames(result) <- result$model
result$model <- NULL
result[] <- lapply(result, function(x) {
  x <- gsub("NO", 0, x)
  x <- gsub("LOH_diff", 1, x)
  x<-gsub("only_late",2,x)
  x<-gsub("common",3,x)
  as.numeric(x)
})

an_col <- as.data.frame(matrix(ncol = 1, nrow = length(colnames(result))))
rownames(an_col) <- colnames(result)
an_col$V1 <- df1_f$chr
names(an_col)[names(an_col)=="V1"] <- "Chr"
#names(an_col)[names(an_col)=="V1"] <- "Enzyme"
an_col$Chr <- factor(an_col$Chr, 
                     levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 
                                'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 
                                'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 
                                'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'))

gappini <- df1_f
gappini2 <- gappini %>% group_by(chr)
gappini2 <- gappini2 %>% summarise(
  segment_id = max(segment_id)
)

righini <- c(1,2,3,4,5,6,7,8)


a<-pheatmap(result, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = an_col, 
         show_colnames = FALSE, legend_breaks = c(0,1,2,3), color = colorRampPalette(c("grey", "black", "red","blue"))(4),
         legend_labels = c("NO","LOH_diff","only_late","common"),gaps_col = gappini2$segment_id, gaps_row = righini) #,
ggsave(output_plot,plot = a,  width = 20,
       height = 20,
       units =  "cm",
       dpi = 300)
