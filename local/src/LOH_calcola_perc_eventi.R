library(tidyr)
library(pheatmap)
library(dplyr)
library(ggplot2)


path_common<-"/scratch/trcanmed/biobanca/dataset/V1/LOH_pdx_pdoearly/check_good"
path_loh<-"/scratch/trcanmed/biobanca/dataset/V1/LOH_pdx_pdoearly/final_good"
chr_len<-"/scratch/trcanmed/biobanca/local/share/data/chr_len.tsv"

files <- list.files(path=path_loh, pattern="*.tsv", full.names=TRUE, recursive=FALSE)
files_common<-list.files(path=path_common,pattern = "*common_loh.tsv",full.names = TRUE,recursive = FALSE)
files_early<-list.files(path=path_common,pattern = "*nosense_loh.tsv",full.names = TRUE,recursive = FALSE)
chr_len_<-read.table(chr_len,quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
total<-sum(chr_len_$bp)
headers<-c('chr',	'start',	'end',	'intersection')
for (i in seq(1,length(files))){
  model<-strsplit(strsplit(files[i],split='/')[[1]][9],split = '_')[[1]][1]
  print(model)
  diff_<-read.table(files[i],quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(diff_)[ncol(diff_)]<-'intersection'
  diff_<-diff_[diff_$V1!='chrX' & diff_$V1!='chrY', ]
  common<-read.table(files_common[i],quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(common)[ncol(common)]<-'intersection'
  common<-common[common$V1!='chrX' & common$V1!='chrY', ]
  early<-read.table(files_early[i],quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(early)[ncol(early)]<-'intersection'
  early<-early[early$V1!='chrX' & early$V1!='chrY', ]
  
  t_diff<-((sum(diff_$intersection))/total)*100
  t_early<-((sum(early$intersection))/total)*100
  t_common<-((sum(common$intersection))/total)*100
  if (i == 1) {
    percentuali <- data.frame(modello=model, error=t_early, gained=t_diff,common=t_common,stringsAsFactors = FALSE)
  } else {
    riga<-c(model,t_early,t_diff,t_common)
    percentuali<-rbind(percentuali,riga,stringsAsFactors = FALSE)
  }
}
percentuali
