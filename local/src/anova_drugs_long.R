library(readxl)
library(stringr)
library(reshape)
library(ggpubr)

# get dir and tsv from snakemake inputs
drugs_tables <- snakemake@input[["drugs_tables"]]
anova_f <- snakemake@output[["anova"]]
out_dir <- snakemake@output[["out_dir"]]
models <- c('CRC0059','CRC0322','CRC0327','CRC1331')

anova_pvals <- data.frame() #col.names=c('MODEL','DRUG','CONDITION','PVALUE')

save.image("anova_next.Rdata")
# do we need to create dir or does snakemake take care of it?
if (dir.exists(out_dir)) {
  unlink(out_dir)
} else {
  dir.create(out_dir)
}

setwd(out_dir)

#i = 1
for (i in seq(1, length(models))) {
  
  Drug_Screening_Tables <- read_excel(drugs_tables, sheet = models[i])
  data <- as.data.frame(Drug_Screening_Tables)
  data$NAME_DRUG <- str_trim(data$NAME_DRUG, "both") 
  data$NAME_DRUG <- str_replace_all(data$NAME_DRUG, " +", "_") 
  
  drugs <- unique(data$NAME_DRUG)
  type <- unique(data$TYPE)
  
  for (j in seq(1, length(drugs))) {
    
    for (k in seq(1, length(type))) {
      
      w_drug <- drugs[j]
      w_type <- type[k]
      
      subset <- data[data$NAME_DRUG == w_drug & data$TYPE == w_type, grepl('CET', colnames(data))]
      
      long <- melt(subset)
      # 
      # avg_sd <- as.data.frame(sapply(unique(long$variable), function(x) { y <- long[long$variable==x,'value']; return(c(mean(y), sd(y)))} ))
      # t_avg_sd <- as.data.frame(t(avg_sd))
      # t_avg_sd$DOSE <- unique(long$variable)
      # colnames(t_avg_sd) <- c('mean','sd','variable')
      # t_avg_sd$upper <- t_avg_sd$mean + t_avg_sd$sd
      # t_avg_sd$lower <- t_avg_sd$mean - t_avg_sd$sd
      # 
      # ggplot(t_avg_sd, aes(x=variable, y=mean)) +  geom_point(stat="identity", shape=1, size=3) +
      #   geom_segment(aes(y=lower, yend=upper, x=variable, xend=variable), size=0.6)+theme_bw()+
      #   geom_point(data=long, aes(x=variable, y=value), stat="identity", size=4, position=position_dodge(0.2))+
      #   theme(axis.text.x = element_text(size=15, angle=90, vjust=0.5, hjust=1), legend.position="none", axis.title.y=element_text(size=15))
      # 
      ggboxplot(long, x = "variable", y = "value",
                fill = "variable",
                ylab = "Viability", xlab = "Variable", add="jitter")
      ggsave(paste0("boxplotjitter_", models[i],"_", drugs[j], "_",type[k],'.png'))
      
      
      long$x <- seq(1, nrow(long))
      ggplot(data=long, aes(y=value, x=x, fill=variable))+geom_col()+theme_bw()+ggtitle(paste0(models[i]," ", drugs[j], " ",type[k]))
      ggsave(paste0("barplot_", models[i],"_", drugs[j], "_",type[k],'.png'))
    
      # Compute the analysis of variance
      res.aov <- aov(value ~ variable, data = long)
      # Summary of the analysis
      sumaov <- summary(res.aov)
      pvalue <- sumaov[[1]]$`Pr(>F)`[1]
      
      res <- c(models[i], drugs[j], type[k], pvalue)
      anova_pvals <- rbind(anova_pvals, res, stringsAsFactors=FALSE)

      #res.dunn <- dunn_test(data = long, formula = value ~ variable, p.adjust.method = "BH")
      #res.dunn.2 <- dunn.test(long$value, long$variable, method = "bh")
    }
  }
  
}

colnames(anova_pvals) <- c('MODEL','DRUG','TYPE','PVALUE')

anova_pvals$padj <- p.adjust(anova_pvals$PVALUE, method="BH")

setwd('..')
write.table(anova_pvals, file = anova_f, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

# output: anova="anova_pvals.tsv", dir=directory('plots')