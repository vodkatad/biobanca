library(readxl)
library(stringr)
library(reshape)
library(ggpubr)
library(tidyverse)

# get dir and tsv from snakemake inputs
drugs_tables <- snakemake@input[["drugs_tables"]]
anova_f <- snakemake@output[["anova"]]
max_value_drug <- snakemake@output[["max_drug"]]
out_dir <- snakemake@output[["out_dir"]]
models <- c('CRC0322','CRC0327','CRC1331','CRC1272','CRC0059')

#save.image("anova_long.Rdata")
anova_pvals <- data.frame() #col.names=c('MODEL','DRUG','CONDITION','PVALUE')
max_value_df <- data.frame()

# do we need to create dir or does snakemake take care of it?
if (dir.exists(out_dir)) {
  unlink(out_dir)
} else {
  dir.create(out_dir)
}

setwd(out_dir)

for (i in seq(1, length(models))) {
  
  Drug_Screening_Tables <- read_excel(drugs_tables, sheet = models[i])
  data <- as.data.frame(Drug_Screening_Tables)
  data$DRUG <- str_trim(data$DRUG, "both") 
  data$DRUG <- str_replace_all(data$DRUG, " +", "_") 
  
  drugs <- unique(data$DRUG)
  conditions <- unique(data$CONDITION)
  
  for (j in seq(1, length(drugs))) {
    
    for (k in seq(1, length(conditions))) {
      
      w_drug <- drugs[j]
      w_condition <- conditions[k]
      
      subset <- data[data$DRUG == w_drug & data$CONDITION == w_condition, grepl('DOSE', colnames(data))]
      
      max_value <- mean((subset$DOSE_1 - subset$DOSE_5) / subset$DOSE_1)
      
      res1 <- c(models[i], drugs[j], conditions[k], max_value)
      
      max_value_df <- rbind(max_value_df, res1, stringsAsFactors=FALSE)
      
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
                ylab = "Viability", xlab = "Dose", add="jitter")
      ggsave(paste0("boxplotjitter_", models[i],"_", drugs[j], "_",conditions[k],'.pdf'))
      
      
      long$x <- seq(1, nrow(long))
      ggplot(data=long, aes(y=value, x=x, fill=variable))+geom_col()+theme_bw()+ggtitle(paste0(models[i]," ", drugs[j], " ",conditions[k]))
      ggsave(paste0("barplot_", models[i],"_", drugs[j], "_",conditions[k],'.pdf'))
      
      # Compute the analysis of variance
      res.aov <- aov(value ~ variable, data = long)
      # Summary of the analysis
      sumaov <- summary(res.aov)
      pvalue <- sumaov[[1]]$`Pr(>F)`[1]
      
      res <- c(models[i], drugs[j], conditions[k], pvalue)
      anova_pvals <- rbind(anova_pvals, res, stringsAsFactors=FALSE)
    }
  }
  
}


colnames(anova_pvals) <- c('MODEL','DRUG','CONDITION','PVALUE')
colnames(max_value_df) <- c('MODEL','DRUG','CONDITION','MAX_VALUE')

save.image("combo.Rdata")
### rimuovo le MONO PERCHé NON SONO MAI STATE FATTE

anova_pvals <- anova_pvals %>% filter(!CONDITION == "Mono")
anova_pvals$padj <- p.adjust(anova_pvals$PVALUE, method="BH")

max_value_df <- max_value_df %>% filter(!CONDITION == "Mono")
max_value_df$MAX_VALUE <- as.numeric(max_value_df$MAX_VALUE)
max_value_df <- max_value_df[order(-max_value_df$MAX_VALUE),]

setwd('..')
write.table(anova_pvals, file = anova_f, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

write.table(max_value_df, file = max_value_drug, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)


# output: anova="anova_pvals.tsv", dir=directory('plots')