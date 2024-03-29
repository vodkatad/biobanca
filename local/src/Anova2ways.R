## two way anova with normalized data

library(readxl)
library(reshape)
library(ggpubr)

drugs_tables <- snakemake@input[["drugs_tables"]]
anova_f <- snakemake@output[["anova"]]
out_dir <- snakemake@output[["out_dir"]]

model <- c("CRC0322", "CRC0059", "CRC1272", "CRC1331", "CRC0327")

res <- data.frame(stringsAsFactors = FALSE)

if (dir.exists(out_dir)) {
   unlink(out_dir)
} else {
   dir.create(out_dir)
}

setwd(out_dir)

for (i in seq(1, length(model))) {
  dnorm <- read_excel(drugs_tables, sheet=model[i])
  dnorm <- as.data.frame(dnorm, stringsAsFactors =FALSE)
  dnorm$DRUG_EXP <- paste0(dnorm$DRUG, dnorm$EXP, "@", dnorm$CONDITION)
  h <- reshape::melt(dnorm, id.vars = "DRUG_EXP", measure.vars = c("DOSE_1", "DOSE_2", "DOSE_3",
                                                            "DOSE_4", "DOSE_5"))
  h$CONDITION <- sapply(strsplit(h$DRUG_EXP, "@", fixed = TRUE), function(x) {x[[2]]})              
  h$DRUG <- sapply(strsplit(h$DRUG_EXP, "@", fixed=TRUE), function(x) {x[[1]]})
  h$DRUG_EXP <- NULL
  h$DRUG <- substr(h$DRUG, 1, nchar(h$DRUG)-1)
  drugs <- unique(h$DRUG)
    for (j in seq(1, length(drugs))) {
      drug <- drugs[j]
      a_f <- subset(h, h$DRUG == drug)
      fit <- aov(value ~ CONDITION*variable, data = a_f)  
      s_fit <- summary(fit)
      sf <- s_fit[[1]][["Pr(>F)"]]
      res <- rbind(res, c(model[i], drug, sf[1:3]), stringsAsFactors = FALSE)
      
      ggline(a_f, x = "variable", y = "value", color = "CONDITION",
             add = c("mean_se", "dotplot"),
             palette = c("#00AFBB", "#E7B800"), title=paste0(model[i], " ", drug))
      ggsave(paste0(model[i], " ", drug,".pdf"))
    }
}

colnames(res) <- c("MODEL", "DRUG", "CONDITION", "VARIABLE", "INTERACTION")  
                   
setwd('..')
write.table(res, file = anova_f, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)                   
                   