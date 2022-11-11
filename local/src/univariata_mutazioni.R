## univariata mutazioni

library(tidyverse)
library(sjPlot)

rf <- snakemake@input[["rdata_targeted_all"]]
buoni_f <- snakemake@input[["pdo"]]
fit_gene <- snakemake@output[["fit_plot_gene"]]
res_genes <- snakemake@output[["res"]]

load(rf)
mut <- as.data.frame(pdobing)
mut <- as.data.frame(t(mut))

mut2 <- mut
for (i in rownames(mut2)) {
  for (j in seq(colnames(mut2))) {
    mut2["numero", j] <- length(mut2[,j][mut2[,j]== TRUE])
  }
}
mut2 <- as.data.frame(t(mut2))

mut2 <- mut2 %>% filter(numero > 4)
mutatipiudicinque <- rownames(mut2)

mut <- as.data.frame(t(mut))
mut$genes <- rownames(mut)

mut <- mut %>% filter(genes %in% mutatipiudicinque)
mut$genes <- NULL
mut <- as.data.frame(t(mut))
mut$smodel <- rownames(mut)


buoni <- read.table(buoni_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

merged <- merge(mut, buoni, by = "smodel")
#setdiff(mut$smodel, buoni$smodel)
#[1] "CRC1448" "CRC1897"

merged$smodel <- NULL
merged$buoni <- as.factor(merged$buoni)

#fit <- glm(buoni ~ KRAS, merged, family = binomial())
#plot_model(fit)
vars <- names(merged[1:26])
formula_strings <- sprintf("buoni ~ %s", vars)
res <- data.frame(matrix(ncol = 2, nrow = 26))
rownames(res) <- vars

plot_list = list()
i = 1
for (f in formula_strings) {
  fit <- glm(formula = as.formula(f), merged, family = binomial())
  sfit <- summary(fit)
  res[i, 1] <- sfit$coefficients[2,4]
  res[i, 2] <- coef(fit)[2]
  #plot_list[[i]] <- plot_model(fit)
  i = i+1
}

colnames(res) <- c("pvalue", "coef(fit)")

write.table(res, res_genes, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

ctnnb1 <- glm(formula = buoni ~ CTNNB1, merged, family = binomial())
pdf(fit_gene)
plot_model(ctnnb1, axis.lim = c(0.01, 1), title = "Validation", axis.labels = "CTNNB1")
dev.off()
##correzione con bh
# pdf("buoni_mut_glm.pdf")
# for (j in 1:26) {
#   print(plot_list[[j]])
# }
# dev.off()
# 
# 
# 
# 
# plot_model(fit)
