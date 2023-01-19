## univariata mutazioni

library(sjPlot)
library(tidyverse)

rf <- snakemake@input[["rdata_targeted_all"]]
buoni_f <- snakemake@input[["pdo"]]
fit_gene <- snakemake@output[["fit_plot_gene"]]
res_genes <- snakemake@output[["res"]]
mut_tab <- snakemake@output[["mut"]]

load(rf)
mut <- as.data.frame(pdobing)
mut <- as.data.frame(t(mut))

mut2 <- mut

save.image("ctnnb1.Rdata")
#for (i in rownames(mut2)) { # solo come nota secondo me il ciclo sulle righe non serve, conti i true totali in ogni colonna no?
  for (j in seq(colnames(mut2))) { # i proprio non lo usi mai dentro il corpo del for, dovrebbe venire esattamente la stessa cosa con il solo for interno
    mut2["numero", j] <- length(mut[,j][mut[,j]== TRUE]) # lavorando su mut2 per la length la riga aggiunta con numero e NA risultava
  }
#}
mut2 <- as.data.frame(t(mut2))

# A fini di esercizio sui famosi apply che a Trento adorano:
# n_muts <- apply(mut, 2, sum)
# n_muts == mut2$numero
# non tornano tutti, in particolare io ho qualche 0 che tu conti come 1 e 1 come 2, mentre quelli piu`grandi sono giusti...capiamo perchè:
# es VHL che per te è 1:
#> mut[,colnames(mut)=='STK11']
#  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [25] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [37] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [49] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [61] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [73] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [85] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [97] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#[109] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#[121] FALSE FALSE FALSE FALSE FALSE

#

mut2 <- mut2 %>% filter(numero > 4)
mutatipiudicinque <- rownames(mut2)

mut <- as.data.frame(t(mut))
mut$genes <- rownames(mut)

mut <- mut %>% filter(genes %in% mutatipiudicinque)
mut$genes <- NULL
mut <- as.data.frame(t(mut))
mut$smodel <- rownames(mut)

buoni <- read.table(buoni_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
buoni <- buoni %>% filter(!type == "Validation not performed")
rownames(buoni) <- buoni$CASE
for (i in rownames(buoni)) {
  if (buoni[i, "type"] == "Validation successful"){
    buoni [i, "type"] <- TRUE
  } else {
    buoni [i, "type"] <- FALSE
  }
}

colnames(buoni) <- c("smodel", "buoni")
merged <- merge(mut, buoni, by = "smodel")
#setdiff(mut$smodel, buoni$smodel)
#[1] "CRC1448" "CRC1897"

merged$smodel <- NULL
merged$buoni <- as.factor(merged$buoni)

merged2 <- merged
merged2$Validation <- NA
for (i in rownames(merged2)) {
  if (merged2[i, "buoni"] == TRUE){
    merged2[i, "Validation"] <- "Validation successful"
  } else {
    merged2[i, "Validation"] <- "Validation not successful"
  }
}
merged2 <- cbind(merged2$Validation, merged2)
merged2$buoni <- NULL
merged2$Validation <- NULL
names(merged2)[names(merged2) == "merged2$Validation"] <- "Validation"
merged2[merged2==TRUE] <- "Mutated"
merged2[merged2==FALSE] <- "WT"

write.table(merged2, mut_tab, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

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
res <- res[order(res$pvalue),]

write.table(res, res_genes, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

##correzione con bh
# pdf("buoni_mut_glm.pdf")
# for (j in 1:26) {
#   print(plot_list[[j]])
# }
# dev.off()
# plot_model(fit)

ctnnb1 <- glm(formula = buoni ~ CTNNB1, merged, family = binomial())
pdf(fit_gene, useDingbats=FALSE)
plot_model(ctnnb1, axis.lim = c(0.01, 1), title = "Validation", axis.labels = "CTNNB1 [True]")
dev.off()

# to put odds ratio in the main text, run 19/01/23 loading ctnnb1.Rdata
#> exp(res[,2])
 #[1] 6.722689e-02 1.939655e-01 1.475410e-01 3.669725e-01 3.559322e-01
 #[6] 3.559322e-01 2.000000e+00 4.102564e-01 4.102564e-01 5.217391e-01
 # [11] 5.789474e-01 1.075472e+00 1.051136e+00 2.502636e+06 4.317079e-08
 # [16] 2.481781e+06 2.481781e+06 2.481781e+06 2.481781e+06 2.481781e+06
 # [21] 7.039492e+06 2.461270e+06 6.860522e+06 6.860522e+06 6.860522e+06
 # [26] 6.860522e+06
