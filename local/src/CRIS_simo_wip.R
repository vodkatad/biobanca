simo <- read.table(simo_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmx <- read.table(lmx_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lmo <- read.table(lmo_in, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#simo[,c(2,3)] <- NULL
lmx <- lmx[,c(1,2)]
names(lmx) <- c("PDX_long_lineage", "PDX")
lmo <- lmo[,c(1,2)]
names(lmo) <- c("PDO_long_lineage", "PDO")
lmx$PDX_lineage <- substr(lmx$PDX_long_lineage, 0, 12)
lmo$PDO_lineage <- substr(lmo$PDO_long_lineage, 0, 12)

lmx_m <- merge(simo, lmx, by="PDX_lineage")
lmx_m[,c(2)] <- NULL
lmo_m <- merge(simo, lmo, by="PDO_lineage")
lmo_m[,c(2)] <- NULL


lmx_m$model <- substr(lmx_m$PDX_lineage, 0, 7)
lmo_m$model <- substr(lmo_m$PDO_lineage, 0, 7)

simo_df <- merge(lmo_m, lmx_m, by = "model")
#model_cris <- simo_df[,c(1,5,3)]

#save.image("pippo.RData")

wrong_pdo <- c()
wrong_pdx <- c()
smodel <- unique(simo_df$model)
right <- data.frame(stringsAsFactors = FALSE)
wrong <- data.frame(stringsAsFactors = FALSE)


for (i in seq(1, length(smodel))) {
  tmp <- simo_df[simo_df$model == smodel[i],]
  if ( nrow(tmp) > 1 ) {
    if ( length(unique(tmp$PDO)) != 1 ) {
      wrong_pdo <- c(wrong_pdo, smodel[i])
      wrong <- rbind(wrong, tmp)
    }
    if ( length(unique(tmp$PDX)) != 1 ) {
      wrong_pdx <- c(wrong_pdx, smodel[i])
      wrong <- rbind(wrong, tmp)
    }
    if ( length(unique(tmp$PDO)) == 1 && length(unique(tmp$PDX)) == 1 ) {
      right <- rbind(right, unique(tmp))
    }
  } else  {
    right <- rbind(right, tmp)
  } 
}

right$PDO_passage <- as.numeric(substr(right$PDO_long_lineage, 13, 14))
right$PDX_passage <- as.numeric(substr(right$PDX_long_lineage, 13, 14))

keep <- data.frame()
single <- data.frame()
i = 2
for (i in seq(2, nrow(right))) {
  if (right[i,'model'] == right[(i-1), "model"]) { 
    if (right[i, "PDO_passage"] <= right[(i-1), "PDO_passage"]) {
      keep <- rbind(keep, right[i,])
    } else {
      keep <- rbind(keep, right[(i-1),])
    }
  } else {
    single <- rbind(single, right[i,])
  }
}

my_data[!duplicated(my_data$Sepal.Width), ]
# if (!(right[i, "model"]) %in% keep[i,"model"]) {
#   keep <- rbind(keep, right[i,])
# } else {
#   double <- rbind(double, right[i,])
# else if (right[i,'model'] != right[(i-1), "model"] && !(right[i, "model"] %in% keep[i,])) {
#   single <- rbind(single, right[i,])

colnames(right) <- c("model", "CRIS_PDX", "CRIS_PDO")
colnames(wrong) <- c("model", "CRIS_PDX", "CRIS_PDO")

sink(snakemake@log[[1]])
print(paste0("wrong pdo: ", length(wrong_pdo)))
print(paste0("wrong pdx: ", length(wrong_pdx)))
sink()

write.table(right, snakemake@output[["right"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(wrong, snakemake@output[["wrong"]], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
