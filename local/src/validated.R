library(stringr)
data <- read.table("~/Downloads/XENTURION_DEF_SML_12-10.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE, comment.char='')

# QC:
# Not passed     Passed 
# 23        247

#          Failed  Successful Successful  
#23          60         186           1 

# fix successful in 

colnames(data) <- c('smodel','QC','Origin','DERIVATION','VALIDATION')
totrim <- c('QC','Origin','DERIVATION','VALIDATION')

for (i in seq(1, length(totrim))) {
  data[,colnames(data)==totrim[i]] <- str_trim(data[,colnames(data)==totrim[i]], side = c("both"))
}

for (i in seq(1, length(totrim))) {
  print(as.data.frame(table(data[,colnames(data)==totrim[i]])))
}

dataWithoutMisteryUniverse <- data[data$VALIDATION != "Not performed",]
dataWithoutMisteryUniverse <- dataWithoutMisteryUniverse[dataWithoutMisteryUniverse$QC  == "Passed",]

dataWithoutMisteryUniverse$buoni <- FALSE
dataWithoutMisteryUniverse$buoni[grepl('Successful', dataWithoutMisteryUniverse$VALIDATION, fixed=FALSE)] <- TRUE
dataWithoutMisteryUniverse$buoni[dataWithoutMisteryUniverse$VALIDATION == "Trusted"] <- TRUE

##
setwd("~/work/biobanca/")
dataWithoutMisteryUniverse <- dataWithoutMisteryUniverse[,c('smodel','buoni')]
write.table(dataWithoutMisteryUniverse, file="biobanca_pdo_buoni.tsv", sep="\t", quote=FALSE, row.names=FALSE)