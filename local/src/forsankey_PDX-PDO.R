#!/opt/R/R-3.6.3/bin/Rscript

library(getopt)

set.seed(42)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'cris', 'c', 1, 'character',
  'meda', 'm', 1, 'character',
  'simo', 's', 1, 'character',
  'output', 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)


simo <- read.table(gzfile(opt$simo), sep='\t', quote="", header=TRUE, stringsAsFactors=FALSE)
cris <- read.table(gzfile(opt$cris), sep='\t', quote="", header=TRUE, stringsAsFactors = FALSE)
meda <- read.table(gzfile(opt$meda), sep='\t', quote="", header=TRUE, stringsAsFactors = FALSE)

df <- merge(cris, meda, by.x=names(cris)[1], by.y=names(meda)[1]) # 14 filtrati
df[is.na(df)] <- ""

df <- df[df$BH.FDR<0.2,] # 140 filtered on reliable classification
df <- df[grep("BASALE", df$type),] # 325 filtered as NON basal
df <- df[!df$RNA_marker=='RNA_marker',] # 3 filtered as lymphomas by RNA_marker
df <- df[!df$RNA_PC=='RNA_PC',] # 0 filtered as lymphomas by RNA_PC
### 743 samples remain

df[,c(3:5,7:15)] <- NULL
df$lineage <- substr(df[,1], 1, 12)
names(df)[c(1,2)] <- c("samples", "prediction")

input_sankey <- as.data.frame(matrix(nrow=nrow(df), ncol=6))
names(input_sankey) <- names(df)

PDO_cris_diff <- 0
PDX_cris_diff <- 0
added <- 0

# DUBBIO era una specificita` del caso che non ci fossero mai due xeno e un solo organoide, vero?
for ( i in 1:nrow(simo) ) {
  tmp_pdx <- df[df[,'lineage']==simo[i,'PDX_lineage'],] 
  tmp_pdo <- df[df[,'lineage']==simo[i,'PDO_lineage'],]
  keep <- TRUE
  if ( nrow(tmp_pdo) > 1 ) {
    if (length(unique(tmp_pdo[,'prediction'])) != 1) {
      keep <- FALSE
      PDO_cris_diff <- PDO_cris_diff + 1    
    }
    if ( nrow(tmp_pdx) > 1 ) {
      if (length(unique(tmp_pdx[,'prediction'])) != 1)  {
        keep <- FALSE
        PDX_cris_diff <- PDX_cris_diff + 1
      }  
    }
    if (keep) { # if I have a unique prediction for PDX and PDO
      input_sankey[added,] <- c(tmp_pdx[1,c("lineage", "prediction","BH.FDR")], tmp_pdo[1,c("lineage", "prediction","BH.FDR")])
      added <- added + 1
    }
  } else if (nrow(tmp_pdx) == 1 & nrow(tmp_pdo) == 1) {
    input_sankey[added,] <- c(tmp_pdx[1,c("lineage", "prediction","BH.FDR")], tmp_pdo[1,c("lineage", "prediction","BH.FDR")])
    added <- added + 1
  }
}

added <-added - 1
names(input_sankey) <- c("LMX_lineage", "prediction_LMX", "BH_FDR_LMX", "LMO_lineage", "prediction_LMO", "BH_FDR_LMO")
input_sankey <- input_sankey[complete.cases(input_sankey),]

for ( r in 1:length(rownames(input_sankey)) ) {
  input_sankey$switched[r] <- ifelse(input_sankey$prediction_LMX[r]!=input_sankey$prediction_LMO[r], "yes", "no")
  input_sankey$switch_type[r] <- ifelse(input_sankey$switched[r]=="yes", paste0(input_sankey$prediction_LMX[r]," > ",input_sankey$prediction_LMO[r]), paste(input_sankey$prediction_LMX[r]))
  input_sankey$pval_switch_sign[r] <- ifelse(input_sankey$BH_FDR_LMO[r]<0.05, "yes", "no")
}

input_sankey <- input_sankey[order(input_sankey$pval_switch_sign),]
# switched_sing <- input_sankey[input_sankey$switched=="yes",]


write.table(input_sankey, opt$output, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)