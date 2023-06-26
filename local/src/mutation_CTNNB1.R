mat <- "/scratch/trcanmed/pdxopedia/dataset/sanger_targeted_v2_genealogy/fl_all_xenos_orig.tsv" 
mat <- read.table(mat, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

mut <- "/scratch/trcanmed/pdxopedia/dataset/sanger_targeted_v2_genealogy/fl_coding_muts_annotated_hg19.tsv"
mut <- read.table(mut, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

mat$type <- substr(mat$V2, 8, 10)
mat <- mat %>% filter(type == "LMX")
mat$type <- NULL
mat <- as.data.frame(mat$V2)
colnames(mat) <- "Sample"
mat$performed <- "P"

mut <- mut %>% filter(Gene == "CTNNB1")
mut$mutated <- "True"
mut$dev <- substr(mut$Sample, 8, 10)
mut <- mut %>% filter(dev == "LMX")

mut_def <- merge(mat, mut, by="Sample", all.x=TRUE)
mut_def <- mut_def[!duplicated(mut_def$Sample), ]
mut_def$mutated <- mut_def$mutated %>% replace_na('False')
mut_def <- mut_def[,c("Sample", "mutated")]
mut_def$CASE <- substr(mut_def$Sample, 1, 7)
names(mut_def)[names(mut_def) == "mutated"] <- "CTNNB1"


write.table(mut_def, file = "/home/mferri/BCAT_mutation.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
