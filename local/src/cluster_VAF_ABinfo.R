library(vcfR)

vcf <- read.vcfR('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/mutect_paired/CRC1599PRX0A02002TUMD03000V2.ABannot.pass.vcf.gz')

VAF <- extract.gt(vcf, element = "AF", as.numeric = TRUE)
A <- extract.info(vcf, element = "A", as.numeric = TRUE)
B <- extract.info(vcf, element = "B", as.numeric = TRUE)

CN <- A+B
CN[is.na(CN)] <- 2
CCF_p1 <- VAF[,2] * CN

CCF_df <- data.frame(row.names=rownames(VAF), CCF=CCF_p1)
hist(CCF_df$CCF)

AD <- extract.gt(vcf, element = "AD", as.numeric = FALSE)
ADs <- strsplit(AD[,2], split=',', fixed=TRUE)
ref <- sapply(ADs, function(x) {as.numeric(x[[1]][1])})
alt <- sapply(ADs, function(x) {as.numeric(x[[2]][1])})

pyclone_in <- data.frame(mutation_id=rownames(VAF), sample_id='CRC1599PRX0A02002TUMD03000V2', ref_counts=ref, alt_counts=alt, major_cn=A,
                         minor_cn=B, normal_cn=2, tumour_content=1, error_rate=0.001)

pyclone_in <- pyclone_in[!is.na(pyclone_in$major_cn) & !is.na(pyclone_in$minor_cn),]
pyclone_in[grepl('chrY', pyclone_in$mutation_id), 'normal_cn'] <- 1

write.table(pyclone_in, file="/mnt/trcanmed/snaketree/prj/snakegatk/dataset/Pri_Mets_godot/pyclone/CRC1599PRX0A02002TUMD03000V2.pyc.tsv", quote=FALSE, sep="\t")
