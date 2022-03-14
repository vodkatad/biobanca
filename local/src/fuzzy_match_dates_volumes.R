tgi_f <- snakemake@input[['tgi']]
lasvol_f <- snakemake@input[['las_vol']]
common_f <- snakemake@input[['common']]
res_f <- snakemake@output[['placebo_vols']]

tgi <- read.table(tgi_f, sep="\t", stringsAsFactors = FALSE, header = TRUE)
lasvol <- read.table(gzfile(lasvol_f), sep="\t", stringsAsFactors = FALSE, header = TRUE)

keep <- read.table(common_f, sep="\t", stringsAsFactors = FALSE, header = FALSE)


tgi$modelarm <- substr(tgi$GenID, 0, 12)
lasvol$modelarm <- substr(lasvol$longen, 0, 12)

tgi <- tgi[tgi$modelarm %in% keep$V1,]
lasvol <- lasvol[lasvol$modelarm %in% keep$V1,]

# quelli in comune sono tutti quelli con un nome di esperimento nel LAS quindi possiamo tenere loro
#length(intersect(tgi$Nome.Gruppo.LAS, lasvol$exp_group))

keep_lasvol <- lasvol[lasvol$exp_group %in% tgi$Nome.Gruppo.LAS,]

write.table(keep_lasvol, file=res_f, sep="\t", quote=FALSE)
