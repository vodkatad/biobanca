snake_bck <- snakemake
rdata_f <- snakemake@input[['rdata']]
load(rdata_f)

xeno_f <- snake_bck@output[['xeno']]
pdo_f <- snake_bck@output[['pdo']]

#pc <- min(xeno_df_bck[xeno_df_bck>0], pdo_df_bck[pdo_df_bck>0])
#pc <- 0.01
pc <- 0.01
xeno_df <- log2(xeno_df_bck+pc)
pdo_df <- log2(pdo_df_bck+pc)

write.table(xeno_df, file=gzfile(xeno_f), quote=FALSE, sep="\t")
write.table(pdo_df, file=gzfile(pdo_f), quote=FALSE, sep="\t")

