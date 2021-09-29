### heatmap for non expressed genes in LMO

library(pheatmap)

## carico i geni non espressi in lmo per avere una lista da cercare in meda (con le altri classi)
## e per la loro espressione utilizzo vsd

genes_f <- "/home/mferri/notexpressed_genes_lmo.tsv"
genes_lmo <- read.table(genes_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meda <- "/mnt/trcanmed/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_f <- filter(meda_f, grepl("LMX_BASALE", type))

vsd <- "/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- rownames(vsd)
rownames(vsd) <- NULL
vsd <- cbind(genes,vsd)
vsd <- vsd %>% mutate(genes = gsub("H_", "", genes))
colnames(vsd)[colnames(vsd) == 'genes'] <- 'symbol'
vsd_subset <- vsd[, colnames(vsd) %in% c(meda_f$sample_id_R, "symbol")]

## vsd_genes continene la traduzione da symbols a ENTREZ ID perchè
## i geni degli LMO sono un ricavato di CMS che lavora solo con gli entrez

vsd_genes <- "/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/VSD_genes.tsv"
vsd_genes <- read.table(vsd_genes, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd_genes <- vsd_genes[!grepl(",", vsd_genes$description),,drop = FALSE]
vsd_genes <- vsd_genes %>% na.exclude

## merge per fare la traduzione dei geni di interesse per gli LMX
vsd_f <- merge(vsd_genes, vsd_subset, by = "symbol")
vsd_f$symbol <- NULL
vsd_f <- vsd_f %>% remove_rownames %>% column_to_rownames(var="description")
genes_lmx <- tibble::rownames_to_column(vsd_f, "GENES")

## stesso lavoro cambiando i nomi per ottenere il subset di geni negli LMH

meda_lmh <- read.table(meda, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meda_lmh <- filter(meda_lmh, grepl("LMH", type))
vsd_subset_lmh <- vsd[, colnames(vsd) %in% c(meda_lmh$sample_id_R, "symbol")]
vsd_f_lmh <- merge(vsd_genes, vsd_subset_lmh, by = "symbol")
vsd_f_lmh$symbol <- NULL
vsd_f_lmh <- vsd_f_lmh %>% remove_rownames %>% column_to_rownames(var="description")
genes_lmh <- tibble::rownames_to_column(vsd_f_lmh, "GENES")

## doppio merge per avere una tabella contenente i geni di interesse 
## e la loro espressione in LMO, LMX e LMH

genes_lmo_lmx <- merge(genes_lmo, genes_lmx, by = "GENES")
genes_all <- merge(genes_lmo_lmx, genes_lmh, by = "GENES")
genes_all <- genes_all %>% remove_rownames %>% column_to_rownames(var="GENES")
genes_all<-as.data.frame(t(genes_all))
genes_all_ph <- data.matrix(genes_all)

## creo un dataframe in cui inserisco i type dei genealogy, necessario a 
## pheatmap per aiutare la classificazione dei casi (la legenda non è comprensibile)
genes_type <- data.frame(row.names = rownames(genes_all), type = substr(rownames(genes_all), 8, 10), stringsAsFactors = FALSE)

#png("heatmap_genes.png")
pheatmap(genes_all_ph, cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE,
         show_colnames = FALSE, annotation_row = genes_type)
#dev.off()
## calcolo media e sds degli lmx per vedere quelli più espressi a che passaggi appartengono

genes_lmx_men <- genes_lmx
genes_lmx_men$GENES <- NULL
sds <- apply(genes_lmx_men, 1, sd)
mean <- colMeans(genes_lmx_men)
zerosd <- as.data.frame(names(sds[sds==0]))

lmx <- data.frame(row.names = colnames(genes_lmx_men), avg_exp = mean)

max <- apply(genes_lmx_men, 2, max)
lmx$max <- max


