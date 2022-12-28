library(ggplot)
cld <- read.table('/scratch/trcanmed/biobanca/local/share/data/clinical_data_done.tsv', sep="\t", header=TRUE, stringsAsFactors = F) # TODO ask Marti

vivod <- read.table('/scratch/trcanmed/biobanca/dataset/V1/cetuximab/cetuxi_perc_w3_buoni.tsv', sep="\t", header=TRUE) 
vitrod <- read.table('/scratch/trcanmed/biobanca/dataset/V1/cetuximab/pdo_cetuxi_buoni.tsv', sep="\t", header=TRUE) 

nrow(cld)
nrow(vivod)
nrow(vitrod)
table(cld$RESPONSE.TO.CETUX.THERAPY)
cld$cet <- cld$RESPONSE.TO.CETUX.THERAPY
cld <- cld[cld$cet != "NA" & !is.na(cld$cet), ]   # we miss all data for some
table(cld$cet, cld$THERAPY..CETUXIMAB..Y.N)
cld <- cld[cld$cet != "N",] # were not treated with ctx

cld[cld$cet == "N.D.", 'cet'] <- 'NA'
cld$cet <- factor(cld$cet, levels=c('PR', 'SD', 'PD','NA'))

m <- merge(cld, vivod, by.x="CASE", by.y="case")

nrow(m)

ggplot(data=m, aes(x=cet, y=perc))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+theme_bw()+ggtitle('vs Vivo 3W')

m <- merge(cld, vitrod, by.x="CASE", by.y="case")

nrow(m)

ggplot(data=m, aes(x=cet, y=CTG_5000))+geom_boxplot(outlier.shape=NULL)+geom_jitter()+theme_bw()+ggtitle('vs CTG5000')
