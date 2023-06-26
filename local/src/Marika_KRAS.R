load('/scratch/trcanmed/biobanca/dataset/V1/targeted/preprocGeneAF_0.05.Rdata')

table(pdobing[rownames(pdobing)=='KRAS',])
table(pdxbing[rownames(pdxbing)=='KRAS',])
table(pdxbing[rownames(pdxbing)=='KRAS',] & pdobing[rownames(pdobing)=='KRAS',])

colnames(pdobing)[pdxbing[rownames(pdxbing)=='KRAS',] != pdobing[rownames(pdobing)=='KRAS',]]

# muts

okras <- gpdo[gpdo$genes=="KRAS",]
xkras <- gpdx[gpdx$genes=="KRAS",]

okras_muts <-  pdo[rownames(okras),]
xkras_muts <-  pdx[rownames(xkras),]

xeno_mutated <- colnames(xkras_muts)[colSums(xkras_muts)>0.05]

# n xenos
length(colnames(xkras_muts))
# n xenos with muts
length(xeno_mutated)

pdo_x_is_mut <- pdo[, xeno_mutated]
pdo_x_is_mut_k <- pdo[rownames(okras), xeno_mutated]

colnames(pdo_x_is_mut_k)[colSums(pdo_x_is_mut_k)<0.05]
length(colnames(pdo_x_is_mut_k)[colSums(pdo_x_is_mut_k)>=0.05])
# always -1 because CRC1241 is still there


#muts
xx <- xkras_muts[,colSums(xkras_muts)>0.05]
gpdx[rownames(xx),]


###

o <- pdo_x_is_mut_k[,colSums(pdo_x_is_mut_k)>=0.05]
x <- xkras_muts[,colSums(xkras_muts)>=0.05]
library(reshape)

o$id <- rownames(o)
x$id <- rownames(x)
ol <- melt(o)
xl <- melt(x)

ol <- ol[ol$value != 0,]

xl <- xl[xl$value != 0,]

ol[,c('id', 'variable')]==xl[,c('id', 'variable')]

###
pdo_mutated <- colnames(okras_muts)[colSums(okras_muts)>0.05]

# n xenos
length(colnames(okras_muts))
# n xenos with muts
length(pdo_mutated)

pdx_o_is_mut <- pdx[, pdo_mutated]
pdx_o_is_mut_k <- pdx[rownames(xkras), pdo_mutated]

colnames(pdx_o_is_mut_k)[colSums(pdx_o_is_mut_k)<0.05]
length(colnames(pdx_o_is_mut_k)[colSums(pdx_o_is_mut_k)>=0.05])

pdo[rownames(okras_muts), colnames(pdx_o_is_mut_k)[colSums(pdx_o_is_mut_k)<0.05]]
gpdo['chr12:25245350:C:A',]
     