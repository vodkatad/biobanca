library(ggplot2)
th <- function() {
  textSize <- 1.5
  theme <- theme_bw() +
  theme(
    strip.text = element_text(size = rel(textSize)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(size = rel(1.8)),
    axis.text.x = element_text(size=rel(1.7)),
    axis.text.y = element_text(angle = 0,
                               size = rel(1.7)),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y.left = unit(3,'mm'),
    
    legend.position = "top",
    legend.justification = "right",
    #legend.margin = margin(unit(0, "cm")),
    legend.title = element_text(size = rel(textSize), face = "bold"),
    legend.text = element_text(size = rel(1.2)),
    legend.background = element_rect(size=0.5, linetype="solid", color="black"),
    plot.title = element_text(
      face = "bold",
      size = rel(2),
      hjust = 0.5
    ),
    panel.border = element_blank(),
    plot.caption = element_text(size=rel(1))
  )
theme
}
current_theme <- th()


exprcor <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/trans_sign/expr/res.tsv', sep='\t', header=T)
cncor <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/diagonal_xo.tsv', sep='\t', header=T)
ecor <- data.frame(smodel=rownames(exprcor), exprPearson=diag(as.matrix(exprcor)))
colnames(cncor)[2] <- 'CNPearson'
m <- merge(ecor, cncor ,by="smodel")
dim(m)
dim(ecor)
dim(cncor)



gglm <- function(x, y, nx, ny, title) {
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y)
  ggplot(d, aes(x=x, y=y)) +geom_point(size=2)+current_theme+xlab(nx)+ylab(ny)+labs(caption=paste0(round(pe$estimate, digits=3), ', pval=', formatC(pe$p.value, format = "e", digits = 3)))+ggtitle(title)+geom_smooth(method=lm)+current_theme
  #ggsave(paste0(title, "_",nx,"_",ny,'.svg'), width=10, height=10, units="in")
}


#ggplot(m, aes(x=CNPearson, y=exprPearson)) +geom_point()+current_theme+geom_smooth(method='lm')
gglm(m$CNPearson, m$exprPearson, nx="CN_pearson", ny="expr_pearson",'Expression/CopyNumber stability PDO/PDX')

mm <- m[!grepl('CRC1888', m$smodel),]
gglm(mm$CNPearson, mm$exprPearson, nx="CN_pearson", ny="expr_pearson",'Expression/CopyNumber stability PDO/PDX')


#### gistic delta scores gsea

library(clusterProfiler)

data <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/dataset/V1/shallowseq/overallgenes_delta_gistic.tsv', header=TRUE, sep="\t")
geneList <- data$delta
names(geneList) <- as.character(data$gene_symbol)
geneList <- sort(geneList, decreasing = TRUE)
library(msigdbr)


m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, human_gene_symbol)

#minSize=15, maxSize=1000, nperm=10000,
em2 <- GSEA(geneList, TERM2GENE = m_t2g, minGSSize=15, maxGSSize=1000, nPerm=10000, pvalueCutoff=0.1, verbose=FALSE)

library(enrichplot)
gseaplot2(em2, geneSetID = 1, title = em2$Description[1])

# TODO decide n based on sign?
#  p.adjust   qvalues 
#n <- sum(em2$p.adjust < 0.05)
n <- 30
ridgeplot(em2, showCategory=20)
dotplot(em2, showCategory=20)
gseaplot2(em2, geneSetID = 10, title = em2$Description[10])

loss <- as.character(data[data$delta < q1, 'gene_symbol'])
gain <- as.character(data[data$delta > q2, 'gene_symbol'])
universe <- as.character(data$gene_symbol)

library('org.Hs.eg.db')
# ggo <- groupGO(gene     = loss,
#                OrgDb    = org.Hs.eg.db,
#                ont      = "CC",
#                level    = 3,
#                readable = TRUE)

ego <- enrichGO(gene          = loss,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.5)

# nothing for BP
barplot(ego, showCategory=20)

ego <- enrichGO(gene          = gain,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.5)

barplot(ego, showCategory=20)

