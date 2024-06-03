setwd('/scratch/trcanmed/DE_RNASeq/dataset/validated_or_not_xeno_biobanca')
dd <- read.table('GO_results_validation_cutoff0.05-successful.vs.notest.failed_down.tsv', sep="\t", quote='', stringsAsFactors = F)
du <- read.table('GO_results_validation_cutoff0.05-successful.vs.notest.failed_up.tsv', sep="\t", quote='',  stringsAsFactors = F)


load('')

dd <- dd[order(dd$pvalue),]
dd$dir <- 'down'
du <- du[order(du$pvalue),]
du$dir <- 'up'
dp <- rbind(head(dd, n=5), head(du, n=5))
dp$count_cor <- ifelse(dp$dir == 'down',  -dp$Count, dp$Count)
dp <- dp[order(-dp$count_cor),]
p1 <- ggplot(d=dp, aes(y=count_cor, x=reorder(Description, count_cor), fill=log10(p.adjust)))+geom_col()+coord_flip()+theme_bw()+
  xlab('GO')+ylab('N. DEG')+#scale_fill_gradient(low="blue", high="red", limits=c(-1,-10))
  scale_fill_continuous(low="blue", high="red")#, limits=c(-10,-1.2))#+scale_fill_continuous(type='viridis')

load('/scratch/trcanmed/git_biobanca/local/data/theme_5.Rdata')

ggsave('/scratch/trcanmed/biobanca/local/share/data/go_doublefailed.pdf', p1 + unmute_theme,  height=2.5, width=5, units='in')
