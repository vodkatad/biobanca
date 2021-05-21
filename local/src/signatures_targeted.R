library(MutationalPatterns)
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
ref_transcriptome <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
library(ref_genome, character.only = TRUE)
library(ref_transcriptome, character.only = TRUE)
library(NMF)
library(gridExtra)
library(ggplot2)
library(reshape)

dx <- '/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdx/mutect/'
do <- '/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_targeted_pdo/mutect/'

vcflistx <- paste0(dx, 'list_pass')
vcflisto <- paste0(do, 'list_pass')

lx <- read.table(vcflistx, sep="\t", header=FALSE)
lo <- read.table(vcflisto, sep="\t", header=FALSE)

lxx <- paste0(dx, lx$V1)
lxo <- paste0(do, lo$V1)

nx <- paste0('xeno',seq(1,length(lxx)))
no <- paste0('pdo',seq(1,length(lxo)))

groups <- c(rep('xeno', length(lxx)), rep('pdo', length(lxo)))

vcfs <- read_vcfs_as_granges(c(lxx, lxo), c(nx, no), ref_genome)

type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
pdf("spectrum.pdf")
plot_spectrum(type_occurrences, CT = TRUE)
graphics.off()
pdf("samples_spectrum.pdf")
plot_spectrum(type_occurrences, CT = TRUE, by=groups)
graphics.off()

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

cosmic <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "") # ???


sp_url <- paste(cosmic, sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe> 
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names>
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
## why NA?

hclust_cosmic = cluster_signatures(cancer_signatures, method = "average") # store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
pdf("cosine_sign_cosmic.pdf")
plot_cosine_heatmap(cos_sim_samples_signatures,col_order = cosmic_order,cluster_rows = TRUE)
graphics.off()

ff <- fit_to_signatures(mut_mat, cancer_signatures)
select <- which(rowSums(ff$contribution) > 10)
# Plot contribution barplot
pdf("cosmic_contribution_topsign.pdf") ## not pretty
plot_contribution(ff$contribution[select,],cancer_signatures[,select],coord_flip = FALSE,mode = "absolute")
graphics.off()
pdf("cosmic_contribution_heatmap.pdf")
plot_contribution_heatmap(ff$contribution,cluster_samples = TRUE,method = "complete")
graphics.off()

rowan <- data.frame(row.names=colnames(ff$contribution), type=groups)
pheatmap(t(as.matrix(ff$contribution)),annotation_row = rowan, cluster_cols=F)


garb <- sapply(sample_names, function(x) { 
  #pdf(paste0("cosmic_",x,".signature_evaluation.pdf"))
  pos <- which(sample_names==x)
  ggsave(file=paste0("cosmic_",x,".signature_evaluation.pdf"), plot=plot_compare_profiles(mut_mat[,pos], ff$reconstructed[,pos],profile_names = c("Original", "Reconstructed"), condensed = TRUE))
  #graphics.off()
})

us <- nmf_res$signatures
colnames(us) <- names_sign
hclust_cosmic_us = cluster_signatures(cbind(us, cancer_signatures), method = "average") # store signatures in new order
pdf("hclust_cosmic_us.pdf")
plot(hclust_cosmic_us)
graphics.off()

cos_sim_ori_rec <- cos_sim_matrix(mut_mat, ff$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
colnames(cos_sim_ori_rec) = "cos_sim"
cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)

ggplot(cos_sim_ori_rec, aes(y=cos_sim, x=sample)) +
  geom_bar(stat="identity", fill = "skyblue4") +
  coord_cartesian(ylim=c(0.8, 1)) +
  coord_flip(ylim=c(0.8,1)) +
  ylab("Cosine similarity\n original VS reconstructed") + 
  xlab("") + xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +theme_bw()+
  theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +geom_hline(aes(yintercept=.95))
ggsave("cosine_samples_cosmic.pdf")