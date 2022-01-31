library(ggplot2)
library(ggsankey)
library(psych)
# installed by hand for me 
# devtools::install_github("davidsjoberg/ggsankey", ref="main")


setwd('/scratch/trcanmed/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected')

#array_f <- "nc_cris_uarray_0.2.tsv"
#us_f <- "cris_fpkm_lmx_nc_arm.tsv"

#us_f <- "cris_tmm_0.2_classes_lmx_basali_ns.tsv"
#array_f <- "cris_uarray_0.2.tsv"


array_f <- "nc_cris_uarray_0.2.tsv"
us_f <- "cris_vsd_lmx_nc_arm.tsv"

array <- read.table(array_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)
us <- read.table(us_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)

array$cris <- gsub('CRIS', 'CRIS-', array$cris, fixed=TRUE)

#https://r-charts.com/flow/sankey-diagram-ggplot2/
m <- merge(array, us, by="genealogy")
rownames(m) <- m$genealogy
m$genealogy <- NULL
colnames(m) <- c("cris_uarray", "cris_vsd_nc_collapsedarm")

df <- m %>%
  make_long("cris_uarray", "cris_vsd_nc_collapsedarm")

ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 16)

# useful metrics
k <- cohen.kappa(m)
k$kappa
dim(m)
table(m[,1] == m[,2])
