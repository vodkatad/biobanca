### Generate a dataset containing short genealogy of a selected class, n of replicates
### and associated long genealogy
### There should be replicates, otherwise a warning is produced and an emtpy file will be put in output_f

library(tidyverse)
#library(janitor)

### Put the input file and load it
#meta_f <- "/Users/Tina/Desktop/Lavoro/bioinfo/ProgettiR/replicati_genealogy/selected_metadata_annot_final_nolinfo_nooutlier"
meta_f <- snakemake@input[['metadata']]
output_f <- snakemake@output[['replicates']]
classes <- snakemake@wildcards[['sclass']]

meta_df <- read.table(meta_f, header = TRUE, sep = "\t", quote = "")

### Select only a specific class
meta_df_lmo_basale<-filter(meta_df, grepl(classes, type))

### Add column with short genealogy and order them
meta_df_lmo_basale$model <- substr(meta_df_lmo_basale$sample_id_R, 0, 7)
meta_df_lmo_basale_ordered <- meta_df_lmo_basale[order(meta_df_lmo_basale$model),] 

### Sort for the model that are not unique
mlmob_ordered_not_unique <- meta_df_lmo_basale_ordered[meta_df_lmo_basale_ordered$model %in% meta_df_lmo_basale_ordered$model[duplicated(meta_df_lmo_basale_ordered$model)],]

# if there arent any duplicates we stop and don't print anything at all
if (nrow(mlmob_ordered_not_unique) == 0) {
  write.table(data.frame(), file=output_f)
  warning('No replicates found')
  quit("no")
}

### For each model obtain the number of replicates
genealogy_replicates <- as.data.frame(table(mlmob_ordered_not_unique$model))
colnames(genealogy_replicates) <- c("model", "number of replicates")
genealogy_replicates$genealogy <- rep("", nrow(genealogy_replicates))

### For each model obtain the long genealogy assoiated
#i=1 ### per vedere cosa fa il for assegnare una variabile così che vedo cosa fa ad un giro
for (i in seq(1, nrow(genealogy_replicates))) {
  shortgen <- genealogy_replicates[i, "model"]
  selection <- mlmob_ordered_not_unique[shortgen == mlmob_ordered_not_unique$model,]
  ###selection <- mlmob_ordered_not_unique[mlmob_ordered_not_unique$model == genealogy_replicate[i, "model"]]
  ###genealogy_replicates[i, "genealogy"] <- paste0(selection[1, "sample_id_R"], ",", selection[2, "sample_id_R"])
  genealogy_replicates[i, "genealogy"] <- paste0(selection[,"sample_id_R"], collapse = ",")
  ###così funziona anche se ci fossero triplicati 
}
### il for poi mi fa vedere solo gli ultimi due valori, giusto 

write.table(genealogy_replicates, file=output_f, sep="\t", quote=FALSE)




