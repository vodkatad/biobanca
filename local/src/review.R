data <- read.table('/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier_replisafe', sep="\t", stringsAsFactors = F, header=T)
head(data)
lmo <- data[grepl('LMO_BASALE', data$type),]
lmo$scong <- as.numeric(substr(lmo$sample_id_R, 13, 14))
lmo$passaggio <-as.numeric(substr(lmo$sample_id_R, 15, 17))

x_axis_labels <- min(lmo[,'scong']):max(lmo[,'scong'])
ggplot(data=lmo, aes(x=scong))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)

ggplot(data=lmo, aes(x=passaggio))+geom_bar()+theme_bw()


targ_buoni <- read.table('/scratch/trcanmed/biobanca/local/share/data/targeted/map_all_validation.tsv',  sep="\t", stringsAsFactors = F, header=T)
lmo_t <- targ_buoni[!is.na(targ_buoni$validation_status) & targ_buoni$validation_status & grepl('LMO', targ_buoni$genalogy),]

lmo_t$scong <- as.numeric(substr(lmo_t$genalogy, 13, 14))
lmo_t$passaggio <- as.numeric(substr(lmo_t$genalogy, 15, 17))
ggplot(data=lmo_t, aes(x=passaggio))+geom_bar()+theme_bw()
x_axis_labels <- min(lmo_t[,'scong']):max(lmo_t[,'scong'])
ggplot(data=lmo_t, aes(x=scong))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)
x_axis_labels <- min(lmo_t[,'passaggio']):max(lmo_t[,'passaggio'])
ggplot(data=lmo_t, aes(x=passaggio))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)



shallow_buoni <- read.table('/scratch/trcanmed/biobanca/local/share/data/shallowseq/map_all_validation.tsv',  sep="\t", stringsAsFactors = F, header=T)
lmo_t <- shallow_buoni[!is.na(shallow_buoni$validation_status) & shallow_buoni$validation_status & grepl('LMO', shallow_buoni$genalogy),]

lmo_t$scong <- as.numeric(substr(lmo_t$genalogy, 13, 14))
lmo_t$passaggio <- as.numeric(substr(lmo_t$genalogy, 15, 17))
ggplot(data=lmo_t, aes(x=passaggio))+geom_bar()+theme_bw()
x_axis_labels <- min(lmo_t[,'scong']):max(lmo_t[,'scong'])
ggplot(data=lmo_t, aes(x=scong))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)


reps <- as.data.frame(table(substr(lmo$sample_id_R, 0,7)))
table(reps$Freq)

has_rep <- reps[reps$Freq==2,'Var1']

lmo$rep <- ifelse(substr(lmo$sample_id_R, 0,7) %in% has_rep, 'yes', 'no')

lmo$model <- substr(lmo$sample_id_R, 0,7)
delta_col <- function(model, data, col) {
  subd <- data[data$model == model,]
  if (nrow(subd)==2) {
    return(abs(subd[1, col]-subd[2, col]))
  } else {
    return(NA)
  }
}

delta_scong <- sapply(has_rep, delta_col, lmo, 'scong')
delta_passaggio <- sapply(has_rep, delta_col, lmo, 'passaggio')

pd <- data.frame(model=has_rep, delta_scong=delta_scong, delta_passaggio=delta_passaggio)
ggplot(data=pd, aes(x=delta_passaggio))+geom_bar()+theme_bw()
ggplot(data=pd, aes(x=delta_scong))+geom_bar()+theme_bw()

lmo <- lmo[lmo$rep == "yes",]
x_axis_labels <- min(lmo[,'passaggio']):max(lmo[,'passaggio'])
ggplot(data=lmo, aes(x=passaggio))+geom_bar()+theme_bw()+scale_x_continuous(labels = x_axis_labels, breaks = x_axis_labels)

early_late <- function(model, data, col, low=3, high=8) {
  subd <- data[data$model == model,]
  subd <- subd[order(subd[,col]),]
  if (nrow(subd)==2) {
    if (subd[1, col] <= low && subd[2, col] >= high) {
      return(TRUE)
    } else{
      return(FALSE)
    }
  } else {
    return(NA)
  }
}

el <- sapply(has_rep, early_late, lmo, 'passaggio')


el2 <- sapply(has_rep, early_late, lmo, 'passaggio', low=10, high=15)
el3 <- sapply(has_rep, early_late, lmo, 'passaggio', low=5, high=15)

