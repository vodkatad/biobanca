lasvol_f <- snakemake@input[['las_vol']]
res_f <- snakemake@output[['slope_exp']]
log_f <- snakemake@log[['log']]


save.image('pippo.Rdata')


###
##prendere le misure -3W +2W attorno a inizio trattamento per fare fit esponenziale senza avere la fase piatta iniziale. 
#length(unique(keep_lasvol$exp_group))
#length(unique(lasvol$exp_group))
#length(unique(keep_lasvol$modelarm))
#length(unique(lasvol$modelarm))

#data <- keep_lasvol
data <- read.table(lasvol_f, sep="\t", header=TRUE, row.names=1)
data$start_date <- as.Date(data$start_date)
data$end_date <- as.Date(data$end_date)
data$measure_date <- as.Date(data$measure_date)


efit_all_mice_expg <- function(exp_group, data) {
  subset_data <- data[data$exp_group == exp_group,]
  start <- unique(subset_data$start_date)
  if (length(start) != 1) {
    print(paste0('Two start dates for ', exp_group, ' wtf?'))
    return(data.frame(expg=exp_group, slope=NA, pval=NA, R2=NA, n_measures=NA, logLik=NA))
  }
  interval <- c(start-21, start+14)
  subset_data <- subset_data[subset_data$measure_date >= interval[1] & subset_data$measure_date <= interval[2],]
  subset_data <- subset_data[order(subset_data$measure_date),]
  vels <- sapply(unique(subset_data$longen), function(x) efit_one_mouse(subset_data[subset_data$longen==x,]) )
  res <- t(vels)
  colnames(res) <- c('slope','pval','R2', 'n_measures', 'logLik')
  res <- as.data.frame(res)
  res$expg <- exp_group
  return(res)
}

lmp <- function (sm) {
  #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  #f <- summary(modelobject)$fstatistic
  f <- sm$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


efit_one_mouse <- function(mouse_d) {
  res = data.frame(vol=rep(0, nrow(mouse_d)), day=rep(0, nrow(mouse_d)))
  res[1,'volume'] <- mouse_d[1,'volume']
  for (i in seq(2, nrow(mouse_d))) {
    res[i,'volume'] <- mouse_d[i, 'volume']
    res[i,'day'] <- mouse_d[i, 'measure_date'] - mouse_d[1,'measure_date']
  }
  
  #ggplot(res,aes(day,volume))+geom_point()+geom_smooth(method='lm')
  #ggsave(paste0(unique(mouse_d$longen),".pdf"))
  res$volume <- log(res$volume)
  model <- lm(data=res, formula="volume~day")
  
  
  sm <- summary(model)
  r2 <- sm$r.squared
  pval <- lmp(sm) 
  ll <- logLik(model)
  vel <- model$coefficients[2]
  attributes(vel) <- NULL
  return(c(vel, pval, r2, nrow(mouse_d),ll))
}

all_expvels <- lapply(unique(data$exp_group), efit_all_mice_expg, data )
dfvel <- do.call(rbind, all_expvels)

#ggplot(data=dfvel, aes(x=R2))+geom_histogram()+theme_bw()
#ggplot(data=dfvel, aes(x=-log10(pval)))+geom_histogram()+theme_bw()
#ggplot(data=dfvel, aes(x=slope))+geom_histogram()+theme_bw()

#dfvel <- dfvel[!is.na(dfvel$pval),] # no NA
dfvel$padj <- p.adjust(dfvel$pval)
dfvel_sign <- dfvel[dfvel$padj < 0.05,]
dfvel_fil <- dfvel[dfvel$pval<0.05 & dfvel$R2 > 0.7,]

sink(log_f)
dim(dfvel)
dim(dfvel_sign)
dim(dfvel_fil)
sink()

summarize_expg <-  function(dfvel_g) {
  avg <- mean(dfvel_g$slope)
  sd <- sd(dfvel_g$slope)
  n <- nrow(dfvel_g)
  ave_days <- mean(dfvel$n_measures)
  return(c(avg, sd, n, ave_days))
}

all_avg <- sapply(unique(dfvel_fil$expg), function(x) summarize_expg(dfvel_fil[dfvel_fil$expg==x,]) )

avg <- as.data.frame(t(all_avg))
colnames(avg) <- c('mean','sd','nmice','ave_measures')
rownames(avg) <- unique(dfvel_fil$expg)
avg$modelarm <- substr(rownames(avg), 0, 12)

write.table(avg, file=res_f, sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

#interval <- read.table('/scratch/trcanmed/biobanca/local/share/data/p1_p2_p3.tsv', sep="\t",header=T)
#
#m <- merge(avg, interval, by.x='modelarm', by.y="Gen_ID")

#ggplot(data=avg, aes(x=mean))+geom_histogram()+theme_bw()


#efit_all_mice_expg('CRC1272LMX0A.2015-08-04.Cetuximab Standard.0', data)


#efit_all_mice_expg('CRC1145LMX0A.2014-11-11.Cetuximab Standard.0', data)


#data[data$exp_group== "CRC1145LMX0A.2014-11-11.Cetuximab Standard.0", ]

#Per i dati in vitro provare sia CTG che Imaging (l’imaging e` vero longitudinale per il ctg abbiamo solo dato a 1w che supponiamo esser partito da == numero cellule piastrate e con stima su media luminescenza (?) di 1 cellula si e` calcolato essere 7500).

save.image('pippo.Rdata')
q(save = "no", status = 0)

ctg <- read.table('/scratch/trcanmed/biobanca/local/share/data/CTG_growth_X_CORR_VALUES.txt', sep="\t",header=T, comment.char="")
colnames(ctg)[1] <- 'smodel'
colnames(ctg)[3] <- 'CTG_75k'

avg$smodel <- substr(avg$modelarm, 0, 7)

get_avg_avg <- function(smodel, data) {
  subset <- data[data$smodel == smodel, ]
  return(mean(subset$mean))
}
all_avg_smodel <- as.data.frame(sapply(unique(avg$smodel), get_avg_avg, avg))
colnames(all_avg_smodel) <- 'eslope'

 #
m <- merge(all_avg_smodel, ctg, by.x='row.names', by.y="smodel")

ggplot(data=m, aes(x=eslope, y=CTG_75k))+geom_point()+geom_smooth(method='lm')+theme_bw()

### volumes when starting fit

volume_all_mice_expg <- function(exp_group, data) {
  subset_data <- data[data$exp_group == exp_group,]
  start <- unique(subset_data$start_date)
  interval <- c(start-21, start+14)
  subset_data <- subset_data[subset_data$measure_date >= interval[1] & subset_data$measure_date <= interval[2],]
  subset_data <- subset_data[order(subset_data$measure_date),]
  vels <- sapply(unique(subset_data$longen), function(x) vol_one_mouse(subset_data[subset_data$longen==x,]) )
  res <- data.frame(volume=vels)
  res$expg <- exp_group
  return(res)
}


vol_one_mouse <- function(mouse_d) {
  res = data.frame(vol=rep(0, nrow(mouse_d)), day=rep(0, nrow(mouse_d)))
  res[1,'volume'] <- mouse_d[1,'volume']
  for (i in seq(2, nrow(mouse_d))) {
    res[i,'volume'] <- mouse_d[i, 'volume']
    res[i,'day'] <- mouse_d[i, 'measure_date'] - mouse_d[1,'measure_date']
  }
  
  return(res[1, 'volume'])
}

all_vols <- lapply(unique(data$exp_group), volume_all_mice_expg, data )
dfvol <- do.call(rbind, all_vols)
ggplot(data=dfvol, aes(x=volume))+geom_histogram(binwidth=20)+theme_bw()
