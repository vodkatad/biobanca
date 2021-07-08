

# We build ourselves two dataframes with some existing pairs (rownames here represents 'CRCxyz', our smodel).
# This is the same format that Martina starts from with Correlation_LMO_LMX-v2.R
# {We assume that we want to compute the correlation for all pairs for which it could be computed - we could add an option for
# the correct pairs (see Marco's script) to get the list of pairs and avoid going down to smodel in reality to have the maximum level of generalization}
data <- matrix(rnorm(100), ncol=20)
pdx <- data.frame(jitter(data), row.names=LETTERS[1:5])
pdo <- data.frame(jitter(data, amount=1), row.names=LETTERS[1:5])
# let's suppose that the last one is not a fair "pair" and that we also have a rogue pdx +1
rownames(pdx)[5] <- 'z'
pdx <- rbind(pdx, data.frame(row.names=c("rogue"), matrix(rnorm(20), ncol=20)))
# and let's plug in a correlation ==1 to check our code
pdo[1,] <- pdx[1,]
pdx <- t(pdx)
pdo <- t(pdo)


# All set!
# Now we get the looping approach (the one that is needed for example to compute lm as Martina did)

merged_df <- merge(pdx, pdo, by="row.names") # here it's not needed because I built the matrix with the same 'genes', but let's keep it
merged_df$Row.names <- NULL

# we define the smodels that are paired:
repeated_smodel <- as.data.frame(table(c(colnames(pdx), colnames(pdo))), stringsAsFactors=FALSE)
has_pair <- repeated_smodel[repeated_smodel$Freq == 2, 'Var1']

res1 <- data.frame(row.names=has_pair, pearson=rep(NA,length(has_pair)))
for (i in seq(1, length(has_pair))) {
  x <- paste0(has_pair[i],'.x')
  y <- paste0(has_pair[i],'.y')
  pearson <- cor.test(merged_df[, colnames(merged_df) == x], merged_df[, colnames(merged_df) == y])
  res1[i,'pearson'] <- pearson$estimate
}

# let's check that pearson for smodel A is 1
res1[rownames(res1)=="A",'pearson'] == 1

# btw in my opinion still easier with sapply:
# we write a function with what we want to do in our loop:
findcol_getcor <- function(smodel, data) {
  x <- paste0(smodel,'.x')
  y <- paste0(smodel,'.y')
  pearson <- cor.test(data[, colnames(data) == x], data[, colnames(data) == y])
  res <- pearson$estimate
  return(unname(res))
}
# use sapply and forget about having to pre-create res1 full of NAs
res2 <- as.data.frame(sapply(has_pair, findcol_getcor, merged_df))
colnames(res2) <- 'pearson'
# check again
all(res1==res2)


# Now if we want some 'wrong' pairs we have two ways, let's follow the loop one.
# We need to extract for each has_pair a 'wrong' match...
res3 <- data.frame(row.names=has_pair, pearson=rep(NA,length(has_pair)))

for (i in indexes) {
  x <- paste0(has_pair[i],'.x')
  # we get a random number that is included in seq(1, length(has_pair)) and is not == to i ...
  random <- sample(indexes[indexes != i], size=1) 
  # check out sample man page and look at why it's not safe to use if the selection on indexes
  # could leave us with a single element from indexes (right now it's not the case, but a corner case with 2 models only...)
  y <- paste0(has_pair[random],'.y')
  pearson <- cor.test(merged_df[, colnames(merged_df) == x], merged_df[, colnames(merged_df) == y])
  res3[i,'pearson'] <- pearson$estimate
}

# check that we have negative correlations

# Here we decided to get a number of random pairs == to the real ones...what if we want to best estimate the random expectations and extract
# all the possible random pairs?
# We can use cor if we are only interested in the pearson estimate (not pvalue or other things)
# We need to remove the not 'matched' pairs before and be sure that the columns are ordered in the same way
cpdo <- pdo[, colnames(pdo) %in% has_pair]
cpdx <- pdx[, colnames(pdx) %in% has_pair]
cpdo <- cpdo[, match(has_pair, colnames(cpdo))]
cpdx <- cpdx[, match(has_pair, colnames(cpdx))]
res4 <- cor(cpdo, cpdx)
pheatmap(res4)
res5 <- data.frame(pearson=diag(res4), row.names=rownames(res4))

all(res5==res1)

