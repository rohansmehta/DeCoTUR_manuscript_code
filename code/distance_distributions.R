setwd('~/GitHub/DeCoTUR_manuscript_code')
source("code/helper_functions.R")

closepairtable <- read.csv('data/closepairtable.csv', stringsAsFactors = F)

distdat <- c()
ks <- read.csv('data/keptsamples.csv', stringsAsFactors = F, header=T)
for(CC in c('CC8', 'CC22', 'CC5', 'CC30', 'CC1', 'CC45', 'CC15', 'CC97', 'CC93')){
  samples <- subset(ks, cc == CC)$sample_id
  fdistmat <- get_matrix(paste('data/', CC, '_distmat.csv', sep=''), issnpmat = F)
  sfdistmat <- fdistmat[rownames(fdistmat) %in% samples, colnames(fdistmat) %in% samples]
  distdat <- rbind(distdat, data.frame(cc = CC, values = as.numeric(sfdistmat[upper.tri(sfdistmat)])))
}
write.csv(distdat, 'data/distdat.csv',row.names = F)
