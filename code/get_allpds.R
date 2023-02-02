setwd(~/Github/DeCoTUR_manuscript_code)
source("code/helper_functions.R")

# Let's get all the close pairs for each cc.
# We skip CC5, CC8, and CC22 for now.
ks <- read.csv('data/keptsamples.csv', stringsAsFactors = F)
ccvec <- unique(ks$cc)
allpd <- c()
for(cc1 in ccvec){
  ssk <- subset(ks, cc == cc1)$sample_id
  fdistmat <- get_matrix(paste0('data/', cc1, '_distmat.csv'), F)
  sfdistmat <- fdistmat[which(rownames(fdistmat) %in% ssk), which(colnames(fdistmat) %in% ssk)]
  mfd <- melt(as.matrix(sfdistmat))
  mfd$cc <- cc1
  allpd <- rbind(allpd, mfd)
  print(cc1)
}
names(allpd) <- c('Sample1', 'Sample2', 'dist', 'cc')
write.csv(x = allpd, file = 'data/all_pds.csv', row.names = F)