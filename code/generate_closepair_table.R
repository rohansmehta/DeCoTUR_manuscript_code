setwd('~/GitHub/DeCoTUR_manuscript_code')
source("code/helper_functions.R")
library(stringr)

gene_matrix <- read.csv('data/gene_presence_absence.csv', stringsAsFactors = F)
samples <- colnames(gene_matrix)
samples <- samples[2:length(samples)]
samples <- str_remove(samples, 'PROKKA_')
samples <- str_remove(samples, '.fa')
ks <- read.csv('data/keptsamples.csv')
numclosepairs <- c()
ccvec_zerodist <- c('CC8', 'CC5', 'CC22')
for(cc in c('CC93', 'CC97', 'CC15', 'CC45', 'CC1', 'CC30')){
  distance_matrix_full <- get_matrix(paste0('data/', cc, '_distmat.csv'), F) # This includes all samples, not just kept samples.
  #gene_matrix_samples <- ks$sample_id[which(ks$cc == cc & ks$sample_id %in% samples)]
  gene_matrix_samples <- ks$sample_id[which(ks$cc == cc)]
  distance_matrix <- 
    distance_matrix_full[which(rownames(distance_matrix_full) %in% gene_matrix_samples),
                         which(colnames(distance_matrix_full) %in% gene_matrix_samples)]
  
  number_close_pairs <- 5000
  number_distances <- dim(distance_matrix)[1]*(dim(distance_matrix)[2]-1)/2
  fraction_distances <- number_close_pairs/number_distances
  distance_cutoff <- quantile(x = as.numeric(as.matrix(distance_matrix)), 
                              probs = fraction_distances, na.rm = T)
  closepairs <- get_closepairs_distance(distance_matrix, distance_cutoff, F, T)
  numclosepairs <- rbind(numclosepairs, c(cc, distance_cutoff, dim(closepairs)[1]))
}
for(cc in ccvec_zerodist){
  distance_matrix_full <- get_matrix(paste0('data/', cc, '_distmat.csv'), F) # This includes all samples, not just kept samples.
  #gene_matrix_samples <- ks$sample_id[which(ks$cc == cc & ks$sample_id %in% samples)]
  gene_matrix_samples <- ks$sample_id[which(ks$cc == cc)]
  distance_matrix <- 
    distance_matrix_full[which(rownames(distance_matrix_full) %in% gene_matrix_samples),
                         which(colnames(distance_matrix_full) %in% gene_matrix_samples)]
  number_close_pairs <- 5000
  number_distances <- dim(distance_matrix)[1]*(dim(distance_matrix)[2]-1)/2
  fraction_distances <- number_close_pairs/number_distances
  distance_cutoff <- quantile(x = as.numeric(as.matrix(distance_matrix)), 
                              probs = fraction_distances, na.rm = T)
  closepairs_many <- get_closepairs_distance(distance_matrix, distance_cutoff, F, T)
  if(distance_cutoff == 0 & dim(closepairs_many)[1] > number_close_pairs){
    set.seed(which(ccvec_zerodist == cc)) # different seed for each CC
    subsample <- sample(1:dim(closepairs_many)[1], number_close_pairs, replace = F)
    closepairs <- closepairs_many[subsample,]
  } else{
    closepairs <- closepairs_many
  }
  numclosepairs <- rbind(numclosepairs, c(cc, distance_cutoff, dim(closepairs)[1]))
}
closepairdat <- as.data.frame(numclosepairs)
names(closepairdat) <- c('Clonal Complex', 'Distance Cutoff', 'No. Close Pairs')
closepairdat$`Distance Cutoff` <- as.numeric(closepairdat$`Distance Cutoff`)
closepairdat <- closepairdat[order(-closepairdat$`Distance Cutoff`),]
write.csv(closepairdat, 'data/closepairtable.csv', row.names = F)
#cptab <- xtable(closepairdat, type = "latex")
#digits(cptab) <- 10
#print(cptab, file = "tables/close_pair_table.tex", include.rownames = F)
