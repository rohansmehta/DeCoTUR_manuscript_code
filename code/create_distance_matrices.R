# User needs to input their own Staphopia %TOKEN
setwd('~/GitHub/DeCoTUR_manuscript_code')
source("code/helper_functions.R")
library(staphopia)   
TOKEN = %TOKEN

# This gets the annotation ids of the core genome
core_index <- read.delim('data/nrd-gene-set.txt', sep = '\t', header=T)
# This is all public samples from staphopia
ks <- read.csv('keptsamples.csv', stringsAsFactors = F)

for(CC in c('CC93', 'CC15', 'CC45', 'CC30', 'CC97', 'CC5','CC22', 'CC8', 'CC1', 'Other')){
  # This code block gets the sequence data in nexus format 
  samples <- ks[which(ks$cc == CC),]$sample_id
  t <- proc.time() # keeps track of elapsed time
  allg <- get_variant_gene_sequence(as.numeric(samples), annotation_ids = core_index$annotation_id) # gets the core genome sequences for the samples in a non-concatenated format
  print(proc.time()-t) # prints elapsed time
  allgmr <- subset(allg, sample_id != 'reference') # gets rid of the reference sequence
  gallgmr <- allgmr %>% group_by(sample_id) %>% mutate(fullseq  = paste0(sequence, collapse = '')) # concatenates the core genome sequences 
  gallgmr <- gallgmr[!duplicated(gallgmr$sample_id),] # I think the previous step has duplicates of samples; this gets rid of any duplicates
  nexusdat <- list() # this is the data structure to store the sequences
  # this loop stores the sequences in the data structure in a way that is understandable by R to write in a nexus format
  for(j in 1:length(samples)){
    nexusdat[[j]] <- strsplit(tolower(gallgmr$fullseq[j]), '')[[1]]
  }
  names(nexusdat) <- c(samples)# adds sample IDs to nexus data
  print(proc.time()-t)
  
  fbin <- as.DNAbin(nexusdat) # converts the nexus format into DNAbin format (compressed losslessly)
  names(fbin) <- samples # names the data properly
  distobj <- dist.dna(fbin)
  fdistmat <- as.matrix(distobj) # turns the "dist" object distobj into a matrix
  write.csv(x = fdistmat, file = paste0('data/', CC, '_distmat.csv')) # writes distance matrix to file
}
