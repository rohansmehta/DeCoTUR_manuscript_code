setwd('~/Github/DeCoTUR_manuscript_code')
library(decotur)
library(staphopia)
library(ggtext)
library(stringr)
library(pbapply)
TOKEN = '683eda5e910b094b15fd5276d0019f03b38200df'
# let us generate data for spydrpick and coinfinder
samplesize <- c(100, 500, 1000, 5000)
data <- read.csv('data/gene_presence_absence_roary.csv', stringsAsFactors = F)
colnames(data) <- str_remove(colnames(data), 'PROKKA_') # gets rid of "PROKKA_"
colnames(data) <- str_remove(colnames(data), '.fa') # gets rid of ".fa"
set.seed(21)
shufflecols <- sample(15:dim(data)[2], length(15:dim(data)[2]), replace  = F)
cshuff <- colnames(data)[shufflecols]
allsamples1 <- cshuff[1:max(samplesize)] # shuffled samples
for(col in 15:dim(data)[2]){
  nas <- which(is.na(data[,col]))
  data[nas,col] <- ''
}
for(numsam in c(100, 500, 1000, 5000)){
  allsamples <- allsamples1[1:numsam] # each subset is nested
  ndat <- cbind(data[,1], data[, which(colnames(data) %in% allsamples)])
  coindat <- cbind(data[,1:14], data[, which(colnames(data) %in% allsamples)])
  gf <- rowSums(ndat[, c(2:dim(ndat)[2])] != '')/(dim(ndat)[2]-1)
  gfs <- order(gf)
  center <- median(which(gf > 0.49 & gf < 0.51))
  centerrank <- gfs[center]
  for(targetno in c(1000, 3000, 5000, 10000)){
    ndat1 <- ndat[which(gfs < (centerrank + targetno/2) & gfs > (centerrank - targetno/2)),]
    names(ndat1)[1] <- 'Gene'
    coindat1 <- coindat[which(gfs < (centerrank + targetno/2) & gfs > (centerrank - targetno/2)),]
    names(coindat1)[1] <- 'Gene'
    for(col in 4:14){
      coindat1[,col] <- as.character(coindat1[,col])
      coindat1[which(is.na(coindat1[,col])), col] <- "NA"
    }
    write.table(ndat1, paste0('data/comparison', numsam, '_', targetno, '_pa.Rtab'), row.names = F, sep = '\t', col.names = T, quote = F)
    write.csv(coindat1,  paste0('data/comparison', numsam, '_', targetno, '_pa.csv'), row.names = F)
  } 
  core_index <- read.delim('data/nrd-gene-set.txt', sep = '\t', header=T)
  t <- proc.time() # keeps track of elapsed time
  allg <- get_variant_gene_sequence(as.numeric(allsamples), annotation_ids = core_index$annotation_id) # gets the core genome sequences for the samples in a non-concatenated format
  print(proc.time()-t) # prints elapsed time
  allgmr <- subset(allg, sample_id != 'reference') # gets rid of the reference sequence
  gallgmr <- allgmr %>% group_by(sample_id) %>% mutate(fullseq  = paste0(sequence, collapse = '')) # concatenates the core genome sequences 
  gallgmr <- gallgmr[!duplicated(gallgmr$sample_id),] # I think the previous step has duplicates of samples; this gets rid of any duplicates
  nexusdat <- list() # this is the data structure to store the sequences
  # this loop stores the sequences in the data structure in a way that is understandable by R to write in a nexus format
  for(j in 1:length(allsamples)){
    nexusdat[[j]] <- strsplit(tolower(gallgmr$fullseq[j]), '')[[1]]
  }
  names(nexusdat) <- c(allsamples)# adds sample IDs to nexus data
  print(proc.time()-t)
  fbin <- as.DNAbin(nexusdat) # converts the nexus format into DNAbin format (compressed losslessly)
  names(fbin) <- allsamples # names the data properly
  write.dna(fbin, paste0('data/comparison', numsam, '.fa'))
  distobj <- dist.dna(fbin)
  fdistmat <- as.matrix(distobj) # turns the "dist" object distobj into a matrix
  write.csv(x = fdistmat, file = paste0('data/comparison', numsam, '_distmat.csv')) # writes distance matrix to file
  fdistmat <- get_matrix('data/comparison', numsam, '_distmat.csv', F)
  pd <- as.phyDat(fbin)
  t <- proc.time()
  njtree <- NJ(fdistmat)
  fit <- pml(tree = njtree, data = pd)
  fitGTR <- optim.pml(fit, model="GTR", optInv=TRUE,
                      rearrangement = "NNI", control = pml.control(trace = 0))
  tree <- fitGTR$tree
  tmin <- min(tree$edge.length[which(tree$edge.length > 0)])
  tree$edge.length[which(tree$edge.length == 0)] <- tmin/100
  #tree$tip.label <- paste0('"', tree$tip.label, '"')
  write.tree(tree, 'data/comparison', numsam, '_tree.newick')
  print(proc.time() - t)
}