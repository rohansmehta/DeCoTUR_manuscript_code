# User has to specify their Staphopia %TOKEN
setwd(~/Github/DeCoTUR_manuscript_code)
source("code/helper_functions.R")
library(staphopia)
TOKEN = %TOKEN

core_index <- read.delim('data/nrd-gene-set.txt', sep = '\t', header=T, stringsAsFactors = F)

ps_map <- read.csv('data/sample_cc_map.csv', stringsAsFactors = F)
 
for(cc in c('CC93', 'CC97', 'CC15', 'CC45', 'CC1', 'CC30', 'CC5', 'CC22', 'CC8', 'Other')){
  ss <- subset(ps_map, clonal_complex == cc)
  samples <- ss$sample_id
  all_genes_cc <- get_genes(as.numeric(ss$sample_id), exclude_sequence = T)
  # So these genes include core genes as well. Well, let's first exclude unnamed genes.
  allgenesn <- subset(all_genes_cc, name != 'none')
  allgenesnc <- subset(allgenesn, name %not in% core_index$gene)
  genetable <- table(allgenesnc$sample_id, allgenesnc$name)
  genemat <- as.matrix(genetable, rownames = Var2)
  genemat[which(genemat > 1)] <- 1
  genemat <- t(genemat)
  accdists <- pdist(genemat, metric = 'hamming')
  rownames(accdists) <- rownames(genemat)
  colnames(accdists) <- rownames(genemat)
  macc <- melt(accdists)
  fdistmat <- get_matrix(paste0('data/', cc, '_distmat.csv'), F)
  fdistmat <- as.matrix(fdistmat)
  rownames(fdistmat) <- colnames(fdistmat)
  coredists <- melt(fdistmat)
  macc1 <- subset(macc, Var1 < Var2)
  core1 <- subset(coredists, Var1 < Var2)
  bothdists <- merge(core1, macc1, by = c('Var1', 'Var2'))
  names(bothdists) <- c('Sample1', 'Sample2', 'coredist', 'accdist')
  write.csv(x = bothdists, file = paste('data/', cc, '_bothdists.csv', sep=''), row.names = F)
  bothzero <- bothdists[which(bothdists$coredist == 0 & bothdists$accdist == 0),]
  if(dim(bothzero)[1] > 0){
    graphdat <- data.frame(v1 = bothzero$Sample1, v2 = bothzero$Sample2, weight = 1)
    graph <- graph_from_data_frame(graphdat, directed = FALSE)
    plot(graph)
    components <- components(graph)
    torem <- c()
    for(i in 1:components$no){
      members <- names(components$membership[which(components$membership == i)])
      torem <- c(torem, members[2:length(members)])
    }
    paredsamples <- samples[which(samples %not in% torem)]
  } else{
    paredsamples <- samples
  keptsamples <- cbind(paredsamples, cc)
  write.table(x = keptsamples, file = 'keptsamples.csv', row.names = F, sep = ',', append = T, col.names = F)
  print(cc)
  }
}