setwd("~/GitHub/DeCoTUR_manuscript_code")
library(data.table)
library(igraph)
library(reshape2)
library(ggplot2)
library(pbapply)
library(stringr)
library(stringi)
source("code/helper_functions.R")
# This code generates graph distances from a panaroo graph and
# creates the spatial-effect figures.
# Read in the graph from panaroo
g <- read_graph('data/final_graph.gml', format = 'gml')
# Read in the consolidated gene matrix (from panaroo, post-processed by 
# get_consolidated_gene_matrix.R)
pa_matrix <- get_matrix('data/consolidated_gene_matrix.csv', F)
# Create data frame of gene frequencies to filter
#gfdat <- data.frame(name = rownames(pa_matrix), freq = rowSums(pa_matrix)/dim(pa_matrix)[2])
#names(gfdat) <- c('name', 'freq')
#infreq <- which(names(V(g)) %in% gfdat$name[which(gfdat$freq < 0.001)])
#g1 <- delete_vertices(g, infreq)
# Compute pairwise graph distances by component; written to be more efficient
# by doing the small components first
cg <- components(g)
dg <- decompose(g)
dgl <- sapply(dg, length)
dglo <- rank(dgl, ties.method = 'first')
dgo <- order(dglo)
i <- dgo[[1]]
d <- dg[[i]]
pd <- distances(d)
mpd <- melt(pd)
ampd <- data.table(mpd)
totsum <- length(d)
for(i in dgo[2:length(dgo)]){
  d <- dg[[i]]
  pd <- distances(d)
  mpd <- melt(pd)
  lmpd <- list(ampd, mpd)
  ampd <- rbindlist(lmpd, use.names = T, fill = T, idcol = F)
  print(length(d))
  totsum <- totsum + length(d)
  print(totsum/length(g))
}
write.csv(ampd, 'data/graph_distances.csv', row.names = F)
ampd <- read.csv('data/graph_distances.csv', stringsAsFactors = F, header=T)
sampd <- c()
blocksize <- 10000000
numblocks <- dim(ampd)[1] %/% blocksize
leftover <-  dim(ampd)[1] %% blocksize
for(i in 1:numblocks){
  start <- (i - 1) * blocksize + 1
  end <- i * blocksize
  sa <-  subset(ampd[start:end,], !str_detect(Var1, 'group') & !str_detect(Var2, 'group'))
  sampd <- rbind(sampd, sa)
}
sa <-  subset(ampd[(numblocks*blocksize+1):dim(ampd)[1],], !str_detect(Var1, 'group') & !str_detect(Var2, 'group'))
sampd <- rbind(sampd, sa)
# and now we have to get gene names in the same format as they are in the score data
ps <- subset(ampd, as.character(Var1) < as.character(Var2))
pss <- ps %>% group_by(Var1, Var2) %>% summarise(valmean =  mean(value, na.rm = T), valmin = min(value, na.rm = T))
v1 <- pbsapply(pss$Var1, function(x){
  x1 <- str_split_fixed(x, '~~~', Inf)
  x2 <- unique(stri_replace(x1, regex = "\\_.*", replacement = ""))
  x3 <- x2[which(x2 != '')]
  if(x3[1] == 'group'){
    return(x1[1])
  } else{
    return(paste0(sort(x3), collapse = '_'))
  }
})
pss$v1 <- v1
v2 <- pbsapply(pss$Var2, function(x){
  x1 <- str_split_fixed(x, '~~~', Inf)
  x2 <- unique(stri_replace(x1, regex = "\\_.*", replacement = ""))
  x3 <- x2[which(x2 != '')]
  if(x3[1] == 'group'){
    return(x1[1])
  } else{
    return(paste0(sort(x3), collapse = '_'))
  }
})
pss$v2 <- v2
pssd <- as.data.table(pss)
write.csv(pssd, 'data/graph_distance_processed.csv', row.names = F)
#test <- read.csv('data/graph_distance_processed.csv', stringsAsFactors = F)
