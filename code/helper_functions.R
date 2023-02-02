# Contains all helper functions and libraries needed to make the figures and data.
library(grid)
library(gtable)
library(cowplot)
library(ggplot2)
library(viridis)
library(dplyr)
library(igraph)
library(ggpubfigs)
library(ggraph)
library(fitdistrplus)
library(pbapply)
library(ape)
library(phangorn)
library(phytools)
library(data.table)
library(ggtree)
library(xtable)
library(rdist)
library(reshape2)
library(scales)
library(Hmisc)
library(PoissonBinomial)

# This function defines "not in"
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# This function reads in a matrix from a csv file
get_matrix <- function(filename, issnpmat){
  snm <- read.csv(filename, header=T)
  # This cleanup is needed because of how R writes and stores matrices as csvs
  row.names(snm) <- snm$X
  snm$X <- NULL
  colnames(snm) <- substring(colnames(snm), 2)
  colnames(snm)[colnames(snm) == 'eference'] <- 'reference' # only matters if we're using the reference.
  if(issnpmat){
    if(max(rowSums(snm)) > dim(snm)[2]/2){
      snm[which(rowSums(snm) > dim(snm)[2]/2),] <- 1-snm[which(rowSums(snm) > dim(snm)[2]/2),]
    }
  }
  # I also need to make sure we label the minor alleles 1, so I have this line of code that does this if so
  return(snm)
}

# This function takes a presence-absence matrix and a set of close pairs and
# keeps only traits that differ more than once across the set of samples
# in the close pairs. It also keeps only samples that are present in the 
# closepairs. it also re-indices the close pairs to match the new matrix
filter_snps_by_closepairs <- function(snm, closepairs){
  closepairinds <- sort(unique(c(closepairs[,1], closepairs[,2])))
  snmsub <- snm[,closepairinds]
  snpsum <- rowSums(snmsub)
  whichsnps <- which(snpsum > 1 & snpsum < length(closepairinds))
  # we need new closepairs which correspond to the new matrix
  cp1 <- match(closepairs[,1], closepairinds)
  cp2 <- match(closepairs[,2], closepairinds)
  newclosepairs <- cbind(cp1, cp2)
  return(list(snmsub[whichsnps,], newclosepairs))
}

