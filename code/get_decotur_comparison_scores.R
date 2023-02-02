setwd('~/GitHub/DeCoTUR_manuscript_code')
source("code/helper_functions.R")
library(decotur)
allpd <- read.csv('data/all_pds.csv', stringsAsFactors = F)
allpd$pair <- paste0(allpd)
t <- proc.time()
for(numsam in c(100, 500, 1000)){
  distmat <- get_matrix(paste0('data/comparison', numsam, '_distmat.csv'), F)
  samps <- rownames(distmat)
  for(numgene in c(500, 1000, 5000, 10000)){
    genemat <- get_matrix(paste0('data/comparison', numsam, '_', numgene, '_pa.csv'), F)
    for(numcp in c(100, 500, 1000)){
      closepairs <- get_closepairs_fixednumber(distmat, numcp, 1, F, F)
      classes <- get_closepair_classes(closepairs)
      classweights <- get_classweights(classes)
      scores <- get_scores_pa_closepairs(genemat, closepairs, classweights, 200, T, 'speed', T)
      ss <- subset(scores, Score > quantile(scores$Score, 0.9))
      topgenes <- unique(c(ss$Trait1, ss$Trait2))
      gmf <- genemat[which(rownames(genemat) %in% topgenes),]
      discdat <- decotur::get_discordances(gmf, closepairs, T)
      dpds <- allpds$dist[match(samps, allpd$pair)]
      dpds[which(dpds <= 0)] <- 1/(2 * 704782)
      disc1 <- discdat$disc[match(ss$Trait1, discdat$trait)]
      disc2 <- discdat$disc[match(ss$Trait2, discdat$trait)]
      sum <- ss$UnweightedPositive + ss$UnweightedNegative
      m <- rep(1, dim(closepairs)[1])
      ndpds <- dpds/sum(dpds)
      blocksize <- 3000
      numblocks <- length(disc1)%/%blocksize
      leftover <- length(disc1)%%blocksize
      res <- c()
      for (i in 1:numblocks) {
        start <- (i - 1) * blocksize + 1
        end <- i * blocksize
        spbmat <- (disc1[start:end] * disc2[start:end]) %*% 
          t(ndpds^2)
        res <- c(res, pbapply(spbmat, 1, function(x) {
          PoissonBinomial::qpbinom(1 - 0.05/length(disc1), 
                                   x, method = "RefinedNormal")
        }))
        print(i)
      }
      if (leftover > 0) {
        start <- i * blocksize + 1
        end <- length(disc1)
        spbmat <- (disc1[start:end] * disc2[start:end]) %*% 
          t(ndpds^2)
        res <- c(res, apply(spbmat, 1, function(x) {
          PoissonBinomial::qpbinom(1 - 0.05/length(disc1), 
                                   x, method = "RefinedNormal")
        }))
      }
      ss$sig <- sum > res
      if(numsam == 500 & numgene == 1000){
        write.csv(ss, paste0('data/decotur_', numcp), row.names = F)
      }
      print(proc.time() - t)
    }
  }
}


