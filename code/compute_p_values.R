setwd('~/GitHub/DeCoTUR_manuscript_code')
source("code/helper_functions.R")
library(decotur)
library(staphopia)
library(ggtext)
library(stringr)
library(pbapply)
library(PoissonBinomial)
verbose <- T
allpd <- read.csv('data/all_pds.csv', stringsAsFactors = F)
gene_matrix_filtered <- get_matrix('data/consolidated_gene_matrix.csv', F)
closepairs <- read.csv('data/closepairs.csv', stringsAsFactors = F)
geneorder <- colnames(gene_matrix_filtered)
samples1 <- geneorder[closepairs[,1]]
samples2 <- geneorder[closepairs[,2]]
samps <- paste0(samples1, samples2)
allpd$pair <- paste0(allpd$Sample1, allpd$Sample2)

ss <- read.csv('data/topscores.csv', stringsAsFactors = F)
#ss <- read.csv('data/lowlinkscores.csv', stringsAsFactors = F)
#ss <- read.csv('data/highlinkscores.csv', stringsAsFactors = F)
topgenes <- unique(c(ss$Trait1, ss$Trait2))
gmf <- gene_matrix_filtered[which(rownames(gene_matrix_filtered) %in% topgenes),]
# This code gets discordance distribution figure
discdat <- decotur::get_discordances(gene_matrix_filtered, closepairs, T)
discplot <- ggplot(discdat, aes(x = discdat)) + geom_histogram() + scale_y_log10() +
  xlab('Number of Discordant Close Pairs') + ylab('Count') + theme_bw()
ggsave('figures/discordance_plot.pdf', discplot)
# and back to p-values
discdat <- decotur::get_discordances(gmf, closepairs, T)
dpds <- allpd$dist[match(samps, allpd$pair)]
dpds[which(dpds <= 0)] <- 1/(2 * 704782)
disc1 <- discdat$disc[match(ss$Trait1, discdat$trait)]
disc2 <- discdat$disc[match(ss$Trait2, discdat$trait)]
sum <- ss$UnweightedPositive + ss$UnweightedNegative
m <- rep(1, dim(closepairs)[1])
ndpds <- dpds/sum(dpds)
ndpds[which(ndpds == 0)] <- min(ndpds[which(ndpds > 0)])/100
blocksize <- min(3000, dim(ss)[1])
numblocks <- length(disc1)%/%blocksize
leftover <- length(disc1)%%blocksize
res <- c()
for (i in 1:numblocks) {
  start <- (i - 1) * blocksize + 1
  end <- i * blocksize
  spbmat <- (disc1[start:end] * disc2[start:end]) %*%
    t(ndpds^2)
  res <- c(res, pbapply(spbmat, 1, function(x) {
    qpbinom_modified(1 - 0.05/length(disc1),
                             x, method = "RefinedNormal")
  }))
  print(i)
}
if (leftover > 0) {
  start <- i * blocksize + 1
  end <- length(disc1)
  spbmat <- (disc1[start:end] * disc2[start:end]) %*%
    t(ndpds^2)
  spbmat[which(is.na(spbmat))] <- 0
  res <- c(res, apply(spbmat, 1, function(x) {
    PoissonBinomial::qpbinom(1 - 0.05/length(disc1),
                             x, method = "RefinedNormal")
  }))
}
ss$sig <- sum > res
write.csv(ss, 'data/topscores_withpval.csv', row.names = F)
#write.csv(ss, 'data/lowlinkscores_withpval.csv', row.names = F)
#write.csv(ss, 'data/highlinkscores_withpval.csv', row.names = F)
