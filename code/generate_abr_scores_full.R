setwd('~/GitHub/DeCoTUR_manuscript_code')
source("code/helper_functions.R")
library(staphopia)
TOKEN = '683eda5e910b094b15fd5276d0019f03b38200df'
library(decotur)


get_scores_pa_closepairs_abr <- function (pa_matrix, close_pairs, classweights, blocksize, verbose, 
          version, withpval) 
{
  if (blocksize > dim(pa_matrix)[1]) {
    stop("Block size is larger than the number of traits")
  }
  if (withpval) {
    uweights <- rep(1, length(classweights))
  }
  gene_names <- rownames(pa_matrix)
  nrows <- dim(pa_matrix)[1]
  ncols <- dim(pa_matrix)[2]
  numblocks <- (nrows%/%blocksize)
  leftover <- nrows%%blocksize
  if (leftover > 0) {
    padding <- blocksize - leftover
    padrows <- matrix(rep(0, padding * (ncols)), nrow = padding, 
                      ncol = dim(pa_matrix)[2])
    rownames(padrows) <- rep("padding", padding)
    colnames(padrows) <- c(colnames(pa_matrix))
    pa_matrix <- rbind(pa_matrix, padrows)
    numblocks <- numblocks + 1
  }
  scoredat <- c()
  if (verbose) {
    pb <- progress_bar$new(format = "  Computing scores [:bar] :percent eta: :eta", 
                           total = numblocks + (leftover > 0), clear = FALSE)
    pb$tick(0)
  }
  for (i in 1:numblocks) {
    print(i)
    starti <- (i - 1) * blocksize + 1
    endi <- i * blocksize
    rangei <- starti:endi
    ind1snp1 <- pa_matrix[rangei, close_pairs[, 1]]
    ind2snp1 <- pa_matrix[rangei, close_pairs[, 2]]
    for (j in i:numblocks) {
      print(j)
      startj <- (j - 1) * blocksize + 1
      endj <- j * blocksize
      rangej <- startj:endj
      if (i == j) {
        pairindices <- combn(1:blocksize, 2)
      }
      else {
        pairindices <- t(expand.grid(1:length(rangei), 
                                     1:length(rangej)))
      }
      ind1snp2 <- pa_matrix[rangej, close_pairs[, 1]]
      ind2snp2 <- pa_matrix[rangej, close_pairs[, 2]]
      snpdiffmat1 <- 1 * (ind1snp1 == 0 & ind2snp1 == 1) + 
        2 * (ind1snp1 == 1 & ind2snp1 == 0)
      snpdiffmat2 <- 1 * (ind1snp2 == 0 & ind2snp2 == 1) + 
        2 * (ind1snp2 == 1 & ind2snp2 == 0)
      if (version == "speed") {
        test1 <- snpdiffmat1[pairindices[1, ], ]
        test2 <- snpdiffmat2[pairindices[2, ], ]
        test12 <- test1 == test2
        test3 <- snpdiffmat1[pairindices[1, ], ] > 0
        test4 <- snpdiffmat2[pairindices[2, ], ] > 0
        pos_test <- test12 * test3 * test4
        snpscoremats <- as.numeric(pos_test %*% classweights)
        if (withpval) {
          us <- as.numeric(pos_test %*% uweights)
        }
        testn12 <- test1 != test2
        neg_test <- testn12 * test3 * test4
        snpscorematd <- as.numeric(neg_test %*% classweights)
        if (withpval) {
          ud <- as.numeric(neg_test %*% uweights)
        }
      }
      else if (version == "memory") {
        snpscoremats <- pbapply(pairindices, 2, function(x) {
          return((as.numeric(snpdiffmat1[x[1], ] == snpdiffmat2[x[2], 
          ] & snpdiffmat1[x[1], ] > 0 & snpdiffmat2[x[2], 
          ] > 0)) %*% classweights)
        })
        snpscorematd <- pbapply(pairindices, 2, function(x) {
          return((as.numeric(snpdiffmat1[x[1], ] != snpdiffmat2[x[2], 
          ] & snpdiffmat1[x[1], ] > 0 & snpdiffmat2[x[2], 
          ] > 0)) %*% classweights)
        })
        if (withpval) {
          us <- pbapply(pairindices, 2, function(x) {
            return((as.numeric(snpdiffmat1[x[1], ] == 
                                 snpdiffmat2[x[2], ] & snpdiffmat1[x[1], 
                                 ] > 0 & snpdiffmat2[x[2], ] > 0)) %*% uweights)
          })
          ud <- pbapply(pairindices, 2, function(x) {
            return((as.numeric(snpdiffmat1[x[1], ] != 
                                 snpdiffmat2[x[2], ] & snpdiffmat1[x[1], 
                                 ] > 0 & snpdiffmat2[x[2], ] > 0)) %*% uweights)
          })
        }
      }
      else {
        stop("Incorrect version. Should be speed or memory")
      }
      snp1 <- rownames(ind1snp1)[pairindices[1, ]]
      snp2 <- rownames(ind1snp2)[pairindices[2, ]]
      score <- pmax(snpscoremats, snpscorematd)
      if (!withpval) {
        us <- NA
        ud <- NA
      }
      scoredat <- rbind(scoredat, data.frame(snp1, snp2, 
                                             snpscoremats, snpscorematd, score, us, ud, stringsAsFactors = F))
    }
    if (verbose) {
      pb$tick()
    }
  }
  names(scoredat) <- c("Trait1", "Trait2", "PositiveAssociation", 
                       "NegativeAssociation", "Score", "UnweightedPositive", 
                       "UnweightedNegative")
  scoredat <- subset(scoredat, Trait1 %in% gene_names & Trait2 %in% 
                       gene_names)
  return(scoredat)
}


get_scores_abr <- function (pa_matrix, distance_matrix, closepair_method, closepair_params, 
          blocksize, downweight = TRUE, withsig = TRUE, verbose = TRUE, 
          version = "speed") 
{
  if (verbose) {
    print("Starting function.")
  }
  if (closepair_method == "distance") {
    if (length(closepair_params) != 2) {
      stop("Incorrect number of closepair_params. (Should be 2).")
    }
    distance_cutoff <- closepair_params[[1]]
    show_hist <- closepair_params[[2]]
    close_pairs <- get_closepairs_distance(distance_matrix, 
                                           distance_cutoff, show_hist, verbose)
  }
  else if (closepair_method == "fraction") {
    if (length(closepair_params) != 2) {
      stop("Incorrect number of closepair_params. (Should be 2).")
    }
    distance_fraction <- closepair_params[[1]]
    show_hist <- closepair_params[[2]]
    close_pairs <- get_closepairs_fraction(distance_matrix, 
                                           distance_fraction, show_hist, verbose)
  }
  else if (closepair_method == "auto") {
    if (length(closepair_params) != 4) {
      stop("Incorrect number of closepair_params. (Should be 4).")
    }
    which_valley <- closepair_params[[1]]
    nbins <- closepair_params[[2]]
    maxvalleyheight <- closepair_params[[3]]
    show_hist <- closepair_params[[4]]
    close_pairs <- get_closepairs_auto(distance_matrix, which_valley, 
                                       nbins, maxvalleyheight, show_hist, verbose)
  }
  else if (closepair_method == "fixednumber") {
    number_close_pairs <- closepair_params[[1]]
    seed <- closepair_params[[2]]
    show_hist <- closepair_params[[3]]
    close_pairs <- get_closepairs_fixednumber(distance_matrix, 
                                              number_close_pairs, seed, show_hist, verbose)
  }
  else {
    stop("Unidentified close pair method (should be one of: distance, fraction, fixednumber, auto).")
  }
  if (verbose) {
    print(paste0("Obtained ", dim(close_pairs)[1], " close pairs."))
  }
  classes <- get_closepair_classes(close_pairs)
  if (verbose) {
    print("Obtained close pair classes.")
  }
  if (downweight) {
    class_weights <- get_classweights(classes)
  }
  else {
    class_weights <- rep(1, length(classes))
  }
  if (verbose) {
    print("Obtained close pair class weights.")
  }
  pa_matrix <- filter_snps_by_closepairs(pa_matrix, close_pairs)
  scores <- get_scores_pa_closepairs_abr(pa_matrix, close_pairs, 
                                     class_weights, blocksize, verbose, version, withsig)
  if (verbose) {
    print("Obtained scores.")
  }
  if (withsig) {
    if (verbose) {
      print("Computing discordance information.")
    }
    discdat <- get_discordances(pa_matrix, close_pairs, verbose)
    if (verbose) {
      print("Obtained discordance information.")
    }
    if (verbose) {
      print("Computing significance")
    }
    dpds <- distance_matrix[close_pairs]
    dpds[which(dpds < 0)] <- 1/(2 * 704782)
    disc1 <- discdat$disc[match(scores$Trait1, discdat$trait)]
    disc2 <- discdat$disc[match(scores$Trait2, discdat$trait)]
    sum <- scores$UnweightedPositive + scores$UnweightedNegative
    m <- rep(1, dim(close_pairs)[1])
    ndpds <- dpds/sum(dpds)
    numblocks <- length(disc1)%/%blocksize
    leftover <- length(disc1)%%blocksize
    res <- c()
    for (i in 1:numblocks) {
      start <- (i - 1) * blocksize + 1
      end <- i * blocksize
      spbmat <- (disc1[start:end] * disc2[start:end]) %*% 
        t(ndpds^2)
      res <- c(res, apply(spbmat, 1, function(x) {
        PoissonBinomial::qpbinom(1 - 0.05/length(disc1), 
                                 x, method = "RefinedNormal")
      }))
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
    scores$sig <- sum > res
    if (verbose) {
      print("Obtained significance")
    }
  }
  return(scores)
}

allpd <- read.csv('data/all_pds.csv', stringsAsFactors = F)
cutoff <- 0.00005
nonz <- subset(allpd, cc %not in% c('CC5', 'CC8', 'CC22') )
nonz <- subset(nonz, !is.na(dist))
nonz <- subset(nonz, Sample1 != Sample2)
cutdists <- subset(nonz, dist <= cutoff)
perrcp <- dim(cutdists)[1]/length(unique(cutdists$cc))
otherpd <- c()
for(cc1 in c('CC5', 'CC8', 'CC22')){
  subpd <- subset(allpd, cc == cc1)
  cutsub <- subset(subpd, dist <= cutoff)
  set.seed(which(c('CC5', 'CC8', 'CC22') == cc1))
  rowsample <- sample(x = 1:dim(cutsub)[1], size = perrcp, replace = F)
  cutsub1 <- cutsub[rowsample,]
  rownames(cutsub1) <- NULL
  otherpd <- rbind(otherpd, cutsub1)
}
cutoq <- rbind(cutdists, otherpd)
cutsamp <- unique(c(cutoq$Sample1, cutoq$Sample2))

abr_res <- staphopia::get_resistance_results(cutsamp, resistance_report = T)
genemat <- t(abr_res)
colnames(genemat) <- genemat[1,]
genemat <- genemat[2:dim(genemat)[1],]

# We need to reformat the closepairs to fit the function
geneorder <- colnames(genemat)

# Some of the samples do not have abr_info, so we must reduce the sample by these
cutq1 <- subset(cutoq, Sample1 %in% geneorder & Sample2 %in% geneorder)

rs1 <- match(cutq1$Sample1, geneorder)
rs2 <- match(cutq1$Sample2, geneorder)

closepairsc <- cbind(rs1, rs2)
set.seed(1)
closepairs <- closepairsc[sample(1:dim(closepairsc)[1], 40000, replace = F),]
classes <- get_closepair_classes(closepairs)
classweights <- get_classweights(classes)
scores <- get_scores_abr(genemat, closepairs, classweights, 22)
scores$polarity <- scores$snpscoremats > scores$snpscorematd
classweights <- rep(1, dim(closepairs)[1])
scores_unweighted <- get_scores_abr(genemat, closepairs, classweights, 22)
gene_matrix_filtered <- filter_snps_by_closepairs(genemat, closepairs)
discdat <- pbapply(gene_matrix_filtered, 1, function(x){get_discordances(x, closepairs)})
ddiscdat <- data.frame(discdat)
ddiscdat$gene <- rownames(ddiscdat)
rownames(ddiscdat) <- NULL
ddiscdat <- ddiscdat[, c(2, 1)]
write.csv(ddiscdat, paste0('data/abr_discordances.csv'), row.names = F)
scores_unweighted$sum <- scores_unweighted$snpscoremats + scores_unweighted$snpscorematd
scores_unweighted$disc1 <- ddiscdat$disc[match(scores_unweighted$snp1, ddiscdat$gene)]
scores_unweighted$disc2 <- ddiscdat$disc[match(scores_unweighted$snp2, ddiscdat$gene)]
# Some of these abrs are unvaried across the whole sample.
scores_unweighted$disc1[which(is.na(scores_unweighted$disc1))] <- 0
scores_unweighted$disc2[which(is.na(scores_unweighted$disc2))] <- 0
numclosepairs <- dim(closepairs)[1]
fex <- data.frame(scores_unweighted$sum, scores_unweighted$disc1 -scores_unweighted$sum, scores_unweighted$disc2 - scores_unweighted$sum, numclosepairs - scores_unweighted$disc1 - scores_unweighted$disc2 + scores_unweighted$sum)
fres <- pbapply(fex,1, function(x) fisher.test(matrix(x,nr=2), alternative = 'greater')$p.value)
scores_unweighted$fpval <- fres
scores_unweighted$fpvalbonferroni <- p.adjust(scores_unweighted$fpval, method = 'bonferroni')
scores_unweighted$fpvalbh <- p.adjust(scores_unweighted$fpval, method = 'BH')
write.csv(x = scores_unweighted, file = paste0('data/abr_scores_unweighted_withnull.csv'))
scores$pval <- scores_unweighted$fpvalbonferroni
scores$polarity <- scores$snpscoremats > scores$snpscorematd
write.csv(scores, paste0('data/abr_scores_withpval.csv'), row.names = F)

scores <- read.csv('data/abr_scores_withpval.csv', stringsAsFactors = F)
scores$signedscore <- rescale(scores$score, to = c(0, 1)) * (-1)^(1-scores$polarity)
# We also need the correlations
correlations <- rcorr(t(genemat))
mcorr <- melt(correlations$r)
mpcorr <- melt(correlations$P)
#mcorr$value <- rescale(mcorr$value, to = c(min(scores$signedscore, na.rm = T), max(scores$signedscore, na.rm = T)))
toadd <- subset(mcorr, as.character(Var1) > as.character(Var2))
scores_toadd <- scores[, c(1, 2, 8)]
#scores1 <- rbind(scores_toadd, data.frame(snp1 = toadd$Var1, snp2 = toadd$Var2, signedscore = toadd$value))
scores1 <- rbind(scores_toadd, data.frame(snp1 = scores_toadd$snp2, snp2 = scores_toadd$snp1, signedscore = scores_toadd$signedscore))
cormat <- acast(scores1, snp1~snp2, value.var = 'signedscore')
cormat[is.na(cormat)] <- 0
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
c1 <- cormat[which(rownames(cormat) %not in% c('Bacitracin', 'Cationic Antimicrobial Peptides',
                                               'Elfamycins', 'Mycobacterium Tuberculosis-Specific Drug',
                                               'Thiostrepton', 'Tunicamycin')),]
c2 <- c1[,which(colnames(c1) %not in% c('Bacitracin', 'Cationic Antimicrobial Peptides',
                                               'Elfamycins', 'Mycobacterium Tuberculosis-Specific Drug',
                                               'Thiostrepton', 'Tunicamycin'))]

reordered_cormat <- reorder_cormat(c2)

mr <- melt(reordered_cormat)
sigs <- subset(scores, pval < 0.05)
sigpairs <- paste0(sigs$snp1, sigs$snp2)
sigpairs1 <- paste0(sigs$snp2, sigs$snp1)
sigpairsall <- c(sigpairs, sigpairs1)
mrs <- subset(mr, paste0(mr$Var1, mr$Var2) %in% sigpairsall)
#OK we're going to do this a bit strangely
for(i in 1:dim(mrs)[1]){
  snp1 <- as.character(mrs$Var1[i])
  snp2 <- as.character(mrs$Var2[i])
  mcorrrow <- which(as.character(mcorr$Var1) == snp1 & as.character(mcorr$Var2) == snp2)
  mrs$value1[i] <- mcorr$value[mcorrrow]
}
# Let's try this
for(i in 1:dim(mrs)[1]){
  if(as.numeric(mrs$Var1[i]) < as.numeric(mrs$Var2[i])){
    whichp <- which(as.character(mpcorr$Var1) == as.character(mrs$Var1[i]) & as.character(mpcorr$Var2) == as.character(mrs$Var2[i]))
    if(mpcorr$value[whichp] < 0.05){
      mrs$value[i] <- mrs$value1[i]
    } else{
      mrs$value[i] <- NA
    }

  }
}
heatplot <- ggplot(mrs, aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = 'white', lwd = 1.5) +
  scale_fill_gradient2(name = 'Scaled Value', na.value = 'white', low = muted('blue'), high = muted('red'), breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1, 1)) + xlab('Coevolution Score') + ylab('Correlation') +
  theme_classic(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_equal()
ggsave('figures/abr_heatplot.pdf', heatplot, scale = 1)
