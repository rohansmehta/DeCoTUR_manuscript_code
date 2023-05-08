setwd("~/GitHub/DeCoTUR_manuscript_code")
library(data.table)
library(igraph)
library(reshape2)
library(ggplot2)
library(pbapply)
source("code/helper_functions.R")
# Now we read in the scores
scores <- read.csv('data/panaroo_decotur.csv', stringsAsFactors = F)
os <- scores
# We now add the graph distances to this dataset, but for interpretability
# if we want, we can remove all interactions involving unidentified gene clusters ("group_XXXXX")
# os <- subset(scores, !str_detect(Trait1, 'group') & !str_detect(Trait2, 'group')) %>% arrange(Trait1, Trait2)
pssd <- read.csv('data/graph_distance_processed.csv', stringsAsFactors = F)
# the order of the genes in pssd within a pair can be different from that in os.
# first, we create two pair labels, one with each order
pssd$pair1 <- paste0(pssd$v1, pssd$v2)
pssd$pair2 <- paste0(pssd$v2, pssd$v1)
os$pair <- paste0(os$Trait1, os$Trait2)
which_os_pair1 <- match(os$pair, pssd$pair1)
which_os_pair2 <- match(os$pair, pssd$pair2)
# So each pair in os has either 0, 1, or 2 distances in this set
# 2 distances is due to a lack of consolidation in the graph distance matrix.
na1 <- is.na(which_os_pair1)
na2 <- is.na(which_os_pair2)
# table(na1+na2)
# there are 56926/5032378 = 1.1% 2, 32440/5032378 = 0.6% 0, and the rest are 1.
# for these, we take the minimum of the distances
dists <- rep(0, length(os$pair))

dists[na1 == 0 & na2 == 1] <- pssd$valmin[which_os_pair1[na1 == 0 & na2 == 1]]
dists[na1 == 1 & na2 == 0] <- pssd$valmin[which_os_pair2[na1 == 1 & na2 == 0]]
dists[na1 == 0 & na2 == 0] <- pmin(
  pssd$valmin[which_os_pair1[na1 == 0 & na2 == 0]], 
  pssd$valmin[which_os_pair2[na1 == 0 & na2 == 0]])
dists[na1 == 1 & na2 == 1] <- NA
os$dist <- dists
write.csv(os, 'data/panaroo_decotur_withdist.csv', row.names = F)
os <- read.csv('data/panaroo_decotur_withdist.csv', stringsAsFactors = F)
# Now we make the figures
os$score <- os$Score * (-1)^(os$PositiveAssociation < os$NegativeAssociation)
oss <- subset(os, !is.na(dist))
#devtools::install_github("EdwinTh/ggoutlier")
library(ggoutlier)
hist <- ggplot(oss) + geom_histogram(aes(x = score), fill = 'gray80') +  theme_bw(base_size = 18) + xlab('Score') + ylab('Count') + 
  scale_y_log10() + geom_vline(xintercept = c(-25, 25), linetype = 'dashed', linewidth=1, color = 'blue') + 
  geom_vline(xintercept = c(-60, 60), linetype = 'dashed', linewidth=1, color = 'orange') + 
  scale_x_continuous(breaks = c(-100, -60, -25, 0, 25, 60, 100), limits = c(-100, 100))
#hist <- ggoutlier_hist(oss, "score", -60, 60, fill = 'gray') + theme_bw() + xlab('Score') + ylab('Count') + 
# scale_y_log10()
ggsave('figures/score_histogram.pdf', hist)
oss$cat <- 1 * (oss$Score > 60 & oss$dist < 30) + 2 * (oss$Score > 25 & oss$dist > 30)
oss1 <- subset(oss, Score > 15)
# OK I think that removing all points that are <= 15 Score removes 97.8 = 98% of the data points while still getting
# the cutoff point across. So let's do that.
spatial_plot <- ggplot(oss1, aes(x = dist, y= score)) + geom_point(alpha = 0.3, aes(color = factor(cat))) + #geom_vline(xintercept = 30) + 
  ylim(-174, 174) + theme_bw() + xlab('Pangenome Graph Distance') + geom_hline(yintercept = c(60, -60, 25, -25), linetype = 'dashed') + 
  ylab('Coevolution Score') + geom_vline(xintercept = c(30), linetype = 'dotted') + 
  scale_color_manual(name = 'Type', labels = c('Other', 'High Linkage, High Score', 
                                               'Low Linkage, High Score'), values = c('gray80', 'orange', 'blue'))
ggsave('figures/spatial_plot.pdf', spatial_plot)
# highest scores
scores_top <- subset(os, Score > 60)
scores_top <- subset(scores_top, !str_detect(Trait1, 'group') & !str_detect(Trait2, 'group'))
write.csv(scores_top, 'data/topscores.csv', row.names = F)
# low linkage high scores
scores_low <- subset(oss, cat == 2)
scores_low <- subset(scores_low, !str_detect(Trait1, 'group') & !str_detect(Trait2, 'group'))
scores_low$cat <- NULL
write.csv(scores_low, 'data/lowlinkscores.csv', row.names = F)
# high linkage high scores
scores_high <- subset(os, Score > 60 & !is.na(dist) & dist < 30)
scores_high <- subset(scores_high, !str_detect(Trait1, 'group') & !str_detect(Trait2, 'group'))
write.csv(scores_high, 'data/highlinkscores.csv', row.names = F)
