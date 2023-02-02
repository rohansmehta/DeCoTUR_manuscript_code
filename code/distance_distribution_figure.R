setwd('~/GitHub/DeCoTUR_manuscript_code')
source("code/helper_functions.R")
library(ggridges)

distdat <- read.csv('data/distdat.csv', stringsAsFactors = F)
cutoffdat <- read.csv('data/closepairtable.csv', stringsAsFactors = F)
names(cutoffdat) <- c('cc', 'cutoff', 'ncp')
distdat$ordered_factor <- factor(distdat$cc, levels = c('CC8', 'CC1', 'CC22', 'CC97', 'CC30', 'CC45', 'CC5', 'CC15', 'CC93'))
cutoffdat$ordered_factor <-  factor(cutoffdat$cc, levels = c('CC8', 'CC1', 'CC22', 'CC97', 'CC30', 'CC45', 'CC5', 'CC15', 'CC93'))
cutoffdat$cutoff[which(cutoffdat$cutoff == 0)] <- NA
p1 <- ggplot() +
  geom_density_ridges(data=subset(distdat, values > 0), aes(x=values, y = ordered_factor), scale = 0.9) + 
  geom_segment(data = cutoffdat, aes(x = cutoff, xend=cutoff, y = as.numeric(ordered_factor), yend = as.numeric(ordered_factor)+0.9), linetype = 'dashed') + 
  xlab('Core Nucleotide Divergence')  + theme_ridges(center_axis_labels = T) + ylab('Clonal Complex')
p2 <- ggplot() +
  geom_density_ridges(data=subset(distdat, values > 0), aes(x=values, y = ordered_factor), scale = 0.9) + 
  geom_segment(data = cutoffdat, aes(x = cutoff, xend=cutoff, y = as.numeric(ordered_factor), yend = as.numeric(ordered_factor)+0.9), linetype = 'dashed') + 
  scale_x_log10() + 
  xlab('Core Nucleotide Divergence')  + theme_ridges(center_axis_labels = T) + ylab('Clonal Complex')


pg <- plot_grid(p1, p2, labels = c('A', 'B'), nrow=1, ncol=2)

save_plot('figures/distance_distributions.pdf', pg, base_width = 12, base_height = 6)
