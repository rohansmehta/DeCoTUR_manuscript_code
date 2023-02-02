setwd('~/GitHub/DeCoTUR_manuscript_code')
source("code/helper_functions.R")

allpd <- read.csv('data/allpds.csv', stringsAsFactors = F)
cutoff <- 0.00005
# what are the samples represented?
nonz <- subset(allpd, cc %not in% c('CC5', 'CC8', 'CC22') )
nonz <- subset(nonz, !is.na(dist))
nonz <- subset(nonz, Sample1 != Sample2)
cutq <- subset(nonz, dist <= cutoff)

quantiles <- quantile(nonz$dist, seq(0.01, 1, 0.01), na.rm = T)
full <- table(nonz$cc)
fullfrac <- full/sum(full)
allq <- c()
for(q in quantiles){
  closepairs <- subset(nonz, dist <= q)
  qtab <- table(closepairs$cc)
  qfrac <- qtab/sum(qtab)
  qtot <- qfrac/fullfrac
  mq <- melt(qtot)
  mq$cutoff <- q
  allq <- rbind(allq, mq)
  print(q)
}
write.csv(allq, 'data/ccfrac_quantiles.csv', row.names = F)

allq <- read.csv('data/ccfrac_quantiles.csv', stringsAsFactors = F)

repres <- ggplot(allq, aes(x = cutoff, y = value, group = Var1, color = Var1)) + geom_point() + geom_line(size=1) +
  geom_hline(yintercept = 1, linetype = 'dashed') + theme_bw() +
  scale_color_manual(values = friendly_pal("bright_seven"), name = 'Clonal Complex') + scale_x_log10() +
  xlab('Close-Pair Cutoff') + ylab('K-Fold Representation in Sample') + geom_vline(xintercept = cutoff, linetype = 'dashed')
save_plot(filename = 'figures/representation.pdf', plot = repres, base_width = 8, base_height = 6)