setwd(~/Github/DeCoTUR_manuscript_code)
# different methods results comparison

for(sp in c(100, 500, 1000)){
  for(gn in c(1000, 3000, 5000, 10000)){
    spytab <- read.table(paste0('data/comparison', sp, '_', gn, '_pa.Rtab'), header=T, stringsAsFactors = F, sep='\t')
    names(spytab) <- str_remove_all(names(spytab), 'X')
    pa_matrix <- pbapply(spytab[,2:(sp+1)], 2, function(x)(return(as.numeric((x != '')))))
    pa_matrix <- cbind(Gene = spytab$Gene, pa_matrix)
    write.table(pa_matrix, paste0('data/comparison', sp, '_', gn, '_pa_new.Rtab'), 
                row.names = F, sep = '\t', quote = F)
  }
}

#let's do 500/1000 with decotur and coinfinder and get spydrpick to work at some point
# let's do the cp-effect plot first.
data100 <- read.csv('data/decotur_100.csv')
data500 <- read.csv('data/decotur_500.csv')
data1000 <- read.csv('data/decotur_1000.csv')
data100$sign <- (-1)^(data100$PositiveAssociation < data100$NegativeAssociation)
data500$sign <- (-1)^(data500$PositiveAssociation < data500$NegativeAssociation)
data1000$sign <- (-1)^(data1000$PositiveAssociation < data1000$NegativeAssociation)
data100$revpair <- paste0(data100$Trait2, data100$Trait1)
data500$revpair <- paste0(data500$Trait2, data500$Trait1)
data1000$revpair <- paste0(data1000$Trait2, data1000$Trait1)
data100$score <- data100$Score * data100$sign
data500$score <- data500$Score * data500$sign
data1000$score <- data1000$Score * data1000$sign
dat1 <- merge(data1000[,c(8, 11, 9, 12)], data500[, c(8, 11, 9, 12)], by = 'pair')
dat2 <- merge(dat1, data100[, c(8, 11, 9, 12)], by = 'pair')
# ok here we go
#write.csv(dat2,  'DECOTUR_COMPARISON.csv')
#dat2 <- read.csv('DECOTUR_COMPARISON.csv', stringsAsFactors = F)
set.seed(100)
ss1 <- subset(dat2, sig.x == T & sig.y == T & sig == T)
ss <- ss1[sample(1:dim(ss1)[1], 10000),]
md <- melt(ss[, c(1, 4, 7, 10)], id = 'pair')
md$x <- 1000*(md$variable == 'score.x') + 500 * (md$variable == 'score.y') + 100 * (md$variable == 'score')
p <- ggplot(md, aes(x = x, y = value, group = pair)) + geom_point() + geom_line(alpha = 0.3) + 
  xlab('Number of Close Pairs') + ylab('Coevolution Score')  + theme_bw() + scale_x_continuous(breaks = c(100, 500, 1000))
ggsave('figures/closepairsplot.pdf', p)

coindata <- read.table('data/co_nodes.tsv', sep = '\t', header=T)
coindatae <- read.table('data/co_edges.tsv', sep = '\t', header=T)
coindatd <- read.table('data/cd_nodes.tsv', sep = '\t', header=T)
coindatde <- read.table('data/cd_edges.tsv', sep = '\t', header=T)
nodedat <- merge(coindata, coindatd, by = 'ID')
nodedat$mean <- (nodedat$Result.x+nodedat$Result.y)/2
coindatde$weight <- -coindatde$weight
coin <- rbind(coindatae, coindatde)
coin$pair <- paste0(coin$Source, coin$Target)
coin1 <- coin %>% group_by(pair) %>% summarise(score = (-1)^(sum(weight) < 0) * pmax(abs(weight)))
write.csv(coin1, 'data/coinfinder_example.csv', row.names = F)
coin <- read.csv('data/coinfinder_example.csv', stringsAsFactors = F)
dec <- dat2[, c(1, 2, 6, 7)]
dec <- subset(dec, sig.y)
dec <- dec[, c(1,2, 4)]
names(dec) <- c('pair', 'revpair', 'decotur')
dc <- merge(dec, coin, by = 'pair')
names(dc) <- c('pair', 'revpair', 'decotur', 'coinfinder')
names(dec) <- c('revpair', 'pair', 'decotur')
dc1 <-  merge(dec, coin, by = 'pair')
names(dc1) <- c('pair', 'revpair', 'decotur', 'coinfinder')
dc2 <- rbind(dc, dc1)
library(scales)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
dc2$coinfindernlp <- dc2$coinfinder
dc2$coinfindernlp[which(dc2$coinfinder > 0)] <- -log10(dc2$coinfinder[which(dc2$coinfinder > 0)])
dc2$coinfindernlp[which(dc2$coinfinder < 0)] <- log10(-dc2$coinfinder[which(dc2$coinfinder < 0)])
dcp <- ggplot(dc2, aes(x = decotur, y = coinfindernlp)) + geom_point(alpha = 0.1) + theme_bw() + 
  xlab('Coevolution Score (DeCoTUR)') + 
  ylab('Negative Log P-value (Coinfinder)') + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
#ggsave('dcplot.pdf',dcp)  
spy <- read.csv('data/gene_pa_spydrpick500_1000_noa.csv')
spy$pair <- paste0(spy$GeneA, spy$GeneB)
spy1 <- spy[, c(5, 3)]
names(dec) <- c('pair', 'revpair', 'decotur')
ds <- merge(dec, spy1, by = 'pair')
names(dec) <- c('revpair', 'pair', 'decotur')
ds1 <-  merge(dec, spy1, by = 'pair')
ds2 <- rbind(ds, ds1)
dsp <- ggplot(ds2, aes(x = decotur, y = MI)) + geom_point(alpha = 0.1) + theme_bw() + 
  xlab('Coevolution Score (DeCoTUR)') + 
  ylab('Mutual Information (SpydrPick)')
#ggsave('dsplot.pdf',dsp)  
ds2$MIs <- ds2$MI
ds2$MIs[which(ds2$decotur < 0)] <- -ds2$MIs[which(ds2$decotur < 0)] 
dspalt <- ggplot(ds2, aes(x = decotur, y = MIs)) + geom_point(alpha = 0.1) + theme_bw() + 
  xlab('Coevolution Score (DeCoTUR)') + 
  ylab('Mutual Information (SpydrPick)') + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
#ggsave('dsplotalt.pdf',dspalt) 

library(cowplot)
both <- plot_grid(dcp, dspalt, labels = c('A', 'B'))
save_plot('figures/methods_comparison.pdf', both, base_asp = 2)
