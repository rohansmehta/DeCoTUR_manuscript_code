setwd("~/GitHub/DeCoTUR_manuscript_code")
library(data.table)
library(igraph)
library(reshape2)
library(ggplot2)
library(pbapply)
library(stringr)
typemap <- read.csv('data/types.csv', stringsAsFactors = F)
scores_low <- read.csv('data/lowlinkscores_withpval.csv', stringsAsFactors = F) 
scores_sig <- subset(scores_low, sig)
genetab <- table(c(scores_sig$Trait1, scores_sig$Trait2))
types <- typemap$type[match(names(genetab), typemap$name)]
typedat <- data.frame(genetab)
typedat$type <- types
names(typedat) <- c('Name', 'Count', 'Type')
write.csv(typedat, 'data/highdist.csv')

scores_high <- read.csv('data/highlinkscores_withpval.csv', stringsAsFactors = F) 
scores_sig <- subset(scores_high, sig)
genetab <- table(c(scores_sig$Trait1, scores_sig$Trait2))
types <- typemap$type[match(names(genetab), typemap$name)]
typedat <- data.frame(genetab)
typedat$type <- types
names(typedat) <- c('Name', 'Count', 'Type')
write.csv(typedat, 'data/lowdist.csv')

highdist <- read.csv('data/highdist.csv', stringsAsFactors = F)
lowdist <- read.csv('data/lowdist.csv', stringsAsFactors = F)
highdist$dist <- 'high'
highdist$frac <- highdist$Count/sum(highdist$Count)
lowdist$frac <- lowdist$Count/sum(lowdist$Count)
lowdist$dist <- 'low'
bardat <- rbind(highdist, lowdist)
bardat <- bardat[order(bardat$dist, bardat$Type, -bardat$frac),]
bardat$ID <- c(1:dim(highdist)[1], 1:dim(lowdist)[1])
manual_scale <- friendly_pal(name=  'muted_nine', n = 9)[c(8, 5, 1, 7, 4, 3)]
labelvec <- c(expression(paste('Antibiotic Resistance (non-SCC', italic('mec'), ')')),
              'Metal Resistance', 'Mobile Genetic Element', 'Other',
              expression(paste('SCC', italic('mec'))), 'Virulence')

bardat <- bardat %>%
  group_by(dist, Type) %>%
  mutate(typesum = sum(frac))
bardat <- bardat %>%
  group_by(dist, Type) %>%
  mutate(cumtypesum = cumsum(frac))
typesums <- bardat %>% group_by(dist) %>% summarise(ut= unique(typesum))
ht <- subset(typesums, dist == 'high')

h <- subset(bardat, dist == 'high')
h$label_y <- rev(cumsum(rev(h$frac))) - 0.5*h$frac
l <- subset(bardat, dist == 'low')
l$label_y <- rev(cumsum(rev(l$frac))) - 0.5*l$frac

bardat <- rbind(h, l)
bardat$label <- NA
whichtolab <- c(1, 3, 5, 6, 7, 10, 11, 12, 13, 23, 24, 25, 27, 28, 29, 30, 31, 
                47, 48, 49, 51, 52, 53, 55, 56, 57, 58, 60, 61, 62, 63, 64, 
                70, 71, 72, 73, 74, 75, 76)
bardat$label[whichtolab] <- bardat$Name[whichtolab]
bardat$label <- str_replace_all(bardat$label, '_', '/')

p <- ggplot(bardat, aes(x = factor(dist), y = frac, group = ID, fill = Type)) + geom_bar(stat = 'identity', color = 'black', size = 0.1 ) + 
  scale_fill_manual(values = manual_scale, labels = labelvec) + geom_text(aes(label = label, y = label_y, color = Type != 'mge')) + 
  xlab('Linkage') + ylab('Fraction') + scale_x_discrete(labels = c('Low', 'High')) + theme_bw() + 
  scale_color_manual(values = c('white', 'black'), guide = NULL)
p

ggsave2(filename = 'figures/linkage_composition.pdf', p, scale = 1.5)
