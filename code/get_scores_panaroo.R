setwd('~/Github/DeCoTUR_manuscript_code')
source("code/helper_functions.R")
library(decotur)
# This script gets the scores for the full dataset analysis. 
# First we need the close pairs
allpd <- read.csv('data/all_pds.csv', stringsAsFactors = F)
cutoff <- 0.00005
# total below cutoff
allpdne <- subset(allpd, Sample1 != Sample2)
numzer <- subset(allpdne, !is.na(dist)) %>% group_by(cc) %>% summarise(fz = sum(dist < cutoff), ld = length(dist))
# First we get all pairs that are not in the CCs with lots below cutoff
nonz <- subset(allpd, cc %not in% c('CC5', 'CC8', 'CC22') )
nonz <- subset(nonz, !is.na(dist))
nonz <- subset(nonz, Sample1 != Sample2)
cutq <- subset(nonz, dist <= cutoff)
cutdists <- subset(nonz, dist <= cutoff)
perrcp <- dim(cutdists)[1]/length(unique(cutdists$cc))
# then we randomly sample from those clonal complexes
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
cutoq <- rbind(cutq, otherpd)
# some of these samples are not in the set we used for the panaroo analysis
cgm <- get_matrix('data/consolidated_gene_matrix.csv', F)
# so we subset based on that
allcp <- subset(cutoq, Sample1 %in% colnames(cgm) & Sample2 %in% colnames(cgm))
# closepairs are identified by their index. We need the column order of cgm
geneorder <- colnames(cgm)
rs1 <- match(allcp$Sample1, geneorder)
rs2 <- match(allcp$Sample2, geneorder)
closepairsc <- cbind(rs1, rs2)
set.seed(1)
closepairs <- closepairsc[sample(1:dim(closepairsc)[1], 10000, replace = F),]
# 10,000 random close pairs were chosen from these 187311 original ones
write.csv(closepairs, 'data/closepairs.csv', row.names = F)
closepairs <- read.csv('data/closepairs.csv', header=T)
gmf <- filter_snps_by_closepairs(cgm, closepairs)
gene_matrix_filtered <- gmf[[1]]
closepairs <- gmf[[2]]
classes <- decotur::get_closepair_classes(closepairs)
classweights <- decotur::get_classweights(classes)
gmf1 <- gene_matrix_filtered[which(rowSums(gene_matrix_filtered) > 500),]
t <- proc.time()
scores <- get_scores_pa_closepairs(gmf1, closepairs, classweights, 200, T, 'speed', T)
write.csv(scores, 'data/panaroo_decotur.csv', row.names = F)
