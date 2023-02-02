# User needs to input their own Staphopia %TOKEN
library(staphopia)   
TOKEN = %TOKEN

ps <- get_public_samples()
cc_map <- read.table('data/bigsdb.txt', header=T, sep='\t')
ps_map <- merge(ps, cc_map, by = "st")
ps_map <- ps_map[,c(1, 2, 3, 14)]
# OK. Let's figure out why there are some samples missing from this
missing_st <- ps$st[which(ps$sample_id %not in% unique(ps_map$sample_id))]
# These are 612 samples with sequence type 0, which is not a sequence type.
# Let's add these guys in "Other"
ps_map_0 <- subset(ps, st == 0)
psmap0 <- ps_map_0[,c(5, 1, 2)]
psmap0$clonal_complex <- 'Other'
ps_map$clonal_complex[which(ps_map$clonal_complex == '')] <- 'Other'
ps_map <- rbind(ps_map, psmap0)
write.csv(x = ps_map, file = 'data/sample_cc_map.csv', row.names = F)