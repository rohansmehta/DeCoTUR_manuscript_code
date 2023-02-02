setwd(~/Github/DeCoTUR_manuscript_code)
source("code/helper_functions.R")

ks <- read.csv('data/keptsamples.csv')
sc <- read.csv('data/sample_cc_map.csv')

ks_table <- table(ks$cc)
sc_table <- table(sc$clonal_complex)
cc_table <- data.frame(sc_table, ks_table)
cc_table <- cc_table[, c(1, 2, 4)]
names(cc_table) <- c('cc', 'Full', 'Distinct')
ccm <- melt(cc_table, id = c('cc'))
p <- ggplot(ccm, aes(x = reorder(cc, -value, 'max'), y = value, group = variable, fill = variable)) + 
  geom_col(position=position_dodge2()) + xlab('Clonal Complex') + ylab('Number of Samples') + 
  theme_cowplot() + scale_fill_manual(values = c("#35274A", "#F2300F"), name = 'Sample', labels = c('Full', "Distinct"))

save_plot(filename = 'figures/cc_sample_hist.pdf', p, base_width = 8 , base_height = 6)

