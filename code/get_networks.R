setwd('~/Github/DeCoTUR_manuscript_code')
source("code/helper_functions.R")

scores <- read.csv('data/topscores_withpval.csv', stringsAsFactors = F)
scores <- read.csv('data/lowlinkscores_withpval.csv', stringsAsFactors = F)
#scores <- read.csv('data/highlinkscores_withpval.csv', stringsAsFactors = F)

# used to get gene names for manually producing "types.csv"
# allgenes <- unique(c(scores_top$Trait1, scores_top$Trait2, scores_high$Trait1,
#                     scores_high$Trait2, scores_low$Trait1, scores_low$Trait2))

typemap <- read.csv('data/types.csv', stringsAsFactors = F)
tograph <- subset(scores, sig)
tograph$polarity <- tograph$PositiveAssociation > tograph$NegativeAssociation
tograph$polarity <- factor(tograph$polarity, levels= c('TRUE', 'FALSE'),
                           labels = c('Association', 'Dissociation'))
graph_data_frame <- tograph[, c(1, 2, 5, 12)]
g <- graph_from_data_frame(graph_data_frame, directed = F)
nv <- names(V(g))
types <- typemap$type[match(nv, typemap$name)]
V(g)$type <- types
manual_scale <- friendly_pal(name=  'muted_nine', n = 9)[c(8, 5, 1, 7, 4, 3)]
labelvec <- c(expression(paste('Antibiotic Resistance (non-SCC', italic('mec'), ')')),
              'Metal Resistance', 'Mobile Genetic Element', 'Other',
              expression(paste('SCC', italic('mec'))), 'Virulence')
gg <- ggraph(g, layout = 'fr') + geom_edge_link(aes(color = Score, linetype = polarity), width = 1) +
  scale_edge_color_gradient(name = 'Interaction', low = 'gray80', high = 'black') +
  scale_edge_linetype(name = 'Polarity') +
  geom_node_label(aes(label = name, fill = type), color = 'black', alpha = 0.7, repel = T) +
  geom_node_point(aes(color = type), size=3) +
  #geom_node_text(aes(label = name), color = 'black') +
  coord_equal() + theme(legend.text.align = 0) +
  guides(
    fill = guide_legend(
      override.aes = aes(label = "")
    )
  ) +
  scale_color_manual(name = 'Type', values = manual_scale, guide = 'none') +
  scale_fill_manual(name = 'Type', values = manual_scale, labels = labelvec) +
  theme_graph(base_family = 'Arial')
ggsave('figures/top_network.pdf', gg, scale = 1.4)
#ggsave('figures/lowlink_network.pdf', gg, scale = 1.4)
