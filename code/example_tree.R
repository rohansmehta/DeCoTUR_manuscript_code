setwd(~/Github/DeCoTUR_manuscript_code)
source("code/helper_functions.R")

cc1 <- 'CC97'
ks <- read.csv('data/keptsamples.csv')
kss <- subset(ks, cc == cc1)$sample_id
distance_matrix <- get_matrix(paste0('data/', cc1, '_distmat.csv'), F)
kept_distance_matrix <- distance_matrix[which(rownames(distance_matrix) %in% kss), which(colnames(distance_matrix) %in% kss)]
distance_matrix <- kept_distance_matrix
# Now we make a NJ tree from this distance matrix
tree <- fastme.bal(as.matrix(distance_matrix))
# what's the distance cutoff?
number_close_pairs <- 5000
number_distances <- dim(distance_matrix)[1]*(dim(distance_matrix)[2]-1)/2
fraction_distances <- number_close_pairs/number_distances
distance_cutoff <- quantile(x = as.numeric(as.matrix(distance_matrix)), 
                            probs = fraction_distances, na.rm = T)
closepairs <- get_closepairs_distance(distance_matrix, distance_cutoff, F, F)
renamed_closepairs <- cbind(rownames(distance_matrix)[closepairs[,1]], rownames(distance_matrix)[closepairs[,2]])
rel <- graph_from_data_frame(renamed_closepairs, directed = F)
classes <- components(rel)
memberships <- classes$membership
names(memberships) <- names(V(rel))
# Let's assign tip labels to closepair classes
tipdat <- data.frame(mem = as.factor(memberships))

colorscale <- c('#72e5ef', '#c82565', '#53f259', '#ba5dd3', '#add465', '#672396', 
                '#699023', '#fe16f4', '#21a708', '#ee0d0e', '#21f0b6', '#823d44', 
                '#239eb3', '#d17778', '#075c62', '#ffc4de', '#270fe2', '#e9c338', 
                '#3d4e92', '#fa7922', '#658bfb', '#957206', '#918ea2', '#115205')
full_colorscale <- rep(colorscale, 10)

p <- ggtree(tree, layout = 'circular')
q <- gheatmap(p, tipdat, width = 0.1, colnames = F) + 
  scale_fill_manual(values = full_colorscale, na.value = 'gray90') + theme(legend.position = 'none')

save_plot('figures/example_tree.pdf', q, base_width = 10, base_height = 10)
