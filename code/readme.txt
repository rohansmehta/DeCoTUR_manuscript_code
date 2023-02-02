%CC = insert a clonal complex name
%S = insert a number of species
%G = insert a number of genes
%CP = insert a number of close pairs

generate_sample_cc_map.R --> takes ST/CC map document from internet and maps public samples to their CCs to create sample_cc_map.csv
remove_redundant_samples.R --> takes nrd-gene-set.txt and sample_cc_map.csv and creates keptsamples.csv
create_distance_matrices.R --> takes nrd-gene-set.txt and keptsamples.csv and creates %CC_distmat.csv
panaroo_code.R --> takes keptsamples.csv and prokka_subsample.txt and obtains gene data, contains procedure to run prokka and panaroo
panaroo results are gene_presence_absence.csv, gene_presence_absence_roary.csv, final_graph.gml, and gene_presence_absence.Rtab
get_consolidated_gene_matrix --> takes gene_presence_absence.csv and creates consolidated_gene_matrix.csv
get_allpds.R --> takes all %CC_distmat.csv and creates all_pds.csv
get_scores_panaroo.R --> takes all_pds.csv and consolidated_gene_matrix.csv and creates closepairs.csv and panaroo_decotur.csv
get_graph_distances.R --> takes final_graph.gml and consolidated_gene_matrix.csv and creates graph_distances.csv and graph_distance_processed.csv
integrate_graph_distances.R --> takes panaroo_decotur.csv and graph_distances_processed.csv and creates panaroo_decotur_withdist.csv, spatial_plot.pdf, topscores.csv, lowlinkscores.csv, and highlinkscores.csv
compute_p_values.R --> takes topscores_withdist.csv, highlinkscores_withdist.csv, and lowlinkscores_withdist.csv and all_pds.csv and creates topscores_withpval.csv and lowlinkscores_withpval.csv and highlinkscores_withpval.csv
get_networks.R --> takes types.csv, topscores_withpval.csv and lowlinkscores_withpval.csv and creates top_network.pdf and lowlinknetwork.pdf
linkage_composition_plot.R --> takes lowlinkscores_withpval.csv, highlinkscores_withpval.csv and creates linkage_composition_plot.pdf, lowdist.csv, and highdist.csv (as intermediates for the plot)
generate_abr_scores_full.R --> takes all_pds.csv and creates abr_discordances.csv, abr_scores_unweighted_withnull.csv, abr_scores_withpval, and abr_heatplot.pdf. Note that this script uses an old version of decotur functions that are specifically defined in the script
create_runtime_comparison.R --> creates timecomp_loglin.pdf
comparison_data.R --> takes gene_presence_absence_roary.csv and creates comparison%S_%G_pa.csv, comparison%S_%G_pa.Rtab, comparison%S.fa, comparison%S_distmat.csv, and comparison%S_tree.newick
get_decotur_comparison_scores.r --> takes comparison_%S_distmat.csv, comparison%S_%G_pa.csv, and all_pds.csv, prints runtimes to be manually recorded in times.R, and creates decotur_%CP.csv.
comparison_methods.R --> takes comparison%S_%G_pa.Rtab and creates comparison%S_%G_pa_new.Rtab, takes decotur_%CP.csv, co_nodes.tsv, co_edges.tsv, cd_nodes.tsv, cd_edges.tsv, and gene_pa_spydrpick500_1000_noa.csv, and creates closepairsplot.pdf, coinfinder_example.csv, and methods_comparison.pdf
sample_sizes.R --> takes keptsamples.csv and sample_cc_map.csv and creates cc_sample_hist.pdf
generate_closepair_table.R --> takes keptsamples.csv, and %CC_distmat.csv and creates closepairtable.csv
distance_distributions.R --> takes %CC_distmat.csv and closepairtable.csv and creates distmat.csv
distance_distribution_figure.R --> takes distmat.csv and closepairtable.csv and creates distance_distributions.pdf
example_tree.R --> takes CC97_distmat.csv and keptsamples.csv and creates example_tree.pdf 
representationfig.R --> takes all_pds.csv and creates representation.pdf
