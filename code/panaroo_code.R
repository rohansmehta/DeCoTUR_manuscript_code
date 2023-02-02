# This code outlines the procedures to get the panaroo results.
setwd(~/Github/DeCoTUR_manuscript_code)
source("~/GitHub/figure_code/code/helper_functions.R")
library(decotur)
library(staphopia)
library(ggtext)
library(stringr)
library(pbapply)
TOKEN = '683eda5e910b094b15fd5276d0019f03b38200df'

# Generates fasta files for all contigs in all samples in our full-dataset analysis
ks <- read.csv('data/keptsamples.csv', stringsAsFactors = F)
samples <- ks$sample_id
blocksize <- 1000
numblocks <- length(samples) %/% blocksize
leftovers <- length(samples) %% blocksize
for(b in 1:numblocks){
  #t <- proc.time()
  contigs <- get_assembly_contigs(samples[((b-1)*blocksize + 1):(b*blocksize)], exclude_sequence = F)
  for(sam in unique(contigs$sample_id)){
    contigseqs <- as.list(subset(contigs, sample_id == sam)$sequence)
    contignames <- paste0(contigs$sample_id, contigs$contig)
    write.fasta(contigseqs, names = contignames, file.out = paste0('fullfa/', sam, '.fa'), as.string = T)   
    #print(sam)
  }
  #print(proc.time() - t)
}
# these files were compiled and are found as fullfa.zip at https://doi.org/10.5061/dryad.2v6wwpzq2

# now we run prokka using the following command

# for k in *.fa; do prokka $k --outdir "$k".prokka.output --prefix PROKKA_$k; echo $k; done

# all resulting gffs were compiled and are found as prokka_gffs.zip at https://doi.org/10.5061/dryad.2v6wwpzq2

# transfer the resulting gff files from all filenames found in prokka_subsample.txt, which was
# obtained using the R function sample into a folder called "gffs", and then run the following line
# for panaroo

# panaroo -i *.gff -o results -t 32 --clean-mode sensitive

# note that clean mode "sensitive" was used in order to maintain forward-compatibility by
# merging with another subset of the prokka files if necessary

# the resulting presence_absence matrix is the input into get_consolidated_gene_matrix.R
