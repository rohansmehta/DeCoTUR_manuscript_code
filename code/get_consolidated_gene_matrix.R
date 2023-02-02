setwd('~/Github/DeCoTUR_manuscript_code')
library(pbapply)
library(stringi)
library(stringr)
# This script computes scores from the panaroo data, using the close pairs
# determined from a full dataset analysis

# First, let us acquire the panaroo data
panaroo <- read.csv('data/gene_presence_absence.csv', stringsAsFactors = F)
pa_matrix <- pbapply(panaroo[,4:10001], 2, function(x)(return(as.numeric((x != '')))))
rn <- panaroo$Non.unique.Gene.name
rn[which(rn == '')] <- panaroo$Gene[which(rn == '')]
rownames(pa_matrix) <- rn
colnames(pa_matrix) <- str_remove(colnames(pa_matrix), 'PROKKA_') # gets rid of "PROKKA_"
colnames(pa_matrix) <- str_remove(colnames(pa_matrix), '.fa') # gets rid of ".fa"
# we want to consolidate annotations
pa1 <- pa_matrix[!str_detect(rownames(pa_matrix), 'group'),]
rn <- rownames(pa1)
rnl <- pblapply(rn, function(x){
  x1 <- str_split_fixed(x, ';', Inf)
  x2 <- unique(stri_replace(x1, regex = "\\_.*", replacement = ""))
  x3 <- x2[which(x2 != '')]
  if(x3[1] == 'group'){
    return(x1[1])
  } else{
    return(paste0(sort(x3), collapse = '_'))
  }
})
newmat <- c()
donotuse <- c()
rownames1 <- c()
for(i in 1:length(rnl)){
  if(i %not in% donotuse){
    eq <- which(rnl[[i]] == rnl)
    donotuse <- c(donotuse, eq)
    rownames1 <- c(rownames1, rnl[[i]])
        submat <- pa1[eq,]
        if(length(eq) > 1){
          cols <- colSums(submat)
          pres <- as.numeric(cols > 0)
        } else{
          pres <- submat
        }
        newmat <- rbind(newmat, pres)
  }
  if(i%%100 == 0){
    print(i/length(rnl)) 
  }
}
# we have to then combine this with the "group" part of pa_matrix
rownames(newmat) <- rownames1
#write.csv(newmat, 'data/intermediate_matrix.csv')
pa2 <- pa_matrix[str_detect(rownames(pa_matrix), 'group'),]
newmat1 <- rbind(newmat,  pa2)
newmat1[which(is.na(newmat1))] <- 0
write.csv(newmat1, 'data/consolidated_gene_matrix.csv', row.names = T)
