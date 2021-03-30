# setwd("")
# 
# tax_tab <- "C:/Users/smajor/Box Sync/Experiments/190501_Annisquam/ARS-bacteria/ARS-bact-97-taxa_table.csv"
# otu_tab <- "C:/Users/smajor/Box Sync/Experiments/190501_Annisquam/ARS-bacteria/ARS-bact-97-otu_table.csv"
# tree_file <- "C:/Users/smajor/Box Sync/Experiments/190501_Annisquam/ARS-bacteria/ARS-bact-97-tree.tre"
# seq_file <- "C:/Users/smajor/Box Sync/Experiments/190501_Annisquam/ARS-bacteria/ARS-bact-97-seqs.csv"

load_Phyles <- function(tax_tab, otu_tab, tree_file, seq_file){
  require(phyloseq)
  
taxa <- as.matrix(read.csv(tax_tab
                           , row.names = 1))
  taxa <- tax_table(taxa)
otu <- as.matrix(read.csv(otu_tab
                          , row.names = 1))
  otu <- otu_table(otu
                   , taxa_are_rows = F)
tree <- read_tree(tree_file)
seq <- as.matrix(read.csv(seq_file
                          , row.names = 1))
  seq <- Biostrings::DNAStringSet(seq)
  names(seq) <- taxa_names(taxa)

test <- phyloseq(otu,
                 taxa,
                 tree,
                 seq)
return(test)
}
