# Running dada2
#
# THIS IS A SIMPLE PIPELINE TO RUN DADA2 SEQUENCE ANALYSIS. 
#
# YOU WILL FIRST READ IN 2 MANDATORY FILES, AND 1 OPTIONAL
# YOU WILL THEN PROCEED TO RUN THE COMMANDS SEQUENTIALLY,
# CHANGING ONLY THE PARAMETERS WITH IN THE FUNCTIONS.
#
# THE CODE IS ALSO GIVEN BELOW.
# I STRONGLY RECOMMEND GOING OVER IT AND FAMILIARIZING YOURSELF WITH 
# FUNCTIONS IN THE DADA2 PACKAGE.
#
# # dataPath = MANDATORY*** YOU SHOULD BE ABLE TO CHANGE YOU WORKING DIRECTORY TO WHERE YOUR RAW DATA IS.
# # database = MANDATORY*** THEN CHOOSE A PATH TO THE DATA BASE TO CHOOSE TAXONOMY (SHOULD BE PROPERLY FORMATED FROM SOURCE)
# # DBSpec = OPTIONAL*** USE THIS IF YOU WANT TO MAKE 100% ID WITH BACTERIA
# #
#####
# Make Simple Quality Plots to determine where you should trim your Data
# RUN AS
# DADA2qualplots(dataPath = dataPath)
#
DADA2qualplots <- function(dataPath = dataPath) { 
  .cran_packages <- c("ggplot2", "gridExtra", "patchwork")
  .bioc_packages <- c("dada2", "msa", "phyloseq")
  # Load packages into session
  sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
  
  fnFs <- sort(list.files(dataPath, pattern = '_R1_001.fastq', full.names = T))
  fnRs <- sort(list.files(dataPath, pattern = '_R2_001.fastq', full.names = T))
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  # Plot Quality scores of a few samples
  print("Plotting the Quality Scores of a few samples")
  fqc <- plotQualityProfile(fnFs[1:6])
  rqc <- plotQualityProfile(fnRs[1:6])
  
  print("Inspect the plots and estimate the trimming length")
  assign("sample.names", sample.names, envir = globalenv())
  assign("fnFs", fnFs, envir = globalenv())
  assign("fnRs", fnRs, envir = globalenv())
  
  print("If all is good, please proceed with quality filtering and trimming with command `DADA2qualfilter()`")
  return(fqc / rqc)
}
#
# Filtering and Trimming with DADA2 
# DECIDE WHERE TO TRIM ON THE READS
# 
# RUN AS
# DADA2QfiltandTaxa(dataPath = dataPath, FwrdStrtTrim = 10, FwrdEndTrim = 220, RevStrtTrim = 10, RevEndTrim = 200)
#
DADA2QfiltandTaxa <- function(dataPath, 
                              FwrdStrtTrim = 10, RevStrtTrim = 10, 
                              FwrdEndTrim = 220, RevEndTrim = 180, 
                              errorRates = F, 
                              database = DB, 
                              databaseSpec = DBspec, 
                              asSpec = F) {
  # Quality filter
  filtFs <- file.path(dataPath, "filtered", 
                      paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(dataPath, "filtered", 
                      paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                       trimLeft = c(FwrdStrtTrim, RevStrtTrim), 
                       truncLen=c(FwrdEndTrim,RevEndTrim),
                       maxN=0, 
                       maxEE=c(2,2), 
                       truncQ=2, 
                       rm.phix=TRUE,
                       compress=TRUE, multithread=F) # On Windows set multithread=FALSE
  # remove the failed files from directory and start over....
  head(out)
  
  # learn error rates
  print("Learning Error RAtes")
  errF <- learnErrors(filtFs, multithread=F)
  errR <- learnErrors(filtRs, multithread=F)
  if (errorRates){
    plotErrors(errF, nominalQ = T)}
  
  else
    
  # apply the core inference algorithm
  dadaFs <- dada(filtFs, err = errF, multithread = F)
  dadaRs <- dada(filtRs, err = errR, multithread = F)
  
  # merge forward and reverse
  print("Merging forward and Reverse reads")
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
  # Inspect the merger data.frame from the first sample
  head(mergers[[1]])
  
  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  dim(seqtab)
  # Inspect distribution of sequence lengths
  table(nchar(getSequences(seqtab)))
  
  # Remove Chimeras
  print("Removing chimeras")
  seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                      method="consensus", 
                                      multithread=F, 
                                      verbose=TRUE) 
  
  assign("seqtab.nochim", seqtab.nochim, envir = globalenv())
  assign("seqtab", seqtab, envir = globalenv())
  
  print(dim(seqtab.nochim))
  print(sum(seqtab.nochim)/sum(seqtab))
  
  print("Assigning taxonomy to the sequences")
  taxa <- assignTaxonomy(seqtab.nochim, database, multithread=F)
  
  if (asSpec == T) {
    print("Hmm... I see you've decided to assign the species... You are a courageous one.")
    taxa <- assignSpecies(seqtab.nochim, databaseSpec)
  }
  else
    
  assign("taxa", taxa, envir = globalenv())
  
  print("If this is all good; please proceed to execute the `phylotree() command...")
}
#
# Building a (unrooted) phylogenetic tree from the DADA2 pipeline file
# https://f1000researchdata.s3.amazonaws.com/manuscripts/9666/d404b7c9-fa49-49ce-a991-df8becc13186_8986_-_susan_holmes.pdf?doi=10.12688/f1000research.8986.1&numberOfBrowsableCollections=17&numberOfBrowsableInstitutionalCollections=4&numberOfBrowsableGateways=22
#
# I STRONGLY RECOMMEND GOING OVER THE HELP SECTION FOR THE LIBRARY AND FUNCTIONS WITHIN THE FOLLOWING FUNCTION.
# RUN AS
# phylotree(seqtab.nochim = seqtab.nochim)
#
phylotree_NJ <- function(seqtab.nochim) {
  .cran_packages <- c("ggplot2")
  .bioc_packages <- c("dada2", "msa", "phyloseq")
  sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
  
  print("This method is generating a non-rooted Neighbor-Joining tree with ClustalW")
  seqs <- getSequences(seqtab.nochim)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  mult <- msa(seqs, method="ClustalW", type="dna", order="input")
  
  # construct a neighbor-joining tree, and then fit a GTR+G+I 
  # maximum likelihood tree using the neighbor-joining tree as a starting point.
  library(phangorn)
  print("Building Phylogenetic Tree")
  phang.align <- as.phyDat(mult, type="DNA", 
                           names = getSequence(seqtab.nochim))
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data=phang.align)
  
  fitGTR.nj <- update(fit, k=4, inv=0.2)
  fitGTR.nj <- optim.pml(fitGTR.nj, 
                         model="GTR", 
                         optInv=TRUE, 
                         optGamma=TRUE,
                         rearrangement = "stochastic", 
                         control = pml.control(trace = 0))
  detach("package:phangorn", unload=TRUE)
  print("Go on! You're almost there! Build your PhyloSeq object!")
  assign("fitGTR.nj", fitGTR.nj, envir = globalenv())
  
}

phylotree_UPGMA <- function(seqtab.nochim) {
  .cran_packages <- c("ggplot2")
  .bioc_packages <- c("dada2", "msa", "phyloseq")
  sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
  
  print("This method is generating a Rooted tree using the UPGMA method with ClustalW")
  seqs <- getSequences(seqtab.nochim)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  mult <- msa(seqs, method="ClustalW", type="dna", order="input")
  
  # construct a rooted UPGMA tree, and then fit a GTR+G+I 
  # maximum likelihood tree using the neighbor-joining tree as a starting point.
  library(phangorn)
  print("Building Phylogenetic Tree")
  phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab.nochim))
  dm <- dist.ml(phang.align)
  treeUPGMA <- upgma(dm) # Note, tip order != sequence order
  fit = pml(treeUPGMA, data=phang.align)
  
  fitGTR.upgma <- update(fit, k=4, inv=0.2)
  fitGTR.upgma <- optim.pml(fitGTR.upgma, 
                            model="GTR", 
                            optInv=TRUE, 
                            optGamma=TRUE,
                            rearrangement = "stochastic", 
                            control = pml.control(trace = 0), 
                            optRooted = T)
  detach("package:phangorn", unload=TRUE)
  print("Go on! You're almost there! Build your PhyloSeq object!")
  assign("fitGTR.upgma", fitGTR.upgma, envir = globalenv())
  
}

#
# RUN and NAME THE PHYLOSEQ OBJECT: THIS IS USED FOR DATA EXPLORATION AND ANALYSIS
# ARS.bact.97 <- buildPhyloseq(seqtab.nochim = seqtab.nochim, taxa = taxa, fitGTR = fitGTR)
#
# Merge the meta data file later
# RUN AS
# OBJECT <- buildPhyloseq(seqtab.nochim = seqtab.nochim, taxa = taxa, fitGTR = fitGTR)
buildPhyloseq <- function(seqtab.nochim,  taxa, fitGTR) {
  .cran_packages <- c("ggplot2", "dplyr", "scales", "grid", "reshape2", "RColorBrewer")
  .bioc_packages <- c("dada2", "msa", "phyloseq", "vegan","Biostrings")
  theme_set(theme_classic())

      phylo<- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F),
                   tax_table(taxa),
                   phy_tree(fitGTR$tree))
  
  
  # save sequences but change the names to ASV
  dna <- Biostrings::DNAStringSet(taxa_names(phylo))
  names(dna) <- taxa_names(phylo)
  phylo <- merge_phyloseq(phylo, dna)
  taxa_names(phylo) <- paste0("ASV", seq(ntaxa(phylo)))
  print(phylo)
  print("This is how many sequences there are in the data set:") 
  print(sum(colSums(otu_table(phylo))))
  return(phylo)
  
}

  