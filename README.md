# GMGI_Microbe_Repository
Identify 16S rRNA gene fragments from raw MiSeq data (fastq)
A script has been written to identify the cultures that were successful in amplifying the 16S rRNA gene. This script/procedure was executed to identify the most likely identity at the genus level (97%; unless species is chosen).

The purpose of this pipeline is to streamline culture classification suing MiSeq data of a fragment of the 16S rRNA gene sequence or the Fungal ITS.

The appropriate databases for taxanomic classification should be downloaded prior to running the command.  It is also suggested that you acquire the most up-to-date version of R and DADA2 which also requires Bioconductor. Attempts to run this on the server with R 3.5.1 and DADA2 1.10 failed (2019-09-03)..

Attached are a a couple of R scripts ("DADA-2 Pipeline.R", "physeq_Object_loading.R", and "SampleID-6.R")  that, when run in sequence (see procedure) will...

Make a simple quality score plot of the first 6 samples (forward and reverse)
Quality filter, trim, learn error rates, merge, and assign the taxonomy.
Build a phylogenetic tree
Build a phyloseq object (R program for data analysis)
Identify all the taxonomy and provide their % identity.
Identify the most commonly found taxonomic classification, how many time this taxonomy was counted, the Representative sequence of the identified taxonomy, the most abundant sequence in each sample, the frequency of this most abundant sequence, and actual sequence. 
All scripts are customizable at the input level. It is highly suggested that you look at the code and identify the function of each command and determine if it needs customization/re-configuration for your dataset.

Required packages.... 

DADA2
Bioconductor
msa
phyloseq
phangorn
Biostrings
vegan
RColorBrewer
reshape2
grid
scales
gridExtra
dplyr
tidyr
patchwork
ggplot2

Databases that have been curated for taxonomic classification by DADA2 are available at the follow website https://benjjneb.github.io/dada2/training.html. Zipped files of the databases are also available in the subdirectory "16S/ITS Reference Databases". These are the most up-to-date databases formatted for the use with DADA2 as of 8/26/2020. Always check if a database has been updated.

Below is an example pipeline of how the pipeline currently works.

This works for both 16S and ITS. However, for fungal ITS identities, the workflow needs to be fiddled with due to the greater variation in size of the ITS sequence, i.e: a one-size-fits-all approach on trimming may be detrimental to the data.

Because our sequences are typically very short (<300 bp) compared to the whole 16S gene (~1600 bp), I advise against using species identification for bacterial cultures.

Please let SRM know if there are any issues if/when this is performed.

The final output is meant to be merged with the "Repository Sequencing Data". The output csv is not perfect and needs to be edited. However, it is effective at identifying the cultures and giving basic statistics. Edits I would like to include:

Remove the "TAXA_TOTAL" field, it is a duplicate of "SAMPLE_TOTAL". 
Add the date the cultures were identified/the code was executed.
Add the database used to identify each culture.

# Running dada2

# THIS IS A SIMPLE PIPELINE TO RUN DADA2 SEQUENCE ANALYSIS. 

# YOU WILL FIRST READ IN 2 MANDATORY FILES, AND 1 OPTIONAL

dataPath = MANDATORY **PATH TO THE RAW DATA**
database = MANDATORY **THE PATH TO THE DATA BASE TO CHOOSE TAXONOMY (SHOULD BE PROPERLY FORMATTED FOR DADA2)**
DBSpec = OPTIONAL **USE THIS IF YOU WANT TO MAKE 100% ID WITH BACTERIA (ANOTHER DATABASE)**

# YOU WILL THEN PROCEED TO RUN THE COMMANDS SEQUENTIALLY, CHANGING ONLY THE HIGHLIGHTED PARAMETERS WITHIN THE FUNCTIONS.

# Set your working directory, identify the location of the raw MiSeq data, and choose the database to reference that data.

'<addr>' library(dada2)
'<addr>' setwd("C:/Users/smajor/Box Sync/Experiments/bacterial_cultures")
'<addr>' dataPath <- "C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/rawdata/"
'<addr>' database <- "C:/databases/rdp_train_set_16.fa.gz"
'<addr>' DBspec <- "C:/databases/rdp_species_assignment_16.fa.gz" # USE THIS FOR 100% ID
# Source the location of the pipeline file.

'<addr>' source("C:/Users/smajor/Box Sync/R-scripts/DADA-2 Pipeline.R")
#####
# Make Simple Quality Plots to determine where you should trim your Data

DADA2qualplots(dataPath = dataPath)
# Decide where to Trim based on these plots. These are very basic plots. Use the program FastQC for more robust QC analysis on the sequences.

#####
# Filter, Trim, Learn the Error rates, Merge, remove chimeras, and assign taxonomy with DADA2

DADA2QfiltandTaxa(dataPath = dataPath
                  , FwrdStrtTrim = 10, FwrdEndTrim = 225
                  , RevStrtTrim = 10, RevEndTrim = 225
                  , errorRates = F, database = database) # ADD database = DBspec to add 100% identification to species
# If any samples are completely removed during filtering and trimming, you will have to remove those samples that failed from the dataPath location and run this code again.

#####
# Build the phylogenetic tree

# Non-rooted neighbor-joining tree with ClustalW

phylotree_NJ(seqtab.nochim = seqtab.nochim)
OR

# Rooted UPGMA method with ClustalW. This is best to use to see consistent relationships; otherwise, non-rooted trees will have the root randomly chosen when doing phylogenetic analyses with phyloseq.

phylotree_UPGMA(seqtab.nochim = seqtab.nochim)
#####

# Build the phyloseq object based on the output from the previous commands.
# The phyloseq object is a convenient way to explore the data in R.

bacterial.cultures <- buildPhyloseq(seqtab.nochim = seqtab.nochim, taxa = taxa, fitGTR = fitGTR.upgma)
# OPTIONAL: Save the individual portions of the phyloseq object to maintain reproducibility

taxa <- as.data.frame(tax_table(bacterial.cultures))
otu <- as.data.frame(otu_table(bacterial.cultures))
tree <- phy_tree(bacterial.cultures)
seqs <- as.data.frame(refseq(bacterial.cultures))

write.csv(taxa, "C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-taxa.csv")
write.csv(otu, "C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-otu.csv")
write.tree(tree, "C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-tree.tre")
write.csv(seqs,"C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-seqs.csv")
# Better to write the sequences as a fasta
Biostrings::writeXStringSet(refseq(bacterial.cultures),"C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-seq.fasta",  format = "fasta") 

# Load the files back into R to create the phyloseq object, if required
source("C:/Users/smajor/Box Sync/Science/SOPs and Protocols/Microbe/CODES for MiSeq Data Analysis/phySeq_Object_loading.R")

taxa <- "C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-taxa.csv"
otu <- "C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-otu.csv"
tree <- "C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-tree.tre"
seq <- "C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-seq.fasta""

bacterial.cultures <- load_Phyles(taxa, otu, tree, seq)
#####

# Identify each individual culture down to the genus level and provide some basic statistics on the representative sequence as well as the most abundant sequence. This is the command that ultimately generates the data fro the "Repository Sequencing Data"

source("C:/Users/smajor/Box Sync/Science/SOPs and Protocols/Microbe/CODES for MiSeq Data Analysis/SampleID-6.R")
bacterial.cultures.genus <- cultureID.2(bacterial.cultures, "Genus")
# Write the object to a csv file.

write.csv(file = "C:/Users/smajor/Box Sync/Experiments/bacterial_cultures/bacterial_cultures-genus.csv", bacterial.cultures.genus)
# Merge the data with the "Repository Sequencing Data " spreadsheet.
 
