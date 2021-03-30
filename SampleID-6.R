# Determining which cultures are mixed and which are not
# I would like to take all the ASV for each sample and name 
# them in a table to determine the purity of each culture

# Use a phyloseq object 

# The output will
# # generate a table with the most identified 
# # taxa along with the most abundant ASV, it's taxonomy, the most abundant ASV sequence, and a representative sequence
#
# SP.bact.ID.97 <- cultureID.2(SP.bact.97, "Genus")
#
# # To find the most abundant OTU
# #
# write.csv(file = "C:/Users/smajor/Box Sync/Experiments/190606_Sponge/Sponge_bacteria/SP-bact-ID-97-genus.csv", SP.bact.ID.97)
#
# 
#
cultureID.2 <- function(physeq, TaxRank){
  
  require(tidyr)
  require(dplyr)
  rank <- tax_glom(physeq, taxrank = TaxRank, NArm = F)
  
  k <- tax_table(rank)
  l <- as.data.frame(otu_table(rank))
  q <- l
  q$TAXA_TOTAL <- rowSums(q)
  q$CULT_NAMES <- rownames(q)
  q <- q[,c(c(ncol(q),(c(ncol(q)-1))),1:(ncol(q)-2))]
  e <- gather(q, key = TAXA_ASV, TAXA_COUNT, 3:(ncol(q)), factor_key  = T)
  dat_max <- e %>% group_by(CULT_NAMES) %>% filter(row_number(TAXA_COUNT) == n())
  dat_max$TAXA_PERC_ID <- round(100*(dat_max$TAXA_COUNT/dat_max$TAXA_TOTAL))
  
  d <- as.data.frame(refseq(physeq))
  d$TAXA_ASV <- rownames(d)
  
  f <- merge(dat_max, d)
  
  g <- f[,c(2,3,4,5,1,6)]
  colnames(g)[6] <- "TAXA_REP_SEQ"
  rownames(g) <-g$CULT_NAMES
  
  b <- as.data.frame(otu_table(physeq))
  j <- b
  j$SAMPLE_TOTAL <- rowSums(j)
  j$CULT_NAMES <- rownames(j)
  
  j <- j[,c(c(ncol(j),(c(ncol(j)-1))),1:(ncol(j)-2))]
  
  w <- gather(j, ABUN_ASV, ABUN_ASV_COUNT, 3:(ncol(j)), factor_key  = T)
  
  max_asv <- w %>% group_by(CULT_NAMES) %>% filter(row_number(ABUN_ASV_COUNT) == n())
  max_asv$ABUN_ASV_PERC_ID <- round(100*(max_asv$ABUN_ASV_COUNT/max_asv$SAMPLE_TOTAL))
  
  d <- as.data.frame(refseq(physeq))
  d$ABUN_ASV <- rownames(d)
  
  v <- merge(max_asv, d)
  
  # g <- f[,c(2,3,4,5,1,6)]
  colnames(v)[6] <- "ABUN_ASV_SEQ"
  rownames(v) <-v$CULT_NAMES
  
  r <- merge(g, v)
  r$CULTURE_QUALITY <- ifelse(r$TAXA_PERC_ID >= 90, "Pure", "Mixed")
  r$CULTURE_QUALITY <- replace_na(r$CULTURE_QUALITY, "Unknown")
  
  tax.tab <- as.data.frame(tax_table(physeq))
  tax.tab$TAXA_ASV <- rownames(tax.tab)
  
  s <- merge(tax.tab, r
             , by = "TAXA_ASV", all.x = T)
  
  t <- s[!is.na(s$CULT_NAMES),]
  t$ID_LEVEL <- TaxRank
  t$Kingdom <- gsub("[[:lower:]]__", "", t$Kingdom)
  t$Phylum <- gsub("[[:lower:]]__", "", t$Phylum)
  t$Class <- gsub("[[:lower:]]__", "", t$Class)
  t$Order <- gsub("[[:lower:]]__", "", t$Order)
  t$Family <- gsub("[[:lower:]]__", "", t$Family)
  t$Genus <- gsub("[[:lower:]]__", "", t$Genus)
  t$Kingdom[is.na(t$Kingdom)] <- "Unknown"
  t$Phylum[is.na(t$Phylum)] <- "Unknown"
  t$Class[is.na(t$Class)] <- "Unknown"
  t$Order[is.na(t$Order)] <- "Unknown"
  t$Family[is.na(t$Family)] <- "Unknown"
  t$Genus[is.na(t$Genus)] <- "Unknown"
  
  names(t) <- sapply(names(t), function(v) {
    if (is.character(v)) 
      return(toupper(v))
    else return(v)})
  if ("SPECIES" %in% names(t) == F){
    t$SPECIES <- NA
    t <- t[c( "CULT_NAMES","CULTURE_QUALITY","ID_LEVEL","KINGDOM","PHYLUM","CLASS",
                          "ORDER","FAMILY","GENUS","SPECIES","SAMPLE_TOTAL","TAXA_TOTAL",
                          "TAXA_COUNT","TAXA_PERC_ID","TAXA_REP_SEQ","TAXA_ASV","ABUN_ASV_COUNT",
                          "ABUN_ASV_PERC_ID","ABUN_ASV_SEQ","ABUN_ASV")]

    return(t)
  }
  else
    t$SPECIES <- gsub("[[:lower:]]__", "", t$SPECIES)
    t$SPECIES[is.na(t$SPECIES)] <- "Unknown"
    t <- t[c( "CULT_NAMES","CULTURE_QUALITY","ID_LEVEL","KINGDOM","PHYLUM","CLASS",
              "ORDER","FAMILY","GENUS","SPECIES","SAMPLE_TOTAL","TAXA_TOTAL",
              "TAXA_COUNT","TAXA_PERC_ID","TAXA_REP_SEQ","TAXA_ASV","ABUN_ASV_COUNT",
              "ABUN_ASV_PERC_ID","ABUN_ASV_SEQ","ABUN_ASV")]
    return(t)
  
  print("You did it! You Identified your bacteria!")
}

