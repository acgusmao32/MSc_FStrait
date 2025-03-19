#ANOSIM
library(phyloseq)
library(vegan)
packageVersion("vegan")



#Loading R DATA ----
load("phy_prok.RData")
load("phy_euk.RData")

#Remove the outlier samples: Fevi36_Up15_prok, Fevi36_Up09_prok
phy_prok_fl <- subset_samples(phy_prok, !(sample_names(phy_prok) %in% c("Fevi36_Up15_prok", "Fevi36_Up09_prok")))
rm(phy_prok)

#Change name variables

to_hab <- function (phy) {
  #extract
  met <- sample_data(phy)
  
  #Change Columns names
  colnames(met)[colnames(met) == "Data_Type"] <- "Habitat"
  
  #Change Variable names
  met$Habitat <- recode(met$Habitat,
                             "RAS_Amplicon" = "Free-living Surface Water",
                             "Sediment_Trap" = "Particle-Attached")
  #Updating the phy object
  sample_data(phy) <- met
  
  return(phy)
}

phy_prok_up <- to_hab(phy_prok_fl)
phy_euk_up <- to_hab(phy_euk)

rm(phy_prok_fl)
rm(phy_euk)

#Create subsets
phy_prok_trap <- subset_samples(phy_prok_up, Habitat == "Particle-Attached")
phy_prok_ras <- subset_samples(phy_prok_up, Habitat == "Free-living Surface Water")

phy_euk_trap <- subset_samples(phy_euk_up, Habitat == "Particle-Attached")
phy_euk_ras <- subset_samples(phy_euk_up, Habitat == "Free-living Surface Water")

#Calculate Dissimilarities ----

hell <- function(phy) {
  
  #Hellinger Transformation
  phyhel <- decostand(as(otu_table(phy), "matrix"), method = "hellinger")
  
  #Recreate phy object
  newphy <- phyloseq(otu_table(phyhel, taxa_are_rows = TRUE), tax_table(phy), sample_data(phy))
  
  return(newphy)
  }

phyhel.p <- hell(phy_prok_up)
phyhel.e <- hell(phy_euk_up)

phyhel.ptrap <- hell(phy_prok_trap)
phyhel.pras <- hell(phy_prok_ras)
phyhel.etrap <- hell(phy_prok_trap)
phyhel.eras <- hell(phy_euk_ras)

#Ordination Bray-Curtis
ord_p_bc <- phyloseq::distance(phyhel.p, method = "bray")
ord_e_bc <- phyloseq::distance(phyhel.e, method = "bray")

ord_ptrap_bc <- phyloseq::distance(phyhel.ptrap, method = "bray")
ord_pras_bc <- phyloseq::distance(phyhel.pras, method = "bray")
ord_etrap_bc <- phyloseq::distance(phyhel.etrap, method = "bray")
ord_eras_bc <- phyloseq::distance(phyhel.eras, method = "bray")

#Run PERMANOVA (Adonis) ----

#Prok RAS
adonis_p_ras <- adonis2(ord_pras_bc ~ sample_data(phyhel.pras)$Season, permutations = 999)
print(adonis_p_ras)

#Prok Trap
adonis_p_trap <- adonis2(ord_ptrap_bc ~ sample_data(phyhel.ptrap)$Season, permutations = 999)
print(adonis_p_trap)

#MicroEuk RAS
adonis_e_ras <- adonis2(ord_eras_bc ~ sample_data(phyhel.eras)$Season, permutations = 999)
print(adonis_e_ras)

#MicroEuk Trap
adonis_e_trap <- adonis2(ord_etrap_bc ~ sample_data(phyhel.etrap)$Season, permutations = 999)
print(adonis_e_trap)

# Perform ANOSIM ----

#Convert categories to factors
to_fac <- function(phy) {
  
  #Month
  sample_data(phy)$month <- as.factor(sample_data(phy)$month)
  
  #Season
  sample_data(phy)$Season <- as.factor(sample_data(phy)$Season)
  
  #Habitat
  sample_data(phy)$Habitat <- as.factor(sample_data(phy)$Habitat)

  return(phy)
}

phyhel.p <- to_fac(phyhel.p)
phyhel.e <- to_fac(phyhel.e)

phyhel.ptrap <- to_fac(phyhel.ptrap)
phyhel.pras <- to_fac(phyhel.pras)
phyhel.etrap <- to_fac(phyhel.etrap)
phyhel.eras <- to_fac(phyhel.eras)

#For the categories tested on the NMDS -> Season/Habitat/Month

#Season
ans.sea.ptrap <- anosim(ord_ptrap_bc, sample_data(phyhel.ptrap)$Season, permutations = 999)
ans.sea.pras <- anosim(ord_pras_bc, sample_data(phyhel.pras)$Season, permutations = 999)
ans.sea.etrap <- anosim(ord_etrap_bc, sample_data(phyhel.etrap)$Season, permutations = 999)
ans.sea.eras <- anosim(ord_eras_bc, sample_data(phyhel.eras)$Season, permutations = 999)

ans.sea.ptrap #sig
ans.sea.pras #sig
ans.sea.etrap #sig
ans.sea.eras #sig

#Month
ans.mon.ptrap <- anosim(ord_ptrap_bc, sample_data(phyhel.ptrap)$month)
ans.mon.pras <- anosim(ord_pras_bc, sample_data(phyhel.pras)$month)
ans.mon.etrap <- anosim(ord_etrap_bc, sample_data(phyhel.etrap)$month)
ans.mon.eras <- anosim(ord_eras_bc, sample_data(phyhel.eras)$month)

ans.mon.ptrap #sig
ans.mon.pras #sig
ans.mon.etrap #sig
ans.mon.eras #sig

#Habitat
ans.hab.p <- anosim(ord_p_bc, sample_data(phyhel.p)$Habitat)
ans.hab.e <- anosim(ord_e_bc, sample_data(phyhel.e)$Habitat)

ans.hab.p #sig
ans.hab.e #sig
