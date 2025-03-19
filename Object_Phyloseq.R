#Phyloseq Objects
#MSc - FStrait - AnaGusmao

library("phyloseq")
library("dplyr")

#---------- LOADING FILES ----
#rownames of TAX needs to be the same as the OTU/ASV table

#Reading files - PROKARYOTES [1]
seq_table_prok <- as.matrix(read.table("Input/ASV_PROK_AfC.txt", header = TRUE, row.names = 1))
tax_table_prok <- as.matrix(read.table("Input/TAX_PROK_AfC.txt", header = TRUE, sep="\t", row.names=1)) 

#Reading files - EUKARYOTES [2]
seq_table_euk <- as.matrix(read.table("Input/ASV_EUK_AfC.txt", header = TRUE, row.names = 1))
tax_table_euk <- as.matrix(read.table("Input/TAX_EUK_AfC.txt", header = TRUE, sep="\t", row.names=1)) 

#Reading files - FOR ALPHA METRICS [3]
seq.alph.prok <- as.matrix(read.table("Input/ASV_PROK_FALPHA_AfC.txt", header = TRUE, row.names = 1))
tax.alph.prok <- as.matrix(read.table("Input/TAX_PROK_FALPHA_AfC.txt", header = TRUE, sep="\t", row.names=1))

seq.alph.euk <- as.matrix(read.table("Input/ASV_EUK_FALPHA_AfC.txt", header = TRUE, row.names = 1))
tax.alph.euk <- as.matrix(read.table("Input/TAX_EUK_FALPHA_AfC.txt", header = TRUE, sep="\t", row.names=1))

#Reading - metadata
metadata <- read.table("Input/MET_fstrait.txt", header=TRUE, sep="\t")

#---------- Subsample the metadata ----

#PROK [1] and [3]
metadata_PROK <- metadata %>%
  filter(!ID_samples %in% c("Fevi40_Up09")) %>%
  select(-sample_names_euk, -ID_samples)

rownames(metadata_PROK) <- metadata_PROK$sample_names_prok

metadata_PROK <- metadata_PROK[, -1]


#EUK [2] and [3]
metadata_EUK <- metadata %>%
  filter(!ID_samples %in% c("Fevi34_Up04",
                            "Fevi34_Up12",
                            "Fevi34_Up15",
                            "Fevi40_Up07",
                            "Fevi40_Up09",
                            "Fevi40_Up10",
                            "Fevi40_Up11",
                            "Fevi40_Up12",
                            "HS2_PS114_E9",
                            "HS2_PS114_E11")) %>%
  select(-ID_samples, -sample_names_prok) 

rownames(metadata_EUK) <- metadata_EUK$sample_names_euk

metadata_EUK <- metadata_EUK[, -1]

#---------- Creating the object phyloseq ----

#[1]
phy_prok <- phyloseq(
  otu_table(seq_table_prok, taxa_are_rows = TRUE),
  tax_table (tax_table_prok),
  sample_data(metadata_PROK))


#[2]
phy_euk <- phyloseq(
  otu_table(seq_table_euk, taxa_are_rows = TRUE),
  tax_table (tax_table_euk),
  sample_data(metadata_EUK))


#[3]
phy_prok_alph <- phyloseq(
  otu_table(seq.alph.prok, taxa_are_rows = TRUE),
  tax_table (tax.alph.prok),
  sample_data(metadata_PROK))

phy_euk_alph <- phyloseq(
  otu_table(seq.alph.euk, taxa_are_rows = TRUE),
  tax_table (tax.alph.euk),
  sample_data(metadata_EUK))


#---------- Change name variables ----

to_hab <- function (phy) {
  #extract
  met <- sample_data(phy)
  
  #Change Columns names
  colnames(met)[colnames(met) == "Data_Type"] <- "Habitat"
  
  #Change Variable names
  met$Habitat <- recode(met$Habitat,
                        "RAS_Amplicon" = "Free-living Surface Water",
                        "Sediment_Trap" = "Particle-Associated")
  #Updating the phy object
  sample_data(phy) <- met
  
  return(phy)
}

phy_prok <- to_hab(phy_prok)
phy_euk <- to_hab(phy_euk)
phy_prok_alph <- to_hab(phy_prok_alph)
phy_euk_alph <- to_hab(phy_euk_alph)

#---------- Filtering ----

#Remove the outlier samples: Fevi36_Up15_prok, Fevi36_Up09_prok
phy_prok_fl <- subset_samples(phy_prok, !(sample_names(phy_prok) %in% c("Fevi36_Up15_prok", "Fevi36_Up09_prok")))
phy_prok_fl_alph <- subset_samples(phy_prok_alph, !(sample_names(phy_prok_alph) %in% c("Fevi36_Up15_prok", "Fevi36_Up09_prok")))

#----------------------- EXPORTING OBJECT PHYLOSEQ ---------------#

save(phy_prok, file = "phy_prok_NOFILTER.RData") #[1]
save(phy_prok_fl, file = "phy_prok.RData") #[1]
save(phy_euk, file = "phy_euk.RData") #[2]
save(phy_prok_fl_alph, file = "phy_prok_alph.RData") #[3]
save(phy_euk_alph, file = "phy_euk_alph.RData") #[3]



