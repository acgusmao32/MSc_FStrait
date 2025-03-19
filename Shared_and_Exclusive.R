#Shared and Exclusive ASVS
#MSc FSTrait - Ana Gusmao

library(phyloseq)
library(dplyr)
library(ggvenn)
library(gridExtra)

# ---------- Load phyloseq objects
load("phy_prok.RData")
load("phy_euk.RData")

#------------ RAREFACTION -----

#Prokaryotes All
head(sort(sample_sums(phy_prok_fl)), 5)  #13614
size_prok <- head(sort(sample_sums(phy_prok_fl)), 5)[2]

# Rarefy
p.rar <- rarefy_even_depth(phy_prok_fl, trimOTUs = T, sample.size = size_prok, rngseed = 42)

#Microbial Eukaryotes All

#Value Used: 11066 
#Sample Removed: Fevi36_Up15_euk (530), HS3_PS121_E35 (1992), Fevi36_Up09_euk (2958)

head(sort(sample_sums(phy_euk)), 5) 
size_euk <- head(sort(sample_sums(phy_euk)), 5)[4]

# Rarefy
e.rar <- rarefy_even_depth(phy_euk, trimOTUs = T, sample.size = size_euk, rngseed = 42)

#---------------------- Venn Diagram - Habitat ----

venndi.hab <- function(phyrar, plot_subname){
  #Define the categories
  categories <- unique(sample_data(phyrar)$Habitat)
  category_names <- c("Particle-Associated", "Free-living Surface Water")
  
  # Subset the phyloseq object for each category
  phyloseq_cat1 <- subset_samples(phyrar, Habitat == "Particle-Associated")
  phyloseq_cat2 <- subset_samples(phyrar, Habitat == "Free-living Surface Water")
  
  # Extract ASVs (OTUs) present in each category
  asvs_cat1 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat1) > 0, phyloseq_cat1))
  asvs_cat2 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat2) > 0, phyloseq_cat2))
  
  # Create a list of ASVs for each category with custom names
  asv_list <- list(
    "Particle-associated" = asvs_cat1,
    "Free-living Surface Water" = asvs_cat2)
  
  #Plot the Venn Diagram
  vennpl <- ggvenn(asv_list, fill_color = c("rosybrown", "cornflowerblue"),
                   stroke_size = 0.5, set_name_size = 4, text_size = 4) +
    labs(subtitle = plot_subname) +
    theme(plot.subtitle = element_text(vjust = 4, hjust = 0.1, size = 13),
          plot.margin = margin(t = 6))
  
  return(vennpl)
}

#run the plots
prok.hab <- venndi.hab(p.rar, "Prokaryotes")
euk.hab <- venndi.hab(e.rar, "Microbial Eukaryotes")

#Grid
hab.all <-grid.arrange(prok.hab, euk.hab, ncol = 1)

#Export
ggsave("VennD_Hab_All.png", plot = hab.all, width = 8.5, height = 7, units = "in")

# --------------------- Venn Diagram - Seasons ----
#Rarefaction for each dataset

#Subset of phyloseq

phy.prok.trap.ir <- subset_samples(phy_prok_fl, Habitat == "Particle-Associated")
phy.prok.ras.ir <- subset_samples(phy_prok_fl, Habitat == "Free-living Surface Water")

phy.euk.trap.ir <- subset_samples(phy_euk, Habitat == "Particle-Associated")
phy.euk.ras.ir <- subset_samples(phy_euk, Habitat == "Free-living Surface Water")

#Rarefaction

#Prok - Trap
#Value Used: 29577 | Sample Removed: None
head(sort(sample_sums(phy.prok.trap.ir)), 5)
size_prok <- head(sort(sample_sums(phy.prok.trap.ir)), 5)[1] # Value used (29577)

p.t.rar.ir <- rarefy_even_depth(phy.prok.trap.ir, trimOTUs = T, sample.size = size_prok, rngseed = 42)

#Prok - RAS
#Value Used: 13614 | Sample Removed: RAS_HS1_P3365
head(sort(sample_sums(phy.prok.ras.ir)), 5)
size_prok2 <- head(sort(sample_sums(phy.prok.ras.ir)), 5)[2] # Value used (13614)

p.r.rar.ir <- rarefy_even_depth(phy.prok.ras.ir, trimOTUs = T, sample.size = size_prok2, rngseed = 42)


#Microbial Eukaryotes - Trap 
#Value Used: 21058 | Sample Removed: Fevi36_Up15_euk, Fevi36_Up09_euk
head(sort(sample_sums(phy.euk.trap.ir)), 5) 
size_euk <- head(sort(sample_sums(phy.euk.trap.ir)), 5)[3]

e.t.rar.ir <- rarefy_even_depth(phy.euk.trap.ir, trimOTUs = T, sample.size = size_euk, rngseed = 42)

#Microbial Eukaryotes - RAS 
#Value Used: 11066 | Sample Removed: 
head(sort(sample_sums(phy.euk.ras.ir)), 5) 
size_euk2 <- head(sort(sample_sums(phy.euk.ras.ir)), 5)[2]

e.r.rar.ir <- rarefy_even_depth(phy.euk.ras.ir, trimOTUs = T, sample.size = size_euk2, rngseed = 42)


#Plots 

venndi.sea.trap <- function(phyrar, plot_subname){
  #Define the categories
  categories <- unique(sample_data(phyrar)$Season)
  category_names <- c("Spring", "Summer", "Autumn", "Winter")
  
  # Subset the phyloseq object for each dataset
  phy_hab <- subset_samples(phyrar, Habitat == "Particle-Associated")
  
  # Subset the phyloseq object for each category
  phyloseq_cat1 <- subset_samples(phy_hab, Season == "Spring")
  phyloseq_cat2 <- subset_samples(phy_hab, Season == "Summer")
  phyloseq_cat3 <- subset_samples(phy_hab, Season == "Autumn")
  phyloseq_cat4 <- subset_samples(phy_hab, Season == "Winter")
  
  # Extract ASVs (OTUs) present in each category
  asvs_cat1 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat1) > 0, phyloseq_cat1))
  asvs_cat2 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat2) > 0, phyloseq_cat2))
  asvs_cat3 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat3) > 0, phyloseq_cat3))
  asvs_cat4 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat4) > 0, phyloseq_cat4))
  
  # Create a list of ASVs for each category with custom names
  asv_list <- list(
    "Spring" = asvs_cat1, "Summer" = asvs_cat2,
    "Autumn" = asvs_cat3, "Winter" = asvs_cat4)
  
  #Plot the Venn Diagram
  vennpl <- ggvenn(asv_list, 
                   fill_color = c("olivedrab", "darkorange", "indianred", "turquoise3"),
                   stroke_size = 0.5, set_name_size = 4, text_size = 3) +
    labs(subtitle = plot_subname) +
    theme(plot.subtitle = element_text(vjust = 6, hjust = 0.1, size = 13))

}
venndi.sea.ras <- function(phyrar, plot_subname){
  #Define the categories
  categories <- unique(sample_data(phyrar)$Season)
  category_names <- c("Spring", "Summer", "Autumn", "Winter")
  
  # Subset the phyloseq object for each dataset
  phy_hab <- subset_samples(phyrar, Habitat == "Free-living Surface Water")
  
  # Subset the phyloseq object for each category
  phyloseq_cat1 <- subset_samples(phy_hab, Season == "Spring")
  phyloseq_cat2 <- subset_samples(phy_hab, Season == "Summer")
  phyloseq_cat3 <- subset_samples(phy_hab, Season == "Autumn")
  phyloseq_cat4 <- subset_samples(phy_hab, Season == "Winter")
  
  # Extract ASVs (OTUs) present in each category
  asvs_cat1 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat1) > 0, phyloseq_cat1))
  asvs_cat2 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat2) > 0, phyloseq_cat2))
  asvs_cat3 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat3) > 0, phyloseq_cat3))
  asvs_cat4 <- taxa_names(prune_taxa(taxa_sums(phyloseq_cat4) > 0, phyloseq_cat4))
  
  # Create a list of ASVs for each category with custom names
  asv_list <- list(
    "Spring" = asvs_cat1, "Summer" = asvs_cat2,
    "Autumn" = asvs_cat3, "Winter" = asvs_cat4)
  
  #Plot the Venn Diagram
  vennpl <- ggvenn(asv_list, 
                   fill_color = c("olivedrab", "darkorange", "indianred", "turquoise3"),
                   stroke_size = 0.5, set_name_size = 4, text_size = 3) +
    labs(subtitle = plot_subname) +
    theme(plot.subtitle = element_text(vjust = 6, hjust = 0.1, size = 13))
  
  return(vennpl)
}

prok.t.sea.ir <- venndi.sea.trap(p.t.rar.ir, "Particle-associated Prokaryotes")
euk.t.sea.ir <- venndi.sea.trap(e.t.rar.ir, "Particle-associated Microbial Eukaryotes")

prok.r.sea.ir <- venndi.sea.ras(p.r.rar.ir, "Free-living Surface Water Prokaryotes")
euk.r.sea.ir <- venndi.sea.ras(e.r.rar.ir, "Free-living Surface Water Microbial Eukaryotes")


#Grid

sea.trap.ir <-grid.arrange(prok.t.sea.ir, euk.t.sea.ir, ncol = 2)
sea.ras.ir <-grid.arrange(prok.r.sea.ir, euk.r.sea.ir, ncol = 2)

#Export
ggsave("VennD_Sea_Trap_IR.png", plot = sea.trap.ir, width = 8.5, height = 7, units = "in")
ggsave("VennD_Sea_RAS_IR.png", plot = sea.ras.ir, width = 8.5, height = 7, units = "in")

#### ---- ####
save.image(file = "R_DATA_Shared_and_Ex.RData")