#Mantel Prok x MicroEuk
#Use the files still without rarefaction

##---------Loading R DATA ----

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
#It needs to be the same number of samples

##---------Removing Samples ----
#It needs to be the same number of samples between both datasets

#Remove samples in one dataset that doesnt have in the other
phy_prok_trap2 <- subset_samples(phy_prok_trap, !(sample_names(phy_prok_trap) %in%
                                                    c("Fevi34_Up04_prok", 
                                                      "Fevi34_Up12_prok",
                                                      "Fevi34_Up15_prok",
                                                      "Fevi40_Up07_prok",
                                                      "Fevi40_Up10_prok",
                                                      "Fevi40_Up11_prok",
                                                      "Fevi40_Up12_prok")))

phy_euk_trap2 <- subset_samples(phy_euk_trap, !(sample_names(phy_euk_trap) %in%
                                                  c("Fevi36_Up09_euk",
                                                    
                                                    "Fevi36_Up15_euk")))

#----------------------- NORMALIZATION - RAREFACTION ----

#ProkTrap
head(sort(sample_sums(phy_prok_trap2)), 5) #29577
size_prok1 <- head(sort(sample_sums(phy_prok_trap2)), 5)[1] #no sample removed
pt.rar <- rarefy_even_depth(phy_prok_trap2, trimOTUs = T, sample.size = size_prok1,
                            rngseed = 42)

#EukTrap
head(sort(sample_sums(phy_euk_trap2)), 5) #21058
size_euk1 <- head(sort(sample_sums(phy_euk_trap2)), 5)[1] #no sample removed
et.rar <- rarefy_even_depth(phy_euk_trap2, trimOTUs = T, sample.size = size_euk1, 
                            rngseed = 42)

#Hellinger transformation and Ordination ----

#Hellinger
proktrap.hel <- decostand(as(otu_table(pt.rar), "matrix"), method = "hellinger")
prokras.hel <- decostand(as(otu_table(pr.rar), "matrix"), method = "hellinger")

euktrap.hel <- decostand(as(otu_table(et.rar), "matrix"), method = "hellinger")
eukras.hel <- decostand(as(otu_table(er.rar), "matrix"), method = "hellinger")

#recreating phy object
newphy.proktrap2 <- phyloseq(otu_table(proktrap.hel, taxa_are_rows = TRUE), tax_table(phy_prok_trap), sample_data(phy_prok_trap))
newphy.euktrap2 <- phyloseq(otu_table(euktrap.hel, taxa_are_rows = TRUE), tax_table(phy_euk_trap), sample_data(phy_euk_trap))

#new distance matrice
ord_ptrap_bc2 <- phyloseq::distance(newphy.proktrap2, method = "bray")
ord_etrap_bc2 <- phyloseq::distance(newphy.euktrap2, method = "bray")

#mantel test
mn_result_all <- mantel(ord_ptrap_bc2, ord_etrap_bc2, method = "spearman", permutations = 999)

sca_mantel <- function(ord, env_ord, mantel_result, plotname) {
  
  #Extract Bray-Curtis ordination distances (as a vector)
  ord_bc_vector <- as.vector(ord)
  
  # Extract Euclidean distances of POC flux (as a vector)
  env_vector <- as.vector(env_ord)
  
  # Create a data frame with the distances
  dist_data <- data.frame(BrayCurtis = ord_bc_vector, EuclideanPOC = env_vector)
  
  # Extract Mantel test result values
  mantel_corr <- mantel_result$statistic  # Mantel correlation coefficient
  mantel_pval <- mantel_result$signif     # Mantel p-value
  
  # Plot using ggplot2
  plot <- ggplot(dist_data, aes(x = BrayCurtis, y = EuclideanPOC)) +
    geom_point(size = 3, alpha = 0.5, colour = "black", shape = 21, fill ="skyblue") + 
    geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
    labs(x = "Bray-Curtis Dissimilarity Prok", y = "Bray-Curtis Dissimilarity MicroEuk",
         subtitle = plotname) +
    theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
          axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
          axis.title= element_text(face = "bold", size = 14, colour = "black"), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, colour = "black")) +
    annotate("text", x = min(dist_data$BrayCurtis) * 0.8, y = max(dist_data$EuclideanPOC) * 0.9, 
             label = paste("Mantel R = ", round(mantel_corr, 3), "\np-value = ", round(mantel_pval, 3)), 
             size = 5, hjust = 0, colour = "black")
}

trap_mnall <- sca_mantel(ord_ptrap_bc2, ord_etrap_bc2, mn_result_all, "Sinking-particle Prokaryotes x Microbial Eukaryotes")
trap_mnall

ggsave("MANTEL_PROK_MICROEUK.png", plot = trap_mnall, width = 7, height = 7, units = "in")
