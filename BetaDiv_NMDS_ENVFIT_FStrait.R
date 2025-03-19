#Beta-diversity - NMDS + ENVFIT
#MSc FStrait - AnaGusmao
#Modified from: my old scripts gerlache strait

library(phyloseq)
library(vegan)
library(ggplot2)
library(patchwork)

#For envfit, it cant have NAs on environmental data, so it will be necessary to filter
#samples with NA in the variables and re-run the NMDS and hellinger transformation

#Load phyloseq objects
load("phy_prok.RData")
load("phy_euk.RData")

#Remove the outlier samples: Fevi36_Up15_prok, Fevi36_Up09_prok
phy_prok_fl <- subset_samples(phy_prok, !(sample_names(phy_prok) %in% c("Fevi36_Up15_prok", "Fevi36_Up09_prok")))

#-------- Organizing ----

#Change name variables
met.prok <- sample_data(phy_prok_fl)
met.euk <- sample_data(phy_euk)


colnames(met.prok)[colnames(met.prok) == "Data_Type"] <- "Habitat"
colnames(met.euk)[colnames(met.euk) == "Data_Type"] <- "Habitat"

met.prok$Habitat <- recode(met.prok$Habitat,
                           "RAS_Amplicon" = "Free-living Surface Water",
                           "Sediment_Trap" = "Particle-Attached")

met.euk$Habitat <- recode(met.euk$Habitat,
                          "RAS_Amplicon" = "Free-living Surface Water",
                          "Sediment_Trap" = "Particle-Attached")

#other variables for envfit
colnames(met.prok)[colnames(met.prok) == "Flux_BothSi"] <- "SIL"
colnames(met.euk)[colnames(met.euk) == "Flux_BothSi"] <- "SIL"
colnames(met.prok)[colnames(met.prok) == "Flux_PoM"] <- "POM"
colnames(met.euk)[colnames(met.euk) == "Flux_PoM"] <- "POM"
colnames(met.prok)[colnames(met.prok) == "Flux_Carbonate"] <- "CARB"
colnames(met.euk)[colnames(met.euk) == "Flux_Carbonate"] <- "CARB"
colnames(met.prok)[colnames(met.prok) == "Flux_POC"] <- "POC"
colnames(met.euk)[colnames(met.euk) == "Flux_POC"] <- "POC"
colnames(met.prok)[colnames(met.prok) == "Flux_PON"] <- "PON"
colnames(met.euk)[colnames(met.euk) == "Flux_PON"] <- "PON"
colnames(met.prok)[colnames(met.prok) == "Oxygen_Conc"] <- "DO"
colnames(met.euk)[colnames(met.euk) == "Oxygen_Conc"] <- "DO"
colnames(met.prok)[colnames(met.prok) == "fraction_AW"] <- "FAW"
colnames(met.euk)[colnames(met.euk) == "fraction_AW"] <- "FAW"
colnames(met.prok)[colnames(met.prok) == "fraction_PW"] <- "FPW"
colnames(met.euk)[colnames(met.euk) == "fraction_PW"] <- "FPW"
colnames(met.prok)[colnames(met.prok) == "Ice_Conc"] <- "IC"
colnames(met.euk)[colnames(met.euk) == "Ice_Conc"] <- "IC"
colnames(met.prok)[colnames(met.prok) == "Ice_Dist"] <- "ID"
colnames(met.euk)[colnames(met.euk) == "Ice_Dist"] <- "ID"
colnames(met.prok)[colnames(met.prok) == "Salinity"] <- "SAL"
colnames(met.euk)[colnames(met.euk) == "Salinity"] <- "SAL"
colnames(met.prok)[colnames(met.prok) == "Temperature"] <- "TEMP"
colnames(met.euk)[colnames(met.euk) == "Temperature"] <- "TEMP"


sample_data(phy_prok_fl) <- met.prok
sample_data(phy_euk) <- met.euk


#Create subsets
phy_prok_trap <- subset_samples(phy_prok_fl, Habitat == "Particle-Attached")
phy_prok_ras <- subset_samples(phy_prok_fl, Habitat == "Free-living Surface Water")

phy_euk_trap <- subset_samples(phy_euk, Habitat == "Particle-Attached")
phy_euk_ras <- subset_samples(phy_euk, Habitat == "Free-living Surface Water")

#Remove samples with NA on Env Parameters ----

checkenv <- sample_data(phy_prok_fl) #check data

#Samples Fevi34_Up01 and Fevi34_Up02 doesnt have data for MLD
for_remove <- c("Fevi34_Up01_prok", "Fevi34_Up02_prok")
phyproktrap_ef <- prune_samples(!(sample_names(phy_prok_trap) %in% for_remove), phy_prok_trap)

for_remove <- c("Fevi34_Up01_euk", "Fevi34_Up02_euk")
phyeuktrap_ef <- prune_samples(!(sample_names(phy_euk_trap) %in% for_remove), phy_euk_trap)

rm(checkenv)
rm(for_remove)

#----------------------- NORMALIZATION - RAREFACTION ----

#ProkTrap
head(sort(sample_sums(phyproktrap_ef)), 5) #29577
size_prok1 <- head(sort(sample_sums(phyproktrap_ef)), 5)[1] #no sample removed
pt.rar <- rarefy_even_depth(phyproktrap_ef, trimOTUs = T, sample.size = size_prok1,
                            rngseed = 42)

#ProkRAS
head(sort(sample_sums(phy_prok_ras)), 5) #13614
size_prok2 <- head(sort(sample_sums(phy_prok_trap)), 5)[2] #RAS_HS1_P3365
pr.rar <- rarefy_even_depth(phy_prok_ras, trimOTUs = T, sample.size = size_prok2,
                            rngseed = 42)

#EukTrap
head(sort(sample_sums(phyeuktrap_ef)), 5) #21058
size_euk1 <- head(sort(sample_sums(phyeuktrap_ef)), 5)[3] #Fevi36 Up15 and Fevi36 Up09
et.rar <- rarefy_even_depth(phyeuktrap_ef, trimOTUs = T, sample.size = size_euk1, 
                            rngseed = 42)

#EukRAS
head(sort(sample_sums(phy_euk_ras)), 5) #11066
size_euk2 <- head(sort(sample_sums(phy_euk_ras)), 5)[2] #HS3 Ps121 E35
er.rar <- rarefy_even_depth(phy_euk_ras, trimOTUs = T, sample.size = size_euk2, 
                            rngseed = 42)


#Hellinger transformation and Ordination ----

#Hellinger
proktrap.hel <- decostand(as(otu_table(pt.rar), "matrix"), method = "hellinger")
prokras.hel <- decostand(as(otu_table(pr.rar), "matrix"), method = "hellinger")

euktrap.hel <- decostand(as(otu_table(et.rar), "matrix"), method = "hellinger")
eukras.hel <- decostand(as(otu_table(er.rar), "matrix"), method = "hellinger")

#recreating phy object
newphy.proktrap <- phyloseq(otu_table(proktrap.hel, taxa_are_rows = TRUE), tax_table(phy_prok_trap), sample_data(phy_prok_trap))
newphy.prokras <- phyloseq(otu_table(prokras.hel, taxa_are_rows = TRUE), tax_table(phy_prok_ras), sample_data(phy_prok_ras))

newphy.euktrap <- phyloseq(otu_table(euktrap.hel, taxa_are_rows = TRUE), tax_table(phy_euk_trap), sample_data(phy_euk_trap))
newphy.eukras <- phyloseq(otu_table(eukras.hel, taxa_are_rows = TRUE), tax_table(phy_euk_ras), sample_data(phy_euk_ras))

#Ordination Bray-Curtis
ord_prok_t_bc <- ordinate(newphy.proktrap, "NMDS", "bray") #Trap
ord_prok_r_bc <- ordinate(newphy.prokras, "NMDS", "bray") #RAS

ord_euk_t_bc <- ordinate(newphy.euktrap, "NMDS", "bray") #Trap
ord_euk_r_bc <- ordinate(newphy.eukras, "NMDS", "bray") #RAS

#Run Envfit - Environmental Parameters ----
ef_envpar <- function(phyhel, ord_nmds){
  #Extract env data from the phyloseq object
  env_data <- as.data.frame(sample_data(phyhel))
  
  #Edit the metadata: Remove variables with NA and biogeochemical part measurements
  #Chl has a lot of NAs
  env_data <- env_data[, !colnames(env_data) %in% c("Cup", "Sampling_days", "Strat",
                                                    "PON", "POC",
                                                    "POM", "Flux_Sil", 
                                                    "Flux_Opal", "CARB",
                                                    "SIL", "Ocean", "Chl", "CO2",
                                                    "Date_begin", "Date_end", "Date_Average",
                                                    "monthFull", "month", "Habitat", "Latitude",
                                                    "Longitude","day", "Equipment", "Year",
                                                    "Depth", "Oxygen_Sat", "Lag_Date1", "Lag_Date2")]
  
  # Run envfit on the NMDS ordination
  env_fit <- envfit(ord_nmds, env_data, permutations = 999)
  
  # Print significant environmental variables
  print(env_fit)
  
  return(env_fit)
} #Without Fluxes
ef_envpar2 <- function(phyhel, ord_nmds){
  #Extract env data from the phyloseq object
  env_data <- as.data.frame(sample_data(phyhel))
  
  #Edit the metadata: Remove variables with NA and biogeochemical part measurements
  #Chl has a lot of NAs
  env_data <- env_data[, !colnames(env_data) %in% c("Cup", "Sampling_days", "Strat",
                                                    "Ocean", "Chl", "CO2","Date_begin",
                                                    "Date_end", "Date_Average","monthFull",
                                                    "month", "Habitat", "Latitude",
                                                    "Longitude", "day", "Equipment", 
                                                    "Year", "Depth", "Oxygen_Sat",
                                                    "POM","Flux_Sil", "Flux_Opal","Lag_Date1", 
                                                    "Lag_Date2")]
  
  # Run envfit on the NMDS ordination
  env_fit <- envfit(ord_nmds, env_data, permutations = 999)
  
  # Print significant environmental variables
  print(env_fit)
  
  return(env_fit)
}

#Run envfit
ef_ptrap <- ef_envpar2 (newphy.proktrap, ord_prok_t_bc)
ef_pras <- ef_envpar (newphy.prokras, ord_prok_r_bc)

ef_etrap <- ef_envpar2 (newphy.euktrap, ord_euk_t_bc)
ef_eras <- ef_envpar (newphy.eukras, ord_euk_r_bc)

#Results were saved on Envfit_Results.xlsx

#Plot NMDS with Envfit ----

# Define custom colors
hab.col <- c("Free-living Surface Water" = "cornflowerblue", 
             "Particle-Attached" = "rosybrown")

season.col <- c("Winter" = "turquoise3","Spring" = "olivedrab",
                "Summer" = "darkorange","Autumn" = "indianred2")

month.col <- c("Jan" = "#0084a8", "Feb" = "#00a2c7", "Mar" = "#00c2d1",
               "Apr" = "#32b48c", "May" = "olivedrab", "Jun" = "#d6a833",
               "Jul" = "darkorange", "Aug" = "indianred", "Sep" = "#b8262e",
               "Oct" = "#8a1707", "Nov" = "#5a4fdf", "Dec" = "#2b74e6")

#Season
#Plot function
plot_nmds_envf <- function(envfit_file, ord, phyhel, plotname){
  # Extract the env scores
  env_scores <- scores(envfit_file, display = "vectors")
  
  #Calculate polygons
  pol <- as.data.frame(scores(ord, display = "sites"))  # Extract NMDS scores
  pol$Season <- sample_data(phyhel)$Season  # Add 'Season' metadata
  
  polygon <- pol %>% group_by(Season) %>% slice(chull(NMDS1, NMDS2))
  
  #Plot
  nmds_plot <- plot_ordination(phyhel, ord, type="samples", color="Season") + 
    geom_point(size=5) +
    geom_polygon(data = polygon, aes(x = NMDS1, y = NMDS2, group = Season, fill = Season), 
                 alpha = 0.3, color = NA) +
    scale_color_manual(values = season.col) +
    scale_fill_manual(values = season.col) +
    theme_bw() +
    labs(title = plotname, subtitle = paste("Stress =", round(ord$stress, 4)))
  
  # Add environmental vectors to the NMDS plot
  finalplot <- nmds_plot + 
    geom_segment(data = as.data.frame(env_scores), 
                 aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                 arrow = arrow(type = "open", length = unit(0.2, "inches")), 
                 color = "black", linewidth = 0.8) + 
    geom_text(data = as.data.frame(env_scores), 
              aes(x = NMDS1*1.1, y = NMDS2*1.1, label = rownames(env_scores)), 
              color = "black", size = 4, hjust = -0.2, vjust = -0.2)

  return(finalplot) }

#Run the Plots
ef.ptrap <- plot_nmds_envf(ef_ptrap, ord_prok_t_bc, newphy.proktrap, "Sinking-Particle Prokaryotes")
ef.pras <- plot_nmds_envf(ef_pras, ord_prok_r_bc, newphy.prokras, "Surface Water Column Prokaryotes")

ef.etrap <- plot_nmds_envf(ef_etrap, ord_euk_t_bc, newphy.euktrap, "Sinking-Particle Microbial Eukaryotes")
ef.eras <- plot_nmds_envf(ef_eras, ord_euk_r_bc, newphy.eukras, "Surface Water Column Microbial Eukaryotes")

#Grid and Display the plots

comb_seap <- (ef.pras / ef.ptrap) +
  plot_layout(guides = "collect") + theme(legend.position = "below") +
  plot_annotation(tag_levels ="A")

comb_seae <- (ef.eras / ef.etrap) +
  plot_layout(guides = "collect") + theme(legend.position = "below") +
  plot_annotation(tag_levels ="A")

ggsave("NMDS_EF_SEA_ProkTrap.png", plot = ef.ptrap, width = 11, height = 9, units = "in")
ggsave("NMDS_EF_SEA_ProkRas.png", plot = ef.pras, width = 11, height = 9, units = "in")
ggsave("NMDS_EF_SEA_EukTrap.png", plot = ef.etrap, width = 11, height = 9, units = "in")
ggsave("NMDS_EF_SEA_EukRas.png", plot = ef.eras, width = 11, height = 9, units = "in")

ggsave("NMDS_EF_SEA_EUK.png", plot = comb_seae, width = 9, height = 11, units = "in")
ggsave("NMDS_EF_SEA_PROK.png", plot = comb_seap, width = 9, height = 11, units = "in")

#Month
plot_nmds_envf_mon <- function(envfit_file, ord, phyhel, plotname){
  
  # Extract the env scores
  env_scores <- scores(envfit_file, display = "vectors")
  
  #Define custom order for the months
  sample_data(phyhel)$month <- factor(sample_data(phyhel)$month, 
                                       levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                                                  "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  #Calculate polygons
  pol <- as.data.frame(scores(ord, display = "sites"))  # Extract NMDS scores
  pol$month <- sample_data(phyhel)$month  # Add 'Month' metadata
  
  polygon <- pol %>% group_by(month) %>% slice(chull(NMDS1, NMDS2))
  
  #Plot
  nmds_plot <- plot_ordination(phyhel, ord, type="samples", color="month") + 
    geom_point(size=5) +
    geom_polygon(data = polygon, aes(x = NMDS1, y = NMDS2, group = month, fill = month), 
                 alpha = 0.3, color = NA) +
    scale_color_manual(values = month.col) +
    scale_fill_manual(values = month.col) +
    theme_bw() +
    labs(title = plotname, subtitle = paste("Stress =", round(ord$stress, 4)))
  
  # Add environmental vectors to the NMDS plot
  finalplot <- nmds_plot + 
    geom_segment(data = as.data.frame(env_scores), 
                 aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                 arrow = arrow(type = "open", length = unit(0.2, "inches")), 
                 color = "black", linewidth = 0.8) + 
    geom_text(data = as.data.frame(env_scores), 
              aes(x = NMDS1*1.1, y = NMDS2*1.1, label = rownames(env_scores)), 
              color = "black", size = 4, hjust = -0.2, vjust = -0.2)
  
  return(finalplot) }
plot_nmds_envf_mon_sig <- function(envfit_file, ord, phyhel, plotname, p_threshold = 0.05) {
  
  # Extract environmental vector scores and p-values
  env_scores <- as.data.frame(scores(envfit_file, display = "vectors"))
  env_pvals <- envfit_file$vectors$pvals
  
  # Filter only significant vectors based on p-value threshold
  env_scores$P_value <- env_pvals  # Add p-values to the table
  env_scores <- env_scores[env_scores$P_value < p_threshold, ]  # Keep only significant vectors
  
  # Define custom order for the months
  sample_data(phyhel)$month <- factor(sample_data(phyhel)$month, 
                                      levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                                                 "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  # Calculate polygons for months
  pol <- as.data.frame(scores(ord, display = "sites"))  # Extract NMDS site scores
  pol$month <- sample_data(phyhel)$month  # Add 'Month' metadata
  
  polygon <- pol %>% group_by(month) %>% slice(chull(NMDS1, NMDS2))  # Convex hull per month
  
  # Base NMDS plot
  nmds_plot <- plot_ordination(phyhel, ord, type = "samples", color = "month") + 
    geom_point(size = 5) +
    geom_polygon(data = polygon, aes(x = NMDS1, y = NMDS2, group = month, fill = month), 
                 alpha = 0.3, color = NA) +
    scale_color_manual(values = month.col) +
    scale_fill_manual(values = month.col) +
    theme_bw() +
    labs(title = plotname, subtitle = paste("Stress =", round(ord$stress, 4)))
  
  # Add environmental vectors (only significant ones)
  finalplot <- nmds_plot
  if (nrow(env_scores) > 0) {  # Only add if there are significant vectors
    finalplot <- finalplot + 
      geom_segment(data = env_scores, 
                   aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                   arrow = arrow(type = "open", length = unit(0.2, "inches")), 
                   color = "black", linewidth = 0.8) + 
      geom_text(data = env_scores, 
                aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = rownames(env_scores)), 
                color = "black", size = 4, hjust = -0.2, vjust = -0.2)
  }
  
  return(finalplot)
}

library(ggrepel)
plot_nmds_envf_mon_sig2 <- function(envfit_file, ord, phyhel, plotname, p_threshold = 0.05) {
  
  # Extract environmental vector scores and p-values
  env_scores <- as.data.frame(scores(envfit_file, display = "vectors"))
  env_pvals <- envfit_file$vectors$pvals
  
  # Filter only significant vectors based on p-value threshold
  env_scores$P_value <- env_pvals  # Add p-values to the table
  env_scores <- env_scores[env_scores$P_value < p_threshold, ]  # Keep only significant vectors
  
  # Define custom order for the months
  sample_data(phyhel)$month <- factor(sample_data(phyhel)$month, 
                                      levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                                                 "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  # Calculate polygons for months
  pol <- as.data.frame(scores(ord, display = "sites"))  # Extract NMDS site scores
  pol$month <- sample_data(phyhel)$month  # Add 'Month' metadata
  
  polygon <- pol %>% group_by(month) %>% slice(chull(NMDS1, NMDS2))  # Convex hull per month
  
  # Base NMDS plot
  nmds_plot <- plot_ordination(phyhel, ord, type = "samples", color = "month") + 
    geom_point(size = 2) +
    geom_polygon(data = polygon, aes(x = NMDS1, y = NMDS2, group = month, fill = month), 
                 alpha = 0.3, color = NA) +
    scale_color_manual(values = month.col) +
    scale_fill_manual(values = month.col) +
    theme_bw() +
    labs(title = plotname, subtitle = paste("Stress =", round(ord$stress, 4)))
  
  # Add environmental vectors (only significant ones)
  finalplot <- nmds_plot
  if (nrow(env_scores) > 0) {  # Only add if there are significant vectors
    finalplot <- finalplot + 
      geom_segment(data = env_scores, 
                   aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                   arrow = arrow(type = "open", length = unit(0.2, "inches")), 
                   color = "black", linewidth = 0.8) + 
      geom_text_repel(data = env_scores, 
                      aes(x = NMDS1 * 1, y = NMDS2 * 1, label = rownames(env_scores)), 
                      color = "black", size = 4, box.padding = 0.5, 
                      point.padding = 0.5, max.overlaps = Inf,
                      segment.color = NA)
  }
  
  return(finalplot)
}

ef.ptrap.mon <- plot_nmds_envf_mon_sig2(ef_ptrap, ord_prok_t_bc, newphy.proktrap, "Sinking-Particle Prokaryotes")
ef.pras.mon <- plot_nmds_envf_mon_sig2(ef_pras, ord_prok_r_bc, newphy.prokras, "Surface Water Column Prokaryotes")
ef.etrap.mon <- plot_nmds_envf_mon_sig2(ef_etrap, ord_euk_t_bc, newphy.euktrap, "Sinking-Particle Microbial Eukaryotes")
ef.eras.mon <- plot_nmds_envf_mon_sig2(ef_eras, ord_euk_r_bc, newphy.eukras, "Surface Water Column Microbial Eukaryotes")

#Grid and Display the plots
comb_mon1 <- (ef.etrap.mon  / ef.ptrap.mon) +
  plot_layout(guides = "collect") + theme(legend.position = "below") +
  plot_annotation(tag_levels ="A")

comb_mon2 <- (ef.eras.mon / ef.pras.mon) +
  plot_layout(guides = "collect") + theme(legend.position = "below") +
  plot_annotation(tag_levels ="A")

comb_mon4 <- (ef.pras.mon/ ef.ptrap.mon | ef.eras.mon / ef.etrap.mon) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right") +
  plot_annotation(tag_levels ="A")

ggsave("NMDS_EF_MON_ProkTrap.png", plot = ef.ptrap.mon, width = 11, height = 9, units = "in")
ggsave("NMDS_EF_MON_ProkRas.png", plot = ef.pras.mon, width = 11, height = 9, units = "in")
ggsave("NMDS_EF_MON_EukTrap.png", plot = ef.etrap.mon, width = 11, height = 9, units = "in")
ggsave("NMDS_EF_MON_EukRas.png", plot = ef.eras.mon, width = 11, height = 9, units = "in")

ggsave("NMDS_EF_MON_TRAP.png", plot = comb_mon1, width = 8.5, height = 10, units = "in")
ggsave("NMDS_EF_MON_RAS.png", plot = comb_mon2, width = 8.5, height = 10, units = "in")
ggsave("NMDS_EF_MON.png", plot = comb_mon4, width = 11, height = 9, units = "in")
# -------------------------------------------
save.image(file = "R_DATA_betadiv_EF.RData")
