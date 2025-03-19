#Beta-diversity - NMDS
#MSc FStrait - AnaGusmao

library("vegan")
library("phyloseq")
library("ggplot2")
library ("dplyr")
library("patchwork")


#Load phyloseq objects
load("nphy_pall.RData") 
load("nphy_eall.RData")

load("nphy_pktrap.RData") 
load("nphy_ektrap.RData")
load("nphy_pkras.RData") 
load("nphy_ekras.RData")

# --------- Calculate Dissimilarities ----

###  Bray Curtis
#Uses the Hellinger transformed data

#Prok
ord_prok_bc <- ordinate(nphy_pall, "NMDS", "bray")
ord_prok_t_bc <- ordinate(nphy_pktrap, "NMDS", "bray")
ord_prok_r_bc <- ordinate(nphy_pkras, "NMDS", "bray")

#MicroEuk
ord_euk_bc <- ordinate(nphy_eall, "NMDS", "bray")
ord_euk_t_bc <- ordinate(nphy_ektrap, "NMDS", "bray")
ord_euk_r_bc <- ordinate(nphy_ekras, "NMDS", "bray")

# --------------- Plot NMDS ----

# Define custom colors
hab.col <- c("Free-living Surface Water" = "cornflowerblue", 
             "Particle-Attached" = "rosybrown")

season.col <- c("Winter" = "turquoise3","Spring" = "olivedrab",
                "Summer" = "darkorange","Autumn" = "indianred2")

month.col <- c("Jan" = "#0084a8", "Feb" = "#00a2c7", "Mar" = "#00c2d1",
               "Apr" = "#32b48c", "May" = "olivedrab", "Jun" = "#d6a833",
               "Jul" = "darkorange", "Aug" = "indianred", "Sep" = "#b8262e",
               "Oct" = "#8a1707", "Nov" = "#5a4fdf", "Dec" = "#2b74e6")



#Plot functions
NMDS.habsea <- function(phyhel, ord, plotName){
  #Calculate polygons
  pol <- as.data.frame(scores(ord, display = "sites"))  # Extract NMDS scores
  pol$Season <- sample_data(phyhel)$Season  # Add 'Season' metadata
  pol$Habitat <- sample_data(phyhel)$Habitat 
  
  polygon <- pol %>% group_by(Habitat) %>% slice(chull(NMDS1, NMDS2))
  
  #Plot
  plot_ordination(phyhel, ord, type="samples", color="Season", shape="Habitat") + 
    geom_point(size=2) +
    #geom_text(aes(label = sample_names(phyhel)), vjust = -1, size = 3) +
    geom_polygon(data = polygon, aes(x = NMDS1, y = NMDS2, group = Habitat, fill = Habitat), 
                 alpha = 0.2, color = NA) +
    scale_color_manual(values = season.col) +
    scale_fill_manual(values = hab.col) +
    theme_bw() +
    labs(title = plotName, subtitle = paste("Stress =", round(ord$stress, 4)))
} #Habitat + Seasons
NMDS.season <- function(phy_hel, ord, plotName){
  
  #Calculate polygons
  pol <- as.data.frame(scores(ord, display = "sites"))  # Extract NMDS scores
  pol$Season <- sample_data(phy_hel)$Season  # Add 'Season' metadata
  
  polygon <- pol %>% group_by(Season) %>% slice(chull(NMDS1, NMDS2))
  
  #Plot
  plot_ordination(phy_hel, ord, type="samples", color="Season") + 
    geom_point(size=2) +
    geom_polygon(data = polygon, aes(x = NMDS1, y = NMDS2, group = Season, fill = Season), 
                 alpha = 0.3, color = NA) +
    scale_color_manual(values = season.col) +
    scale_fill_manual(values = season.col) +
    theme_bw() +
    labs(title = plotName, subtitle = paste("Stress =", round(ord$stress, 4)))
    
} #By Season
NMDS.months <- function(phy_hel, ord, plotName){
  
  #Define custom order for the months
  sample_data(phy_hel)$month <- factor(sample_data(phy_hel)$month, 
                                       levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                                                  "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  
  #Calculate polygons
  pol <- as.data.frame(scores(ord, display = "sites"))  # Extract NMDS scores
  pol$month <- sample_data(phy_hel)$month  # Add 'Month' metadata
  
  polygon <- pol %>% group_by(month) %>% slice(chull(NMDS1, NMDS2))
  
  
  #Plot
    plot_ordination(phy_hel, ord, type="samples", color="month") + 
    geom_point(size=2) +
    geom_polygon(data = polygon, aes(x = NMDS1, y = NMDS2, group = month, fill = month), 
                 alpha = 0.3, color = NA) +
    scale_color_manual(values = month.col) +
    scale_fill_manual(values = month.col) +
    theme_bw() +
    labs(title = plotName, subtitle = paste("Stress =", round(ord$stress, 4)))
} #By months

#Habitat + Season
plot_prokbc <- NMDS.habsea(nphy_pall,ord_prok_bc,"Prokaryotes") 
plot_euk_bc <- NMDS.habsea(nphy_eall,ord_euk_bc, "Microbial Eukaryotes")

#Season
sea.trapprok_bc <- NMDS.season(nphy_pktrap, ord_prok_t_bc, "Particle-Attached Prokaryotes")
sea.rasprok_bc <- NMDS.season(nphy_pkras, ord_prok_r_bc, "Surface Water Column Prokaryotes")
sea.trapeuk_bc <- NMDS.season(nphy_ektrap, ord = ord_euk_t_bc, "Particle-Attached Microbial Eukaryotes")
sea.raseuk_bc <- NMDS.season(nphy_ekras, ord = ord_euk_r_bc, "Surface Water Column Microbial Eukaryotes")

#Months
mon.trapprok_bc <- NMDS.months(nphy_pktrap, ord_prok_t_bc, "Particle-Attached Prokaryotes")
mon.rasprok_bc <- NMDS.months(nphy_pkras, ord_prok_r_bc, "Surface Water Column Prokaryotes")
mon.trapeuk_bc <- NMDS.months(nphy_ektrap, ord = ord_euk_t_bc, "Particle-Attached Microbial Eukaryotes")
mon.raseuk_bc <- NMDS.months(nphy_ekras, ord = ord_euk_r_bc, "Surface Water Column Microbial Eukaryotes")

#Grid and Display the plots

#Habitat + Season
comb_hab <- (plot_prokbc / plot_euk_bc) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right") +
  plot_annotation(tag_levels ="A")

comb_hab

ggsave("NMDS_HAB_BC.png", plot = comb_hab, width = 8.5, height = 8, units = "in")

#Season
comb_sea <- (sea.rasprok_bc / sea.trapprok_bc | sea.raseuk_bc / sea.trapeuk_bc) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right") +
  plot_annotation(tag_levels ="A")

comb_sea

ggsave("NMDS_SEASON_BC.png", plot = comb_sea, width = 11, height = 9, units = "in")
ggsave("NMDS_SEASON_BC_Klein.png", plot = comb_sea, width = 9, height = 7, units = "in")

#Month
comb_mon <- (mon.rasprok_bc / mon.trapprok_bc | mon.raseuk_bc / mon.trapeuk_bc) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right") +
  plot_annotation(tag_levels ="A")

comb_mon

ggsave("NMDS_MONTH_BC.png", plot = comb_mon, width = 11, height = 9, units = "in")
ggsave("NMDS_MONTH_BC_Klein.png", plot = comb_mon, width = 9, height = 7, units = "in")

save.image(file = "R_DATA_betadiv.RData")


# --------------- Testing Other categories for Sinking Part Prok ----

env <- sample_data(nphy_pktrap)
env2 <- sample_data(nphy_ektrap)

env$Flux_Mag <- cut(env$Flux_PoM, 
                   breaks = c(-Inf, 10, 30, Inf), 
                   labels = c("Weak", "Moderate", "Strong"))

env2$Flux_Mag <- cut(env2$Flux_PoM, 
                    breaks = c(-Inf, 10, 30, Inf), 
                    labels = c("Weak", "Moderate", "Strong"))

sample_data(nphy_pktrap) <- env
sample_data(nphy_ektrap) <- env2


NMDS.magni <- function(phy_hel, ord, plotName){
  
  #Calculate polygons
  pol <- as.data.frame(scores(ord, display = "sites"))  # Extract NMDS scores
  pol$Flux_Mag <- sample_data(phy_hel)$Flux_Mag  # Add 'Season' metadata
  
  polygon <- pol %>% group_by(Flux_Mag) %>% slice(chull(NMDS1, NMDS2))
  
  #Plot
  plot_ordination(phy_hel, ord, type="samples", color="Flux_Mag") + 
    geom_point(size=2) +
    geom_polygon(data = polygon, aes(x = NMDS1, y = NMDS2, group = Flux_Mag, fill = Flux_Mag), 
                 alpha = 0.3, color = NA) +
    scale_color_manual(values = c("purple", "goldenrod", "forestgreen")) +  # Directly specifying 3 colors
    scale_fill_manual(values = c("purple", "goldenrod", "forestgreen")) + 
    theme_bw() +
    labs(title = plotName, subtitle = paste("Stress =", round(ord$stress, 4)))
  
} #By Flux Magnitude

magni.pt <- NMDS.magni(nphy_pktrap, ord_prok_t_bc, "Flux Magnitude of Particle-Attached Prokaryotes")
magni.pt

magni.et <- NMDS.magni(nphy_ektrap, ord_euk_t_bc, "Flux Magnitude of Particle-Attached MicroEukaryotes")
magni.et

#Habitat + Season + MAGNI PT
comb3 <- (comb_hab / magni.pt) +
  plot_annotation(tag_levels ="A")

comb3

library(gridExtra)

grid.arrange(comb_hab, magni.pt, ncol = 1)

ggsave("NMDS_HAB_MAGNI_PT_BC.png", plot = comb3, width = 8.5, height = 8, units = "in")

# --------------- Jaccard ----

#Transform in Binary Data (Presence/Absence) 

bin_prok <- transform_sample_counts(phy_prok_fl, function(x) ifelse(x > 0, 1, 0))
bin_proktrap <- transform_sample_counts(phy_prok_trap, function(x) ifelse(x > 0, 1, 0))
bin_prokras <- transform_sample_counts(phy_prok_ras, function(x) ifelse(x > 0, 1, 0))

bin_euk <- transform_sample_counts(phy_euk, function(x) ifelse(x > 0, 1, 0))
bin_euktrap <- transform_sample_counts(phy_euk_trap, function(x) ifelse(x > 0, 1, 0))
bin_eukras <- transform_sample_counts(phy_euk_ras, function(x) ifelse(x > 0, 1, 0))

# --------- Calculate Dissimilarities

### Jaccard
#Uses RAW data

#Prok
ord_prok_ja <- ordinate(bin_prok, "NMDS", "jaccard")
ord_prok_t_ja <- ordinate(bin_proktrap, "NMDS", "jaccard") #Trap
ord_prok_r_ja <- ordinate(bin_prokras, "NMDS", "jaccard") #RAS

#MicroEuk
ord_euk_ja <- ordinate(bin_euk, "NMDS", "jaccard")
ord_euk_t_ja <- ordinate(bin_euktrap, "NMDS", "jaccard") #Trap
ord_euk_r_ja <- ordinate(bin_eukras, "NMDS", "jaccard") #RAS

#Plots
#Habitat + Season
plot_prok_ja <- NMDS.habsea(nphy_pall,ord_prok_ja,"Prokaryotes") 
plot_euk_ja <- NMDS.habsea(nphy_eall,ord_euk_ja, "Microbial Eukaryotes")

#Season
sea.trapprok_ja <- NMDS.season(nphy_pktrap, ord_prok_t_ja, "Particle-Attached Prokaryotes")
sea.rasprok_ja <- NMDS.season(nphy_pkras, ord_prok_r_ja, "Surface Water Column Prokaryotes")
sea.trapeuk_ja <- NMDS.season(nphy_ektrap, ord = ord_euk_t_ja, "Particle-Attached Microbial Eukaryotes")
sea.raseuk_ja <- NMDS.season(nphy_ekras, ord = ord_euk_r_ja, "Surface Water Column Microbial Eukaryotes")

#Months
mon.trapprok_ja <- NMDS.months(nphy_pktrap, ord_prok_t_bc, "Particle-Attached Prokaryotes")
mon.rasprok_ja <- NMDS.months(nphy_pkras, ord_prok_r_bc, "Surface Water Column Prokaryotes")
mon.trapeuk_ja <- NMDS.months(nphy_ektrap, ord = ord_euk_t_bc, "Particle-Attached Microbial Eukaryotes")
mon.raseuk_ja <- NMDS.months(nphy_ekras, ord = ord_euk_r_bc, "Surface Water Column Microbial Eukaryotes")


#PERMANOVA(ADONIS)

#RAS Prok
otu3 <- as.matrix(otu_table(nphy_pkras))
if (taxa_are_rows(nphy_pkras)) {
  otu3 <- t(otu3)
}

env3 <- sample_data(nphy_pkras)
env3$Season <- as.factor(env3$Season)

perm2pr <- adonis2(otu3 ~ env3$Season, 
                   data = otu3, 
                   method = "bray", 
                   permutations = 999)

perm3pr <- adonis2(otu3 ~ env3$WaterFrac, 
                   data = otu3, 
                   method = "bray", 
                   permutations = 999)

print(perm2pr)
print(perm3pr)

#RAS Euk
otu2 <- as.matrix(otu_table(nphy_ekras))
if (taxa_are_rows(nphy_ekras)) {
  otu2 <- t(otu2)
}

env2 <- sample_data(nphy_ekras)
env2$Season <- as.factor(env2$Season)

perm2er <- adonis2(otu2 ~ env2$Season, 
                  data = otu2, 
                  method = "bray", 
                  permutations = 999)

perm3er <- adonis2(otu2 ~ env2$WaterFrac, 
                  data = otu2, 
                  method = "bray", 
                  permutations = 999)

print(perm2er)
print(perm3er)

#Sinking Particle Euk
otu1 <- as.matrix(otu_table(nphy_ektrap))
if (taxa_are_rows(nphy_ektrap)) {
  otu1 <- t(otu1)
}

env1 <- sample_data(nphy_ektrap)
env1$Flux_Mag <- as.factor(env1$Flux_Mag)
env1$Season <- as.factor(env1$Season)

perm1e <- adonis2(otu1 ~ env1$Flux_Mag, 
                 data = otu1, 
                 method = "bray", 
                 permutations = 999)

perm2e <- adonis2(otu1 ~ env1$Season, 
                 data = otu1, 
                 method = "bray", 
                 permutations = 999)

perm3e <- adonis2(otu1 ~ env1$WaterFrac, 
                 data = otu1, 
                 method = "bray", 
                 permutations = 999)
print(perm1e)
print(perm2e)
print(perm3e)

#Sinking Particle Prok
otu <- as.matrix(otu_table(nphy_pktrap))
if (taxa_are_rows(nphy_pktrap)) {
  otu <- t(otu)
}

env <- sample_data(nphy_pktrap)
env$Flux_Mag <- as.factor(env$Flux_Mag)
env$Season <- as.factor(env$Season)

perm1 <- adonis2(otu ~ env$Flux_Mag, 
                            data = otu, 
                            method = "bray", 
                            permutations = 999)

perm2 <- adonis2(otu ~ env$Season, 
                            data = otu, 
                            method = "bray", 
                            permutations = 999)

perm3 <- adonis2(otu ~ env$WaterFrac, 
                 data = otu, 
                 method = "bray", 
                 permutations = 999)
print(perm1)
print(perm2)
print(perm3)

# Run betadisper to check dispersion within each group
betadisp1 <- betadisper(vegdist(otu, method = "bray"), env$Flux_Mag)
betadisp2 <- betadisper(vegdist(otu, method = "bray"), env$Season)

anova(betadisp1)
anova(betadisp2)