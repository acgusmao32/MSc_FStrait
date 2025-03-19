#Mantel Test Correlogram

library(dplyr)
library(linkET)
library(ggplot2)
library(FD)

#Loading R DATA ----
load("nphy_pktrap.RData")
load("nphy_ektrap.RData")

load("nphy_pkras.RData")
load("nphy_ekras.RData")


#Change Env names
to_env_codes <- function (phy) {
  #extract
  met <- sample_data(phy)
  
  #Change Columns names
  colnames(met)[colnames(met) == "Flux_BothSi"] <- "SIL"
  colnames(met)[colnames(met) == "Flux_Carbonate"] <- "CARB"
  colnames(met)[colnames(met) == "Flux_POC"] <- "POC"
  colnames(met)[colnames(met) == "Flux_PON"] <- "PON"
  colnames(met)[colnames(met) == "Oxygen_Conc"] <- "DO"
  colnames(met)[colnames(met) == "fraction_AW"] <- "FAW"
  colnames(met)[colnames(met) == "fraction_PW"] <- "FPW"
  colnames(met)[colnames(met) == "Ice_Conc"] <- "SIC"
  colnames(met)[colnames(met) == "Ice_Dist"] <- "SID"
  colnames(met)[colnames(met) == "Salinity"] <- "SAL"
  colnames(met)[colnames(met) == "Temperature"] <- "TEMP"
  
  #Updating the phy object
  sample_data(phy) <- met
  
  return(phy)
}

nphy_pktrap <- to_env_codes(nphy_pktrap)
nphy_ektrap <- to_env_codes(nphy_ektrap)
nphy_pkras <- to_env_codes(nphy_pkras)
nphy_ekras <- to_env_codes(nphy_ekras)

#Mantel Correlogram ----
#Extract files

#Extract env for correlogram
env <- sample_data(nphy_pktrap)
env <- env[, !names(env) %in% c("Habitat", "Equipment", "Cup", "monthFull", "month", "day",
                                  "Latitude", "Longitude", "Ocean", "Season", "Year", 
                                  "Sampling_days", "Date_begin", "Date_end", "Lag_Date1",
                                  "Lag_Date2", "Date_Average", "IceCat", "WaterFrac",
                                  "CO2", "Chl", "MLD", "Depth", "Oxygen_Sat", "Flux_Sil",
                                  "Flux_PoM", "Flux_Opal")]

env1 <- sample_data(nphy_pkras)
env1 <- env1[, !names(env1) %in% c("Habitat", "Equipment", "Cup", "monthFull", "month", "day",
                                "Latitude", "Longitude", "Ocean", "Season", "Year", 
                                "Sampling_days", "Date_begin", "Date_end", "Lag_Date1",
                                "Lag_Date2", "Date_Average", "IceCat", "WaterFrac",
                                "CO2", "Chl", "MLD", "Depth", "Oxygen_Sat", "Flux_Sil",
                                "Flux_PoM", "Flux_Opal", "PON", "POC","CARB", "SIL")]

#Mn test
mntest <- function(phy, who) {
  
  #Extract OTU table
  otu <- as.data.frame(otu_table(phy))
  otu <- t(otu)
  
  #Extract env
  env <- sample_data(phy)
  env <- env[, !names(env) %in% c("Habitat", "Equipment", "Cup", "monthFull", "month", "day",
                                  "Latitude", "Longitude", "Ocean", "Season", "Year", 
                                  "Sampling_days", "Date_begin", "Date_end", "Lag_Date1",
                                  "Lag_Date2", "Date_Average", "IceCat", "WaterFrac",
                                  "CO2", "Chl", "MLD", "Depth", "Oxygen_Sat", "Flux_Sil",
                                  "Flux_PoM", "Flux_Opal")]
  
  #Run Mantel test
  mantel <- mantel_test(otu, env, method = "spearman", permutations = 999) %>%
    mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
           pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
  
  #> `mantel_test()` using 'bray' dist method for 'spec'.
  #> `mantel_test()` using 'euclidean' dist method for 'env'.
  
  #Change Columns names
  mantel$spec <- recode(mantel$spec,
                        "spec" = who)
  
  return(mantel)
}
mntest.ras <- function(phy, who) {
  
  #Extract OTU table
  otu <- as.data.frame(otu_table(phy))
  otu <- t(otu)
  
  #Extract env
  env <- sample_data(phy)
  env <- env[, !names(env) %in% c("Habitat", "Equipment", "Cup", "monthFull", "month", "day",
                                  "Latitude", "Longitude", "Ocean", "Season", "Year", 
                                  "Sampling_days", "Date_begin", "Date_end", "Lag_Date1",
                                  "Lag_Date2", "Date_Average", "IceCat", "WaterFrac",
                                  "CO2", "Chl", "MLD", "Depth", "Oxygen_Sat", "Flux_Sil",
                                  "Flux_PoM", "Flux_Opal", "PON", "POC","CARB", "SIL")]
  
  #Run Mantel test
  mantel <- mantel_test(otu, env, method = "spearman", permutations = 999) %>%
    mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
           pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
  
  #> `mantel_test()` using 'bray' dist method for 'spec'.
  #> `mantel_test()` using 'euclidean' dist method for 'env'.
  
  #Change Columns names
  mantel$spec <- recode(mantel$spec,
                        "spec" = who)
  
  return(mantel)
}

mnprok <- mntest(nphy_pktrap, "Prokaryotes")
mneuk <- mntest(nphy_ektrap, "MicroEukaryotes") 

mnprok.ras <- mntest.ras(nphy_pkras, "Prokaryotes")
mneuk.ras <- mntest.ras(nphy_ekras, "MicroEukaryotes") 

#Merge Both Mantel Results
merged_mn <- rbind(mnprok, mneuk)
merged_mn2 <- rbind(mnprok.ras, mneuk.ras)

#Plot - Particle-associated
mn_corr <- qcorrplot(correlate(env), type = "upper", diag = FALSE) +
  geom_tile() +
  geom_couple(aes(colour = pd, size = rd), 
              data = merged_mn, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "PiYG")) +
  scale_size_manual(values = c(0.5, 1.5, 2)) +
  scale_colour_manual(values = c("< 0.01" = "deepskyblue", 
                                 "0.01 - 0.05" = "orange", 
                                 ">= 0.05" = "gray")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

#Plot - Surface Water
mn_corr2 <- qcorrplot(correlate(env1), type = "upper", diag = FALSE) +
  geom_tile() +
  geom_couple(aes(colour = pd, size = rd), 
              data = merged_mn2, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "PiYG")) +
  scale_size_manual(values = c(0.5, 1.5, 2)) +
  scale_colour_manual(values = c("< 0.01" = "deepskyblue", 
                                 "0.01 - 0.05" = "orange", 
                                 ">= 0.05" = "gray")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

#Export image
ggsave("MANTEL_CORR_TRAP.png", plot = mn_corr, width = 16, height = 9, units = "in")
ggsave("MANTEL_CORR_TRAP_KL.png", plot = mn_corr, width = 8.5, height = 10, units = "in")


ggsave("MANTEL_CORR_RAS.png", plot = mn_corr2, width = 16, height = 9, units = "in")
ggsave("MANTEL_CORR_RAS_KL.png", plot = mn_corr2, width = 10, height = 7, units = "in")

