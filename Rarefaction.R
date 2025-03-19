### RAREFACTION 

#MSc FStrait - AnaGusmao
#Modified from RarefacDiversity/DataProcessing/DataLoad - MWietz

library(iNEXT)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(vegan)
library(patchwork)

#------------- Loading files and transforming in OTU TABLE-----------#
#It can also use directly a phyloseq object

# -- PROK
ASV.bac <- read.table("ASV_PROK_FALPHA_AfC.txt", header = TRUE, sep = "\t", row.names = 1)
iNEXT.bac <- otu_table(ASV.bac, taxa_are_rows=T)

# -- EUK
ASV.euk <- read.table("ASV_EUK_FALPHA_AfC.txt", header = TRUE, sep = "\t", row.names = 1)
iNEXT.euk <- otu_table(ASV.euk, taxa_are_rows=T)

#------------- SETTING VARIABLES OF METADATA -----------#

MET.all <- read.table("MET_fstrait.txt", header = TRUE, sep = "\t", row.names = 1)
  
#Subset metadata for each dataset

#For the prokaryotes ALL
met.prok <- MET.all %>%
  select(ID_samples,sample_names_prok, Season)

met.prok <- met.prok[!(met.prok$ID_samples == "PS121_26-3_Fevi40_Up09"), ]
met.prok <- met.prok[, !colnames(met.prok) %in% c("ID_samples")]
rownames(met.prok) <- met.prok[, 1]

#For the prokaryotes trap
met.prok.t <- MET.all %>%
  filter(Data_Type == "Sediment_Trap") %>%
  select(ID_samples,sample_names_prok, Season)

met.prok.t <- met.prok.t[!(met.prok.t$ID_samples == "PS121_26-3_Fevi40_Up09"), ]
met.prok.t <- met.prok.t[, !colnames(met.prok.t) %in% c("ID_samples")]
rownames(met.prok.t) <- met.prok.t[, 1]

#For the prokaryotes ras
met.prok.r <- MET.all  %>%
  filter(Data_Type == "RAS_Amplicon") %>%
  select(ID_samples,sample_names_prok, Season)

met.prok.r <- met.prok.r[, !colnames(met.prok.r) %in% c("ID_samples")]
rownames(met.prok.r) <- met.prok.r[, 1]


#For the eukaryotes ALL
met.euk <- MET.all %>%
  select(ID_samples,sample_names_euk, Season)

met.euk <- met.euk[!(met.euk$ID_samples %in% c("PS99_72-1_Fevi34_Up04",
                                               "PS99_72-1_Fevi34_Up12",
                                               "PS99_72-1_Fevi34_Up15",
                                               "PS121_26-3_Fevi40_Up07",
                                               "PS121_26-3_Fevi40_Up09",
                                               "PS121_26-3_Fevi40_Up10",
                                               "PS121_26-3_Fevi40_Up11",
                                               "PS121_26-3_Fevi40_Up12", 
                                               "RAS_HS2_02_Y18_1",
                                               "RAS_HS2_03_Y18_1")), ]

met.euk <- met.euk[, !colnames(met.euk) %in% c("ID_samples")]
rownames(met.euk) <- met.euk[, 1]

#For the eukaryotes trap
met.euk.t <- MET.all %>%
  filter(Data_Type == "Sediment_Trap") %>%
  select(ID_samples,sample_names_euk, Season)

met.euk.t <- met.euk.t[!(met.euk.t$ID_samples %in% c("PS99_72-1_Fevi34_Up04",
                                                     "PS99_72-1_Fevi34_Up12",
                                                     "PS99_72-1_Fevi34_Up15",
                                                     "PS121_26-3_Fevi40_Up07",
                                                     "PS121_26-3_Fevi40_Up09",
                                                     "PS121_26-3_Fevi40_Up10",
                                                     "PS121_26-3_Fevi40_Up11",
                                                     "PS121_26-3_Fevi40_Up12")), ]

met.euk.t <- met.euk.t[, !colnames(met.euk.t) %in% c("ID_samples")]
rownames(met.euk.t) <- met.euk.t[, 1]

#For the eukaryotes ras
met.euk.r <- MET.all %>%
  filter(Data_Type == "RAS_Amplicon") %>%
  select(ID_samples,sample_names_euk, Season)

met.euk.r <- met.euk.r[!(met.euk.r$ID_samples %in% c("RAS_HS2_02_Y18_1",
                                                     "RAS_HS2_03_Y18_1")), ]

met.euk.r <- met.euk.r[, !colnames(met.euk.r) %in% c("ID_samples")]
rownames(met.euk.r) <- met.euk.r[, 1]

#------------- SUBSET ASV FILE IN TRAP X RAS -----------#

# Get the sample names that match
trap_samples_names_prok <- met.prok.t$sample_names_prok 
ras_samples_names_prok <- met.prok.r$sample_names_prok 

trap_samples_names_euk <- met.euk.t$sample_names_euk 
ras_samples_names_euk <- met.euk.r$sample_names_euk 

# Subset ASV table (keeping only selected samples)
ASV_PROK_TRAP <- ASV.bac %>% select(all_of(trap_samples_names_prok))
ASV_PROK_RAS <- ASV.bac %>% select(all_of(ras_samples_names_prok))

ASV_EUK_TRAP <- ASV.euk %>% select(all_of(trap_samples_names_euk))
ASV_EUK_RAS <- ASV.euk %>% select(all_of(ras_samples_names_euk))

#Creating OTU Table 
iNEXT.bac.trap <- otu_table(ASV_PROK_TRAP, taxa_are_rows=T)
iNEXT.bac.ras <- otu_table(ASV_PROK_RAS, taxa_are_rows=T)

iNEXT.euk.trap <- otu_table(ASV_EUK_TRAP, taxa_are_rows=T)
iNEXT.euk.ras <- otu_table(ASV_EUK_RAS, taxa_are_rows=T)

#------------- GROUP SAMPLES BY SEASONS -----------#

inext_obj_by_seasons <- function (asv_table, met_file)
{
  #Merge ASV and metadata
  ASV.seasons <- as.data.frame(t(asv_table)) # Transpose ASV table so samples are rows
  merged_data <- merge(met_file, as.data.frame(ASV.seasons), by = "row.names", all.x = TRUE)
  merged_data <- merged_data[, -1]
  
  # Aggregate by season
  seasonal_data <- merged_data %>%
    group_by(Season) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE))
  
  # Convert back to ASV table format for iNEXT
  seasonal.asv <- as.data.frame(t(seasonal_data))
  colnames(seasonal.asv) <- seasonal_data$Season
  seasonal.asv <- seasonal.asv[-1,]
  
  # Convert to numeric/integer for iNEXT
  seasonal.asv$Autumn <- as.integer(as.numeric(trimws(seasonal.asv$Autumn)))
  seasonal.asv$Spring <- as.integer(as.numeric(trimws(seasonal.asv$Spring)))
  seasonal.asv$Summer <- as.integer(as.numeric(trimws(seasonal.asv$Summer)))
  seasonal.asv$Winter <- as.integer(as.numeric(trimws(seasonal.asv$Winter)))
  
  #creating the OTU table
  iNEXT.otu <- otu_table(seasonal.asv, taxa_are_rows=T)
  
  return(iNEXT.otu)
}

iNEXT.prok.season <- inext_obj_by_seasons(ASV.bac, met.prok)
iNEXT.prok.trap.season <- inext_obj_by_seasons(ASV_PROK_TRAP, met.prok.t)
iNEXT.prok.ras.season <- inext_obj_by_seasons(ASV_PROK_RAS, met.prok.r)

iNEXT.euk.season <- inext_obj_by_seasons(ASV.euk, met.euk)
iNEXT.euk.trap.season <- inext_obj_by_seasons(ASV_EUK_TRAP, met.euk.t)
iNEXT.euk.ras.season <- inext_obj_by_seasons(ASV_EUK_RAS, met.euk.r)

save.image(file = "iNEXT_PART1.RData")

#------------- SPECIES RICHNESS ESTIMATION -----------#
#It can takes some days to run

#The species richness estimation, with rarefaction and bootstrapping on the abundance data

#Parameters: 
#q = diversity order (0 = species richness)
#nboot = boostrap / how many times will repeat the counting process

### Prokaryotes --- [1]
iNEXT.bac.result <- iNEXT(as.data.frame(iNEXT.bac), q=c(0),
                        datatype="abundance", conf = 0.95, nboot = 100)

iNEXT.bac.trap.result <- iNEXT(as.data.frame(iNEXT.bac.trap), q=c(0),
                         datatype="abundance", conf = 0.95, nboot = 100)

iNEXT.bac.ras.result <- iNEXT(as.data.frame(iNEXT.bac.ras), q=c(0),
                         datatype="abundance", conf = 0.95, nboot = 100)

save.image(file = "iNEXT_100B_PART1.RData")

### Eukaryotes --- [2]
iNEXT.euk.result <- iNEXT(as.data.frame(iNEXT.euk), q=c(0),
                          datatype="abundance", conf = 0.95, nboot = 100)

iNEXT.euk.trap.result <- iNEXT(as.data.frame(iNEXT.euk.trap), q=c(0),
                          datatype="abundance", conf = 0.95, nboot = 100)

iNEXT.euk.ras.result <- iNEXT(as.data.frame(iNEXT.euk.ras), q=c(0),
                          datatype="abundance", conf = 0.95, nboot = 100)

save.image(file = "iNEXT_100B_PART2.RData")

#[by seasons]

iNEXT.season.prok.t.result <- iNEXT(as.data.frame(iNEXT.prok.trap.season), q=c(0),
                                  datatype="abundance", conf = 0.95, nboot = 100)

iNEXT.season.prok.r.result <- iNEXT(as.data.frame(iNEXT.prok.ras.season), q=c(0),
                                  datatype="abundance", conf = 0.95, nboot = 100)


iNEXT.season.euk.t.result <- iNEXT(as.data.frame(iNEXT.euk.trap.season), q=c(0),
                                  datatype="abundance", conf = 0.95, nboot = 100)

iNEXT.season.euk.r.result <- iNEXT(as.data.frame(iNEXT.euk.ras.season), q=c(0),
                                 datatype="abundance", conf = 0.95, nboot = 100)

save.image(file = "iNEXT_100B_PART3.RData")


#[if with a phyloseq]
#iNEXT.example.result <- iNEXT(as.data.frame(otu_table(iNEXT.euk)), q=c(0),datatype="abundance", conf = 0.95, nboot = 100)


#------------- PROCESSING THE RESULTS -----------#
#Type: type of output. 
#Type 1 is extrapolation/rarefaction curves (Sample-sized based)
#Type 3 is coverage-based

#All ----

#Rarefaction Curves
rich.curves <- function (iNextresult, title) {
  
  covcurve <- ggiNEXT(iNextresult, type=1) +
    scale_color_manual(values = rep("black",148)) +
    scale_fill_manual(values = rep("grey", 148)) + 
    theme(legend.position = "none") +
    geom_line(linewidth = 0.3) +
    ggtitle(title)
}

prokrich <- rich.curves (iNEXT.bac.result, "Prokaryotes")
eukrich <- rich.curves (iNEXT.euk.result, "Microbial Eukaryotes")

#Coverage
cov.curves <- function (iNextresult, title) {
  
  covcurve <- ggiNEXT(iNextresult, type=2) +
    scale_color_manual(values = rep("black",148)) +
    scale_fill_manual(values = rep("grey", 148)) + 
    theme(legend.position = "none") +
    geom_line(linewidth = 0.3) +
    ggtitle(title)
}

prokcov <- cov.curves (iNEXT.bac.result, "Prokaryotes")
eukcov <- cov.curves (iNEXT.euk.result, "Microbial Eukaryotes")

#Grid
comb_rar <- (prokrich / prokcov | eukrich / eukcov) +
  plot_annotation(tag_levels ="A") 

comb_rar

ggsave("RAR_ALL.png", plot = comb_rar, width = 16, height = 9, units = "in")

#By datasets ----

bydf.rar.curves <- function (iNextresult, title) {
  
  bydf.curve <- ggiNEXT(iNextresult, type=1) + 
  scale_color_manual(values = rep("black",148)) +
  scale_fill_manual(values = rep("grey", 148)) + 
  theme(legend.position = "none") +
  geom_line(linewidth = 0.3) +
  ggtitle(title)
  
}

ptrap.curve <- bydf.rar.curves(iNEXT.bac.trap.result, "Sinking-Particle Prokaryotes")
pras.curve <- bydf.rar.curves(iNEXT.bac.ras.result, "Surface Water Prokaryotes")
etrap.curve <- bydf.rar.curves(iNEXT.euk.trap.result, "Sinking-Particle MicroEukaryotes")
eras.curve <- bydf.rar.curves(iNEXT.euk.ras.result, "Surface Water MicroEukaryotes")

comb4 <- (pras.curve / ptrap.curve | eras.curve / etrap.curve) +
  plot_annotation(tag_levels ="A") 

comb4

ggsave("RAR_4DATASETS.png", plot = comb4, width = 16, height = 9, units = "in")

#By Seasons ----

sea.rar.curves <- function (iNextresult, title){
  
  seacurve <- ggiNEXT(iNextresult, type=1, se = T) +
    scale_color_manual(values = c("Spring" = "olivedrab", 
                                  "Summer" = "darkorange", 
                                  "Autumn" = "indianred2", 
                                  "Winter" = "turquoise3")) +
    ggtitle(title)
  
}

seaprok.t <- sea.rar.curves(iNEXT.season.prok.t.result, "Sinking-Particle Prokaryotes")
seaprok.r <- sea.rar.curves(iNEXT.season.prok.r.result, "Surface Water Prokaryotes")
seaeuk.t <- sea.rar.curves(iNEXT.season.euk.t.result, "Sinking-Particle Microbial Eukaryotes")
seaeuk.r <- sea.rar.curves(iNEXT.season.euk.r.result, "Surface Water Microbial Eukaryotes")


comb.sea <- (seaprok.r / seaprok.t | seaeuk.r / seaeuk.t) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels ="A") & 
  theme(legend.position = "bottom")

comb.sea

ggsave("RAR_SEASONS.png", plot = comb.sea, width = 16, height = 9, units = "in")

