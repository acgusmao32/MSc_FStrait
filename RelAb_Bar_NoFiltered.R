#Relative Abundance - Bar Plots
#MSc FStrait - AGUSMAO

library(phyloseq)
library(ggplot2)
library(dplyr)

# ---------- Load phyloseq objects ----

load("phy_prok_NOFILTER.RData")

# ---------- Subset the phyloseq object
# Prok
prok_trap2 <- subset_samples(phy_prok, Habitat == "Particle-Associated") #No filter

# ---------- Aggregate at taxon level ----
taxon_level <- "Family"

gloom_pt_fam2 <- tax_glom(prok_trap2, taxon_level)

# ---------- Transform to relative abundance ----

pt_relfam2 <- transform_sample_counts(gloom_pt_fam2, function(x) x / sum(x)* 100) #in %

#Functions ----

# Function to merge rare taxa into Others < threshold
merge_rare_taxa1 <- function(taxa, threshold, taxa_name) { 
  
  # Find mean relative abundance for each taxon across all samples
  taxa_mean_abundance <- taxa_sums(taxa) / nsamples(taxa)
  
  # Identify taxa below the threshold
  rare_taxa <- taxa_names(taxa)[taxa_mean_abundance < threshold]
  
  # Merge rare taxa into "Others"
  physeq_merged <- merge_taxa(taxa, rare_taxa)
  
  # Update the tax table to name the merged taxa as "Others"
  tax_table(physeq_merged)[taxa_names(physeq_merged) %in% rare_taxa, taxa_name] <- paste("Others < ",as.character(threshold),"%")
  
  return(physeq_merged)
}


# Function that organizes Season and Years into a chronological order
organize_data_SeasonYear <- function(dates, seasons){
  month <- as.numeric(format(dates, "%m"))
  year <- as.numeric(format(dates, "%y"))
  
  if(seasons == "Summer")
  {
    return(paste("Sum",as.character(year)))
  }
  else if(seasons == "Autumn")
  {
    return(paste("Aut",as.character(year)))
  }
  else if(seasons == "Spring")
  {
    return(paste("Spr",as.character(year)))
  }
  else
  {
    if(month<5)
    {
      return(paste("Win",as.character(year-1)))
    }
    else
    {
      return(paste("Win",as.character(year)))
    }
  }
}

#Plot relative abundance in chrono order
ab_plot_prok <- function(taxa, threshold, taxa_name, graphName){
  # Organize the data into Season+Year
  sample_data(taxa)$Season_Year <- mapply(organize_data_SeasonYear, 
                                          as.Date(sample_data(taxa)$Date_Average), sample_data(taxa)$Season)
  
  #Order the Season and Year
  sample_data(taxa)$Season_Year <- factor(sample_data(taxa)$Season_Year, 
                                          levels = c("Sum 16","Aut 16","Win 16","Spr 17","Sum 17","Aut 17",
                                                     "Win 17","Spr 18","Sum 18","Aut 18","Win 18","Spr 19",
                                                     "Sum 19","Aut 19","Win 19","Spr 20","Sum 20","Aut 20",
                                                     "Win 20","Spr 21","Sum 21","Aut 21","Win 21"))
  
  #First merge all rare samples using the function merge_rare_taxa1
  taxa_merged <- merge_rare_taxa1(taxa, threshold, taxa_name)
  
  # Add row names as a new column in the sample data (metadata) of the phyloseq object
  sample_data(taxa_merged)$Row_Names <- rownames(sample_data(taxa_merged))
  
  # Modify the Row_Names column to remove '_prok'
  sample_data(taxa_merged)$Row_Names <- gsub("_prok$", "", sample_data(taxa_merged)$Row_Names)
  
  #Now plot the merged taxa into a nice bar plot
  plot_bar(taxa_merged, x = "Row_Names",fill = taxa_name) +
    geom_area(stat = "identity") +
    facet_grid(~ Season_Year, scales = "free_x", space = "free_x") +
    theme_minimal() +
    labs(title = "Relative Abundance", subtitle = graphName , x = "", y = "Relative Abundance (%)") +
    scale_fill_brewer(palette = "Spectral") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.5),
          strip.text = element_text(size = 8),
          legend.position = "bottom",
          legend.box.spacing = unit(-18, "pt"),
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill = "lightblue", color = "black"),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.2, "lines"),
          plot.margin = margin(b = 1, t = 5, r =5, l = 5)) +
    guides(fill = guide_legend(nrow = 2))
  
}
  
# Set threshold for rare taxon and organization of data in percentage 2%
thr_abund <- 2

#Plots ----

# Prok trap

#By Samples
plot_prok_trap_NOFILT <- ab_plot_prok(pt_relfam2,thr_abund,taxon_level,"Prokaryotes Particle-Associated - No Filtered")
plot_prok_trap_NOFILT

ggsave("BAR_RelAb_Prok_Trap_NO_FILTERED.png", plot = plot_prok_trap_NOFILT, width = 16, height = 9, units = "in")
