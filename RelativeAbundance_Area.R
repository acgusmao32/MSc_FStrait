#Relative Abundance
#MSc FStrait - AnaGusmao
library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

# ---------- Load phyloseq objects ----
load("phy_prok.RData")
load("phy_euk.RData")

#---------- Number of ASVs ----

#Number of ASVs per sample
phy_prok_fl #13658 ASVs - 148 samples
phy_euk #10635 ASVs - 139 samples

#Number of reads per sample
mean(sample_sums(phy_prok_fl))
sd(sample_sums(phy_prok_fl))

mean(sample_sums(phy_euk))
sd(sample_sums(phy_euk))

#[Large SD = large variation between samples]

#How many 
proktax <- tax_table(phy_prok_fl)
sum(proktax[, "Domain"] == "Bacteria", na.rm = TRUE)
sum(proktax[, "Domain"] == "Archaea", na.rm = TRUE)

#How many
euktax <- tax_table(phy_euk)
unique(euktax[, "Kingdom"]) #10

sum(euktax[, "Kingdom"] == "TSAR", na.rm = TRUE)
sum(euktax[, "Kingdom"] == "Haptista", na.rm = TRUE)
sum(euktax[, "Kingdom"] == "Obazoa", na.rm = TRUE)
sum(euktax[, "Kingdom"] == "Archaeplastida", na.rm = TRUE)
sum(euktax[, "Kingdom"] == "Cryptista", na.rm = TRUE)
sum(euktax[, "Kingdom"] == "Cryptista:nucl", na.rm = TRUE)
sum(euktax[, "Kingdom"] == "Amoebozoa", na.rm = TRUE)

#---------- Subset the phyloseq object ----
phy_ptrap <- subset_samples(phy_prok_fl, Habitat == "Particle-Associated")
phy_pras <- subset_samples(phy_prok_fl, Habitat == "Free-living Surface Water")

# Micro Euk
phy_etrap <- subset_samples(phy_euk, Habitat == "Particle-Associated")
phy_eras <- subset_samples(phy_euk, Habitat == "Free-living Surface Water")

# ---------- Transform to relative abundance ----
pt_rel <- transform_sample_counts(phy_ptrap, function(x) x / sum(x)* 100) #in %
pr_rel <- transform_sample_counts(phy_pras, function(x) x / sum(x)* 100) #in %

et_rel <- transform_sample_counts(phy_etrap, function(x) x / sum(x)* 100) #in %
er_rel <- transform_sample_counts(phy_eras, function(x) x / sum(x)* 100) #in %

# ---------- Aggregate at taxon level ----

taxon_level <- "Family"

gloom_pt_fam <- tax_glom(pt_rel, taxon_level)
gloom_pr_fam <- tax_glom(pr_rel, taxon_level)

gloom_et_fam <- tax_glom(et_rel, taxon_level)
gloom_er_fam <- tax_glom(er_rel, taxon_level)

#taxon_level2 <- "Class"

#gloom_et_cla <- tax_glom(et_rel, taxon_level2)
#gloom_er_cla <- tax_glom(er_rel, taxon_level2)

taxon_level3 <- "Genus"

gloom_pt_gen <- tax_glom(phy_ptrap, taxon_level3)
gloom_pr_gen <- tax_glom(phy_pras, taxon_level3)

gloom_et_gen <- tax_glom(phy_etrap, taxon_level3)
gloom_er_gen <- tax_glom(phy_eras, taxon_level3)

# ---------- The 5 most abundant of a taxon level -----
the5most <- function(phyrel, taxon_level) {
  
  # Extract OTU table and sum across samples
  otu_table_df <- as.data.frame(otu_table(phyrel))
  tax_sums <- rowSums(otu_table_df)  # Summing across all samples
  
  # Create a data frame with taxon names and their total relative abundance
  tax_abundance <- data.frame(Taxon = tax_table(phyrel)[, taxon_level],
                              Relative_Abundance = tax_sums / sum(tax_sums) * 100)
  
  # Get the top 5 most abundant taxa
  top_taxa <- tax_abundance %>%
    arrange(desc(Relative_Abundance)) %>%
    head(5)
  
  return(top_taxa)
}

the5most(gloom_pr_fam, "Family")
the5most(gloom_pt_fam, "Family") 

#the5most(gloom_er_cla, "Class")
#the5most(gloom_et_cla, "Class")

the5most(gloom_er_fam, "Family")
the5most(gloom_et_fam, "Family")

the5most(gloom_er_gen, "Genus")
the5most(gloom_et_gen, "Genus")

test <- as.data.frame(tax_table(gloom_er_gen))

# ---------- Bar Plot ----

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

# Set threshold for rare taxon (%)
thr_abund <- 2

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
ab_plot_euk <- function(taxa, threshold, taxa_name, graphName){
  
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
  sample_data(taxa_merged)$Row_Names <- gsub("_euk$", "", sample_data(taxa_merged)$Row_Names)
  sample_data(taxa_merged)$Row_Names <- gsub("^RAS_", "", sample_data(taxa_merged)$Row_Names)
  
  # Custom colors (13)
  custom_palette <- c(brewer.pal(12, "Paired"),
    "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", 
    "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", 
    "#5E4FA2", "#FF8C94", "#9F5F52")
  
  #Now plot the merged taxa into a nice bar plot
  plot_bar(taxa_merged, x = "Row_Names",fill = taxa_name) +
    geom_area(stat = "identity") +
    facet_grid(~ Season_Year, scales = "free_x", space = "free_x") +
    theme_minimal() +
    labs(title = "Relative Abundance", subtitle = graphName , x = "", y = "Relative Abundance (%)") +
    scale_fill_manual(values = custom_palette) +
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

#Run the plots

# Micro Euk - TRAP
plot_euk_trap_fam <- ab_plot_euk(gloom_pr_fam,thr_abund,
                                  taxon_level,
                                  "Microbial Eukaryotes Particle-Attached")
plot_euk_trap_fam



#Prok
plot_prok_trap_NOFILT <- abundance_plot(pt_relfam2,thr_abund,taxon_level,"Prokaryotes Particle-Attached - No Filtered")
ggsave("BAR_RelAb_Prok_Trap_NO_FILTERED.png", plot = plot_prok_trap_NOFILT, width = 16, height = 9, units = "in")

# Micro Euk
plot_euk_trap_fam <- ab_plot_euk2(gloom_pr_fam,thr_abund,taxon_level,"Microbial Eukaryotes Particle-Attached")
plot_euk_trap_fam



plot_euk_trap_gen <- ab_plot_euk(et_relgen,thr_abund,taxon_level3,"Microbial Eukaryotes Particle-Attached")
ggsave("BAR_RelAb_MicroEuk_Trap_GEN.png", plot = plot_euk_trap_gen, width = 16, height = 9, units = "in")

plot_euk_ras_gen <- ab_plot_euk(er_relgen,thr_abund,taxon_level3,"Microbial Eukaryotes Particle-Attached")
ggsave("BAR_RelAb_MicroEuk_RAS_GEN.png", plot = plot_euk_ras_gen, width = 16, height = 9, units = "in")


# ---------- Area Plot ----

# Remove rare taxa below a threshold 
remove_rare_taxa <- function(taxa) { 
  
  # Find mean relative abundance for each taxon across all samples
  taxa_mean_abundance <- taxa_sums(taxa) / nsamples(taxa)
  
  # Identify taxa above the threshold
  abundant_taxa <- taxa_names(taxa)[taxa_mean_abundance >= 2]
  
  # Prune dataset to keep only abundant taxa
  physeq_filtered <- prune_taxa(abundant_taxa, taxa)
  
  return(physeq_filtered)
}

phy_pt_abund <- remove_rare_taxa(gloom_pt_fam)
phy_pr_abund <- remove_rare_taxa(gloom_pr_fam)

phy_et_abund <- remove_rare_taxa(gloom_et_fam)
phy_er_abund <- remove_rare_taxa(gloom_er_fam)

#Custom Palletes ----

#for prokaryotes
tax_colors1 <- c(
  "Flavobacteriaceae" = "#5A9BD5", 
  "Coxiellaceae" = "#4CAF50",       
  "Beijerinckiaceae" = "#FDCB82",   
  "Mycoplasmataceae" = "#D36B6D",  
  "Colwelliaceae" = "#FF8C94",     
  "Spongiibacteraceae" = "#9F5F52", 
  "Pirellulaceae" = "#7A5CFd",      
  "Rubritaleaceae" = "#E9B44C",    
  "Paracoccaceae" = "#2C6E49")

tax_colors2 <- c(
  "AEGEAN-169 marine group" = "#FDCB82",
  "Flavobacteriaceae" = "#5A9BD5",
  "HOC36 UC" = "#FF8C94",
  "Magnetospiraceae" = "#4CAF50",
  "Methylophilaceae" = "#7E428B",
  "Nitrincolaceae" = "#E9B44C",
  "Nitrosopumilaceae" = "#5B9AA0",
  "Paracoccaceae" = "#C98C5E",
  "Porticoccaceae" = "#A7C7E7",
  "SAR11 Clade I" = "#9B2D30",
  "SAR11 Clade II" = "#8A9A5B",
  "SAR86 clade" = "#D5A6BD",
  "Thioglobaceae" = "#1F6F70")

#For MicroEuk
custom_palette <- c(brewer.pal(8, "Dark2"),"#9E0142", "#D53E4F", "#F46D43", 
                    "#3288BD","#FEE08B", "#66C2A5", 
                    "#FFFFBF", "#E6F598", "#ABDDA4", 
                    "#5E4FA2", "#FF8C94", "#9F5F52")


#Functions Plots ----
#no colors or order
area_plot1 <- function (phy_gloom, graphName, pallete){
  
  # Transform phyloseq in dataframe
  df <- psmelt (phy_gloom)
  
  # Group by relevant variables and calculate mean abundance
  df_summary <- df %>%
    group_by(Family, Sample, Date_Average) %>%
    summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>%
    ungroup() %>%
    rename(Family = Family)
  
  #Date format
  df_summary <- df_summary %>% mutate(Date_Average = as.Date(Date_Average, format = "%Y-%m-%d"))
  
  #Create labels
  x2lns <- data.frame(
    xstart = as.Date(c("05-Jul-2016", "05-Jan-2017", "05-Jan-2018", "05-Jan-2019",
                       "05-Jan-2020", "05-Jan-2021"), format="%d-%b-%Y"),
    xend = as.Date(c("20-Dec-2016", "20-Dec-2017", "20-Dec-2018", "20-Dec-2019",
                     "20-Dec-2020", "12-Jun-2021"), format="%d-%b-%Y")
  )
  
  x2lbls <- x2lns |>
    rowwise() |>
    summarize(x=mean(c(xstart,xend))) |>
    ungroup() |>
    mutate(lbl=c("2016", "2017", "2018", "2019", "2020", "2021"))
  
  #plot
  plot_area <- ggplot(df_summary, aes(x = Date_Average, y = Abundance)) +
    geom_area(aes(fill = Family, color = Family), alpha = 0.7) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_fill_manual(values = custom_palette) +
    scale_color_manual(values = custom_palette) +
    guides(fill = guide_legend(keyheight = 1, keywidth = 1)) +
    coord_cartesian(ylim=c(0,100),clip="off") +
    labs(subtitle = graphName , x = "", y = "Relative Abundance (%)") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = as.Date(c("2016-07-01", "2021-06-15")), expand = c(0, 0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.margin=margin(t=1,b=1,unit="lines"),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.margin = margin(l = -10)) +
    geom_segment(data=x2lns,mapping=aes(x=xstart,xend=xend),y=-18,yend=-18,
                 linewidth=0.6) +
    geom_text(data=x2lbls,mapping=aes(x=x,label=lbl),y=-22,
              size=10/.pt,fontface="bold")
} 

#For each plot
area_plot_fam_RAS <- function (phy_gloom, graphName, cord, pallete){
  
  # Transform phyloseq in dataframe
  df <- psmelt (phy_gloom)
  
  # Group by relevant variables and calculate mean abundance
  df_summary <- df %>%
    group_by(Family, Sample, Date_Average) %>%
    summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>%
    ungroup() %>%
    rename(Family = Family)
  
  #Date format
  df_summary <- df_summary %>% mutate(Date_Average = as.Date(Date_Average, format = "%Y-%m-%d"))
  
  #Create labels
  x2lns <- data.frame(
    xstart = as.Date(c("05-Jul-2016", "05-Jan-2017", "05-Jan-2018", "05-Jan-2019",
                       "05-Jan-2020"), format="%d-%b-%Y"),
    xend = as.Date(c("20-Dec-2016", "20-Dec-2017", "20-Dec-2018", "20-Dec-2019",
                     "20-Sep-2020"), format="%d-%b-%Y")
  )
  
  x2lbls <- x2lns |>
    rowwise() |>
    summarize(x=mean(c(xstart,xend))) |>
    ungroup() |>
    mutate(lbl=c("2016", "2017", "2018", "2019", "2020"))
  
  # Custom order
  df_summary$Family <- factor(df_summary$Family, levels = cord)
  
  #plot
  plot_area <- ggplot(df_summary, aes(x = Date_Average, y = Abundance)) +
    geom_area(aes(fill = Family, color = Family), alpha = 0.7) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_fill_manual(values = pallete) +
    scale_color_manual(values = pallete) +
    guides(fill = guide_legend(keyheight = 1, keywidth = 1)) +
    coord_cartesian(ylim=c(0,100),clip="off") +
    labs(subtitle = graphName , x = "", y = "Relative Abundance (%)") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = as.Date(c("2016-07-01", "2020-09-31")), expand = c(0, 0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.margin=margin(t=1,b=1,unit="lines"),
          legend.text = element_text(size = 9),
          legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.margin = margin(l = -1)) +
    geom_segment(data=x2lns,mapping=aes(x=xstart,xend=xend),y=-18,yend=-18,
                 linewidth=0.6) +
    geom_text(data=x2lbls,mapping=aes(x=x,label=lbl),y=-22,
              size=10/.pt,fontface="bold")
  
}

#Custom orders ----
corder_er <- c("Chaetocerotaceae", "Thalassiosiraceae",
               "Ephelotidae", "Phaeocystaceae", "Strombidiidae", 
               "RAD-C ucb", "Suessiaceae", "Dinophyceae UC", "Gymnodiniaceae",
               "Bacillariaceae", "Dino-I-1")

corder_pr <- c("AEGEAN-169 marine group", "Magnetospiraceae","Paracoccaceae",
               "SAR11 Clade I","SAR11 Clade II", 
               "HOC36 UC", "Methylophilaceae", "Nitrincolaceae", 
               "Porticoccaceae", "SAR86 clade", "Thioglobaceae", 
               "Flavobacteriaceae", "Nitrosopumilaceae")

#Run the Plots ----

#Surface Water
area_er <- area_plot_fam_RAS(phy_er_abund, "Microbial Eukaryotes",
                          corder_er, custom_palette)

area_pr <- area_plot_fam_RAS(phy_pr_abund, "Prokaryotes",
                              corder_pr, tax_colors2)

comb_arear <- (area_er / area_pr)
comb_arear

ggsave("AREA_RelAb_RAS.png", plot = comb_arear, width = 8.5, height = 10, units = "in")

#Particle-associated
area_plot_fam <- function (phy_gloom, graphName, cord, pallete){
  
  # Transform phyloseq in dataframe
  df <- psmelt (phy_gloom)
  
  # Group by relevant variables and calculate mean abundance
  df_summary <- df %>%
    group_by(Family, Sample, Date_Average) %>%
    summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>%
    ungroup() %>%
    rename(Family = Family)
  
  #Date format
  df_summary <- df_summary %>% mutate(Date_Average = as.Date(Date_Average, format = "%Y-%m-%d"))
  
  #Create labels
  x2lns <- data.frame(
    xstart = as.Date(c("05-Jul-2016", "05-Jan-2017", "05-Jan-2018", "05-Jan-2019",
                       "05-Jan-2020", "05-Jan-2021"), format="%d-%b-%Y"),
    xend = as.Date(c("20-Dec-2016", "20-Dec-2017", "20-Dec-2018", "20-Dec-2019",
                     "20-Dec-2020", "12-Jun-2021"), format="%d-%b-%Y")
  )
  
  x2lbls <- x2lns |>
    rowwise() |>
    summarize(x=mean(c(xstart,xend))) |>
    ungroup() |>
    mutate(lbl=c("2016", "2017", "2018", "2019", "2020", "2021"))
  
  # Custom order
  df_summary$Family <- factor(df_summary$Family, levels = cord)
  
  #plot
  plot_area <- ggplot(df_summary, aes(x = Date_Average, y = Abundance)) +
    geom_area(aes(fill = Family, color = Family), alpha = 0.7) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_fill_manual(values = pallete) +
    scale_color_manual(values = pallete) +
    guides(fill = guide_legend(keyheight = 1, keywidth = 1)) +
    coord_cartesian(ylim=c(0,100),clip="off") +
    labs(subtitle = graphName , x = "", y = "Relative Abundance (%)") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = as.Date(c("2016-07-01", "2021-06-15")), expand = c(0, 0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.margin=margin(t=1,b=1,unit="lines"),
          legend.text = element_text(size = 9),
          legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.margin = margin(l = -10)) +
    geom_segment(data=x2lns,mapping=aes(x=xstart,xend=xend),y=-18,yend=-18,
                 linewidth=0.6) +
    geom_text(data=x2lbls,mapping=aes(x=x,label=lbl),y=-22,
              size=10/.pt,fontface="bold")
}

corder_et <- c("Acantharea_C3", "Chaetocerotaceae",
               "Cryothecomonas-lineage", "Protaspa-lineage", 
               "MAST-9", "Chloropicaceae",
               "Microsporida uc","Pirsoniales uc", "Suessiaceae", 
               "Dinophyceae UC", "Gymnodiniales UC","Dino-I-1", 
               "Dino-I-2", "Dino-I-3")

corder_pt <- c("Beijerinckiaceae", "Paracoccaceae", "Rubritaleaceae",
               "Pirellulaceae",
               "Colwelliaceae", "Coxiellaceae", "Spongiibacteraceae",
               "Mycoplasmataceae","Flavobacteriaceae")

#Plots
area_et <- area_plot_fam(phy_et_abund, "Microbial Eukaryotes",
                         corder_et,custom_palette)
area_pt <- area_plot_fam(phy_pt_abund, "Prokaryotes",
                      corder_pt, tax_colors1)

#Grid, 2 plots together
library(patchwork)

comb_areat <- (area_et / area_pt) 
comb_areat

ggsave("AREA_RelAb_TRAP.png", plot = comb_areat, width = 8.5, height = 10, units = "in")


##### For Katja

comb_eukr <- (area_er2 / area_er) +
  theme(legend.position = "right")

comb_eukt <- (area_et2 / area_et) +
  theme(legend.position = "right")

area_er2

ggsave("AREA_RelAb_MicroEukRAS.png", plot = comb_eukr, width = 9, height = 7, units = "in")
ggsave("AREA_RelAb_MicroEukTrap.png", plot = comb_eukt, width = 9, height = 7, units = "in")

# ---------- Inside the Family ----

#Flavobacteriaceae
flavo1 <- subset_taxa(phy_pras, Family == "Flavobacteriaceae")
flavo1.gen <- tax_glom(flavo1, taxrank = "Genus")

the5most(flavo1.gen, "Genus")

flavo2 <- subset_taxa(phy_ptrap, Family == "Flavobacteriaceae")
flavo2.gen <- tax_glom(flavo2, taxrank = "Genus")

the5most(flavo2.gen, "Genus")


#Paracoccaceae
para1 <- subset_taxa(phy_pras, Family == "Paracoccaceae")
para1.gen <- tax_glom(para1, taxrank = "Genus")

the5most(para1.gen, "Genus")

para2 <- subset_taxa(phy_ptrap, Family == "Paracoccaceae")
para2.gen <- tax_glom(para2, taxrank = "Genus")

the5most(para2.gen, "Genus")

#Bacillariaceae
baci1 <- subset_taxa(phy_eras, Family == "Bacillariaceae")
baci1.gen <- tax_glom(baci1, taxrank = "Genus")

the5most(baci1.gen, "Genus")










