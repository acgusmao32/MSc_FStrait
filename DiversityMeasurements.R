#Diversity Measurements
#MSc - FStrait - AnaGusmao

library("ggplot2")
library("gridExtra")
library("patchwork")
library("phyloseq")
library("tidyr")
library("dplyr")
library("ggpubr")

#----------------------- LOADING PHYLOSEQ OBJECTS

load("p_rar_alph.RData")
load("e_rar_alph.RData")

#-------- ALPHA DIVERSITY METRICS ----

alpha.p <- estimate_richness(p.rar.alph, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
alpha.e <- estimate_richness(e.rar.alph, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))

#export as a table
write.csv(alpha.p,"alpha_prok.csv")
write.csv(alpha.e,"alpha_euk.csv")

#Add the metrics to the metadata
sample_data(p.rar.alph) <- cbind(sample_data(p.rar.alph), alpha.p)
sample_data(e.rar.alph) <- cbind(sample_data(e.rar.alph), alpha.e)

#Summary
summary_fun <- function (alpha_metrics_file) {
  
  the_sum <- alpha_metrics_file %>%
    summarise(
      Observed_min = min(Observed, na.rm = TRUE),
      Observed_max = max(Observed, na.rm = TRUE),
      Observed_mean = mean(Observed, na.rm = TRUE),
      Observed_sd = sd(Observed, na.rm = TRUE),
      
      Chao1_min = min(Chao1, na.rm = TRUE),
      Chao1_max = max(Chao1, na.rm = TRUE),
      Chao1_mean = mean(Chao1, na.rm = TRUE),
      Chao1_sd = sd(Chao1, na.rm = TRUE),
      
      Shannon_min = min(Shannon, na.rm = TRUE),
      Shannon_max = max(Shannon, na.rm = TRUE),
      Shannon_mean = mean(Shannon, na.rm = TRUE),
      Shannon_sd = sd(Shannon, na.rm = TRUE),
      
      InvSimpson_min = min(InvSimpson, na.rm = TRUE),
      InvSimpson_max = max(InvSimpson, na.rm = TRUE),
      InvSimpson_mean = mean(InvSimpson, na.rm = TRUE),
      InvSimpson_sd = sd(InvSimpson, na.rm = TRUE)
    )
  
  return(the_sum)
  
}

alphap_summary <- summary_fun (alpha.p) 
alphap_summary

alphae_summary <- summary_fun (alpha.e) 
alphae_summary

#-------- EDIT METADATA ---- 

met.prok <- sample_data(p.rar.alph) 
met.prok <- subset(met.prok, select = 
                     c(Equipment, Habitat, Date_Average, 
                       Year, Season))
alpha.p.merged <- merge(alpha.p, met.prok, by = "row.names")


met.euk <- sample_data(e.rar.alph) 
met.euk <- subset(met.euk, select = c(Equipment, Habitat, 
                                      Date_Average, Year, Season))
alpha.e.merged <- merge(alpha.e, met.euk, by = "row.names")

#Defining colors Habitat
col.hab <- c("Particle-Attached" = "rosybrown",
             "Surface Water Column" = "cornflowerblue")

#-------- Collect the Season and Years in the correct way
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

# Update the metadata with Season and Year in the correct order
organize_taxa_SeasonYear <- function(taxa){
  
  # Organize the data into Season+Year
  taxa$Season_Year <- mapply(organize_data_SeasonYear,as.Date(taxa$Date_Average), taxa$Season)
  
  #Order the Season and Year
  taxa$Season_Year <- factor(taxa$Season_Year, 
                            levels = c("Sum 16","Aut 16","Win 16","Spr 17","Sum 17","Aut 17",
                                       "Win 17","Spr 18","Sum 18","Aut 18","Win 18","Spr 19",
                                       "Sum 19","Aut 19","Win 19","Spr 20","Sum 20","Aut 20",
                                       "Win 20","Spr 21","Sum 21","Aut 21","Win 21"))
  # Return the organized data
  return(taxa)
}

#-------- Plot functions - TimeS ----

#TEST ALL COMPARISONS
plot_allcomp <- function(taxa, metric, plotname){
  
  # First organize the data into Season and Year
  orgTaxa <- organize_taxa_SeasonYear(taxa)
  
  # Now the Water part
  orgTaxa$Habitat <- factor(orgTaxa$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
  
  #Things to compare
  comp.new.surf.all <- list( c("Sum 17", "Sum 18"), c("Sum 17", "Sum 19"), c("Sum 17", "Sum 20"), #summer
                             c("Sum 18", "Sum 19"), c("Sum 18", "Sum 20"), c("Sum 19", "Sum 20"), #summer
                             c("Aut 16", "Aut 17"), c("Aut 16", "Aut 18"), c("Aut 16", "Aut 19"), #autumn
                             c("Aut 16", "Aut 20"), c("Aut 17", "Aut 18"), c("Aut 17", "Aut 19"), #autumn
                             c("Aut 17", "Aut 20"), c("Aut 18", "Aut 19"), c("Aut 18", "Aut 20"), #autumn
                             c("Aut 19", "Aut 20"), #autumn
                             c("Win 16", "Win 17"), c("Win 16", "Win 18"), c("Win 16", "Win 19"), #winter
                             c("Win 17", "Win 18"), c("Win 17", "Win 19"), c("Win 18", "Win 19"), #winter
                             c("Spr 17", "Spr 18"), c("Spr 17", "Spr 19"), c("Spr 17", "Spr 20"), #spring
                             c("Spr 18", "Spr 19"), c("Spr 18", "Spr 20"), c("Spr 19", "Spr 20")) #spring
  
  comp.new.part.all <- list( c("Sum 16", "Sum 17"), c("Sum 16", "Sum 20"), c("Sum 17", "Sum 20"), #summer
                             c("Aut 16", "Aut 17"), c("Aut 16", "Aut 18"), c("Aut 16", "Aut 19"), #autumn
                             c("Aut 16", "Aut 20"), c("Aut 17", "Aut 18"), c("Aut 17", "Aut 19"), #autumn
                             c("Aut 17", "Aut 20"), c("Aut 18", "Aut 19"), c("Aut 18", "Aut 20"), #autumn
                             c("Aut 19", "Aut 20"), #autumn
                             c("Win 16", "Win 17"), c("Win 16", "Win 18"), c("Win 16", "Win 19"), #winter
                             c("Win 16", "Win 20"), c("Win 17", "Win 18"), c("Win 17", "Win 19"), #winter
                             c("Win 17", "Win 20"), c("Win 18", "Win 19"), c("Win 18", "Win 20"), #winter
                             c("Win 19", "Win 20"), #winter
                             c("Spr 17", "Spr 18"), c("Spr 17", "Spr 20"), c("Spr 17", "Spr 21"), #spring
                             c("Spr 18", "Spr 20"), c("Spr 18", "Spr 21"), c("Spr 20", "Spr 21")) #spring)
  
  
  # Now plot the organized data
  plotDone <- ggplot(orgTaxa, aes(x = Season_Year, y = .data[[metric]], fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.new.surf.all, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.new.part.all, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Particle-Associated")) +
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = plotname, y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  # Return the plot
  return(plotDone)
}

combinePlots1 <- function(taxa,graphName){
  # Do the Obs and Chao plots
  plot1 <- plot_allcomp(taxa, "Observed", "Observed ASVs")
  plot2 <- plot_allcomp(taxa, "Chao1", "Chao1")
  
  # Grid
  combPlot <- (plot1 | plot2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add nice things
  finalPlot <- combPlot + plot_annotation(
    subtitle = graphName) &
    theme(plot.tag = element_text (size = 10))
  
  # Return the nice plot
  return(finalPlot)
} #Obs and Chao1
combinePlots2 <- function(taxa,graphName){
  # Do the Obs and Chao plots
  plot1 <- plot_allcomp(taxa, "Shannon", "Shannon")
  plot2 <- plot_allcomp(taxa, "InvSimpson", "Inverse Simpson")
  
  # Grid
  combPlot <- (plot1 | plot2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add nice things
  finalPlot <- combPlot + plot_annotation(
    title = 'Alpha Diversity Metrics',
    subtitle = graphName) &
    theme(plot.tag = element_text (size = 10))

  # Return the nice plot
  return(finalPlot)
} #Sha and InvSimp

combinePlots1(alpha.p.merged,"Prokaryotes")
combinePlots2(alpha.p.merged,"Prokaryotes")

combinePlots1(alpha.e.merged,"MicroEukaryotes")
combinePlots2(alpha.e.merged,"MicroEukaryotes")


#ONLY SIGNIFICATIVE 

#Prok
plotp_Observed <- function(taxa){
  # First organize the data into Season and Year
  orgTaxa <- organize_taxa_SeasonYear(taxa)
  
  # Now the Water part
  orgTaxa$Habitat <- factor(orgTaxa$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
  
  #Things to compare
  comp.new.surf.prok <- list( c("Win 17", "Win 19"))
  comp.new.part.prok <- list( c("Win 16", "Win 18"), c("Win 18", "Win 20"))
  
  # Now plot the organized data
  plotDone <- ggplot(orgTaxa, aes(x = Season_Year, y = Observed, fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.new.surf.prok, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.new.part.prok, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Particle-Associated")) +
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = "Observed ASVs", y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  # Return the plot
  return(plotDone)
}
plotp_Chao1 <- function(taxa){
  # First organize the data into Season and Year
  orgTaxa <- organize_taxa_SeasonYear(taxa)
  
  # Now the Water part
  orgTaxa$Habitat <- factor(orgTaxa$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
  
  #Things to compare
  comp.new.surf.prok <- list( c("Win 17", "Win 19"), c("Sum 17", "Sum 19"))
  comp.new.part.prok <- list( c("Win 16", "Win 18"))
  
  # Now plot the organized data
  plotDone <- ggplot(orgTaxa, aes(x = Season_Year, y = Chao1, fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.new.surf.prok, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.new.part.prok, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Particle-Associated")) +
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = "Chao1", y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  # Return the plot
  return(plotDone)
}
plotp_Shannon <- function(taxa){
  # First organize the data into Season and Year
  orgTaxa <- organize_taxa_SeasonYear(taxa)
  
  # Now the Water part
  orgTaxa$Habitat <- factor(orgTaxa$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
  
  #Things to compare
  comp.new.surf.prok <- list( c("Aut 16", "Aut 20"), c("Win 17", "Win 19"), c("Win 18", "Win 19"))
  comp.new.part.prok <- list( c ("Win 16", "Win 18"))
  
  # Now plot the organized data
  plotDone <- ggplot(orgTaxa, aes(x = Season_Year, y = Shannon, fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.new.surf.prok, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.new.part.prok, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Particle-Associated")) +
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = "Shannon", y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  # Return the plot
  return(plotDone)
}
plotp_InvSimp <- function(taxa){
  # First organize the data into Season and Year
  orgTaxa <- organize_taxa_SeasonYear(taxa)
  
  # Now the Water part
  orgTaxa$Habitat <- factor(orgTaxa$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
  
  #Things to compare
  comp.new.surf.prok <- list( c("Win 17", "Win 19"), c("Win 18", "Win 19"))
  comp.new.part.prok <- list( c("Win 16", "Win 18"), c("Win 16", "Win 19"))
  
  # Now plot the organized data
  plotDone <- ggplot(orgTaxa, aes(x = Season_Year, y = InvSimpson, fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.new.surf.prok, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.new.part.prok, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Particle-Associated")) +
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = "Inverse Simpson Diversity", y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  # Return the plot
  return(plotDone)
}

combinePlotsP1 <- function(taxa,graphName){
  # Do the Obs and Chao plots
  plot1 <- plotp_Observed(taxa)
  plot2 <- plotp_Chao1(taxa)
  
  # Grid
  combPlot <- (plot1 | plot2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add nice things
  finalPlot <- combPlot + plot_annotation(
    subtitle = graphName) &
    theme(plot.tag = element_text (size = 10))
  
  # Return the nice plot
  return(finalPlot)
} #Obs and Chao1
combinePlotsP2 <- function(taxa,graphName){
  # Do the Obs and Chao plots
  plot1 <- plotp_Shannon(taxa)
  plot2 <- plotp_InvSimp(taxa)
  
  # Grid
  combPlot <- (plot1 | plot2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add nice things
  finalPlot <- combPlot + plot_annotation(
    subtitle = graphName) &
    theme(plot.tag = element_text (size = 10))
  
  # Export the plot (Big PNG)
  fileName1 <- paste0("AlphaDiv_",graphName,"_Sha_Inv_Gross.png")
  ggsave(fileName1, plot = finalPlot, width = 16, height = 9, units = "in", dpi = 300)
  
  # Export the plot (Small PNG)
  fileName3 <- paste0("AlphaDiv_",graphName,"_Sha_Inv_Klein.png")
  ggsave(fileName3, plot = finalPlot, width = 8.5, height = 7, units = "in", dpi = 300)
  
  # Return the nice plot
  return(finalPlot)
} #Sha and InvSimp

#MicroEuk
plote_Observed <- function(taxa){
  # First organize the data into Season and Year
  orgTaxa <- organize_taxa_SeasonYear(taxa)
  
  # Now the Water part
  orgTaxa$Habitat <- factor(orgTaxa$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
  
  #Things to compare
  comp.new.surf.euk <- list( c("Sum 17", "Sum 19"))
  comp.new.part.euk <- list( c("Win 18", "Win 19"))
  
  # Now plot the organized data
  plotDone <- ggplot(orgTaxa, aes(x = Season_Year, y = Observed, fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.new.surf.euk, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.new.part.euk, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Particle-Associated")) +
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = "Observed ASVs", y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  # Return the plot
  return(plotDone)
}
plote_Chao1 <- function(taxa){
  # First organize the data into Season and Year
  orgTaxa <- organize_taxa_SeasonYear(taxa)
  
  # Now the Water part
  orgTaxa$Habitat <- factor(orgTaxa$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
  
  #Things to compare
  comp.new.surf.euk <- list( c("Spr 17", "Spr 20"), c("Spr 17", "Spr 18"))
  comp.new.part.euk <- list( c("Win 18", "Win 19"))
  
  # Now plot the organized data
  plotDone <- ggplot(orgTaxa, aes(x = Season_Year, y = Chao1, fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.new.surf.euk, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.new.part.euk, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Particle-Associated")) +
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = "Chao1", y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  # Return the plot
  return(plotDone)
}
plote_Shannon <- function(taxa){
  # First organize the data into Season and Year
  orgTaxa <- organize_taxa_SeasonYear(taxa)
  
  # Now the Water part
  orgTaxa$Habitat <- factor(orgTaxa$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
  
  #Things to compare
  comp.new.surf.euk <- list( c("Win 18", "Win 19"), c("Win 17", "Win 18"), c("Win 16", "Win 19"),
                             c("Win 16", "Win 17"), c("Sum 17", "Sum 19"))
  comp.new.part.euk <- list()
  
  # Now plot the organized data
  plotDone <- ggplot(orgTaxa, aes(x = Season_Year, y = Shannon, fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.new.surf.euk, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.new.part.euk, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Particle-Associated")) +
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = "Shannon", y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  # Return the plot
  return(plotDone)
}
plote_InvSimp <- function(taxa){
  # First organize the data into Season and Year
  orgTaxa <- organize_taxa_SeasonYear(taxa)
  
  # Now the Water part
  orgTaxa$Habitat <- factor(orgTaxa$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
  
  #Things to compare
  comp.new.surf.euk <- list( c("Win 18", "Win 19"), c("Win 17", "Win 18"), c("Sum 17", "Sum 19"))
  comp.new.part.euk <- list( c("Win 18", "Win 20"))
  
  # Now plot the organized data
  plotDone <- ggplot(orgTaxa, aes(x = Season_Year, y = InvSimpson, fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.new.surf.euk, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.new.part.euk, method = "wilcox.test",
                       label = "p.signif", hide.ns = FALSE, data = subset(orgTaxa, Habitat == "Particle-Associated")) +
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = "Inverse Simpson Diversity", y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  # Return the plot
  return(plotDone)
}

combinePlotsE1 <- function(taxa,graphName){
  # Do the Obs and Chao plots
  plot1 <- plote_Observed(taxa)
  plot2 <- plote_Chao1(taxa)
  
  # Grid
  combPlot <- (plot1 | plot2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add nice things
  finalPlot <- combPlot + plot_annotation(
    subtitle = graphName) &
    theme(plot.tag = element_text (size = 10))
  
  # Return the nice plot
  return(finalPlot)
} #Obs and Chao1
combinePlotsE2 <- function(taxa,graphName){
  # Do the Obs and Chao plots
  plot1 <- plote_Shannon(taxa)
  plot2 <- plote_InvSimp(taxa)
  
  # Grid
  combPlot <- (plot1 | plot2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add nice things
  finalPlot <- combPlot + plot_annotation(
    subtitle = graphName) &
    theme(plot.tag = element_text (size = 10))
  
  # Export the plot (Big PNG)
  fileName1 <- paste0("AlphaDiv_",graphName,"_Sha_Inv_Gross.png")
  ggsave(fileName1, plot = finalPlot, width = 16, height = 9, units = "in", dpi = 300)
  
  # Export the plot (Small PNG)
  fileName3 <- paste0("AlphaDiv_",graphName,"_Sha_Inv_Klein.png")
  ggsave(fileName3, plot = finalPlot, width = 8.5, height = 7, units = "in", dpi = 300)
  
  # Return the nice plot
  return(finalPlot)
} #Sha and InvSimp

#Run the plots

#All plots
prokOBS_CHA <- combinePlotsP1(alpha.p.merged,"Prokaryotes")
prokSHA_INVS <- combinePlotsP2(alpha.p.merged,"Prokaryotes")
eukOBS_CHA <- combinePlotsE1(alpha.e.merged,"Microbial Eukaryotes")
eukSHA_INVS <- combinePlotsE2(alpha.e.merged,"Microbial Eukaryotes")

#Export
ggsave("AlphaDiv_Prok_TimeS1.png", plot = prokOBS_CHA, width = 16, height = 9, units = "in", dpi = 300)
ggsave("AlphaDiv_Prok_TimeS2.png", plot = prokSHA_INVS, width = 16, height = 9, units = "in", dpi = 300)
ggsave("AlphaDiv_MicroEuk_TimeS1.png", plot = eukOBS_CHA, width = 16, height = 9, units = "in", dpi = 300)
ggsave("AlphaDiv_MicroEuk_TimeS2.png", plot = eukSHA_INVS, width = 16, height = 9, units = "in", dpi = 300)

#Obs and Sha
combinePlotsP3 <- function(taxa,graphName){
  # Do the Obs and Chao plots
  plot1 <- plotp_Observed(taxa)
  plot2 <- plotp_Shannon(taxa)
  
  # Grid
  combPlot <- (plot1 | plot2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add nice things
  finalPlot <- combPlot + plot_annotation(
    subtitle = graphName) &
    theme(plot.tag = element_text (size = 10))
  
  # Return the nice plot
  return(finalPlot)
} 
combinePlotsE3 <- function(taxa,graphName){
  # Do the Obs and Chao plots
  plot1 <- plote_Observed(taxa)
  plot2 <- plote_Shannon(taxa)
  
  # Grid
  combPlot <- (plot1 | plot2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add nice things
  finalPlot <- combPlot + plot_annotation(
    subtitle = graphName) &
    theme(plot.tag = element_text (size = 10))
  
  
  
  # Return the nice plot
  return(finalPlot)
} 

prok3 <- combinePlotsP3(alpha.p.merged,"Prokaryotes")
euk3 <- combinePlotsE3(alpha.e.merged,"Microbial Eukaryotes")

ggsave("AlphaDiv_Prok_TimeS3_Klein.png", plot = prok3, width = 8.5, height = 7, units = "in", dpi = 300)
ggsave("AlphaDiv_MicroEuk_TimeS3_Klein.png", plot = euk3, width = 8.5, height = 7, units = "in", dpi = 300)


#-------- NORMALITY TESTS ----
#Shapiro Test
#Assumes the data follow a normal distribution
#p < 0.05 = diverges from normal
#W close to 1 indicate is close of a normal distribution

#QQplot -> straight line = theoretical normal distribution

#[1] Prokaryotes
shapiro.test(alpha.p$Observed) 
ggplot(alpha.p, aes(sample = Observed)) + stat_qq() + stat_qq_line() + #QQplot 
  labs(title = "Q-Q Plot for Normality", x = "Theoretical Quantiles", y = "Sample Quantiles")
#Result: not normal

shapiro.test(alpha.p$Chao1) 
ggplot(alpha.p, aes(sample = Chao1)) + stat_qq() + stat_qq_line() + #QQplot
  labs(title = "Q-Q Plot for Normality", x = "Theoretical Quantiles", y = "Sample Quantiles")
#Result: not normal

shapiro.test(alpha.p$Shannon) 
ggplot(alpha.p, aes(sample = Shannon)) + stat_qq() + stat_qq_line() + #QQplot
  labs(title = "Q-Q Plot for Normality", x = "Theoretical Quantiles", y = "Sample Quantiles")
#Result: not normal

shapiro.test(alpha.p$InvSimpson)
ggplot(alpha.p, aes(sample = InvSimpson)) + stat_qq() + stat_qq_line() + #QQplot
  labs(title = "Q-Q Plot for Normality", x = "Theoretical Quantiles", y = "Sample Quantiles")
#Result: not normal

#[2] Eukaryotes
shapiro.test(alpha.e$Observed) 
ggplot(alpha.e, aes(sample = Observed)) + stat_qq() + stat_qq_line() + #QQplot
  labs(title = "Q-Q Plot for Normality", x = "Theoretical Quantiles", y = "Sample Quantiles")
#Result: not normal

shapiro.test(alpha.e$Chao1) 
ggplot(alpha.e, aes(sample = Chao1)) + stat_qq() + stat_qq_line() + #QQplot
  labs(title = "Q-Q Plot for Normality", x = "Theoretical Quantiles", y = "Sample Quantiles")
#Result: not normal

shapiro.test(alpha.e$Shannon)
ggplot(alpha.e, aes(sample = Shannon)) + stat_qq() + stat_qq_line() + #QQplot
  labs(title = "Q-Q Plot for Normality", x = "Theoretical Quantiles", y = "Sample Quantiles")
#Result: not normal

shapiro.test(alpha.e$InvSimpson) 
ggplot(alpha.e, aes(sample = InvSimpson)) + stat_qq() + stat_qq_line() + #QQplot
  labs(title = "Q-Q Plot for Normality", x = "Theoretical Quantiles", y = "Sample Quantiles")
#Result: not normal

###---------- SEASONS COLLAPSED ----

#All pairwise comparisons using Wilcox Test
alpha.p.merged$Habitat <- factor(alpha.p.merged$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))
alpha.e.merged$Habitat <- factor(alpha.e.merged$Habitat, levels = c("Free-living Surface Water", "Particle-Associated"))

# Plot function
# Correction for multiple wilcox test -> Benjamini-Hochberg (BH)

alpha_sea_collap.compall <- function(df, variable, plot_name){
  
  #Things to compare
  comp.all <- list( c("Winter", "Spring"), c("Winter", "Summer"), 
                    c("Winter", "Autumn"), c("Spring", "Autumn"),
                    c("Spring", "Summer"), c("Summer", "Autumn"))
  
  #plot
  plot_byseason <- ggplot(df, aes(x = Season, y = .data[[variable]], fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.all, method = "wilcox.test", p.adjust.method = "BH", 
                       label = "p.signif", hide.ns = FALSE, data = subset(df, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.all, method = "wilcox.test", p.adjust.method = "BH",
                       label = "p.signif", hide.ns = FALSE, data = subset(df, Habitat == "Particle-Associated")) +
    stat_compare_means(label.y = 0) + #global standard = Kruskal-Wallis
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = plot_name, y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position ="none",
          strip.text = element_text(size = 10))
  
  return(plot_byseason)
}

#Prok
obs_allp <- alpha_sea_collap.compall(alpha.p.merged, "Observed", "Observed ASVs")
cha_allp <- alpha_sea_collap.compall(alpha.p.merged, "Chao1", "Chao1")
sha_allp <- alpha_sea_collap.compall(alpha.p.merged, "Shannon", "Shannon")
ivs_allp <- alpha_sea_collap.compall(alpha.p.merged, "InvSimpson", "Inverse Simpson")

#Grid
combined.allp <- (obs_allp | cha_allp | sha_allp | ivs_allp) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

byseasons.allp <- combined.allp + plot_annotation(
  title = 'Alpha Diversity Metrics',
  subtitle = 'Prokaryotes') &
  theme(plot.tag = element_text (size = 10))

byseasons.allp

#MicroEuk
obs_alle <- alpha_sea_collap.compall(alpha.e.merged, "Observed", "Observed ASVs")
cha_alle <- alpha_sea_collap.compall(alpha.e.merged, "Chao1", "Chao1")
sha_alle <- alpha_sea_collap.compall(alpha.e.merged, "Shannon", "Shannon")
ivs_alle <- alpha_sea_collap.compall(alpha.e.merged, "InvSimpson", "Inverse Simpson")

#Grid
combined.alle <- (obs_alle | cha_alle | sha_alle | ivs_alle) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

byseasons.alle <- combined.alle + plot_annotation(
  title = 'Alpha Diversity Metrics',
  subtitle = 'Microbial Eukaryotes') &
  theme(plot.tag = element_text (size = 10))

byseasons.alle


#-------- Only Significant Comparisons----

#Prokaryotes
comp.obs.p.surf <- list( c("Autumn", "Winter"), c("Autumn", "Summer"), c("Spring", "Summer"),
                         c("Spring", "Winter"), c("Summer", "Winter"))
comp.obs.p.part <- list( c("Summer", "Winter"), c("Spring", "Winter"), c ("Autumn", "Winter"))

comp.cha.p.surf <- list( c("Autumn", "Summer"), c("Autumn", "Winter"), c("Spring", "Summer"), 
                         c("Spring", "Winter"), c("Summer", "Winter"))
comp.cha.p.part <- list( c("Spring", "Winter"), c("Summer", "Winter"), c("Autumn", "Winter"))

comp.sha.p.surf <- list( c("Autumn", "Summer"), c("Autumn", "Winter"), c("Spring", "Winter"), 
                         c("Summer", "Winter"))
comp.sha.p.part <- list()

comp.ivs.p.surf <- list( c("Autumn", "Summer"), c("Autumn", "Winter"), c("Spring", "Winter"),
                          c("Summer", "Winter"))
comp.ivs.p.part <- list()

#Microbial Eukaryotes
comp.obs.e.surf <- list( c("Autumn", "Winter"), c("Spring", "Summer"), c("Autumn", "Spring"))
comp.obs.e.part <- list( c("Autumn", "Summer"), c("Autumn", "Spring"), c("Spring", "Winter"),
                         c("Summer", "Winter"))

comp.cha.e.surf <- list( c("Autumn", "Spring"), c("Autumn", "Winter"), c("Spring", "Summer"), 
                         c("Summer", "Winter"))
comp.cha.e.part <- list( c("Autumn", "Spring"), c("Spring", "Winter"), c("Summer", "Winter"))

comp.sha.e.surf <- list() 
comp.sha.e.part <- list( c("Summer", "Winter"))


comp.ivs.e.surf <- list()
comp.ivs.e.part <- list()

#-------- Plots Season Collapsed ---- 
#p-value label size = 3.2 (small) and 3.8 (big)

alpha_sea_collap <- function(df, variable, comp.surf, comp.part, plot_name){
  
  plot_byseason <- ggplot(df, aes(x = Season, y = .data[[variable]], fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.surf, method = "wilcox.test", p.adjust.method = "BH", 
                       label = "p.signif", hide.ns = FALSE, data = subset(df, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.part, method = "wilcox.test", p.adjust.method = "BH", 
                       label = "p.signif", hide.ns = FALSE, data = subset(df, Habitat == "Particle-Associated")) +
    stat_compare_means(method = "kruskal.test", label = "..p..", label.y = 0, size = 3.2) + #global standard = Kruskal-Wallis
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = plot_name, y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    #customization
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position ="none",
          panel.grid = element_blank(),  
          strip.text = element_text(size = 9))
  
  return(plot_byseason)
}
alpha_sea_collap_withoutlabel <- function(df, variable, comp.surf, comp.part, plot_name){
  
  plot_byseason <- ggplot(df, aes(x = Season, y = .data[[variable]], fill = Season)) +
    geom_boxplot() +
    stat_compare_means(comparisons = comp.surf, method = "wilcox.test", p.adjust.method = "BH", 
                       label = "p.signif", hide.ns = FALSE, data = subset(df, Habitat == "Free-living Surface Water")) +
    stat_compare_means(comparisons = comp.part, method = "wilcox.test", p.adjust.method = "BH", 
                       label = "p.signif", hide.ns = FALSE, data = subset(df, Habitat == "Particle-Associated")) +
    stat_compare_means(method = "kruskal.test", label = "..p..", label.y = 0, size = 3.2) + #global standard = Kruskal-Wallis
    facet_grid(Habitat ~., scales = "free", space = "free_x") +
    theme_bw() +
    labs(title = plot_name, y = "") +
    scale_fill_manual(values = c("Spring" = "olivedrab",
                                 "Summer" = "darkorange", 
                                 "Autumn" = "indianred", 
                                 "Winter" = "turquoise3")) +
    #customization
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position ="none",
          panel.grid = element_blank(),
          strip.text = element_blank())
  
  return(plot_byseason)
}

#Prok
obs.p.season <- alpha_sea_collap_withoutlabel(alpha.p.merged, "Observed", comp.obs.p.surf, comp.obs.p.part,
                                 "Observed ASVs")

cha.p.season <- alpha_sea_collap_withoutlabel(alpha.p.merged, "Chao1", comp.cha.p.surf, comp.cha.p.part, 
                                 "Chao1")


sha.p.season <- alpha_sea_collap_withoutlabel(alpha.p.merged, "Shannon", comp.sha.p.surf, comp.sha.p.part, 
                                 "Shannon")

ivs.p.season <- alpha_sea_collap(alpha.p.merged, "InvSimpson",comp.ivs.p.surf, comp.ivs.p.part,
                                 "Inverse Simpson")
#Grid Prok
combined.p <- (obs.p.season | cha.p.season | sha.p.season | ivs.p.season) +
  plot_layout(guides = "collect", heights = c(1, 1, 1 ,1), widths = c(1, 1, 1 ,1), ncol = 4) &
  theme(legend.position = "bottom")

byseasons.collap.p <- combined.p + plot_annotation(
  subtitle = 'Prokaryotes')

byseasons.collap.p 

#Export the plot
ggsave("AlphaDiv_Prok_Season.png", plot = byseasons.collap.p, width = 16, height = 9, units = "in", dpi = 300)
ggsave("AlphaDiv_Prok_Season_Klein.png", plot = byseasons.collap.p, width = 8.5, height = 7.5, units = "in", dpi = 300)

#MicroEuk
obs.e.season <- alpha_sea_collap_withoutlabel(alpha.e.merged, "Observed", comp.obs.e.surf, comp.obs.e.part,
                                 "Observed ASVs")

cha.e.season <- alpha_sea_collap_withoutlabel(alpha.e.merged, "Chao1", comp.cha.e.surf, comp.cha.e.part,
                                 "Chao1")

sha.e.season <- alpha_sea_collap_withoutlabel(alpha.e.merged, "Shannon", comp.sha.e.surf, comp.sha.e.part,
                                 "Shannon")

ivs.e.season <- alpha_sea_collap(alpha.e.merged, "InvSimpson", comp.ivs.e.surf, comp.ivs.e.part, 
                                 "Inverse Simpson")


#Grid MicroEuk
combined.e <- (obs.e.season | cha.e.season | sha.e.season | ivs.e.season) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

byseasons.collap.e <- combined.e + plot_annotation(
  subtitle = 'Microbial Eukaryotes')

byseasons.collap.e

#Export the plot
ggsave("AlphaDiv_MicroEuk_Season.png", plot = byseasons.collap.e, width = 16, height = 9, units = "in", dpi = 300)
ggsave("AlphaDiv_MicroEuk_Season_Klein.png", plot = byseasons.collap.e, width = 8.5, height = 7, units = "in", dpi = 300)

#Save Image RDATA ----
save.image(file = "alphadiv.RData")
