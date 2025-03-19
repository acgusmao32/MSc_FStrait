#Particle Fluxes
#MSc FStrait - agusmao

library (ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Loading data ----
biogeo <- read.table("MET_fstrait.txt", header = TRUE, sep = "\t")

#Remove RAS samples
biogeo_up <- biogeo[biogeo$Data_Type == "Sediment_Trap", ]

# Ensure Date_begin and Date_end columns are in Date format
biogeo_up$Date_Average <- as.Date(biogeo_up$Date_Average, format = "%Y-%m-%d")
biogeo_up$Date_begin <- as.Date(biogeo_up$Date_begin, format = "%Y-%m-%d")

# Gather all flux columns to long format for easier plotting
biogeo_long <- biogeo_up %>%
  pivot_longer(cols = starts_with("Flux_"), 
               names_to = "FluxType", 
               values_to = "FluxValue")

# ---- Plot - TimeSeries Bar ----

#All
ggplot(biogeo_long, aes(x = Date_Average, y = FluxValue, color = FluxType)) +
  geom_line(linewidth = 1) +  # Line thickness
  labs(y = "Particle Flux (mg/m²/day)", color = "Flux Type") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y",
               limits = as.Date(c("2016-07-01", "2021-06-15")), expand = c(0, 0)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.position = "top")

#Per Flux
indline_upper <- function (biofile, flux){
  
  #Calculate Mean of the flux
  meanflux <- mean(biofile[[flux]], na.rm = TRUE)  
  cat("The mean is:", (meanflux), "\n")
  
  #custom colors
  season.col <- c("Winter" = "turquoise3","Spring" = "olivedrab",
                  "Summer" = "darkorange","Autumn" = "indianred2")
  
  #Plot
  plot <- ggplot(biofile, aes(x = Date_Average, y = .data[[flux]])) +
    geom_rect(aes(xmin = as.Date("2020-01-01"), xmax = as.Date("2020-09-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01, color = NA) + #cold anomaly
    geom_col(aes(fill = Season), width = 13, alpha = 1) +
    geom_hline(yintercept = meanflux, linetype = "dashed", color = "black", linewidth = 0.6) +
    labs(y = "Particle Flux (mg/m²/day)") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y",
                 limits = as.Date(c("2016-07-01", "2021-06-15")), expand = c(0, 0)) +
    scale_y_continuous(breaks = pretty(biofile[[flux]], n = 6)) +
    scale_color_manual(values = season.col) +
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          legend.position = "none")
  
  return (plot)
}
indline <- function (biofile, flux){
  
  #Calculate Mean of the flux
  meanflux <- mean(biofile[[flux]], na.rm = TRUE)  
  cat("The mean is:", (meanflux), "\n")
  
  #custom colors
  season.col <- c("Winter" = "turquoise3","Spring" = "olivedrab",
                  "Summer" = "darkorange","Autumn" = "indianred2")
  
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
  
  #Plot
  plot <- ggplot(biofile, aes(x = Date_Average, y = .data[[flux]])) +
    geom_rect(aes(xmin = as.Date("2020-01-01"), xmax = as.Date("2020-09-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01, color = NA) + #cold anomaly
    geom_col(aes(fill = Season), width = 13, alpha = 1) +
    geom_hline(yintercept = meanflux, linetype = "dashed", color = "black", linewidth = 0.6) +
    labs(y = "Particle Flux (mg/m²/day)") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b",
                 limits = as.Date(c("2016-07-01", "2021-06-15")), expand = c(0, 0)) +
    scale_y_continuous(breaks = pretty(biofile[[flux]], n = 6)) +
    scale_color_manual(values = season.col) +
    coord_cartesian(ylim=c(0,75),clip="off") +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.margin=margin(t=1,b=1.8,unit="lines"),
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          legend.position = "none")+
    geom_segment(data=x2lns,mapping=aes(x=xstart,xend=xend), y=-26,yend=-26,
                 linewidth=0.6) +
    geom_text(data=x2lbls,mapping=aes(x=x,label=lbl),y=-32,
              size=10/.pt,fontface="bold")
  
  return (plot)
}

poc <-indline_upper(biogeo_up, "Flux_POC") #Mean: 6.71
pon <- indline_upper(biogeo_up, "Flux_PON") #Mean: 0.92
carb <- indline_upper(biogeo_up, "Flux_Carbonate") #Mean: 12.66
si <- indline(biogeo_up, "Flux_BothSi") #Mean: 5.27
pom <- indline(biogeo_up, "Flux_PoM") #Mean: 18.09

comb_bio <- (poc / pon / carb / si) +
  plot_layout(guides = "collect", axis_titles = "collect") + 
  plot_annotation(tag_levels = list(c("POC", "PON", "CARB", "SIL"))) & 
  theme(plot.tag.position = c(0.98,0.98))

comb_bio

ggsave("TimeS_FluxParticles.png", plot = comb_bio, width = 16, height = 9, units = "in")
ggsave("TimeS_FluxParticles_Klein2.png", plot = comb_bio, width = 8.5, height = 8, units = "in")

# ---- Comparison between years ----

#Check distribution
#p < 0.05 = not normal
shapiro.test(biogeo_up$Flux_POC[biogeo_up$Year == 2020]) #normal
shapiro.test(biogeo_up$Flux_POC[biogeo_up$Year == 2018]) 
shapiro.test(biogeo_up$Flux_POC[biogeo_up$Year == 2017]) 

shapiro.test(biogeo_up$Flux_Carbonate[biogeo_up$Year == 2020]) #normal
shapiro.test(biogeo_up$Flux_Carbonate[biogeo_up$Year == 2018]) #normal
shapiro.test(biogeo_up$Flux_Carbonate[biogeo_up$Year == 2017]) #normal

shapiro.test(biogeo_up$Flux_BothSi[biogeo_up$Year == 2020])
shapiro.test(biogeo_up$Flux_BothSi[biogeo_up$Year == 2018])
shapiro.test(biogeo_up$Flux_BothSi[biogeo_up$Year == 2016])

shapiro.test(biogeo_up$Flux_PON[biogeo_up$Year == 2020])
shapiro.test(biogeo_up$Flux_PON[biogeo_up$Year == 2018])
shapiro.test(biogeo_up$Flux_PON[biogeo_up$Year == 2017])

#Wilcoxon
wilcox.test(biogeo_up$Flux_POC[biogeo_up$Year == 2020], 
            biogeo_up$Flux_POC[biogeo_up$Year == 2018])

wilcox.test(biogeo_up$Flux_POC[biogeo_up$Year == 2020], 
            biogeo_up$Flux_POC[biogeo_up$Year == 2017])

wilcox.test(biogeo_up$Flux_BothSi[biogeo_up$Year == 2020], 
            biogeo_up$Flux_BothSi[biogeo_up$Year == 2018])

wilcox.test(biogeo_up$Flux_BothSi[biogeo_up$Year == 2020], 
            biogeo_up$Flux_BothSi[biogeo_up$Year == 2016])

wilcox.test(biogeo_up$Flux_PON[biogeo_up$Year == 2020], 
            biogeo_up$Flux_PON[biogeo_up$Year == 2018])

wilcox.test(biogeo_up$Flux_PON[biogeo_up$Year == 2020], 
            biogeo_up$Flux_PON[biogeo_up$Year == 2017])


#T Test
t.test(Flux_Carbonate ~ Year, data = biogeo_up, subset = Year %in% c(2020, 2018),
       var.equal = TRUE)

t.test(Flux_Carbonate ~ Year, data = biogeo_up, subset = Year %in% c(2020, 2017),
       var.equal = TRUE)

t.test(Flux_POC ~ Year, data = biogeo_up, subset = Year %in% c(2020, 2018),
       var.equal = TRUE)

t.test(Flux_POC ~ Year, data = biogeo_up, subset = Year %in% c(2020, 2017),
       var.equal = TRUE)







