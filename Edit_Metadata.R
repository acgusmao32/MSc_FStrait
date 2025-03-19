#Editing the metadata
#MSc FStrait - AnaGusmao

#Modified from DataLoad - M.Wietz et al (2021)

library("dplyr")
library("tidyr")

#----------- LOADING FILES ----

metadata <- read.table("metadata1.txt", header=TRUE, sep="\t")

metadata$Date_begin <- as.Date(metadata$Date_begin, format="%Y-%m-%d")
metadata$Date_end <- as.Date(metadata$Date_end, format="%Y-%m-%d")

#loading Env data with daily resolution
env.daily <- read.table("Daily_Res_Env_RAS_FramStrait.txt", header=TRUE, sep="\t") %>%
  mutate(Date = as.Date(Date, format="%m/%d/%Y"))

#----------- Assign environmental data to samples ----

#New columns with Lag_Date1 (begin) and Lag_Date2 (end)
metadata <- metadata %>% mutate( Lag_Date1 = Date_begin - 4, Lag_Date2 = Date_end - 4)

#Extract the information necessary
lag_met <- metadata %>%
  filter(Data_Type == "Sediment_Trap") %>%
  select(ID_samples, Lag_Date1, Lag_Date2)

ras_met <- metadata %>%
  filter(Data_Type == "RAS_Amplicon") %>%
  select(ID_samples, Date_begin, Date_end)

# Compute average environmental parameters for the period
avg_lag_data <- lag_met %>%
  rowwise() %>%
  mutate(avg_env1 = list(env.daily %>% # Filter daily data within the date range for each sample
        filter(Date >= Lag_Date1 & Date <= Lag_Date2) %>%
          summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))))) %>%
  unnest_wider(avg_env1)  # Expand the list column into separate columns


avg_ras_data <- ras_met %>%
  rowwise() %>%
  mutate(avg_env2 = list(env.daily %>% # Filter daily data within the date range for each sample
                          filter(Date >= Date_begin & Date <= Date_end) %>%
                          summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))))) %>%
  unnest_wider(avg_env2)  # Expand the list column into separate columns


#Delete columns
avg_lag_data <- avg_lag_data %>%
  select(-Depth, -Sigma, -pH)

avg_ras_data <- avg_ras_data %>%
  select(-Sigma, -pH, -Date_begin, -Date_end)

# Merge values of env_parameters
metadata_new <- metadata %>%
  left_join(avg_lag_data, by = "ID_samples") %>%
  select(-Lag_Date1.x, -Lag_Date2.x) 

colnames(metadata_new)[colnames(metadata_new) == "Lag_Date1.y"] <- "Lag_Date1"
colnames(metadata_new)[colnames(metadata_new) == "Lag_Date2.y"] <- "Lag_Date2"

metadata_new2 <- metadata_new %>%
  left_join(avg_ras_data, by = "ID_samples") %>%
  mutate(Temperature = coalesce(Temperature.x, Temperature.y),
    Depth = coalesce(Depth.x, Depth.y), Salinity = coalesce(Salinity.x, Salinity.y),
    Ice_Dist = coalesce(Ice_Dist.x, Ice_Dist.y), Ice_Conc = coalesce(Ice_Conc.x, Ice_Conc.y),
    MLD = coalesce(MLD.x, MLD.y), fraction_AW = coalesce(fraction_AW.x, fraction_AW.y),
    fraction_PW = coalesce(fraction_PW.x, fraction_PW.y), CO2 = coalesce(CO2.x, CO2.y),
    Chl = coalesce(Chl.x, Chl.y), Oxygen_Sat = coalesce(Oxygen_Sat.x, Oxygen_Sat.y), 
    Oxygen_Conc = coalesce(Oxygen_Conc.x, Oxygen_Conc.y)) %>%
  select(-Temperature.x, -Temperature.y, -Depth.x, -Depth.y, -Salinity.x, -Salinity.y,
         -Ice_Dist.x, -Ice_Dist.y, -Ice_Conc.x, -Ice_Conc.y, -MLD.x, -MLD.y,
         -fraction_AW.x, -fraction_AW.y, -fraction_PW.x, -fraction_PW.y, 
         -CO2.x, -CO2.y, -Chl.x, -Chl.y, -Oxygen_Sat.x, -Oxygen_Sat.y,
         -Oxygen_Conc.x, -Oxygen_Conc.y)

metadata <- metadata_new2

#----------- New column with Date Average ----

metadata <- metadata %>%
  rowwise %>%
  mutate(Date_Average = mean.Date(c(Date_begin, Date_end), na.rm = TRUE))

#-----------Add day and month information (date = Date Average) ----

metadata$Date_Average <- as.Date(metadata$Date_Average) #convert to date format

metadata$monthFull <- format(metadata$Date_Average, "%b-%y")  # Month + Year
metadata$month <- format(metadata$Date_Average, "%b")
metadata$day <- format(metadata$Date_Average, "%d") 

#----------- Assign Seasons based on Date Average ---- 

assignSeason <- function(Date_Average) {
  if (is.na(Date_Average)) {
    return(NA)
  }
  month <- format(as.Date(Date_Average), "%b")
  day <- as.numeric(format(as.Date(Date_Average), "%d"))
  
    if ((month == "Apr" & day > 15) || (month == "May") || (month == "Jun" & day <= 15)) {
      return("Spring")
    } else if ((month == "Jun" & day > 15) || (month == "Jul")) {
      return("Summer")
    } else if ((month == "Aug") || (month == "Sep") ||(month == "Oct")) {
      return("Autumn")
    } else {
      return("Winter")
    }
  
  
  return("Unknown")
}

metadata$Season <- mapply(assignSeason, metadata$Date_Average)

#----------- Add Ice categories ----
metadata$IceCat <- case_when(
  metadata$Ice_Conc > 70 ~ ">70%",
  metadata$Ice_Conc > 50 & metadata$Ice_Conc < 70 ~ "30-70%",
  metadata$Ice_Conc > 30 & metadata$Ice_Conc <= 50 ~ "30-50%",
  metadata$Ice_Conc > 10 & metadata$Ice_Conc <= 30 ~ "10-30%",
  metadata$Ice_Conc < 10 ~ "<10%")


#----------- Add Water Fraction categories ----
metadata$WaterFrac <- case_when(
  metadata$fraction_AW > 0.8 ~ "AW",
  metadata$fraction_PW > 0.8 ~ "PW", 
  metadata$fraction_AW <= 0.8 & metadata$fraction_PW <= 0.8 ~ "MIX")

#----------- Organize Columns ----
metadata <- metadata[, c(1, 2, 3, 4, 5, 6, 12, 7, 8, 9, 10, 11, 20, 21, 34, 38,
                         13:19, 22:33,35:37, 39, 40)]

#-----------  HANDLING THE NA VALUES ---- 

metadata[] <- lapply(metadata, function(x) if(is.numeric(x)) ifelse(is.nan(x), NA, x) else x)

#---------------- EXPORT
write.table(metadata, file="MET_fstrait.txt", sep="\t", row.names=T, quote=F)
save.image(file = "R_DATA_edmet.RData")
