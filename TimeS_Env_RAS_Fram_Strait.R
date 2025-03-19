#Environmental Data Moorings RAS
#Project: Fram Strait

#HG-IV-S1: From 12-July-2016 to 27-July-2017
#HG-IV-S2: From 11-August-2017 to 16-July-2018
#HG-IV-S3: From 17-July-2018 to 24-August-2019
#HG-IV-S4: From 26-August-2019 to 01-June-2021

library (dplyr)
library (ggplot2)
library (lubridate)
library (patchwork)
library (gridExtra)


#----------- Loading the data ----
HS1 <- read.table("HG_IV_S1_RAS_Env.txt", header = TRUE, sep = "\t")
HS2 <- read.table("HG_IV_S2_RAS_Env.txt", header = TRUE, sep = "\t")
HS3 <- read.table("HG_IV_S3_RAS_Env.txt", header = TRUE, sep = "\t")
HS4 <- read.table("HG_IV_S4_RAS_Env.txt", header = TRUE, sep = "\t")

HG_Sea_Ice <- read.table("HG_IV_SeaIce.txt", header = TRUE, sep = "\t")

#Joining all the tables
All_moorings <- rbind(HS1, HS2, HS3, HS4)

# Removing time from Date Column
All_moorings$DateTime <- as.POSIXct(All_moorings$DateTime, format="%m/%d/%Y %H:%M", tz = "UTC") #defining the Date Column as Date-time object
HG_Sea_Ice$Date <- as.Date(HG_Sea_Ice$Date, format="%m/%d/%Y")

All_moorings$Date <- as.Date(All_moorings$DateTime) #Creating the Date column and extracting only the Date

#Remove Column DateTime
All_moorings <- All_moorings %>%
  select(-DateTime)

#Changing position of Date Column
All_moorings <- All_moorings[, c(1, 14, 2, 3, 4, 5, 6, 7 ,8 ,9 ,10, 11, 12, 13)]

#----------- Averages - From hour to daily resolution ----

# Group by Day and calculate the average for each parameter 
All_moorings_average <- All_moorings %>%
  group_by(Date) %>%
  summarize(
    Mooring = Mooring[1],
    across(where(is.numeric), mean, na.rm = TRUE))

All_moorings_average <- All_moorings_average[, c(2, 1, 3, 4, 5, 6, 7 ,8 ,9 ,10, 11, 12, 13, 14)]

#Merge with Sea Ice data
All_dailyres_moorings <- full_join(All_moorings_average, HG_Sea_Ice, by = "Date")

All_dailyres_moorings <- All_dailyres_moorings %>%
  select(-Mooring.y)

#Arrange by Date
All_dailyres_moorings <- All_dailyres_moorings %>% arrange(Date)

#Remove the rows that only have sea-ice data in the beginning (no mooring values)
#From  to 2016-05-01 to 2016-07-11

All_dailyres_moorings <- All_dailyres_moorings %>%
  filter(!(Date >= as.Date("2016-05-01") & Date <= as.Date("2016-07-11")))


#Adding the column Season
All_dailyres_moorings_Season <- All_dailyres_moorings %>%
  mutate(Season = case_when(
    #Spring: April 16 to June 15 = mid-April to mid-June
    (month(Date) == 4 & day(Date) > 15) | month(Date) %in% c(5) | (month(Date) == 6 & day(Date) <= 15) ~ "Spring", 
    #Summer: June 16 to July 31 = mid-June to late-July
    (month(Date) == 6 & day(Date) > 15) | month(Date) %in% c(7) ~ "Summer",   
    #Autumn: August to October
    (month(Date) == 8 & day(Date) >= 1) | month(Date) %in% c(9, 10) ~ "Autumn",
    #Winter: November to April 15 = November to mid-April
    (month(Date) %in% c(11, 12, 1, 2, 3)) | (month(Date) == 4 & day(Date) <= 15) ~ "Winter",
    # In case of missing or unexpected values                                           
    TRUE ~ NA_character_ ))

All_dailyres_moorings_Season <- All_dailyres_moorings_Season[, c(1, 2, 17, 3, 4, 5, 6, 7 ,8 ,9 ,10, 11, 12, 13, 14, 15, 16)]

#---------- Exporting the File ----

write.csv(All_dailyres_moorings_Season, "Daily_Res_Env_RAS_FramStrait.csv", row.names = FALSE) #to save as .csv

#---------- Plots in Daily Resolution ----

#Check if Date is as Date type
str(All_dailyres_moorings_Season$Date)
All_dailyres_moorings_Season$Year <- format(All_dailyres_moorings_Season$Date, "%Y") #Adding the Year column


# TS Daily - With Event + Seasons
#Event: Cold Anomaly

daily_env_upper_notrend <- function(env_variable, plot_name, variable_name) {
  plot <- ggplot(All_dailyres_moorings_Season, aes(x = Date, y = .data[[env_variable]])) +
    geom_rect(aes(xmin = as.Date("2017-04-15"), xmax = as.Date("2017-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-04-15"), xmax = as.Date("2018-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-04-15"), xmax = as.Date("2019-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-04-15"), xmax = as.Date("2020-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2021-04-15"), xmax = as.Date("2021-06-01"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-07-01"), xmax = as.Date("2016-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-06-15"), xmax = as.Date("2017-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-06-15"), xmax = as.Date("2018-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-06-15"), xmax = as.Date("2019-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-06-15"), xmax = as.Date("2020-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-08-01"), xmax = as.Date("2016-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-08-01"), xmax = as.Date("2017-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-08-01"), xmax = as.Date("2018-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-08-01"), xmax = as.Date("2019-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-08-01"), xmax = as.Date("2020-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-11-01"), xmax = as.Date("2017-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-11-01"), xmax = as.Date("2018-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-11-01"), xmax = as.Date("2019-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-11-01"), xmax = as.Date("2020-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-11-01"), xmax = as.Date("2021-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_line(linewidth = 0.5, color = "black") +
    geom_vline(xintercept = as.Date("2020-01-01"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = as.Date("2020-09-15"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    labs(title = plot_name, x = "", y = variable_name, tag = "Cold Anomaly") +
    scale_x_date(limits = as.Date(c("2016-07-01", "2021-06-01")),date_breaks = "1 month", date_labels = "%b %Y",expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          plot.tag.position = c(.759,.98),
          plot.tag = element_text(hjust = 0, color = "blue", size = 13))
  
  return(plot)
}
daily_env_upper_trend <- function(env_variable, plot_name, variable_name) {
  plot <- ggplot(All_dailyres_moorings_Season, aes(x = Date, y = .data[[env_variable]])) +
    geom_rect(aes(xmin = as.Date("2017-04-15"), xmax = as.Date("2017-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-04-15"), xmax = as.Date("2018-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-04-15"), xmax = as.Date("2019-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-04-15"), xmax = as.Date("2020-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2021-04-15"), xmax = as.Date("2021-06-01"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-07-01"), xmax = as.Date("2016-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-06-15"), xmax = as.Date("2017-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-06-15"), xmax = as.Date("2018-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-06-15"), xmax = as.Date("2019-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-06-15"), xmax = as.Date("2020-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-08-01"), xmax = as.Date("2016-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-08-01"), xmax = as.Date("2017-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-08-01"), xmax = as.Date("2018-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-08-01"), xmax = as.Date("2019-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-08-01"), xmax = as.Date("2020-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-11-01"), xmax = as.Date("2017-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-11-01"), xmax = as.Date("2018-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-11-01"), xmax = as.Date("2019-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-11-01"), xmax = as.Date("2020-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-11-01"), xmax = as.Date("2021-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_line(linewidth = 0.5, color = "black") +
    geom_vline(xintercept = as.Date("2020-01-01"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = as.Date("2020-09-15"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_smooth(method="gam", formula=y ~ s(x, k=20), color ="#DC143C", fill ="indianred", se = TRUE, linewidth = 0.2)+
    labs(title = plot_name, x = "", y = variable_name, tag = "Cold Anomaly") +
    scale_x_date(limits = as.Date(c("2016-07-01", "2021-06-01")),date_breaks = "1 month", date_labels = "%b %Y",expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          plot.tag.position = c(.759,.98),
          plot.tag = element_text(hjust = 0, color = "blue", size = 13))
  
  return(plot)
}

#TAG specifications
#Big: .759, .98 | size = 13
#Small: .742, .98 | size = 11

daily_env_plot_trend <- function(env_variable, plot_name, variable_name) {
  plot <- ggplot(All_dailyres_moorings_Season, aes(x = Date, y = .data[[env_variable]])) +
    geom_rect(aes(xmin = as.Date("2017-04-15"), xmax = as.Date("2017-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-04-15"), xmax = as.Date("2018-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-04-15"), xmax = as.Date("2019-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-04-15"), xmax = as.Date("2020-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2021-04-15"), xmax = as.Date("2021-06-01"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-07-01"), xmax = as.Date("2016-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-06-15"), xmax = as.Date("2017-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-06-15"), xmax = as.Date("2018-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-06-15"), xmax = as.Date("2019-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-06-15"), xmax = as.Date("2020-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-08-01"), xmax = as.Date("2016-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-08-01"), xmax = as.Date("2017-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-08-01"), xmax = as.Date("2018-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-08-01"), xmax = as.Date("2019-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-08-01"), xmax = as.Date("2020-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-11-01"), xmax = as.Date("2017-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-11-01"), xmax = as.Date("2018-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-11-01"), xmax = as.Date("2019-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-11-01"), xmax = as.Date("2020-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-11-01"), xmax = as.Date("2021-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_line(linewidth = 0.5, color = "black") +
    geom_smooth(method="gam", formula=y ~ s(x, k=20), color ="#DC143C", fill ="indianred", se = TRUE, linewidth = 0.2)+
    geom_vline(xintercept = as.Date("2020-01-01"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = as.Date("2020-09-15"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    labs(title = plot_name, x = "", y = variable_name) +
    scale_x_date(limits = as.Date(c("2016-07-01", "2021-06-01")),date_breaks = "1 month", date_labels = "%b %Y",expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_blank())
  
  return(plot)
}

daily_env_below_notrend <- function(env_variable, plot_name, variable_name) {
  plot <- ggplot(All_dailyres_moorings_Season, aes(x = Date, y = .data[[env_variable]])) +
    geom_rect(aes(xmin = as.Date("2017-04-15"), xmax = as.Date("2017-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-04-15"), xmax = as.Date("2018-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-04-15"), xmax = as.Date("2019-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-04-15"), xmax = as.Date("2020-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2021-04-15"), xmax = as.Date("2021-06-01"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-07-01"), xmax = as.Date("2016-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-06-15"), xmax = as.Date("2017-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-06-15"), xmax = as.Date("2018-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-06-15"), xmax = as.Date("2019-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-06-15"), xmax = as.Date("2020-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-08-01"), xmax = as.Date("2016-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-08-01"), xmax = as.Date("2017-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-08-01"), xmax = as.Date("2018-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-08-01"), xmax = as.Date("2019-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-08-01"), xmax = as.Date("2020-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-11-01"), xmax = as.Date("2017-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-11-01"), xmax = as.Date("2018-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-11-01"), xmax = as.Date("2019-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-11-01"), xmax = as.Date("2020-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-11-01"), xmax = as.Date("2021-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    labs(title = plot_name, x = "", y = variable_name) +
    geom_line(linewidth = 0.5, color = "black") +
    geom_vline(xintercept = as.Date("2020-01-01"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = as.Date("2020-09-15"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    scale_x_date(limits = as.Date(c("2016-07-01", "2021-06-01")),date_breaks = "1 month", date_labels = "%b %Y",expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  return(plot)
}
daily_env_below_trend <- function(env_variable, plot_name, variable_name) {
  plot <- ggplot(All_dailyres_moorings_Season, aes(x = Date, y = .data[[env_variable]])) +
    geom_rect(aes(xmin = as.Date("2017-04-15"), xmax = as.Date("2017-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-04-15"), xmax = as.Date("2018-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-04-15"), xmax = as.Date("2019-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-04-15"), xmax = as.Date("2020-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2021-04-15"), xmax = as.Date("2021-06-01"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-07-01"), xmax = as.Date("2016-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-06-15"), xmax = as.Date("2017-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-06-15"), xmax = as.Date("2018-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-06-15"), xmax = as.Date("2019-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-06-15"), xmax = as.Date("2020-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-08-01"), xmax = as.Date("2016-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-08-01"), xmax = as.Date("2017-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-08-01"), xmax = as.Date("2018-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-08-01"), xmax = as.Date("2019-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-08-01"), xmax = as.Date("2020-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-11-01"), xmax = as.Date("2017-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-11-01"), xmax = as.Date("2018-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-11-01"), xmax = as.Date("2019-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-11-01"), xmax = as.Date("2020-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-11-01"), xmax = as.Date("2021-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    labs(title = plot_name, x = "", y = variable_name) +
    geom_smooth(method="gam", formula=y ~ s(x, k=20), color ="#DC143C", fill ="indianred", se = TRUE, linewidth = 0.2)+
    geom_line(linewidth = 0.5, color = "black") +
    geom_vline(xintercept = as.Date("2020-01-01"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = as.Date("2020-09-15"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    scale_x_date(limits = as.Date(c("2016-07-01", "2021-06-01")),date_breaks = "1 month", date_labels = "%b %Y",expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  return(plot)
}

#Temp, Sal, Sea-ice Conc and MLD
daily_icecon <- daily_env_upper_notrend("Ice_Conc", "Sea-ice Concentration","Sea-ice [%]")
daily_temp <- daily_env_plot_trend("Temperature", "Temperature","Temperature [°C]")
daily_sal <- daily_env_plot_trend("Salinity", "Salinity","Salinity")
daily_mld <- daily_env_below_notrend("MLD", "Mixed Layer Depth","MLD [m]")

plot4 <- (daily_icecon/ daily_temp / daily_sal / daily_mld) + plot_layout(heights = c(2, 2, 2, 2.25))
plot4

#To export as PNG image: Klein(w850 x H700), Gross(W1600xH900)
#To export as PDF (inches)

ggsave("TimeS_Temp_Sal_IceCon_MLD_Gross.png", plot = plot4, width = 16, height = 9, units = "in")
ggsave("TimeS_Temp_Sal_IceCon_MLD_Klein.png", plot = plot4, width = 8.5, height = 8, units = "in")


#Plot Chl, Sea-ice distance, Oxy Con, Oxy Sat
daily_chl <- daily_env_upper_notrend("Chl", "Chlorophyll Concentration (Uncalibrated)","Chl Concentration [~mg/m³]")
daily_icedist <- daily_env_plot_trend("Ice_Dist", "Distance of sea-ice edge","distance [km]")
daily_oxycon <- daily_env_below_trend("Oxygen_Conc", "Dissolved Oxygen Concentration","O2 Concentration [µmol/l]")

plot5 <- (daily_chl / daily_icedist / daily_oxycon) + plot_layout(heights = c(2, 2, 2.25))
plot5

#To export as PNG image: Klein(w850 x H700), Gross(W1600xH900)
#To export as PDF (inches)

ggsave("TimeS_Oxy_Chl_IceDist_Gross.png", plot = plot5, width = 16, height = 9, units = "in")
ggsave("TimeS_Oxy_Chl_IceDist_Klein.png", plot = plot5, width = 8.5, height = 7, units = "in")


# ---- Time-Series Daily - With seasons - Water Fraction

#Transform in %
All_dailyres_moorings_Season <- All_dailyres_moorings_Season %>%
  mutate(fraction_AW = fraction_AW * 100,
    fraction_PW = fraction_PW * 100
  )


daily_env_below_trendblue <- function(env_variable, plot_name, variable_name){
  plot <- ggplot(All_dailyres_moorings_Season, aes(x = Date, y = .data[[env_variable]])) +
    geom_rect(aes(xmin = as.Date("2017-04-15"), xmax = as.Date("2017-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-04-15"), xmax = as.Date("2018-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-04-15"), xmax = as.Date("2019-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-04-15"), xmax = as.Date("2020-06-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2021-04-15"), xmax = as.Date("2021-06-01"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.6, 0.98, 0.6), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-07-01"), xmax = as.Date("2016-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-06-15"), xmax = as.Date("2017-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-06-15"), xmax = as.Date("2018-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-06-15"), xmax = as.Date("2019-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-06-15"), xmax = as.Date("2020-07-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.85, 0.73), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-08-01"), xmax = as.Date("2016-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-08-01"), xmax = as.Date("2017-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-08-01"), xmax = as.Date("2018-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-08-01"), xmax = as.Date("2019-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-08-01"), xmax = as.Date("2020-10-31"), ymin = -Inf, ymax = Inf), 
              fill = rgb(1, 0.71, 0.76), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2016-11-01"), xmax = as.Date("2017-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2017-11-01"), xmax = as.Date("2018-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2018-11-01"), xmax = as.Date("2019-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2019-11-01"), xmax = as.Date("2020-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    geom_rect(aes(xmin = as.Date("2020-11-01"), xmax = as.Date("2021-04-15"), ymin = -Inf, ymax = Inf), 
              fill = rgb(0.68, 0.85, 1), alpha = 0.01)+
    labs(title = plot_name, x = "", y = variable_name) +
    geom_vline(xintercept = as.Date("2020-01-01"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = as.Date("2020-09-15"), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_line(linewidth = 0.5, color = "black") +
    geom_smooth(method="gam", formula=y ~ s(x, k=20), color ="blue", fill ="steelblue", se = TRUE, linewidth = 0.2)+
    scale_x_date(limits = as.Date(c("2016-07-01", "2021-06-01")),date_breaks = "1 month", date_labels = "%b %Y",expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  return(plot)
}

daily_aw <- daily_env_upper_trend("fraction_AW","Proportion of Atlantic Water ", "Fraction of Water [%]")
daily_pw <- daily_env_below_trendblue("fraction_PW","Proportion of Polar Water", "Fraction of Water [%]")

plot6 <- (daily_aw / daily_pw) + plot_layout(heights = c(2, 2.25))
plot6

#To export as PNG image: Klein(w850 x H700), Gross(W1600xH900)
#To export as PDF (inches)

ggsave("TimeS_WaterFrac.png", plot = plot6, width = 16, height = 9, units = "in")
ggsave("TimeS_WaterFrac_Klein.png", plot = plot6, width = 8.5, height = 6, units = "in")

#---------- Extra info ----
#Min, Max and Average + SD of each parameter
#Date = y/m/d
#Obs: No Chl values, since it is not calibrated

library(writexl)

general_info <- function (df, variablecol) {
  min <- df %>% 
    filter(df[[variablecol]] == min(df[[variablecol]], na.rm = TRUE)) %>%
    select(all_of(c("Date", variablecol, "Season")))
  
  max <- df %>%
    filter(df[[variablecol]] == max(df[[variablecol]], na.rm = TRUE)) %>%
    select(all_of(c("Date", variablecol, "Season")))
  
  #average
  avg <- mean(df[[variablecol]], na.rm = TRUE)
  sd <- sd(df[[variablecol]], na.rm = TRUE)
  
  return(list(
    Minimum = min, Maximum = max,
    Mean = avg,
    Standard_Deviation = sd))
}

info_temp <- general_info(All_dailyres_moorings_Season, "Temperature")
info_temp

info_sal <- general_info(All_dailyres_moorings_Season, "Salinity")
info_sal

info_oxyc <- general_info(All_dailyres_moorings_Season, "Oxygen_Conc")
info_oxyc

info_oxys <- general_info(All_dailyres_moorings_Season, "Oxygen_Sat")
info_oxys

info_mld <- general_info(All_dailyres_moorings_Season, "MLD")
info_mld

info_iced <- general_info(All_dailyres_moorings_Season, "Ice_Dist")
info_iced

info_aw <- general_info(All_dailyres_moorings_Season, "fraction_AW")
info_aw

info_pw <- general_info(All_dailyres_moorings_Season, "fraction_PW")
info_pw

info_icec <- general_info(All_dailyres_moorings_Season, "Ice_Conc")
info_icec 

#For Seasons
general_info_sea <- function(df, variablecol) {
  # Filter data for Winter season
  cold_data <- df %>% 
    filter(Season == "Winter" | Season == "Spring")
  
  warm_data <- df %>% 
    filter(Season == "Summer" | Season == "Autumn")

  # Average and Standard Deviation 
  
  avg1 <- mean(winter_data[[variablecol]], na.rm = TRUE)
  sd_val1 <- sd(winter_data[[variablecol]], na.rm = TRUE)
  
  avg2 <- mean(summer_data[[variablecol]], na.rm = TRUE)
  sd_val2 <- sd(summer_data[[variablecol]], na.rm = TRUE)
  
  return(list(
    Mean.Win = avg1, Standard_Deviation.W = sd_val1, 
    Mean.Sum = avg2, Standard_Deviation.S = sd_val2 ))
  }

info_MLD2 <- general_info_temp(All_dailyres_moorings_Season, "MLD")
info_MLD2

#I saved the results manually to an Excel file: Env_Info

save.image(file = "R_DATA_env.RData")


