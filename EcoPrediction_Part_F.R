library(ncdf4)
library(ncdf4.helpers)
library(readr)
library(dplyr)
library(tidyr)
library(plyr)
library(purrr)

setwd("/Volumes/GoogleDrive/My Drive/VaughnData/PhenoPrediction")    

# Read in a list of all files with a .nc extension/file-type within the directory
files <- list.files(path="noaa/NOAAGEFS_1hr/BART/2021-01-27/00/", pattern="*.nc", full.names=TRUE)

# Create a list to populate with extracted values
dat.list <- list()

# Create a for loop which will iterate through the .nc files and do some math/processing
for(i in 1:length(files)){
  
  # Read in an individual file from the list we created and extract air_temperature variable
  # Convert it to a data.frame object
  mydata <- nc_open(files[i])
  airtemp <- ncvar_get(mydata, varid="air_temperature")
  airtemp2=data.frame(airtemp)
  
  # Convert the temperature to Celsius from Kelvin.
  # Calculate the mean daily temperature by aggregating the extracted air_temperature variable
  # by every 24 hours.
  airtemp2 <- airtemp2-273.15
  dailyTemp <- aggregate(airtemp2, by = list(gl(ceiling(nrow(airtemp2)/24), 24)[1:nrow(airtemp2)]), FUN = mean)
  
  # Rename the variable columns for clarity
  names(dailyTemp)[1] <- "DOY"
  names(dailyTemp)[2] <- "GDD"
  
  # Use tidyverse mutate functions to create a new variable which adjusts the GDD accumulation by +5
  # considering how trees accumulate GDD. This method saves time over having two nested for-loops
  dailyTemp <- dailyTemp %>% dplyr::mutate(GDD2=ifelse(GDD >= -5, GDD+5, 0)) 
  
  # Add cumulative corrected GDD value
  dailyTemp <- dailyTemp %>% dplyr::mutate(cumGDD2=cumsum(GDD2))
  
  dat.list[[i]] <- dailyTemp
}

# This returns the row index of the first value greater than or equal to a target threshold
# Naresh, you can replace these with the target GDD value for the 50% green-up.
min(which(dat.list[[1]]$cumGDD2 >= 10))

# This returns the actual DOY of the above index search
as.numeric(dat.list[[i]][min(which(dat.list[[1]]$cumGDD2 >= 10)),]$DOY)

# The most straightforward way to automate for many .nc models from multiple locations would be
# to ensure the model is processing files in the order you want and to create a vector object
# in which the GDD threshold number is set for EACH .nc file (so there will be many repeating
# thresholds).

# Set an arbitrary threshold at 10 GDD as an example 10 times (for the 10 .nc files)
# In real life this would vary per site.
# GDD accumulated from 1 Jan-30 Jan 2021, values in () indicate the threshold for 50% Peak
# CLBJ 2130 (1640) 
# STEI 739 (1152)
# SCBI 1239 (1597)
# GRSM 1538 (1551)
# DELA 2416 (1635)
# BART 641 (1004)
# HARV 931 (1277)
# UKFS 1316 (1387)

thres <- c(rep(10, 10)) 

# Create a vector of DOY outcomes based on GDD threshold to save DOYs to
DOYs <- c()

for(i in 1:length(dat.list)){
  DOYs[i] <- as.numeric(dat.list[[i]][min(which(dat.list[[i]]$cumGDD2 >= 10)),]$DOY)
}

# Calculate statistics for this vector of DOYs (you would want to make sure you subset future
# vectors so that you don't calculate statistics across all models/sites)
sd(DOYs)
mean(DOYs)



