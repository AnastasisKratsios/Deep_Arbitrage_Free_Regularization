#-############ -----------------------------------################-#
# Load required packages if available...if not install - BEGIN
#-############ -----------------------------------################-#
# Some packages should install this package manually and separately
### The Script is in the file: Installations.R

  # Check for missing packages
  list.of.packages <- c("beepr","foreach","FKF","ggplot2","dplyr","magrittr","neuralnet","mvtnorm","ggplot2","YieldCurve","Meucci","xtable","gridExtra","mvtnorm","ggplot2","YieldCurve","Meucci","xtable","Matrix","MASS","reshape2")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  # Download missing packages
  if(length(new.packages)) install.packages(new.packages)
  # Load packages requried for code...
  lapply(list.of.packages, require, character.only = TRUE)
  #
  # Check for If Auxliary Functions Are Loaded
  list.of.required.functions<-c("activation","gradientR") # List of required functions
  Missing.Q<-lapply(list.of.required.functions, exists) # Check if each function is loaded
  Missing.Q<-(min(unlist(Missing.Q))==0) # Check if any function is missing
  # If Missing.Q is TRUE then load Auxilair functions from external file (should be kept in same working directory!!!)
  if(Missing.Q){
    # Read Working Directory
    # Move to Greater Stored Functions Directory
    setwd(paste(Default.directory,"Auxiliary_Code_and_Functions",sep = "/"))
    source("SGD.R") # Loads SGD Code NOTE:(written by Dr. Dennis Prangle of NewCastle University)
    # Load Data
    # Move Working Directory
    setwd(paste(Default.directory,"Datasets",sep = "/"))
    #bond.data<-spliceddata<-USrates.52to04
    bond.data<-spliceddata<-spliceddata<-YC.DEU#read.csv(paste(getwd(),"YC-DEU.csv",sep="/"), header=FALSE)
    # Reset/Move back to Primary Working Directory
    setwd(Default.directory)
  }
  # Ensure that this bit is not re-run until next time code it loaded
  NO.CHECK_Q<-TRUE
  # Signal that packags have been loaded and installed
  beep()

#
#-############ -----------------------------------################-#
# Load required packages if available...if not install - END
#-############ -----------------------------------################-#