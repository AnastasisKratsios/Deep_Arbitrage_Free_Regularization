################################################-#
#         Deep Arbitrage-Free Learning           #
################################################-#
# Implementation for the paper of:
# Anastasis Kratsios - ETH, Zurich, Switzerland
# Cody Hyndman       - Concordia, Montreal, Canada
# 2019


# Notes: Both Experiments Should be run separately.

# Specify Put Data and Maturities
maturity  <- Maturities.Grid<-c(.5,seq(from=1,to=30,by=1)) #<-c(0.25,0.5,1,2,5,10)
spliceddata<-YC.DEU




# Run Files
#-----------#
#-----------------------------------------------#
### Predictions for Nelson-Siegel Experiments ###
#-----------------------------------------------#
# Initializations
#setwd("/scratch/users/kratsioa/Dropbox/Numerics/Deep_Arbitrage_Free_Regularization/Full_Version/dPCA/")
wd<-getwd() 
source(paste(wd,"Initializations.R",sep="/")) # Initialize

# Nelson-Siegel Models
source(paste(wd,"PCA_Kalman_full.R",sep="/")) # Nelson-Siegel Model
source(paste(wd,"SGD.R",sep="/")) # Loads SGD Code
source(paste(wd,"Deep_Arbitrage_Free_Regularization_full_PCA.R",sep="/")) # AF-Regularization of AFNS
source(paste(wd,"AFREG_dPCA_Kalman_full.R",sep="/")) # Compiles Predictive Algorithm


# Make & Report Predictions
source(paste(wd,"Predictions_dPCA.R",sep="/"))
source(paste(wd,"Comparions_as_a_Function_of_Lambda_dPCA.R",sep="/"))

