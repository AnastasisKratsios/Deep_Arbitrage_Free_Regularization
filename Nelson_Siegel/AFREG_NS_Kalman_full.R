# Note:
# Must run first: 
# 1) NS_Kalman 
# 2) Deep_Arbitrage_Free_Regularization

##########------######################---#
# AF-Reg: Kalman Filter for dNS Model
###########--###--######-##-##########---#
# OUTLINE
# A) Read Data
# B) Load Packages
# C) Estimate dNS Parameters
# D) Write dNS Kalman Filtering Function
# E) Predict (external)


#---------#-############-#-------#
#--------------------------------#
#         A: Read-Data           #
#--------------------------------#
#---------#-############-#-------#
# Run NS_Kalman_Full file!




#        #-############-#        #
#--------------------------------#
#        B: Load Packages        #
#--------------------------------#
#        #-############-#        #

## Load/Install Packages
if(!exists("NO.CHECK.Vas.Kalm")){
  # Check for missing packages
  list.of.packages <- c("beepr","foreach","FKF","ggplot2","dplyr","magrittr","gridExtra","mvtnorm","ggplot2","YieldCurve","Meucci","xtable")
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
  # Ensure that this bit is not re-run until next time code it loaded
  NO.CHECK.Vas.Kalm<-TRUE
  # Signal that packags have been loaded and installed
  beep()
}
#
#-############ -----------------------------------################-#
# Load required packages if available...if not install - END
#-############ -----------------------------------################-#


#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#                  C) Initialize Parameters
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#


##### INITIALIZATION Parameters
log.bond.dat<-matrix(as.vector(unlist(bond.data)),ncol=ncol(bond.data))
log.bond.dat<-t(diff(t(cbind(rep(0,nrow(log.bond.dat)),log.bond.dat))))
log.bond.dat<-t(ginv(AF.Reg.Factors)%*%t(log.bond.dat))
ginv(AF.Reg.Factors)
Estimated.OU.Parameters.Betas.AFREG<-FitOrnsteinUhlenbeck(log.bond.dat,1)




#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#                  D) Write Filter
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
# Random initialization of parameters
lambda_init <- runif(m, min=-1.0, max=1.0)
measurement_err_init.AFReg <- diag(runif(d, 
                                   min=0.0, 
                                   max=sqrt(sum((eigen(Estimated.OU.Parameters.Betas$Sigma)$values)^2))
))


# Kalman NS
AFReg_NS_KF <- function(observations)
{ 
  # initial state variable (a0: m x 1)
  r_initAF.Reg.NS <- as.vector(Estimated.OU.Parameters.Betas.AFREG$Mu)
  
  # variance of state variable (P0: m x m)
  P_initAF.Reg.NS <- t(Estimated.OU.Parameters.Betas.AFREG$Sigma)%*%Estimated.OU.Parameters.Betas.AFREG$Sigma
  P_initAF.Reg.NS <- .5*P_initAF.Reg.NS
  
  # intercept of state transition equation (dt: m x 1)
  CAF.Reg.NS <- (Estimated.OU.Parameters.Betas.AFREG$Theta)%*%(Estimated.OU.Parameters.Betas.AFREG$Mu)
  
  # factor of transition equation (Tt: m x m x 1)
  TtAF.Reg.NS <- -Estimated.OU.Parameters.Betas.AFREG$Theta
  
  # factor of measurement equation (Zt: d x m x 1)
  BAF.Reg.NS <- NS.samples[,-1]
  
  # intercept of measurement equation (ct: d x 1)
  AAF.Reg.NS <- rep(0,d)
  
  
  # variance of innovations of transition (HHt: m x m x 1)
  QAF.Reg.NS <- t(Estimated.OU.Parameters.Betas.AFREG$Sigma)%*%Estimated.OU.Parameters.Betas.AFREG$Sigma
  
  # variance of measurement error (GGt: d x d x 1)
  RAF.Reg.NS <-  array(measurement_err_init.AFReg,dim=c(d,d,1))
  
  filtered_process <- fkf(a0=r_initAF.Reg.NS, P0=P_initAF.Reg.NS, dt=CAF.Reg.NS, ct=AAF.Reg.NS, Tt=TtAF.Reg.NS, Zt=BAF.Reg.NS, HHt=QAF.Reg.NS, GGt=RAF.Reg.NS, yt=matrix(observations))
  return(filtered_process)
}





#---------#-############-#-------#
#--------------------------------#
#         E: Prediction          #
#--------------------------------#
#---------#-############-#-------#



# SEE PREDICTION FILE!!