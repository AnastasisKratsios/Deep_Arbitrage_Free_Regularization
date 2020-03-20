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
# C) Estimate A-Reg(dPCA) Parameters
# D) Write A-Reg(dPCA) Kalman Filtering Function
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
log.bond.dat<-t(ginv(AF.Reg.Factors.PCA)%*%t(log.bond.dat))
ginv(AF.Reg.Factors)
Estimated.OU.Parameters.Betas.PCA<-FitOrnsteinUhlenbeck(log.bond.dat,1)




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
                                   max=sqrt(sum((eigen(Estimated.OU.Parameters.Betas.PCA$Sigma)$values)^2))
))


# Kalman NS
AFReg_dPCA_KF <- function(observations)
{ 
  # initial state variable (a0: m x 1)
  r_initAF.Reg.AFReg <- as.vector(Estimated.OU.Parameters.Betas.PCA$Mu)
  
  # variance of state variable (P0: m x m)
  P_initAF.Reg.AFReg <- t(Estimated.OU.Parameters.Betas.PCA$Sigma)%*%Estimated.OU.Parameters.Betas.PCA$Sigma
  P_initAF.Reg.AFReg <- .5*P_initAF.Reg.AFReg
  
  # intercept of state transition equation (dt: m x 1)
  CAF.Reg.AFReg <- (Estimated.OU.Parameters.Betas.PCA$Theta)%*%(Estimated.OU.Parameters.Betas.PCA$Mu)
  
  # factor of transition equation (Tt: m x m x 1)
  TtAF.Reg.AFReg <- -Estimated.OU.Parameters.Betas.PCA$Theta
  
  # factor of measurement equation (Zt: d x m x 1)
  BAF.Reg.AFReg <- dPCA.samples
  
  # intercept of measurement equation (ct: d x 1)
  AAF.Reg.AFReg <- rep(0,d)
  
  
  # variance of innovations of transition (HHt: m x m x 1)
  QAF.Reg.AFReg <- t(Estimated.OU.Parameters.Betas.PCA$Sigma)%*%Estimated.OU.Parameters.Betas.PCA$Sigma
  
  # variance of measurement error (GGt: d x d x 1)
  RAF.Reg.AFReg <-  array(measurement_err_init.AFReg,dim=c(d,d,1))
  
  
  filtered_process <- fkf(a0=r_initAF.Reg.AFReg, P0=P_initAF.Reg.AFReg, dt=CAF.Reg.AFReg, ct=AAF.Reg.AFReg, Tt=TtAF.Reg.AFReg, Zt=BAF.Reg.AFReg, HHt=QAF.Reg.AFReg, GGt=RAF.Reg.AFReg, yt=matrix(observations))
  return(filtered_process)
}





#---------#-############-#-------#
#--------------------------------#
#         E: Prediction          #
#--------------------------------#
#---------#-############-#-------#



# SEE PREDICTION FILE!!