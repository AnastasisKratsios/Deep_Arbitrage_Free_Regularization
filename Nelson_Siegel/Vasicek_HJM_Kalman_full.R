#################################################################################
##                                                                             ##
##                 Code For Estimating Vasicek-HJM Model                       ##
##                                                                             ##
#################################################################################
# dX = kappa(theta -X)dt + sigma dW
# dX/dt = a + bX + epsilon

### OUTLINE
# A) Read and Format Data 
# B) Load/Install Packages and write Auxiliary Functions
# C) Estimate auxiliry Variables for OU process
# D) Write Kalman Filtering Algorithm
# E) Perform Kalman Filtering and Bond price Prediction for 1-day ahead bond prices





#---------#-############-#-------#
#--------------------------------#
#         A: Read-Data           #
#--------------------------------#
#---------#-############-#-------#

# Initialize Auxliliary Quanities
maturity  <- Maturities.Grid<-c(.5,seq(from=1,to=30,by=1)) #<-c(0.25,0.5,1,2,5,10)
spliceddata<-YC.DEU#USrates.52to04
spliceddata<-spliceddata[-1,-1]
bond.data<-matrix(unlist(spliceddata),ncol=length(maturity))
bond.data<-matrix(as.numeric(bond.data),ncol=length(maturity))


# Format and Initialize Related Quantities
obs<-bond.data[,2]
dates<-as.Date.ts(YC.DEU[,1])
data<-bond.data[,1] # Rates for Shortest observed bond

# Initialize Auxliliary Quanities
d         <- length(maturity)   # dimension of observations
N         <- length(data)-1 # number of observations
dt<-1/12
m<-1
lagdata<-obs[1:N]
delta_t<-1/365



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


#        #-############-#        #
#--------------------------------#
#     C: Estimate Parameters     #
#--------------------------------#
#        #-############-#        #


# Vasicek Model
bhat.Vas<-(sum(data[-1]*lagdata) - sum(data[-1])*sum(lagdata)/N)/(sum(lagdata*lagdata) - sum(lagdata)*sum(lagdata)/N)
kappahat.Vas.Vas<--log(bhat.Vas)/dt
ahat.Vas<-sum(data)/N-bhat.Vas*sum(lagdata)/N
thetahat.Vas<-ahat.Vas/(1-bhat.Vas)

s2hat<-sum((data[-1]-lagdata*bhat.Vas-ahat.Vas)^2)/N
sigmahat.Vas<-sqrt(2*kappahat.Vas.Vas*s2hat/(1-bhat.Vas^2))


#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#               D) Write Auxiliary Functions
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#

# SIMULATION 
# simulate short rate process and zero-coupon rates under Vasicek model
# ----------------------------------------------------------------------------------------
# Model setup
# PARAMETER ESTIMATION
# ----------------------------------------------------------
# Random initialization of parameters
lambda_init <- runif(m, min=-1.0, max=1.0)
measurement_err_init <- runif(d, min=0.0, max=sigmahat.Vas)
# ----------------------------------------------------------

### HJM-ification (extension)

# r is current estimated short rate
# x is Musiela parameterization 
FRC.Vasicek<-function(r,x){
  as.vector(r*exp(-x*bhat.Vas) + (ahat.Vas/bhat.Vas)*(1-exp(-x*bhat.Vas)) - ((sigmahat.Vas^2)/(2*bhat.Vas^2))*(1-exp(-x*bhat.Vas))^2)
}


vasicek_KF <- function(kappa=kappahat.Vas.Vas, theta=thetahat.Vas, sigma=sigmahat.Vas, lambda=lambda_init, measurement_errs=measurement_err_init, observations)
{ 
  # initial state variable (a0: m x 1)
  r_init <- as.vector(thetahat.Vas)
  
  # variance of state variable (P0: m x m)
  P_init <- (sigma^2/(2*kappa))*diag(1,m,m) # unconditional variance of state variable
  
  # intercept of state transition equation (dt: m x 1)
  C <- matrix((theta)*(1-exp(-kappa*delta_t)))
  
  # factor of transition equation (Tt: m x m x 1)
  F_ <- array(exp(-kappa*delta_t)*diag(m),dim=c(m,m,1))
  
  # factor of measurement equation (Zt: d x m x 1)
  B <- array(1/matrix(rep(kappa,d),d,byrow=TRUE)*(1-exp(-matrix(rep(kappa,d),d,byrow=TRUE) * matrix(rep(maturity,m),d))),dim=c(d,m,1))
  
  # intercept of measurement equation (ct: d x 1)
  gamma <- kappa^2*(theta-sigma*lambda/kappa) - 0.5*sigma^2
  A <- t(gamma/(kappa^2)) %*% t(B[,,1]-as.vector(matrix(rep(maturity,m),d))) + t(-sigma^2/(4*kappa)) %*% t(B[,,1]^2)
  A <- matrix(-A/maturity,nrow=length(maturity))
  
  B <- array(B[,,1]/matrix(rep(maturity,m),d),dim=c(d,m,1))
  
  # variance of innovations of transition (HHt: m x m x 1)
  Q <- array(sigma^2/(2*kappa)*(1-exp(-2*kappa*delta_t))*diag(m),dim=c(m,m,1))
  
  # variance of measurement error (GGt: d x d x 1)
  R <- array(diag(d)*measurement_errs^2,dim=c(d,d,1))
  
  filtered_process <- fkf(a0=r_init, P0=P_init, dt=C, ct=A, Tt=F_, Zt=B, HHt=Q, GGt=R, yt=matrix(observations))
  return(filtered_process)
}
#END B









#---------#-############-#-------#
#--------------------------------#
#         E: Prediction          #
#--------------------------------#
#---------#-############-#-------#



# SEE PREDICTION FILE!!