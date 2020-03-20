################################---#
# Kalman Filter for dNS Model
################################---#
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

# Initialize Auxliliary Quanities
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
dt<-1/365
m<-1
lagdata<-obs[1:N]


#        #-############-#        #
#--------------------------------#
#        B: Load Packages        #
#--------------------------------#
#        #-############-#        #

## Load/Install Packages
if(!exists("NO.CHECK.Vas.Kalm")){
  # Check for missing packages
  list.of.packages <- c("beepr","foreach","FKF","ggplot2","dplyr","magrittr","gridExtra","mvtnorm","ggplot2","YieldCurve","Meucci","xtable","Matrix","MASS")
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

# Auxiliary Functions
#----------------------------#
# Nelson-Siegel Curves
NS.curves<-function(T,lambda=1){
  c(1,exp(-T*lambda),T*exp(-T*lambda))
}


# Estimate Nelson-Siegel Shape Parameter
if(!exists("shape.estimate")){#Run only if needed (since slow)
# Fit Nelson-Siegel on each window
NSParameters <- Nelson.Siegel(rate=bond.data,maturity=Maturities.Grid)
# Extracts Shape Estimate on Training Set
shape.estimate<-colMeans(NSParameters);shape.estimate<-shape.estimate[4]
}

# Initialize Sample from NS Manifold
NS.samples<-matrix(data=NA,ncol=4,nrow=length(Maturities.Grid))
colnames(NS.samples)<-c("T.Mat","b1","b2","b3")
rownames(NS.samples)<-Maturities.Grid

# Populate Matrix
# Initialize Progress Bar
pb <- txtProgressBar(0, nrow(NS.samples), style = 3)
for(i in 1:nrow(NS.samples)){
  # Update Progress Bar
  setTxtProgressBar(pb, i)
  # Populate Next Entry in matrix of FRC Manifold Samples
  NS.samples[i,]<-c(Maturities.Grid[i],NS.curves(T=Maturities.Grid[i],lambda = shape.estimate))
}
# Close Progress Bar
close(pb)

# Visualize Head
head(NS.samples)
# Visualize Tail
tail(NS.samples)
# Test Dimensions
dim(NS.samples)


##### INITIALIZATION Parameters
log.bond.dat<-matrix(as.vector(unlist(bond.data)),ncol=ncol(bond.data))
log.bond.dat<-t(diff(t(cbind(rep(0,nrow(log.bond.dat)),log.bond.dat))))
log.bond.dat<-t(ginv(NS.samples[,-1])%*%t(log.bond.dat))

Estimated.OU.Parameters.Betas<-FitOrnsteinUhlenbeck(log.bond.dat,1)




#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#                  D) Write Filter
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
# Random initialization of parameters
lambda_init <- runif(m, min=-1.0, max=1.0)
measurement_err_init.NS <- diag(runif(d, 
                              min=0.0, 
                              max=sqrt(sum((eigen(Estimated.OU.Parameters.Betas$Sigma)$values)^2))
                              ))


# Kalman NS
NS_KF <- function(observations)
{ 
  # initial state variable (a0: m x 1)
  r_init.NS <- as.vector(Estimated.OU.Parameters.Betas$Mu)
  
  # variance of state variable (P0: m x m)
  P_init.NS <- t(Estimated.OU.Parameters.Betas$Sigma)%*%Estimated.OU.Parameters.Betas$Sigma
  P_init.NS <- .5*P_init.NS
  
  # intercept of state transition equation (dt: m x 1)
  C.NS <- (Estimated.OU.Parameters.Betas$Theta)%*%(Estimated.OU.Parameters.Betas$Mu)
  
  # factor of transition equation (Tt: m x m x 1)
  Tt.NS <- -Estimated.OU.Parameters.Betas$Theta
  
  # factor of measurement equation (Zt: d x m x 1)
  B.NS <- NS.samples[,-1]
  
  # intercept of measurement equation (ct: d x 1)
  A.NS <- rep(0,d)
  
  
  # variance of innovations of transition (HHt: m x m x 1)
  Q.NS <- t(Estimated.OU.Parameters.Betas$Sigma)%*%Estimated.OU.Parameters.Betas$Sigma
  
  # variance of measurement error (GGt: d x d x 1)
  R.NS <- array(measurement_err_init.NS,dim=c(d,d,1))
  
  filtered_process <- fkf(a0=r_init.NS, P0=P_init.NS, dt=C.NS, ct=as.matrix(A.NS), Tt=Tt.NS, 
                          Zt=B.NS, HHt=Q.NS, GGt=R.NS, yt=as.matrix(observations))
  return(filtered_process)
}




#---------#-############-#-------#
#--------------------------------#
#         E: Prediction          #
#--------------------------------#
#---------#-############-#-------#



# SEE PREDICTION FILE!!