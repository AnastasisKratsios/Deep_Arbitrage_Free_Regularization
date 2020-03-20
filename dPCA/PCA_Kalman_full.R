N.factors<-3
################################-------#
# Kalman Filter for AF-Reg(dynamic PCA)
################################-------#
# OUTLINE
# A) Read Data
# B) Load Packages
# C) Estimate dPCA Parameters
# D) Write dPCA Kalman Filtering Function
# E) Predict (external)


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
  list.of.packages <- c("beepr","foreach","FKF","ggplot2","dplyr","magrittr","gridExtra","mvtnorm","ggplot2","YieldCurve","Meucci","xtable","Matrix","MASS","glmnet")
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



# Initialize Sample from NS Manifold
prin.comps<-princomp(bond.data[c(1:100),])
dPCA.samples<-cbind(NS.samples[,1],prin.comps$loadings[,c(1:N.factors)])
rownames(dPCA.samples)<-Maturities.Grid



# Visualize Head
head(dPCA.samples)
# Visualize Tail
tail(dPCA.samples)
# Test Dimensions
dim(dPCA.samples)

# Create list of betas for factor model
require(glmnet)
# Initialize betas
betas.mat<-matrix(NA,nrow=(-1+nrow(bond.data)),ncol=(N.factors+1))

# Initialize Progress Bar
pb <- txtProgressBar(0, (-1+nrow(bond.data)), style = 3)

# Populate Betas Matrix
for(t.j in 1:(-1+nrow(bond.data))){
test<-cv.glmnet(x=dPCA.samples,y=bond.data[t.j,])
betas.mat[t.j,]<-as.vector(glmnet(x=dPCA.samples,y=bond.data[t.j,],lambda = test$lambda.min)$beta)
# Update Progress Bar
setTxtProgressBar(pb, t.j)
}
# Close Progress Bar
close(pb)


##### Fit OU Process to Betas

Estimated.OU.Parameters.Betas.PCA<-FitOrnsteinUhlenbeck(betas.mat,1)
Estimated.OU.Parameters.Betas.PCA

# Add intercept & Remove Dates Column
dPCA.samples<-cbind(rep(1,nrow(dPCA.samples)),dPCA.samples[,-1])

#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#                  D) Write Filter
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
#-############ -----------------------------------################-#
# Random initialization of parameters
lambda_init <- runif(m, min=-1.0, max=1.0)
measurement_err_init.dPCA <- diag(runif(d, 
                              min=0.0, 
                              max=sqrt(sum((eigen(Estimated.OU.Parameters.Betas.PCA$Sigma)$values)^2))
                              ))


# Kalman NS
dPCA_KF <- function(observations)
{ 
  # initial state variable (a0: m x 1)
  r_init.dPCA <- as.vector(Estimated.OU.Parameters.Betas.PCA$Mu)
  
  # variance of state variable (P0: m x m)
  P_init.dPCA <- t(Estimated.OU.Parameters.Betas.PCA$Sigma)%*%Estimated.OU.Parameters.Betas.PCA$Sigma
  P_init.dPCA <- .5*P_init.dPCA
  
  # intercept of state transition equation (dt: m x 1)
  C.dPCA <- (Estimated.OU.Parameters.Betas.PCA$Theta)%*%(Estimated.OU.Parameters.Betas.PCA$Mu)
  
  # factor of transition equation (Tt: m x m x 1)
  Tt.dPCA <- -Estimated.OU.Parameters.Betas.PCA$Theta
  
  # factor of measurement equation (Zt: d x m x 1)
  B.dPCA <- dPCA.samples
  
  # intercept of measurement equation (ct: d x 1)
  A.dPCA <- rep(0,d)
  
  
  # variance of innovations of transition (HHt: m x m x 1)
  Q.dPCA <- t(Estimated.OU.Parameters.Betas.PCA$Sigma)%*%Estimated.OU.Parameters.Betas.PCA$Sigma
  
  # variance of measurement error (GGt: d x d x 1)
  R.dPCA <- array(measurement_err_init.dPCA,dim=c(d,d,1))
  
  filtered_process <- fkf(a0=r_init.dPCA, P0=P_init.dPCA, dt=C.dPCA, ct=as.matrix(A.dPCA), Tt=Tt.dPCA, 
                          Zt=B.dPCA, HHt=Q.dPCA, GGt=R.dPCA, yt=as.matrix(observations))
  return(filtered_process)
}




#---------#-############-#-------#
#--------------------------------#
#         E: Prediction          #
#--------------------------------#
#---------#-############-#-------#



# SEE PREDICTION FILE!!