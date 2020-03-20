################################################-#
#         Deep Arbitrage-Free Learning           #
################################################-#
# Implementation for the paper of:
# Anastasis Kratsios - ETH, Zurich, Switzerland
# Cody Hyndman       - Concordia, Montreal, Canada
# 2019


#--------------------------------------------------------#
# What this Algorithm Does:
#--------------------------------------------------------#
# This implementation is for the Nelson-Siegel family
# However, code can be easliy modified to suit other 
# factor models for the forward-rate curve.
#--------------------------------------------------------#


#-------------------------------------------------#
# How the Algorithm's Workflow is Designed
#-------------------------------------------------#
# -1) Initializes Packages

# A) Estimating FRC Model Parameters


# B) Nelson-Siegel Curve Generator
# 1) Write Grid of Maturities (important for integration in last step)
# 2) Simulates points on the NS manifold (for training in phase B)
# 3) Build Graph

# C) Pre-Processing Data
# 1) Standardize sample points from NS manifold
# 2) Partition training and testing samples (to evalulate goodness of fit)

# D)
# Overview of deep-learning workflow
# 1) Initialize and Load Auxiliary Functions
# 2) Pre-process Data
# -) (Optional) Visualize Network  - Not in use atm
# 3) Build Model by reading in and then looping over generic depth and height of network
# 4) Train Network
# 5) Evaluate Network

# E) Define Arbitrage-Penalty
# 1) Define Reservoir
# 2) Define Integrated Reservoir
# 3) Write Arbitrage-Penalty Function

# F) Perform Arbitrage-Free Regularization
# 1) Minimize new redout using AF Regularization


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#--------------------------------------------------------#
#                  ! ! ! Read Me ! ! !
#--------------------------------------------------------#
# If used on a different machine, certain paths should be 
# Change accordingly so that things work.  
#--------------------------------------------------------#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#















#---------#       #-#        #---------#
#--------------------------------------#
# Keras based Deep Learning Parameters #
#--------------------------------------#
#---------#       #-#        #---------#

# Note:
# These parameters design and fit FFNN to Nelson-Siegel Model
# Then truncate the last layer to create a high-dimensional reservoir.
# This reservoir then serves as a "smart" random initialization for the 
# Arbitrage-Free Regularization Procedure in step E.

#  Network and Reservoir Parameters 
Depth<-5
Height<-10^1
Sparsity.percentrage<-.5 #(~~Dropout/layer)
# SGD Parameters
N.epochs<-10^2
validation.proportion<-.1
# Nelson-Siegel Sampling Parameters
maturity  <- Maturities.Grid<-c(.5,seq(from=1,to=30,by=1)) #<-c(0.25,0.5,1,2,5,10)
#maturity  <- Maturities.Grid<-c(0.25,0.5,1,2,5,10)
N.mat<-length(Maturities.Grid)
# Training/Testing Ratio
  Train.to.Test.ratio<-.1
# Regularization Parameters
l1.pen<-.01
l2.pen<-.02



#---------#         #-#          #---------#
#------------------------------------------#
# Arbitrage-Free Regularization Parameters #
#------------------------------------------#
#---------#         #-#          #---------#
# Penalization Amount
lambda.AFreg.param.homotopic<-0.9999 # Between 0 and 1; 0 no change 1 (AF)
# Finite-Difference Step for AF-Loss Function
FD.stepsize<-(10^-3)
# Stochastic Gradient Descent Parameters 
learning.rate.SGD<-5
N.iterations.SGD<-10^1
# Decay Tuning Parameter
kappa<- 1
p.Lp<-2
# Parameter (Re)Calibration?
RECALIBRATEPARAMETERS.Q<-TRUE

#---------------------#
#  General Parameters #
#---------------------#
# Set working directory
getwd()
Default.directory<- "/scratch/users/kratsioa/Dropbox/Numerics/Deep_Arbitrage_Free_Regularization/Deep_Arbitrage_Free_Regularization" #If unknown, use: getwd()


# Notes to user:
# Too sparse = poor results
# Too deep = poor convergence rate
# Too tall = unstable predictive performance
# Too much regularization = linear
# Regularization should be low as the goal is to fit to the NS model and not to predict it









print("-1: Loading Packages, Auxiliary Functions, and Datasets.")
#########
#########
#########
#########
#########
#########
###############################################################
#                                                             #
# -1  Loading Packages, Auxiliary Functions, and Datasets     #
#                                                             #
###############################################################
#########
#########
#########
#########
#########
#########
# Load packages
### INITIALIZATION PORTION
#-############ -----------------------------------################-#
# Load required packages if available...if not install - BEGIN
#-############ -----------------------------------################-#
# Some packages should install this package manually and separately
### The Script is in the file: Installations.R

if(!exists("NO.CHECK_Q")){
  # Check for missing packages
  list.of.packages <- c("beepr","keras","dplyr","magrittr","neuralnet","mvtnorm","ggplot2","YieldCurve","Meucci","xtable")
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
}
#
#-############ -----------------------------------################-#
# Load required packages if available...if not install - END
#-############ -----------------------------------################-#

if((!exists("long.term.level") | RECALIBRATEPARAMETERS.Q)){
  #################################################
  #################################################
  #################################################
  ###    A) Estimation of Dynamic Parameters    ###
  #################################################
  #################################################
  #################################################
  # Note this is performed externally within the NS_Kalman file!!!! Load that first
  
  # (Estimated) Long Term Mean
  long.term.level<-Estimated.OU.Parameters.Betas.PCA$Mu
  # (Estimated) Mean-Reversion Speed
  mean.reversion.speed<-Estimated.OU.Parameters.Betas.PCA$Theta
  # (Estimated)  Volatility Matrix
  Vol.mat.estim<-Estimated.OU.Parameters.Betas.PCA$Sigma
  
  # Returns Estimated Parameters of Multivariate OU Process (in notation of paper)
  # (Affine) Drift Process Parameters
  gam<-as.vector(mean.reversion.speed%*%long.term.level)
  Gam<-as.matrix(mean.reversion.speed)
  # (Affine) Variance Process Parameters
  alpha<-as.matrix(t(Vol.mat.estim)%*%Vol.mat.estim)
  Alpha<-0 #Note!!!! Alpha matrix is zero for implementation because of OU dynamics!

}# END PARAMETER ESTIMATIONS







print("B: Sampling from FRC Manifold.")
#########
#########
#########
#########
#########
#########
#####################################
#                                   #
#  B  Sampling from FRC Manifold    #
#                                   #
#####################################
#########
#########
#########
#########
#########
#########

## PARAMETERS
###
# Evaludate Distribution of means and loadings
###

# Packages
#----------------------------#

# Auxiliary Functions
#----------------------------#
# Nelson-Siegel Curves
NS.curves<-function(T,lambda=1){
  c(1,exp(-T*lambda),T*exp(-T*lambda))
}

# Transforms homotopic parameterization of Lambda parameter into standard parameterization
lambda.AFreg.param<- lambda.AFreg.param.homotopic/(1-lambda.AFreg.param.homotopic)

# Build data-set to train Network on
# ---------------------------------------# 

# Sample from NS Manifold
### --------------- ----- --- ###

# Initialize Sample from PCA Manifold
dPCA.samples<-cbind(Maturities.Grid,dPCA.samples)

# Visualize Head
head(dPCA.samples)
# Visualize Tail
tail(dPCA.samples)
# Test Dimensions
dim(dPCA.samples)








print("C: Pre-Processing Data.")
#########
#########
#########
#########
#########
#########
#####################################
#                                   #
#  C)    Pre-Processing Data        #
#                                   #
#####################################
#########
#########
#########
#########
#########
#########



# Read Data and Format Data
#-------------------------------#
# Read Data
data<-dPCA.samples
# Format Data
data<-as.matrix(data)
# Pre-Processing Data
########
#Partition Data
set.seed(123)
ind<-sample(2,nrow(data),replace = T,prob=c(.7,.3))
# Training
training<-data[,1]
trainingtarget<-data[,-1]
# Normalize Data
m<-mean(training)#colMeans(training)
s<-sd(training)#apply(training,2,sd)
training<-scale(training,center=m,scale=s)









print("D: Learning FRC Manifold Structure: ie: Learning Good Initialization for AF-Regularization Problem!")
#########
#########
#########
#########
#########
#########
  ##################################
 # #                              # #
# D)     Learning FRC Manifold       #
 # #                              # #
  ##################################
#########
#########
#########
#########
#########
#########



# Building Deep FF-Network
#--------------------------------#
model<-keras_model_sequential()
# Define bulk of the network
model %>% layer_dense(units=Height,activation = "relu",input_shape = c(1))

for(i in 1:Depth){
    model %>% layer_dense(units=Height,activation = "relu",input_shape = c(4),kernel_regularizer = regularizer_l1_l2(l1 = l1.pen, l2 = l2.pen)) 
}
# Readout Layer
model %>% layer_dense(units=4)

# Compile
model %>% keras::compile(loss="mse",
                  optimizer="adam",
                  metrics="mse")

# Fit Model
fittedmodel<- model %>%
  keras::fit(training,
      trainingtarget,
      epochs=N.epochs,
      batch_size=35,
      validation_split=validation.proportion)

# Notify: Training Complete
beep()



print("E: Defining The Arbitrage-Free Regularization Problem's Objective-Function")
#################################################
#################################################
#################################################
###   E) Arbitrage-Free Regularization        ###
#################################################
#################################################
#################################################

# Create Reservoir (remove final layer)
pop_layer(model)
Reservoir.BU<- model %>% predict(data[,1])
# Normalize Reservoir (mathematically makes no difference but numerically adds so much stability)
mr<-colMeans(Reservoir.BU)
sr<-apply(Reservoir.BU,2,sd)
Reservoir<-Reservoir.BU
Reservoir[,!sr==0]<-scale(Reservoir.BU[,!sr==0],center = mr[!sr==0],scale=sr[!sr==0])
Reservoir<- cbind(dPCA.samples[,-1],Reservoir) # Add intercept and NS basis functions for numerical stability
colnames(Reservoir)<-c()


## Build Integrated Reservoir
# Initialized Integrated Columns Matrix
Integrated.Reservoir<-apply(Reservoir,2,cumsum)


##########################################
#        Define Objective Function       #
##########################################
# Note: the multiplicative factor of lambda.AFreg will be incorporated directly into the SGD step.
AFReg_Loss<-function(W.dummy,lambda.AFreg=2){
  #---------------------------------------------#  
  #---------------------------------------------#
  #             Goodness of Fit                 #
  #---------------------------------------------#
  #---------------------------------------------#
  # Evaluates 0th Order Fit
  goodness.of.fit<-(dPCA.samples[,-1]-Reservoir%*%W.dummy)
  goodness.of.fit<-as.vector(goodness.of.fit)
  goodness.of.fit<-sqrt(sum(goodness.of.fit^p.Lp))
  #---------------------------------------------#  
  #---------------------------------------------#
  #               Arb.Penalty                   #
  #---------------------------------------------#
  #---------------------------------------------#
  # Write Arbitrage-Penlaty
  
  # Initialize penalty value
  pen.val<-0
  
  # Define AF Penalty Function
  #----------------------------#
  # 0th term separately
  #--------------------------
  zer.val<-0
  for(u in 1:nrow(Reservoir)){
    # compute first term
    zer.val.iter<-as.numeric(Reservoir[u,]%*%W.dummy[,1]-Reservoir[1,]%*%W.dummy[,1]-gam%*%as.vector(Reservoir[u,]%*%W.dummy))
    # compute last term
    aFWWF0th<-0
    for(i in 1:length(gam)){
      for(j in 1:length(gam)){
        aFWWF0th<- aFWWF0th + alpha[i,j]*as.numeric(Reservoir[j,]%*%W.dummy%*%t(Reservoir[j,]%*%W.dummy))
      }
    }
    zer.val<- zer.val + exp(-abs(u)^kappa)*(zer.val.iter+ .5*aFWWF0th)^2 #Note: weighting helps numerical stability but doesn't affect validity of penalty 
    
    #----------------------------#
    # kth term separately
    #--------------------------
    k.th.val.total<-0
    for(k in 2:length(gam)){
      # compute first term 
      k.th.val.iter<-(as.numeric(Reservoir[u,]%*%W.dummy[,k]-Reservoir[1,]%*%W.dummy[,k]-Gam[k,]%*%as.vector(Reservoir[u,]%*%W.dummy)))^2
      k.th.val.total<-k.th.val.total +exp(-abs(u)^kappa)*k.th.val.iter
    }
    # Update penalth value at uth maturity
    pen.val<-((zer.val+k.th.val.total)^(1/lambda.AFreg))
  }
  
  # Adds constants
  pen.val<-pen.val * (1/((gamma(1+(1/kappa)))^(1/lambda.AFreg)))
  
  # Combine Fit Loss and Arbitrage Penalty Loss Values
  loss.value.output<-as.numeric(pen.val+goodness.of.fit)
  # Return Value of Loss-Function
  return(loss.value.output)
}




print("F: Performing SGD to Optimize Loss Function; initialized using Pre-trained Deep Network")
#################################################
#-----------------------------------------------#
#################################################
# F)       Perform SGD on Objective Function    #
print("F: Performing SGD on Objective Function")
#################################################
#-----------------------------------------------#
#################################################
#################################################
## Externally-based/Auxiliary SGD Parameters
N.basis<-(ncol(dPCA.samples)-1)
set.seed(0)
#-------------#
####  SGD  ####
#-------------#
# SGD Iterations
n = N.iterations.SGD

#
# Initialize near original model
pars<-pars.ini<- rbind(diag(4),matrix(0,ncol=N.basis,nrow=(ncol(Reservoir)-N.basis)))
m.pars.ini<-colMeans(pars)
s.pars.ini<-apply(pars,2,sd)
# Initialize Target Data
Target.AF.Reg<-dPCA.samples[,-1];colnames(Target.AF.Reg)<-c();colnames(Reservoir)<-c()


# Initialize SGD History
output = matrix(nrow=n, ncol=length(pars))
output[1,] = pars

# Prepare SGD
sgd = SGD$new(initialLearningRate=learning.rate.SGD)
# Initialize Progress Bar
pb <- txtProgressBar(0, n, style = 3)

# Save Scale and Center of Original Data
m.Af<-colMeans(Target.AF.Reg)
s.Af<-apply(Target.AF.Reg,2,sd)

# Performs SGD
for (i in 2:n) {
  #---------------------------------------------------------#
  # Compute Gradient at current position of parameters
  #---------------------------------------------------------#
  # Initialize Current W matrix
  W.dummy<-matrix(data=pars,ncol=N.basis,nrow=ncol(Reservoir))
  
  # Compute Goodness of Fit Gradient
  grad.Fit.loss<- 2*t(Reservoir)%*%(Target.AF.Reg-Reservoir%*%W.dummy)
  
  # Approximate AF-Loss's Gradient
  grad.AF.loss<-pars
  ## Gradient
  for(k in 1:length(pars)){
    # Perturbs Parameter
    grad.AF.loss[k]<-grad.AF.loss[k]+FD.stepsize
    # Forms into Matrix
    W.dummy.perturbed<-matrix(data=grad.AF.loss,ncol=N.basis,nrow=ncol(Reservoir))
    W.dummy<-matrix(data=pars,ncol=N.basis,nrow=ncol(Reservoir))
    # Performs Derivative Computation
    ith.derivative<-AFReg_Loss(W.dummy.perturbed,lambda.AFreg=lambda.AFreg.param)- AFReg_Loss(W.dummy,lambda.AFreg=(lambda.AFreg.param))
    # Note: Adding 10 doesn't make an asymptotic difference but it forbodes explosions near lambda~~0 (which really isn't the goal of AF since lambda should be large)
    ith.derivative<-ith.derivative*FD.stepsize
    # Updates ith gradient entry
    grad.AF.loss[k]<-ith.derivative
  }# END:END Gradient Computation
  

  
  # Compute total Gradient
  gradest = lambda.AFreg.param.homotopic*as.vector(grad.AF.loss)+ (1-lambda.AFreg.param.homotopic)*as.vector((pars - pars.ini))
  pars = pars + sgd$step(gradest) 
  
  
  
  # Correct numerical problems (if any: rare)
  W.dummy.temp<-matrix(data=pars,ncol=N.basis,nrow=ncol(Reservoir))
  W.dummy.temp[is.infinite(W.dummy.temp)]<-max(W.dummy.temp[!is.infinite(W.dummy.temp)])
  W.dummy.temp[is.nan(W.dummy.temp)]<-0
  W.dummy.temp[is.na(W.dummy.temp)]<-0
  ## Renormalize pars
  for(i.scale in 1:ncol(W.dummy.temp)){
    W.dummy.temp[,i.scale]<-W.dummy.temp[,i.scale]/(max(W.dummy.temp[,i.scale])-min(W.dummy.temp[,i.scale]))
  }
  # Vectorize
  pars<-as.vector(W.dummy.temp)
    
  # Vectorize
  pars<-as.vector(W.dummy.temp)
  #
  output[i,] = pars 
  # %%% #
  # Update Progress bar (inloop)
  setTxtProgressBar(pb, i)
}
# Close Progress Bar
close(pb)
W.optimal<-W.dummy.temp
AF.Reg.Factors.PCA<-Reservoir%*%W.optimal


