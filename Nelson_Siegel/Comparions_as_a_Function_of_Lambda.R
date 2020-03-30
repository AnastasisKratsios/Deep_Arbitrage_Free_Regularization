######################---#
# Learning Stage
######################---#
# Initialize Lambda Sequence
Lambda_seq = seq(from=(10^(-7)),to=(1- 10^(-7)),length.out = (10^2))


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

# Set Working Directory
wd<-"/scratch/users/kratsioa/Dropbox/Numerics/Deep_Arbitrage_Free_Regularization/Full_Version/Datasets"
setwd(wd)



#---------#         #-#          #---------#
#------------------------------------------#
# Arbitrage-Free Regularization Parameters #
#------------------------------------------#
#---------#         #-#          #---------#
# Penalization Amount
# Finite-Difference Step for AF-Loss Function
FD.stepsize<-(10^-3)
# Stochastic Gradient Descent Parameters 
learning.rate.SGD<-5
N.iterations.SGD<-10^2
# Decay Tuning Parameter
kappa<- 1
p.Lp<-2
# Parameter (Re)Calibration?
RECALIBRATEPARAMETERS.Q<-TRUE


# -------------------- #
# Initialization
# -------------------- #
# Load initializations from other source files...
if(!exists("NO.CHECK_Q")){
  # Nelson-Siegel Models
  source("Initializations.R") # Initializations
  source("Vasicek_HJM_Kalman_full.R") # Nelson-Siegel Model
  source("NS_Kalman_full.R") # Nelson-Siegel Model
  source("AFNS_full.R") # AFNS Correction Term (Christiensen et al.)
  source("SGD.R") # Loads SGD Code
  Missing.Q = TRUE
  source("Deep_Arbitrage_Free_Regularization_full.R") # AF-Regularization of AFNS
  source(paste(wd,"AFREG_NS_Kalman_full.R",sep="/")) # Compiles Predictive Algorithm
}


#---------------------#
#  General Parameters #
#---------------------#
# Set working directory
getwd()
Default.directory<- "/scratch/users/kratsioa/Dropbox/Numerics/Deep_Arbitrage_Free_Regularization/Deep_Arbitrage_Free_Regularization" #If unknown, use: getwd()


# Initialize Readings
Lambda_comparisons<-c() # INITIALIZE MSE Comparisons
Lambda_comparisons_Bond_Prices<-c() # INITIALIZE Bond-Price Comparisons
Lambda_to_date<-c()


for(i.lambda in 1:length(Lambda_seq)){
## Set Working Lambda
  lambda.AFreg.param.homotopic<-Lambda_seq[i.lambda]


#################################################
#################################################
#################################################
###   E) Arbitrage-Free Regularization        ###
#################################################
#################################################
#################################################



##########################################
#        Define Objective Function       #
##########################################
# Note: the multiplicative factor of lambda.AFreg will be incorporated directly into the SGD step.
AFReg_Loss<-function(W.dummy,lambda.AFreg=(lambda.AFreg.param.homotopic/(1-lambda.AFreg.param.homotopic))){
  #---------------------------------------------#  
  #---------------------------------------------#
  #             Goodness of Fit                 #
  #---------------------------------------------#
  #---------------------------------------------#
  # Evaluates 0th Order Fit
  goodness.of.fit<-(NS.samples[,-1]-Reservoir%*%W.dummy)
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
N.basis<-(ncol(NS.samples)-1)
set.seed(0)
#-------------#
####  SGD  ####
#-------------#
# SGD Iterations
n = N.iterations.SGD

#
# Initialize near original model
pars<-pars.ini<- rbind(diag(3),matrix(0,ncol=N.basis,nrow=(ncol(Reservoir)-N.basis)))
m.pars.ini<-colMeans(pars)
s.pars.ini<-apply(pars,2,sd)
# Initialize Target Data
Target.AF.Reg<-NS.samples[,-1];colnames(Target.AF.Reg)<-c();colnames(Reservoir)<-c()


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
AF.Reg.Factors<-Reservoir%*%W.optimal













######################------------#
# Compiling Predictive Algorithm
######################------------#
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


# Make Prediction








###############----------#
# Update Bond Prices
###############----------#




# Initialize one-day ahead predicted bond prices!
bond.prices.NS.HJM.Kalman.errors<-matrix(data=NA,nrow=N,ncol=d)
colnames(bond.prices.NS.HJM.Kalman.errors)<-Maturities.Grid
bond.prices.NS.AF.Reg.HJM.Kalman<-bond.prices.NS.HJM.Kalman<-bond.prices.NS.HJM.Kalman.errors
# Initialize Bond Pricing Errors
bond.prices.NS.AF.Reg.HJM.Kalman.errors<-bond.prices.NS.HJM.Kalman.errors


# Populate Matrix of one-day ahead predicted bond prices
for(t.i in 1:(nrow(bond.prices.NS.AF.Reg.HJM.Kalman.errors))){
  # Write realised next day bond price
  day.ahead.bond.price<-unlist(as.vector(exp(
    cumsum((Maturities.Grid)*(bond.data[t.i,])))
    ))
  
  # Filter betas 
  KF_NS_loop<-NS_KF(observations = bond.data[t.i,])
  AF.Reg.KF_NS_loop<-AFReg_NS_KF(observations = bond.data[t.i,])

  # Write Predicted FRC for NS Models
  FRC.NS.loop<-as.vector(NS.samples[,-1]%*%(KF_NS_loop$att))
  FRC.NS.AF.Reg.loop<-as.vector(AF.Reg.Factors%*%(AF.Reg.KF_NS_loop$att))
  
  
  # Write Bond prices
  prices.NS.loop<-exp(cumsum(FRC.NS.loop))
  prices.NS.AF.Reg.loop<-exp(cumsum(FRC.NS.AF.Reg.loop))
  bond.prices.NS.HJM.Kalman[t.i,]<- prices.NS.loop
  bond.prices.NS.AF.Reg.HJM.Kalman[t.i,]<- prices.NS.AF.Reg.loop
  # Write Errors
  bond.prices.NS.HJM.Kalman.errors[t.i,]<- (prices.NS.loop-day.ahead.bond.price)
  bond.prices.NS.AF.Reg.HJM.Kalman.errors[t.i,]<- (prices.NS.AF.Reg.loop-day.ahead.bond.price)
  
}# END Prediction Step








# Save Readings
#--------------#
# MSEs
MSE.NS.errors<-colMeans(bond.prices.NS.AF.Reg.HJM.Kalman.errors^2)
MSE.NS.errors_bond<-colMeans(bond.prices.NS.AF.Reg.HJM.Kalman)


# Update Readings
Lambda_comparisons<-rbind(Lambda_comparisons,MSE.NS.errors)
Lambda_comparisons_Bond_Prices<-rbind(Lambda_comparisons_Bond_Prices,MSE.NS.errors_bond)
Lambda_to_date<-sort(c(Lambda_to_date,lambda.AFreg.param.homotopic))

# Update User
print(MSE.NS.errors)
print((i.lambda/length(Lambda_seq)))
}






# Format Data
#--------------#
rownames(Lambda_comparisons)<-Lambda_to_date
rownames(Lambda_comparisons_Bond_Prices)<-Lambda_seq

# Write Tables
#--------------#
#cat("\014")



# Generate Plots
#----------------#
# Write Data_Frame
Lambda_comparisons = as.data.frame(Lambda_comparisons)
head(Lambda_comparisons)
Lambda_comparisons_Bond_Prices = as.data.frame(Lambda_comparisons_Bond_Prices)
head(Lambda_comparisons_Bond_Prices)

#----------------------------------------------------------------#
# Compare Change in Day-Ahead MSE
#----------------------------------------------------------------#
Colorss1 <- scales::seq_gradient_pal("red", "blue", "Lab")(seq(0,1,length.out=nrow(Lambda_comparisons_dPCA)))

data_plot = as.data.frame(cbind(Maturities.Grid,t(Lambda_comparisons)))
colnames(data_plot) = c("Maturities",round(c(Lambda_seq),4))
colnames(data_plot)[length(colnames(data_plot))] = Lambda_seq[length(Lambda_seq)]
test_data_long <- melt(data_plot, id="Maturities")  # convert to long format

# Generate (gg)plot(2)
ggplot(data=test_data_long,
       aes(x=Maturities, y=value, colour=variable)) +
      geom_line() +
  labs(x ="Maturities", y = "Day-Ahead Prediction MSE", color="Lambdas") +
  #scale_color_hue(l=40, c=35)+
  scale_colour_manual(values=Colorss1)




#----------------------------------------------------------------#
# Compare Change in Predeicted Prices
#----------------------------------------------------------------#
Colorss <- scales::seq_gradient_pal("red", "blue", "Lab")(seq(0,1,length.out=nrow(Lambda_comparisons)))


# Write Data_Frame
Lambda_comparisons_Bond_Prices = as.data.frame(Lambda_comparisons_Bond_Prices)
data_plot2 = as.data.frame(cbind(Maturities.Grid,Lambda_comparisons_Bond_Prices))
colnames(data_plot2)[1] = "Maturities"
test_data_long2 <- melt(data_plot2, id="Maturities")  # convert to long format

# Generate (gg)plot(2)
ggplot(data=test_data_long2,
       aes(x=Maturities, y=value, colour=variable)) +
      geom_line() +
  labs(x ="Maturities", y = "Average Day-Ahead Bond Yield",color="Lambda") +
  theme(legend.position="none")+
  scale_colour_manual(values=Colorss)


# Save Output Data
#-----------------#
# Save Table
save(Lambda_comparisons, file = "Lambda_Comparisons_File.RData")
save(Lambda_comparisons_Bond_Prices, file = "Lambda_Comparisons_File_Bond_Prices.RData")