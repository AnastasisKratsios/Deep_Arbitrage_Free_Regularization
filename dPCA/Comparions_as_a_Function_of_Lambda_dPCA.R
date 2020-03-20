######################---#
# Learning Stage
######################---#
# Initialize Lambda Sequence
Lambda_seq = seq(from=(10^(-7)),to=(1- 10^(-7)),length.out = 3)


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
  source(paste(wd,"Initializations.R",sep="/")) # Initialize
  source(paste(wd,"PCA_Kalman_full.R",sep="/")) # Nelson-Siegel Model
  source(paste(wd,"SGD.R",sep="/")) # Loads SGD Code
  source(paste(wd,"Deep_Arbitrage_Free_Regularization_full_PCA.R",sep="/")) # AF-Regularization of AFNS
  source(paste(wd,"AFREG_dPCA_Kalman_full.R",sep="/")) # Compiles Predictive Algorithm
}


#---------------------#
#  General Parameters #
#---------------------#
# Set working directory
getwd()
Default.directory<- "/scratch/users/kratsioa/Dropbox/Numerics/Deep_Arbitrage_Free_Regularization/Deep_Arbitrage_Free_Regularization/dPCA" #If unknown, use: getwd()


# Initialize Readings
Lambda_comparisons_dPCA<-c() # INITIALIZE MSE Comparisons
Lambda_comparisons_Bond_Prices_dPCA<-c() # INITIALIZE Bond-Price Comparisons
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
AFReg_Loss<-function(W.dummy,lambda.AFreg=2){
    #---------------------------------------------#  
    #---------------------------------------------#
    #             Goodness of Fit                 #
    #---------------------------------------------#
    #---------------------------------------------#
    # Evaluates 0th Order Fit
    goodness.of.fit<-(dPCA.samples-Reservoir%*%W.dummy)
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
  N.basis<-ncol(dPCA.samples)
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
  Target.AF.Reg<-dPCA.samples;colnames(Target.AF.Reg)<-c();colnames(Reservoir)<-c()
  
  
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


# Make Prediction








###############----------#
# Update Bond Prices
###############----------#
  # Initialize one-day ahead predicted bond prices!
  bond.prices.dPCA_AF.high.Kalman<-matrix(data=NA,nrow=N,ncol=d)
  colnames(bond.prices.dPCA_AF.high.Kalman)<-Maturities.Grid
  
  bond.prices.dPCA_AF.high.Kalman<-bond.prices.dPCA.Kalman<-bond.prices.NS.AF.Reg.HJM.Kalman<-bond.prices.AFNS.HJM.Kalman<-bond.prices.NS.HJM.Kalman<-bond.prices.Vasicek.HJM.Kalman
  # Initialize Bond Pricing Errors
  bond.prices.dPCA_AF.high.Kalman.error<-bond.prices.dPCA.Kalman.error<-bond.prices.NS.AF.Reg.HJM.Kalman.errors<-bond.prices.AFNS.HJM.Kalman.errors<-bond.prices.NS.HJM.Kalman.errors<-bond.prices.Vasicek.HJM.Kalman.errors<-bond.prices.Vasicek.HJM.Kalman
  
  
  # Populate Matrix of one-day ahead predicted bond prices
  for(t.i in 1:(nrow(bond.prices.Vasicek.HJM.Kalman.errors))){
    
    # Filter betas for dPCA
    dPCA_loop<-dPCA_KF(observations = bond.data[t.i,])
    dPCA.high_loop<-dPCA_KF(observations = bond.data[t.i,])
    AF_dPCA.high_loop<-AFReg_dPCA_KF(observations = bond.data[t.i,])
    
    
    
    # Write for dPCA Models
    FRC.dPCA.loop<-as.vector(dPCA.samples%*%(dPCA_loop$att))
    FRC.dPCA.AF.Reg.loop<-as.vector(AF.Reg.Factors.PCA%*%(AF_dPCA.high_loop$att))
    
    
    # Write Bond prices
    prices.dPCA.loop<-exp(-cumsum(FRC.dPCA.loop))
    prices.dPCA_REG.loop<-exp(-cumsum(FRC.NS.loop))
    
    
    # Write Predicted Bond Prices
    bond.prices.dPCA.Kalman[t.i,]<-prices.dPCA.loop
    bond.prices.dPCA_AF.high.Kalman[t.i,]<-prices.dPCA_REG.loop
    
    # Evaluate Prediction Errors (1-day ahead)
    # Write realised next day bond price
    day.ahead.bond.price<-unlist(as.vector(exp(-(Maturities.Grid)*(bond.data[t.i,]))))
    # Evaluate Error
    bond.prices.dPCA.Kalman.error[t.i,]<- (prices.dPCA.loop-day.ahead.bond.price)
    bond.prices.dPCA_AF.high.Kalman.error[t.i,]<- (prices.dPCA_REG.loop-day.ahead.bond.price)
    
    
  }# END Prediction Step








# Save Readings
#--------------#
# MSEs
MSE.dPCA.errors<-colMeans(bond.prices.dPCA_AF.high.Kalman.error^2)
MSE.dPCA.errors_bond<-colMeans(bond.prices.dPCA_AF.high.Kalman)


# Update Readings
Lambda_comparisons_dPCA<-rbind(Lambda_comparisons_dPCA,MSE.NS.errors)
Lambda_comparisons_Bond_Prices_dPCA<-rbind(Lambda_comparisons_Bond_Prices_dPCA,MSE.NS.errors_bond)
Lambda_to_date<-sort(c(Lambda_to_date,lambda.AFreg.param.homotopic))

# Update User
print(MSE.NS.errors)
print((i.lambda/length(Lambda_seq)))
}






# Format Data
#--------------#
rownames(Lambda_comparisons_dPCA)<-Lambda_to_date
rownames(Lambda_comparisons_Bond_Prices_dPCA)<-Lambda_to_date


# Write Tables
#--------------#
#cat("\014")



# Generate Plots
#----------------#
# Write Data_Frame
Lambda_comparisons_dPCA = as.data.frame(Lambda_comparisons_dPCA)
head(Lambda_comparisons_dPCA)
Lambda_comparisons_Bond_Prices_dPCA = as.data.frame(Lambda_comparisons_Bond_Prices_dPCA)
head(Lambda_comparisons_Bond_Prices_dPCA)

#----------------------------------------------------------------#
# Compare Change in Day-Ahead MSE
#----------------------------------------------------------------#
data_plot_dPCA = as.data.frame(cbind(Lambda_to_date,Lambda_comparisons_dPCA))
colnames(data_plot_dPCA) = c("Lambdas",Maturities.Grid)
test_data_long_dPCA <- melt(test_data_long_dPCA, id="Lambdas")  # convert to long format

# Generate (gg)plot(2)
ggplot(data=test_data_long_dPCA,
       aes(x=Lambdas, y=value, colour=variable)) +
      geom_line() +
  labs(x ="Lambdas", y = "Prediction MSE", color="Maturities") +
  scale_color_hue(l=40, c=35)




#----------------------------------------------------------------#
# Compare Change in Predeicted Prices
#----------------------------------------------------------------#
# Write Data_Frame
Lambda_comparisons_Bond_Prices_dPCA = as.data.frame(Lambda_comparisons_Bond_Prices_dPCA)
head(Lambda_comparisons)
data_plot2_dPCA = as.data.frame(cbind(Lambda_to_date,Lambda_comparisons_Bond_Prices_dPCA))
colnames(data_plot2_dPCA) = c("Lambdas",Maturities.Grid)
test_data_long2 <- melt(data_plot2, id="Lambdas")  # convert to long format

# Generate (gg)plot(2)
ggplot(data=test_data_long2_dPCA,
       aes(x=Lambdas, y=value, colour=variable)) +
      geom_line() +
  labs(x ="Lambdas", y = "Day-Ahead Bond Price", color="Maturities") +
  scale_color_hue(l=40, c=35)


# Save Output Data
#-----------------#
# Save Table
save(Lambda_comparisons_dPCA, file = "Lambda_Comparisons_File_dPCA.RData")
save(Lambda_comparisons_Bond_Prices_dPCA, file = "Lambda_Comparisons_File_Bond_Prices_dPCA.RData")