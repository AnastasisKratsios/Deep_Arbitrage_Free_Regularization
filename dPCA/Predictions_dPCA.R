# Initialize one-day ahead predicted bond prices!
bond.prices.Vasicek.HJM.Kalman<-matrix(data=NA,nrow=N,ncol=d)
colnames(bond.prices.Vasicek.HJM.Kalman)<-Maturities.Grid
#rownames(bond.prices.Vasicek.HJM.Kalman)<-round(dates[-1],1)
bond.prices.dPCA_AF.high.Kalman<-bond.prices.dPCA.Kalman<-bond.prices.NS.AF.Reg.HJM.Kalman<-bond.prices.AFNS.HJM.Kalman<-bond.prices.NS.HJM.Kalman<-bond.prices.Vasicek.HJM.Kalman
# Initialize Bond Pricing Errors
bond.prices.dPCA_AF.high.Kalman.error<-bond.prices.dPCA.Kalman.error<-bond.prices.NS.AF.Reg.HJM.Kalman.errors<-bond.prices.AFNS.HJM.Kalman.errors<-bond.prices.NS.HJM.Kalman.errors<-bond.prices.Vasicek.HJM.Kalman.errors<-bond.prices.Vasicek.HJM.Kalman


# Populate Matrix of one-day ahead predicted bond prices
for(t.i in 1:(nrow(bond.prices.Vasicek.HJM.Kalman.errors))){
  # Filter Short Rate for Vasicek
  #KF_Vasicek_loop<-vasicek_KF(observations = bond.data[t.i,])
  r.predict<-as.vector(KF_Vasicek_loop$att)
  
  # Filter betas for dPCA
  dPCA_loop<-dPCA_KF(observations = bond.data[t.i,])
  dPCA.high_loop<-dPCA_KF(observations = bond.data[t.i,])
  AF_dPCA.high_loop<-AFReg_dPCA_KF(observations = bond.data[t.i,])
  
  
  # Write FRCs
  # Write Predicted FRC for Vasicek
  FRC.Vas.HJM.loop<-rep(NA,d)
  for(u in 1:d){
    FRC.Vas.HJM.loop[u]<-FRC.Vasicek(r.predict,1)
  }
  # Write for dPCA Models
  FRC.dPCA.loop<-as.vector(dPCA.samples[,-1]%*%(dPCA_loop$att))
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


# MSEs
MSE.dPCA.errors<-rbind(colMeans(bond.prices.Vasicek.HJM.Kalman.errors^2),
                     colMeans(bond.prices.dPCA.Kalman.error^2),
                     colMeans(bond.prices.dPCA_AF.high.Kalman.error^2))
rownames(MSE.dPCA.errors)<-c("Vsk","dPCA","AF-Reg(dPCA)")
round(MSE.dPCA.errors,3)

# Compute Error Statistics
mean.errors<-colMeans(bond.prices.Vasicek.HJM.Kalman.errors)
conf.bound<-(qnorm(0.95)/sqrt(N))*apply(bond.prices.Vasicek.HJM.Kalman.errors,2,sd)
error.statistics<-rbind(mean.errors-conf.bound,mean.errors,mean.errors+conf.bound)
rownames(error.statistics)<-c("L 95","Mean Err.","U 95")

# Report Findings
head(round(bond.prices.Vasicek.HJM.Kalman,3))
round(error.statistics,3)

if(FALSE){
  cat("\014")
  ## Write Tables
  xtable::xtable(x=rbind(MSE.dPCA.errors[,c(1:5)],MSE.dPCA.errors[,c(6:10)]),
                 caption="(Short): MSE Comparisons for 1-day ahead predictions.",
                 label="tab_AFReg_dPCA_Compare_short",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  #
  xtable::xtable(x=rbind(MSE.dPCA.errors[,c(11:15)],MSE.dPCA.errors[,c(16:20)]),
                 caption="(Mid): MSE Comparisons for 1-day ahead predictions.",
                 label="tab_AFReg_dPCA_Compare_mid",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  #
  xtable::xtable(x=rbind(MSE.dPCA.errors[,c(21:25)],MSE.dPCA.errors[,c(26:30)]),
                 caption="(Long): MSE Comparisons for 1-day ahead predictions.",
                 label="tab_AFReg_dPCA_Compare_long",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  xtable::xtable(x=as.matrix(MSE.dPCA.errors[,30]),
                 caption="(30 Year): MSE Comparisons for 1-day ahead predictions.",
                 label="tab_AFReg_dPCA_Compare_30year_long",
                 display=c("s",rep("e",1)),
                 digits=c(0,rep(3,1))
  )
}
###
# Difference Between NS, AFNS, AFReg(NS)
differences<-rbind(colMeans((bond.prices.NS.HJM.Kalman-bond.prices.AFNS.HJM.Kalman)/bond.prices.NS.HJM.Kalman),
                   colMeans((bond.prices.NS.HJM.Kalman-bond.prices.NS.AF.Reg.HJM.Kalman)/bond.prices.NS.HJM.Kalman))
colnames(differences)<-Maturities.Grid
rownames(differences)<-c("(NS-AFNS)/NS","(NS-AFReg(NS))/NS")

if(FALSE){
  cat("\014")
  ## Write Tables
  xtable::xtable(x=rbind(differences[,c(1:5)],differences[,c(6:10)]),
                 caption="(Short): Relative differences between NS and Abitrage-Free Versions for 1 day-ahead predictions.",
                 label="tab_AFReg_dPCA_Compare_difference_short",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  #
  xtable::xtable(x=rbind(differences[,c(11:15)],differences[,c(16:20)]),
                 caption="(Mid): Relative differences between NS and Abitrage-Free Versions for 1 day-ahead predictions.",
                 label="tab_AFReg_dPCA_Compare_differences_mid",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  #
  xtable::xtable(x=rbind(differences[-1,c(21:25)],differences[-1,c(26:30)]),
                 caption="(Long): Relative differences between NS and Abitrage-Free Versions for 1 day-ahead predictions.",
                 label="tab_AFReg_dPCA_Compare_differences_long",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  xtable::xtable(x=as.matrix(differences[-1,30]),
                 caption="(30 Year): Relative differences between NS and Abitrage-Free Versions for 1 day-ahead predictions.",
                 label="tab_AFReg_dPCA_Compare_30year_differences_long",
                 display=c("s",rep("e",1)),
                 digits=c(0,rep(3,1))
  )
}

############-#
############-#
############-#
############-#





# Reports
round(MSE.dPCA.errors[-3,]*10^4,10)
sum(apply(MSE.dPCA.errors[-3,],2,which.max)==2)/ncol(MSE.dPCA.errors)
sum(apply(MSE.dPCA.errors[-3,],2,which.min)==2)/ncol(MSE.dPCA.errors)
