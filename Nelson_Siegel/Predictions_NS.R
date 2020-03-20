# Initialize one-day ahead predicted bond prices!
bond.prices.Vasicek.HJM.Kalman<-matrix(data=NA,nrow=N,ncol=d)
colnames(bond.prices.Vasicek.HJM.Kalman)<-Maturities.Grid
#rownames(bond.prices.Vasicek.HJM.Kalman)<-round(dates[-1],1)
bond.prices.dPCA.high.Kalman<-bond.prices.dPCA.Kalman<-bond.prices.NS.AF.Reg.HJM.Kalman<-bond.prices.AFNS.HJM.Kalman<-bond.prices.NS.HJM.Kalman<-bond.prices.Vasicek.HJM.Kalman
# Initialize Bond Pricing Errors
bond.prices.dPCA.high.Kalman.error<-bond.prices.dPCA.Kalman.error<-bond.prices.NS.AF.Reg.HJM.Kalman.errors<-bond.prices.AFNS.HJM.Kalman.errors<-bond.prices.NS.HJM.Kalman.errors<-bond.prices.Vasicek.HJM.Kalman.errors<-bond.prices.Vasicek.HJM.Kalman


# Populate Matrix of one-day ahead predicted bond prices
for(t.i in 1:(nrow(bond.prices.Vasicek.HJM.Kalman.errors))){
  # Filter Short Rate for Vasicek
  #KF_Vasicek_loop<-vasicek_KF(observations = bond.data[t.i,])
  r.predict<-as.vector(KF_Vasicek_loop$att)
  
  # Filter betas for dPCA
  #dPCA_loop<-dPCA_KF(observations = bond.data[t.i,])
  #dPCA.high_loop<-dPCA_KF(observations = bond.data[t.i,])
  
  # Filter betas for NS
  KF_NS_loop<-NS_KF(observations = bond.data[t.i,])
  # Filter betas for AF.Reg(NS)
  AF.Reg.KF_NS_loop<-AFReg_NS_KF(observations = bond.data[t.i,])
  
  
  # Write FRCs
  # Write Predicted FRC for Vasicek
  FRC.Vas.HJM.loop<-rep(NA,d)
  for(u in 1:d){
    FRC.Vas.HJM.loop[u]<-FRC.Vasicek(r.predict,1)
  }
  # Write Predicted FRC for NS Models
  FRC.NS.loop<-as.vector(NS.samples[,-1]%*%(KF_NS_loop$att))
  FRC.AFNS.loop<-FRC.NS.loop
  FRC.NS.AF.Reg.loop<-as.vector(AF.Reg.Factors%*%(AF.Reg.KF_NS_loop$att))
  
  
  # Write Bond prices
  prices.Vas.loop<-exp(-cumsum(FRC.Vas.HJM.loop))
  prices.NS.loop<-exp(-cumsum(FRC.NS.loop))
  prices.AFNS.loop<-rep(NA,length(prices.NS.loop))
    # Write-in Analytical Adjustment Term for AFNS
    for(adj.T.i in 1:length(Maturities.Grid)){
      prices.AFNS.loop[adj.T.i]<-AFNS.yield.adjustment.term(time = Maturities.Grid[adj.T.i])
    }
  prices.AFNS.loop<-prices.NS.loop*exp(cumsum(prices.AFNS.loop))# correct
  prices.NS.AF.Reg.loop<-exp(-cumsum(FRC.NS.AF.Reg.loop))
  
  # Write Predicted Bond Prices
  bond.prices.Vasicek.HJM.Kalman[t.i,]<-prices.Vas.loop
  bond.prices.NS.HJM.Kalman[t.i,]<-prices.NS.loop
  bond.prices.AFNS.HJM.Kalman[t.i,]<-prices.AFNS.loop
  bond.prices.NS.AF.Reg.HJM.Kalman[t.i,]<-prices.NS.AF.Reg.loop
  
  # Evaluate Prediction Errors (1-day ahead)
  # Write realised next day bond price
  day.ahead.bond.price<-unlist(as.vector(exp(-(Maturities.Grid)*(bond.data[t.i,]))))
  # Evaluate Error
  bond.prices.Vasicek.HJM.Kalman.errors[t.i,]<- (prices.Vas.loop-day.ahead.bond.price)
  bond.prices.NS.HJM.Kalman.errors[t.i,]<- (prices.NS.loop-day.ahead.bond.price)
  bond.prices.AFNS.HJM.Kalman.errors[t.i,]<- (prices.AFNS.loop-day.ahead.bond.price)
  bond.prices.NS.AF.Reg.HJM.Kalman.errors[t.i,]<- (prices.NS.AF.Reg.loop-day.ahead.bond.price)
  
}# END Prediction Step


# MSEs
MSE.NS.errors<-rbind(colMeans(bond.prices.Vasicek.HJM.Kalman.errors^2),
                     colMeans(bond.prices.NS.HJM.Kalman.errors^2),
                     colMeans(bond.prices.AFNS.HJM.Kalman.errors^2),
                     colMeans(bond.prices.NS.AF.Reg.HJM.Kalman.errors^2))
rownames(MSE.NS.errors)<-c("Vsk","NS","AFNS","AF-Reg(NS)")
round(MSE.NS.errors,3)

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
  xtable::xtable(x=rbind(MSE.NS.errors[,c(1:5)],MSE.NS.errors[,c(6:10)]),
                 caption="(Short): MSE Comparisons for 1-day ahead predictions.",
                 label="tab_AFReg_NS_Compare_short",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  #
  xtable::xtable(x=rbind(MSE.NS.errors[,c(11:15)],MSE.NS.errors[,c(16:20)]),
                 caption="(Mid): MSE Comparisons for 1-day ahead predictions.",
                 label="tab_AFReg_NS_Compare_mid",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  #
  xtable::xtable(x=rbind(MSE.NS.errors[,c(21:25)],MSE.NS.errors[,c(26:30)]),
                 caption="(Long): MSE Comparisons for 1-day ahead predictions.",
                 label="tab_AFReg_NS_Compare_long",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  xtable::xtable(x=as.matrix(MSE.NS.errors[,30]),
                 caption="(30 Year): MSE Comparisons for 1-day ahead predictions.",
                 label="tab_AFReg_NS_Compare_30year_long",
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
                 label="tab_AFReg_NS_Compare_difference_short",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  #
  xtable::xtable(x=rbind(differences[,c(11:15)],differences[,c(16:20)]),
                 caption="(Mid): Relative differences between NS and Abitrage-Free Versions for 1 day-ahead predictions.",
                 label="tab_AFReg_NS_Compare_differences_mid",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  #
  xtable::xtable(x=rbind(differences[-1,c(21:25)],differences[-1,c(26:30)]),
                 caption="(Long): Relative differences between NS and Abitrage-Free Versions for 1 day-ahead predictions.",
                 label="tab_AFReg_NS_Compare_differences_long",
                 display=c("s",rep("e",5)),
                 digits=c(0,rep(3,5))
  )
  xtable::xtable(x=as.matrix(differences[-1,30]),
                 caption="(30 Year): Relative differences between NS and Abitrage-Free Versions for 1 day-ahead predictions.",
                 label="tab_AFReg_NS_Compare_30year_differences_long",
                 display=c("s",rep("e",1)),
                 digits=c(0,rep(3,1))
  )
}

############-#
############-#
############-#
############-#





# Reports
round(MSE.NS.errors[-3,]*10^4,10)
sum(apply(MSE.NS.errors[-3,],2,which.max)==2)/ncol(MSE.NS.errors)
sum(apply(MSE.NS.errors[-3,],2,which.min)==2)/ncol(MSE.NS.errors)
