### AFNS of Christiansen et Al

AFNS.yield.adjustment.term<-function(sigma.mat=(Estimated.OU.Parameters.Betas$Sigma),NS.shape.param=lambda_init,time){
  
# Write components in Christiansen et Al.'s notation:  
 Atilde<-sigma.mat[1,1]^2+sigma.mat[1,2]^2+sigma.mat[1,3]^2

 Btilde<-sigma.mat[2,1]^2+sigma.mat[2,2]^2+sigma.mat[2,3]^2

 Ctilde<-sigma.mat[3,1]^2+sigma.mat[3,2]^2+sigma.mat[3,3]^2

 Dtilde<-sigma.mat[1,1]*sigma.mat[2,1]+sigma.mat[1,2]*sigma.mat[2,2]+sigma.mat[1,3]*sigma.mat[2,3]

 Etilde<-sigma.mat[1,1]*sigma.mat[3,1]+sigma.mat[1,2]*sigma.mat[3,2]+sigma.mat[1,3]*sigma.mat[3,3]

 Ftilde<-sigma.mat[2,1]*sigma.mat[3,1]+sigma.mat[2,2]*sigma.mat[3,2]+sigma.mat[2,3]*sigma.mat[3,3]

 # Compute each component of the adjustment term:
 adj.term1<--Atilde*time^2/6

 adj.term2<--Btilde*(1/(2*NS.shape.param^2)-(1-exp(-NS.shape.param*time))/(NS.shape.param^3*time)+(1-exp(-2*NS.shape.param*time))/(4*NS.shape.param^3*time))

 adj.term3<--Ctilde*(1/(2*NS.shape.param^2)+exp(-NS.shape.param*time)/(NS.shape.param^2)-time*exp(-2*NS.shape.param*time)/(4*NS.shape.param)-3*exp(-2*NS.shape.param*time)/(4*NS.shape.param^2)-2*(1-exp(-NS.shape.param*time))/(NS.shape.param^3*time)+5*(1-exp(-2*NS.shape.param*time))/(8*NS.shape.param^3*time))

 adj.term4<--Dtilde*(time/(2*NS.shape.param)+exp(-NS.shape.param*time)/(NS.shape.param^2)-(1-exp(-NS.shape.param*time))/(NS.shape.param^3*time))

 adj.term5<--Etilde*(3*exp(-NS.shape.param*time)/(NS.shape.param^2)+time/(2*NS.shape.param)+time*exp(-NS.shape.param*time)/(NS.shape.param)-3*(1-exp(-NS.shape.param*time))/(NS.shape.param^3*time))

 adj.term6<--Ftilde*(1/(NS.shape.param^2)+exp(-NS.shape.param*time)/(NS.shape.param^2)-exp(-2*NS.shape.param*time)/(2*NS.shape.param^2)-3*(1-exp(-NS.shape.param*time))/(NS.shape.param^3*time)+3*(1-exp(-2*NS.shape.param*time))/(4*NS.shape.param^3*time))

 # Sum terms making up yield Curve
 adjustment.out<-adj.term1+adj.term2+adj.term3+adj.term4+adj.term5+adj.term6
 # Return sum
 return(as.numeric(adjustment.out))
}
