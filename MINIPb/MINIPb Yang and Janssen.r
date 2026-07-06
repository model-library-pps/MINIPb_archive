#Implementation of MINIPb, the Yang & Janssen model for SOC breakdown
#Sources:
#Yang, H. S. and B. H. Janssen (2001). "A mono-component model of carbon mineralization with a dynamic rate constant." European Journal of Soil Science 51(3): 517-529.
#Yang, H. S. and B. H. Janssen (2002). "Relationship between substrate initial reactivity and residues ageing speed in carbon mineralization." Plant and Soil 239(2): 215-224.
#
#Examples are using data from Table 6.1, PhD thesis H.S. Yang, 1996
#Yang, H. S. (1996). Modelling organic matter mineralization and exploring options for organic matter management in arable farming in Northern China Proefschrift Wageningen, Yang.
#
#Algebraic solution to equilibrium values are derived from De Willigen et al., 1997 Appendix2:
# de Willigen, P., B. H. Janssen, H. I. M. Heesmans, J. G. Conijn, G. J. Velthof and W. J. Chardon (2008). 
# Decomposition and accumulation of organic matter in soil; comparison of some models. 
# Alterra-rapport 1726, http://edepot.wur.nl/15401 Wageningen, Alterra.
#
#The solution implemented here provides mean residence time and relative decay values for fresh and humified material ("_lit")
#
#AGT Schut, Wageningen University, May 2018


#Yang and Janssen model. 
#Assume that after 10000 yrs, all material is gone.
fCremaining<-function(t=time, Cin=Cin, R=R, S=S, Temp=Temp){
  
  #Determine amounts left after t years
  Ct<-YangJanssen(t=t, Cin=Cin, R=R,S=S,Temp=Temp)

  #Cin, Annual addition of C unit/ha
  #Algebraic solution to equilibrium values, from Appendix2, De Willigen et al., 1997
  #Is only valid when temperature correction is not needed!!
  
  #When using R=1.39 and S=0.64 (green manure), the SEQ=1.831, 
  #This matches with values in table 6 De Willigen et al. 1997 for Temp=9oC and R9==R
  f<-2^((Temp-9)/9)
  #Temperature is constant, so weighted average fm = f
  R9<-R / f^(S-1)
  
  CEQ<-Cin*(R9^(1/(S-1))*gamma(1/(1-S)))/(1-S)
  #Amount, in this case also the amount left after 1 year
  C1<-YangJanssen(t=1, Cin=Cin,R=R,S=S,Temp=Temp)$Ct[1]
  
  #MRT is time it takes to replace the pool or:
  #MRT=Ceq/Cin, here with annually 1 unit  additon, MRT=Ceq
  #In equilibrium, amount in equals amount out, so Cin=km*Ceq and:
  #Ceq=Cin/km and MRT=Ceq/Cin=1/km
  km<-Cin/CEQ
  MRT=1/km
  #SOC measured in the field is typically before application and hence represence the amount present after one year. 
  #Or with an Cin of 1, and the C(1) the amount left after one year i.e. the annual addition of SOC1a:
  #SOC1=Ceq-C1
  #Amount that is in the soil is sum of material left from previous years
  SOC1<-CEQ-C1
  
  #In equilibrium,:
  #SOC1=C1/km1, and km1=C1/SOC1
  km1=C1/SOC1
  MRT1=1/km1
  Ct[,"CEQ"]<-CEQ
  Ct[,"km"]<-km
  Ct[,"MRT"]<-MRT
  Ct[,"C1"]<-C1
  Ct[,"CEQ_lit"]<-SOC1
  Ct[,"km1"]<-km1
  Ct[,"MRT_lit"]<-MRT1
  return(Ct)
}
YangJanssen<-function(t=time, Cin=Cin, R=R,S=S,Temp=Temp){
  f<-2^((Temp-9)/9)
  #Temperature is constant, so weighted average fm = f
  R9<-R / f^(S-1)
  logK9<-log10(R9)-S*log10(f*t)
  K9<-10^logK9
  Ct=Cin*exp(-R9*t^(1-S))
  k<-(1-S)*R9*t^(-S)
  Ct<-data.frame(t=t,f=f,k=k,Ct=Ct)
  return(Ct)
}
#SOC buildup, with specified annual inputs
SoilSOC<-function(SOCini=0,Rsoc=R,Ssoc=S, CinAnn=CinAnn, R=R,S=S,Temp=Temp){
  ncomp<-ncol(CinAnn)
  ny<-nrow(CinAnn)
  components<-seq(1,ncomp,1)
  years<-seq(1,ny,1)
  SoilC<-fCremaining(t=years, Cin=SOCini, R=Rsoc,S=Ssoc,Temp=Temp)$Ct
  for(y in years){
    soilc<-0*years
    for(comp in components){
      soilc[y:ny]<-soilc[y:ny]+fCremaining(t=seq(1,ny+1-y,1), Cin=CinAnn[y,comp], R=R[comp],S=S[comp],Temp=Temp)$Ct
    }
    SoilC<-rbind(SoilC, soilc)
  }
    
  SC<-colSums(SoilC)
  return(list(YEARS=years, SOC=SC, CONTRIBUTIONS=SoilC))
}
SoilSOCFuture<-function(SOCini=0,Rsoc=R,Ssoc=S, CinAnn=CinAnn, R=R,S=S,Temp=Temp,Years=Years){
  #Determines SOC for a moment in the future after annually applying same inputs 
  nscenarios<-nrow(CinAnn)
  eSOC=NULL
  for(ns in seq(1,nscenarios,1)){
    CinAnnR<-NULL
    for(y in seq(1,Years,1)){
     CinAnnR<-rbind(CinAnnR,CinAnn)
    }
    tSOC<-SoilSOC(SOCini=SOCini,Rsoc=Rsoc,Ssoc=Ssoc, CinAnn=CinAnnR, R=R,S=S,Temp=Temp)
    eSOC<-rbind(eSOC,tSOC$SOC[length(tSOC$SOC)])
  }
  return(eSOC)
}


#=========================================================================================================================
#SOME EXAMPLES
#Annual C input into soil with straw and roots
#=========================================================================================================================

#SOM + annual addition of straw and roots. All in t/ha
CinAnn<-data.frame(Straw=rep(5,1000),
                   Roots=rep(1,1000))
RD<-0.2
#t topsoil / ha
SoilWeight<-(RD*1.25*1E4)
tSOC9<-SoilSOC(SOCini=0.01*SoilWeight,Rsoc=0.057,Ssoc=0.46, CinAnn=CinAnn,R=c(1.11,0.8),S=c(0.66,0.67),Temp=9)
tSOC21<-SoilSOC(SOCini=0.01*SoilWeight,Rsoc=0.057,Ssoc=0.46, CinAnn=CinAnn,R=c(1.11,0.8),S=c(0.66,0.67),Temp=21)
layout(matrix(c(1,2,3,4)))
par(mfrow=c(2,2),mar=c(0.5,4,2,1))
plot(tSOC9$YEARS,tSOC9$SOC,type="l",main="1% SOC with annual additions of 1 t root + 5 t straw-C/ha", 
     xlim=c(1,1000),ylim=c(0,75),xlab="Time, y",ylab="SOC, t/ha")
lines(tSOC21$YEARS,tSOC21$SOC,col="red")
plot(tSOC9$YEARS,100*tSOC9$SOC/SoilWeight,type="l",main="Annual additions of 1 t root + 5 t straw-C/ha", 
     xlim=c(1,1000),ylim=c(0,2.5),xlab="Time, y",ylab="SOC in 0.2 m topsoil, %")
lines(tSOC21$YEARS,100*tSOC21$SOC/SoilWeight,col="red")



#=========================================================================================================================
#SOME EXAMPLES
#Values for k (relative decay rate), mean residence times of C in a range of scenarios varying in time, input type and temperature
#=========================================================================================================================

#Data from Table 6.1, PhD thesis H.S. Yang, 1996
#Straw, MRT fresh=4.10 y ,humific=0.33, MRT of humified material = 11.45 y
fCremaining(t=c(1,5,10,20,100), Cin=1, R=1.11,S=0.66,Temp=9)
#Straw, MRT of humified material = 11.45, 6.84, 4.33, 2.92 years
fCremaining(t=1, Cin=1,R=1.11,S=0.66,Temp=c(9,15,20,24))
#Green manure, MRT fresh=1.83,humific=0.25, MRT of humified material = 6.35 
fCremaining(t=c(1,5,10,20,100), Cin=1,R=1.39,S=0.64,Temp=9)
#Roots, MRT fresh=12.26,humific=0.45, MRT of humified material = 26.28
fCremaining(t=c(1,5,10,20,100), Cin=1,R=0.8,S=0.67,Temp=9)
#FYM, MRT  fresh=2.85,humific=0.44, MRT of humified material = 5.46
fCremaining(t=c(1,5,10,20,100), Cin=1,R=0.82,S=0.49,Temp=9)
#SOM, MRT  = 352.8 years, k=0.0307, annual addition of 1 t/ha gives 352.8 tC/ha in equill.
fCremaining(t=c(1,5,10,20,100), Cin=1,R=0.057,S=0.46,Temp=9)
fCremaining(t=c(1,5,10,20,100), Cin=1,R=0.057,S=0.46,Temp=16)
fCremaining(t=c(1,5,10,20,100), Cin=1,R=0.057,S=0.46,Temp=21)


#=========================================================================================================================
#SOME EXAMPLES
#Comparing differences between materials
#=========================================================================================================================

layout(matrix(c(1,2,3,4)))
par(mfrow=c(2,2),mar=c(0.5,4,2,1))

x<-seq(0.3,50,0.1)
plot(x,fCremaining(t=x, Cin=1,R=2.04,S=0.86,Temp=9)$k,ylim=c(0,0.2),
     xlab="Time, y", ylab="k",type="l",lty=1,lwd=2.0,col="blue",main="T=9oC")
lines(x,fCremaining(t=x, Cin=1,R=1.5,S=0.83,Temp=9)$k,col="red",lty=2,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=1.2,S=0.8,Temp=9)$k,col="brown",lty=3,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=1.15,S=0.75,Temp=9)$k,col="green",lty=5,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=0.2,S=0.7,Temp=9)$k,col="black",lty=5,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=0.4,S=0.7,Temp=9)$k,col="black",lty=5,lwd=2.0)
legend("topright",legend=c("Glucose","cellulose","Stalks/straw","Ryegrass","FYM"),
       lty=c(1,2,3,5),col=c("blue","red","brown","green","black"),lwd=c(2.0,2.0,2.0,2.0))

par(mfg=c(2,1),mar=c(4,4,0.5,1))
plot(x,fCremaining(t=x, Cin=1,R=2.04,S=0.86,Temp=9)$Ct,ylim=c(0,0.8),
     xlab="Time, y", ylab="SOC remaining,kg/kg",lwd=2.0,type="l",lty=1,col="blue")
lines(x,fCremaining(t=x, Cin=1,R=1.5,S=0.83,Temp=9)$Ct,col="red",lty=2,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=1.2,S=0.8,Temp=9)$Ct,col="brown",lty=3,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=1.15,S=0.75,Temp=9)$Ct,col="green",lty=4,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=0.2,S=0.7,Temp=9)$Ct,col="black",lty=5,lwd=2.0)

par(mfg=c(1,2),mar=c(0.5,4,2,1))
plot(x,fCremaining(t=x, Cin=1,R=2.04,S=0.86,Temp=18)$k,ylim=c(0,0.2),
     xlab="Time, y", ylab="k",type="l",lty=1,lwd=2.0,col="blue",main="T=20oC")
lines(x,fCremaining(t=x, Cin=1,R=1.5,S=0.83,Temp=18)$k,col="red",lty=2,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=1.2,S=0.8,Temp=18)$k,col="brown",lty=3,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=1.15,S=0.75,Temp=18)$k,col="green",lty=5,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=0.2,S=0.7,Temp=18)$k,col="black",lty=5,lwd=2.0)

par(mfg=c(2,2),mar=c(4,4,0.5,1))
plot(x,fCremaining(t=x, Cin=1,R=2.04,S=0.86,Temp=18)$Ct,ylim=c(0,0.8),xlab="Time, y", ylab="SOC remaining,kg/kg",lty=1,lwd=2.0,type="l",col="blue")
lines(x,fCremaining(t=x, Cin=1,R=1.5,S=0.83,Temp=18)$Ct,col="red",lty=2,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=1.2,S=0.8,Temp=18)$Ct,col="brown",lty=3,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=1.15,S=0.75,Temp=18)$Ct,col="green",lty=5,lwd=2.0)
lines(x,fCremaining(t=x, Cin=1,R=0.2,S=0.7,Temp=18)$Ct,col="black",lty=5,lwd=2.0)

