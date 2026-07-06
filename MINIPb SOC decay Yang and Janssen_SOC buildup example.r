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
  #Ceq=Cin/km
  km<-Cin/CEQ
  MRT=Cin/km
  #SOC measured in the field is typically before application and hence represence the amount present after one year. 
  #Or with an Cin of 1, and the C(1) the amount left after one year i.e. the annual addition of SOC1a:
  #SOC1=Ceq-C1
  #Amount that is in the soil is sum of material left from previous years
  SOC1<-CEQ-C1
  
  #In equilibrium,:
  #SOC1=C1/km1, and km1=C1/SOC1
  km1=C1/SOC1
  MRT1=Cin/km1
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
  #R9 is the temperature corrected (relative) decomposition rate in the first year
  #S is the speed of aging of the material
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
  SoilC<-fCremaining(t=c(0,years), Cin=SOCini, R=Rsoc,S=Ssoc,Temp=Temp)$Ct
  for(y in years){
    soilc<-0*years
    for(comp in components){
      soilc[y:ny]<-soilc[y:ny]+fCremaining(t=seq(1,ny+1-y,1), Cin=CinAnn[y,comp], R=R[comp],S=S[comp],Temp=Temp)$Ct
    }
    SoilC<-rbind(SoilC, c(0,soilc))
  }
    
  SC<-colSums(SoilC)
  return(list(YEARS=c(0,years), SOC=SC, CONTRIBUTIONS=SoilC))
}

#SOM + annual addition of straw and roots
  CinAnn<-data.frame(Straw=rep(5,1000),
                     Roots=rep(2,1000))
  Straw9<-fCremaining(t=1, Cin=5, R=1.11,S=0.66,Temp=9)
  Roots9<-fCremaining(t=1, Cin=2,R=0.8,S=0.67,Temp=9)
  Straw21<-fCremaining(t=1, Cin=5, R=1.11,S=0.66,Temp=21)
  Roots21<-fCremaining(t=1, Cin=2,R=0.8,S=0.67,Temp=21)
  
  fSOC9<-SoilSOC(SOCini=0.01*(0.2*1.25*1E4),Rsoc=0.057,Ssoc=0.46, CinAnn=data.frame(Straw=rep(0,1000),Roots=rep(0,1000)),R=c(1.11,0.8),S=c(0.66,0.67),Temp=9)
  tSOC9<-SoilSOC(SOCini=0.01*(0.2*1.25*1E4),Rsoc=0.057,Ssoc=0.46, CinAnn=CinAnn,R=c(1.11,0.8),S=c(0.66,0.67),Temp=9)
  fSOC21<-SoilSOC(SOCini=0.01*(0.2*1.25*1E4),Rsoc=0.057,Ssoc=0.46, CinAnn=data.frame(Straw=rep(0,1000),Roots=rep(0,1000)),R=c(1.11,0.8),S=c(0.66,0.67),Temp=21)
  tSOC21<-SoilSOC(SOCini=0.01*(0.2*1.25*1E4),Rsoc=0.057,Ssoc=0.46, CinAnn=CinAnn,R=c(1.11,0.8),S=c(0.66,0.67),Temp=21)
  t9<-SoilSOC(SOCini=0,Rsoc=0.057,Ssoc=0.46, CinAnn=CinAnn,R=c(1.11,0.8),S=c(0.66,0.67),Temp=9)
  t21<-SoilSOC(SOCini=0,Rsoc=0.057,Ssoc=0.46, CinAnn=CinAnn,R=c(1.11,0.8),S=c(0.66,0.67),Temp=21)

#Determine amount in equilibrium for different inputs of straw and 
  StrawEQ9<-data.frame(AnnIn=seq(0,10,0.5),SOCeq=0)
  RootsEQ9<-data.frame(AnnIn=seq(0,10,0.5),SOCeq=0)
  StrawEQ21<-StrawEQ9
  RootsEQ21<-RootsEQ9
  StrawEQ9[,"SOCeq"]<-fCremaining(t=1, Cin=StrawEQ9[,"AnnIn"], R=1.11,S=0.66,Temp=9)["CEQ"]
  RootsEQ9[,"SOCeq"]<-fCremaining(t=1, Cin=RootsEQ9[,"AnnIn"], R=0.8,S=0.67,Temp=9)["CEQ"]
  StrawEQ21[,"SOCeq"]<-fCremaining(t=1, Cin=StrawEQ21[,"AnnIn"], R=1.11,S=0.66,Temp=21)["CEQ"]
  RootsEQ21[,"SOCeq"]<-fCremaining(t=1, Cin=RootsEQ21[,"AnnIn"], R=0.8,S=0.67,Temp=21)["CEQ"]

  
    
    
pdf("..\\Results\\MINIPb_predicted SOC for annual additions of 2tCroots.ha roots and 5tCstover.ha.pdf")
  layout(matrix(c(1,2,3,4)))
  par(mfrow=c(2,2),mar=c(0.5,4,2,0.5))
  plot(tSOC9$YEARS,tSOC9$SOC,type="l",col="black",main="Annually 2 tC/ha root and 5 tC/ha straw", 
       xlim=c(1,1000),ylim=c(0,60),xlab="Time, y",ylab="SOC, t/ha")
  lines(fSOC9$YEARS,fSOC9$SOC,col="black",lty=5)
  lines(tSOC9$YEARS,rep(Roots9$CEQ_lit+Straw9$CEQ_lit,length(tSOC9$YEARS)),col="black",lty=3)
  lines(tSOC21$YEARS,tSOC21$SOC,col="red")
  lines(fSOC21$YEARS,fSOC21$SOC,col="red",lty=5)
  lines(tSOC21$YEARS,rep(Roots21$CEQ_lit+Straw21$CEQ_lit,length(tSOC21$YEARS)),col="red",lty=3)
  
  par(mfg=c(1,2),mar=c(0.5,4,2,0.5))
  plot(tSOC9$YEARS,100*tSOC9$SOC/(0.2*1.25*1E4),type="l",col="black",main="Annually 2 tC/ha root and 5 tC/ha straw", 
       xlim=c(1,1000),ylim=c(0,4),xlab="Time, y",ylab="SOC in 0.2 m topsoil, %")
  lines(tSOC9$YEARS,rep(100*(Roots9$CEQ_lit+Straw9$CEQ_lit)/(0.2*1.25*1E4),length(tSOC9$YEARS)),col="black",lty=3)
  lines(fSOC9$YEARS,100*fSOC9$SOC/(0.2*1.25*1E4),col="black",lty=5)
  lines(tSOC21$YEARS,100*tSOC21$SOC/(0.2*1.25*1E4),col="red")
  lines(fSOC21$YEARS,100*fSOC21$SOC/(0.2*1.25*1E4),col="red",lty=5)
  lines(tSOC21$YEARS,rep(100*(Roots21$CEQ_lit+Straw21$CEQ_lit)/(0.2*1.25*1E4),length(tSOC21$YEARS)),col="red",lty=3)
  legend("topleft",legend=c("9dC SOC+7tC/ha/y","21dC SOC+7tC/ha/y","9dC Equil.","21dC Equil.",
                            "9dC SOC+0tC/ha/y","21dC SOC+0tC/ha/y"),col=c("black","red","black","red","black","red"),lty=c(1,1,3,3,5,5))
  
  par(mfg=c(2,1),mar=c(4,4,0.5,0.5))
  plot(t9$YEARS,t9$SOC,type="l",col="black", 
       xlim=c(1,1000),ylim=c(0,50),xlab="Time, y",ylab="SOC, t/ha")
  lines(t9$YEARS,rep(Roots9$CEQ_lit+Straw9$CEQ_lit,length(t9$YEARS)),col="black",lty=3)
  lines(t21$YEARS,t21$SOC,col="red")
  lines(t21$YEARS,rep(Roots21$CEQ_lit+Straw21$CEQ_lit,length(t21$YEARS)),col="red",lty=3)
  
  par(mfg=c(2,2),mar=c(4,4,0.5,0.5))
  plot(t9$YEARS,100*t9$SOC/(0.2*1.25*1E4),type="l",col="black", 
       xlim=c(1,1000),ylim=c(0,2),xlab="Time, y",ylab="SOC in 0.2 m topsoil, %")
  lines(t9$YEARS,rep(100*(Roots9$CEQ_lit+Straw9$CEQ_lit)/(0.2*1.25*1E4),length(t9$YEARS)),col="black",lty=3)
  lines(t21$YEARS,100*t21$SOC/(0.2*1.25*1E4),col="red")
  lines(t21$YEARS,rep(100*(Roots21$CEQ_lit+Straw21$CEQ_lit)/(0.2*1.25*1E4),length(t21$YEARS)),col="red",lty=3)
  
  
  layout(matrix(c(1,2,3,4)))
  par(mfrow=c(2,2),mar=c(0.5,4,2,0.5))
  
  plot(StrawEQ9[,"AnnIn"],StrawEQ9[,"SOCeq"],col="black",type="l",lty=2,
       xlim=c(0,10),ylim=c(0,200),xlab="Annual input, t C/ha/y",ylab="SOC in equilibrium, t C/ha")  
  lines(RootsEQ9[,"AnnIn"],RootsEQ9[,"SOCeq"],col="black",lty=3)  
  lines(RootsEQ9[,"AnnIn"],RootsEQ9[,"SOCeq"]*1/3+StrawEQ9[,"SOCeq"]*2/3,col="black")  
  lines(StrawEQ21[,"AnnIn"],StrawEQ21[,"SOCeq"],col="red",lty=2)  
  lines(RootsEQ21[,"AnnIn"],RootsEQ21[,"SOCeq"],col="red",lty=3)  
  lines(RootsEQ21[,"AnnIn"],RootsEQ21[,"SOCeq"]/3 +StrawEQ21[,"SOCeq"]*2/3,col="red")  
  legend("topleft",legend=c("9oC roots","9oC stover","9oC 1/3 roots + 2/3 stover",
                            "21oC roots","21oC stover","21oC 1/3 roots + 2/3 stover"),
         col=c("black","black","black","red","red","red"),lty=c(2,3,1,2,3,1))
  
  par(mfg=c(1,2),mar=c(0.5,4,2,0.5))
  plot(StrawEQ9[,"AnnIn"],100*StrawEQ9[,"SOCeq"]/(0.2*1.25*1E4),col="black",type="l",lty=2,
       xlim=c(0,10),ylim=c(0,5),xlab="Annual input, t C/ha/y",ylab="SOC in equilibrium, %")  
  lines(RootsEQ9[,"AnnIn"],100*RootsEQ9[,"SOCeq"]/(0.2*1.25*1E4),col="black",lty=3)  
  lines(RootsEQ9[,"AnnIn"],100*(RootsEQ9[,"SOCeq"]*1/3+StrawEQ9[,"SOCeq"]*2/3)/(0.2*1.25*1E4),col="black")  
  lines(StrawEQ21[,"AnnIn"],100*StrawEQ21[,"SOCeq"]/(0.2*1.25*1E4),col="red",lty=2)  
  lines(RootsEQ21[,"AnnIn"],100*RootsEQ21[,"SOCeq"]/(0.2*1.25*1E4),col="red",lty=3)  
  lines(RootsEQ21[,"AnnIn"],100*(RootsEQ21[,"SOCeq"]/3 +StrawEQ21[,"SOCeq"]*2/3)/(0.2*1.25*1E4),col="red")  
  legend("topleft",legend=c("9oC roots","9oC stover","9oC 1/3 roots + 2/3 stover",
                            "21oC roots","21oC stover","21oC 1/3 roots + 2/3 stover"),
         col=c("black","black","black","red","red","red"),lty=c(2,3,1,2,3,1))
  
  
dev.off()
  

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
fCremaining(t=c(1,5,10,20,100), Cin=1,R=0.057,S=0.46,Temp=21)
